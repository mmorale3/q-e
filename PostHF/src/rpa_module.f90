!--------------------------------------------------------------------
! Written by Miguel A. Morales, LLNL, 2020 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE rpa_module
  !----------------------------------------------------------------------
  ! 
  ! Implements RPA + SOSEX in the auxiliary-basis formulation. 
  ! Assumes that Cholesky factorization of the 2-electron integrals has 
  ! has already been calculated. See the cholesky routine in twobody_hamiltonian.f90  
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY: omega, alat, tpiba, tpiba2, at, bg
  USE posthf_mod, ONLY: nksym,norb,numspin,e0,efc,e2Ha,&
                        nelec_tot,nup,ndown,ke_factorization,&
                        xksym,igksym,ngksym,QKtoK2,kminus,Qpts,DM_mf,nmax_DM, &
                        nQuniq,wQ,xQ
  USE read_orbitals_from_file, ONLY: open_esh5_read,close_esh5_read,h5file_type
  USE gvecw, ONLY : ecutwfc
  USE wavefunctions, ONLY : evc, psic
  USE control_flags, ONLY : gamma_only
  USE gvect, ONLY: ngm, ngm_g, g, gstart, gg
  USE io_files, ONLY: nwordwfc, iunwfc, tmp_dir, prefix
  USE io_global, ONLY: stdout, ionode,  ionode_id
  USE wvfct, ONLY: nbnd, npwx, nbnd, g2kin
  USE klist,  ONLY : nkstot, wk, nks, xk, ngk, igk_k
  USE mp,           ONLY: mp_sum, mp_max, mp_bcast, mp_barrier
  USE mp_images, ONLY: intra_image_comm, me_image, root_image, nproc_image
  USE mp_pools,     ONLY: inter_pool_comm, intra_pool_comm, npool, nproc_pool, &
                          me_pool,root_pool,my_pool_id
  USE noncollin_module,     ONLY : noncolin, npol
  USE lsda_mod, ONLY: lsda, nspin
  USE read_orbitals_from_file, ONLY: get_orbitals_set
  use fft_interfaces,       ONLY : invfft, fwfft
  USE fft_types, ONLY: fft_type_descriptor
  USE orbital_generators, ONLY: get_noccK
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: verbose
  CHARACTER(6), PARAMETER  :: orbsG = 'OrbsG'
  CHARACTER(6), PARAMETER  :: orbsR = 'OrbsR'
  !
  CONTAINS
  !
  ! hamil_esh5 must contain Cholesky matrix, eigenvalues and mf density matrix 
  ! right now assumes a HF reference. Generalize to DFT reference with singles
  !    term later
  ! 
  SUBROUTINE rpa(edrpa,dfft,hamil_esh5,esosex) 
    !
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE gvect, ONLY : ecutrho
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_esh5
    COMPLEX(DP), INTENT(OUT) :: edrpa
    COMPLEX(DP), INTENT(IN), OPTIONAL :: esosex
    !
    TYPE(h5file_type) :: h5id_hamil
    !
    LOGICAL :: get_sosex
    INTEGER :: iq, Q, ka, ki, ia, ii, ispin, ik, ikk, iw, iwQ 
    INTEGER :: P_beg, P_end, Q_beg, Q_end, Naux, nP, nQ
    INTEGER :: k_beg, k_end, nkloc, nia_max
    INTEGER :: nfreq, nQuniq, maxnorbs
    !
    INTEGER :: nel(2), maxocc, maxvir
    INTEGER, ALLOCATABLE :: noccK(:,:), nvirK(:,:)
    INTEGER, ALLOCATABLE :: norbK(:), ncholQ(:)
    REAL(DP), ALLOCATABLE :: weight(:,:,:),eigval(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Qpq(:,:,:) ! Q matrix 
    COMPLEX(DP), ALLOCATABLE :: Ypq(:,:,:) ! Y matrix 
    COMPLEX(DP), ALLOCATABLE :: MOs(:,:,:) ! Matrix of MO coefficients  
    COMPLEX(DP), ALLOCATABLE :: Chol(:,:) ! Cholesky matrix in SPO basis 
    COMPLEX(DP), ALLOCATABLE :: CholT(:,:) ! Cholesky matrix in SPO basis 
    COMPLEX(DP), ALLOCATABLE :: T1(:) ! Work space
    COMPLEX(DP), ALLOCATABLE :: T2(:,:) ! Work space 
    COMPLEX(DP), ALLOCATABLE :: Lia(:,:) ! left hand side of MO Cholesky matrix 
    COMPLEX(DP), ALLOCATABLE :: Ria(:,:) ! right hand side of MO Cholesky matrix 
    COMPLEX(DP), ALLOCATABLE :: xfreq(:), wfreq(:)

    get_sosex = present(esosex)
    CALL esh5_posthf_open_file_read(h5id_hamil%id,TRIM(hamil_esh5),  &
                                                  LEN_TRIM(hamil_esh5),error)
    if(error .ne. 0 ) &
      call errore('rpa','error opening hamil file',1)

    maxnorbs = h5id_hamil%maxnorbs 
    allocate( MOs(maxnorbs, maxnorbs, 2), ncholQ(nksym), norbK(nksym) )

    call esh5_posthf_read_nchol(h5id_hamil%id,ncholQ, error)
    if(error .ne. 0) call errore('rpa','Error: esh5_posthf_read_nchol',error)
    nchol_max = maxval(ncholQ(:))
    call esh5_posthf_read_norbk(h5id_hamil%id,norbK, error)
    if(error .ne. 0) call errore('rpa','Error: esh5_posthf_read_norbk',error)

    !  Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    ! Partition M*M 
    call find_2d_partition(nchol_max,me_pool+1,nproc_pool,P_beg,P_end,Q_beg,Q_end)
    nP = P_end-P_beg+1
    nQ = Q_end-Q_beg+1

    ! 1. Read eigenvalues, fermi weights and determinte occupied/virtual partitioning
    allocate( noccK(nksym,nspin), nvirK(nksym,nspin) )
    allocate(eigval(maxnorb,nksym*nspin),weight(maxnorb,nksym*nspin))
    weight(:,:) = 0.d0
    eigval(:,:) = 0.d0
    noccK(:,:)=0
    nvirK(:,:)=0
    do ispin=1,nspin
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        call esh5_posthf_read_et(h5id_hamil%id,ikk-1,eigval(1,ikk), &
                               weight(1,ikk),error)
        if(error .ne. 0 ) &
          call errore('mp2_g','error reading weights',1)
      enddo
    enddo
    ! find number of electrons
    call get_noccK(noccK,nel,maxnorb,nksym,nspin,weight,maxnorb)
    write(*,*) ' Number of electrons per spin channel: ',(nel(i),i=1,nspin)
    do ispin=1,nspin
      nvirK(:,ispin) = norbK(:) - noccK(:,ispin)
    enddo
    maxocc = maxval(noccK(1:nksym,:))
    maxvir = maxval(nvirK(1:nksym,:))
    nia_max = maxocc*maxvir

    ! limiting to insulators for now, need changes in case of metals
    do ispin=1,nspin
      do ik=2,nksym
        if(noccK(ik,ispin) .ne. noccK(1,ispin)) &
          call errore('rpa','Error: Only insulators for now!!!',1)
      enddo
    enddo

    ! reading all cholesky vectors for now, until you have a routine that reads
    ! a section
    allocate( MOs(maxnorb,maxnorb,2) )
    allocate( Lia(nia_max,nP), Ria(nia_max,nQ), Qpq(nchol_max,nchol_max,2) )
    allocate( T1(maxocc*maxnorb*max(nP,nQ)) )
    if( get_sosex ) then
      allocate( Ypq(nchol_max,nchol_max,2) ) 
    endif

    ! 2. Generate frequency grid and weights
    ! for simplicity now, use transform Gauss-Legendre grid

    ! 3. Loop over Q and frequency 
    edrpa = (0.d0,0.d0)
    if(get_sosex) esosex=(0.d0,0.d0)
    Qpq(:,:,:) = (0.d0,0.d0)
    iwQ=0
    do iw = 1, nfreq
      do iq=1,nQuniq

        Q = xQ(iq)
        Qpq(:,:,1) = (0.d0,0.d0)
        allocate( Chol(ncholQ(Q),maxnorb*maxnorb) )

        do ik=1,nkloc

          ki = k_beg+ik-1
          ka = QKtoK2(Q,ki)  ! ka = ki-Q+G

          !   3.a Read Cholesky matrix 
          call esh5_posthf_read_cholesky(h5id_hamil%id,Q-1,ki-1,ncholQ(Q),&
                maxnorb*maxnorb,Chol(1,1),error)
          if(error .ne. 0 ) &
            call errore('rpa','error reading Chol',1)
          ! should I call blas?

          allocate( CholT(norbK(ki)*norbK(ka),max(nP,nQ)) )

          do ispin=1,nspin

            !  3.b read MO matrix
            call esh5_posthf_read_dm(h5id_hamil%id,"MO_Mat",6,ki-1,ispin-1,MOs(1,1,1),error)
            if(error .ne. 0 ) &
              call errore('rpa','error reading MOs',1)
            call esh5_posthf_read_dm(h5id_hamil%id,"MO_Mat",6,ka-1,ispin-1,MOs(1,1,2),error)
            if(error .ne. 0 ) &
              call errore('rpa','error reading MOs',1)

            !  3.c rotate to spin-orbital basis 
            CholT(1:norbK(ki)*norbK(ka),1:nP) = TRANSPOSE(Chol(P_beg:P_end,1:norbK(ki)*norbK(ka)))
        
            ! Lia,P = sum_uv conj(M(u,i)) CholT(u,v,P) M(v,a)  
            Lia(:,:) = (0.d0,0.d0)
            ! T1(i,v,P) = sum_u conj(M(u,i)) CholT(u,v,P)
            call zgemm('C','N',noccK(ki,ispin),norbK(ka)*nP,norbK(ki), &
                       (1.d0,0.d0),MOs(1,1,1),maxnorb,CholT(1,1),norbK(ki), &
                       (0.d0,0.d0),T1(1),noccK(ki,ispin)) 
            do n=1,nP
              ! Lia(i,a,P) = sum_v T1(i,v,P) * M(v,a)
              call zgemm('N','N',noccK(ki,ispin),nvirK(ka,ispin),norbK(ka), &
                       (1.d0,0.d0),T1((n-1)*noccK(ki,ispin)*norbK(ka)+1),noccK(ki,ispin), &
                       MOs(1,noccK(ka,ispin)+1,2),norbK(ka), &
                       (0.d0,0.d0),Lia(1,n),noccK(ki,ispin)) 
            enddo

This doesn't make sense, is it (ia|ia) or (ia|ai)???
If it is (ia|ia), then Ria is the rotated Cholesky from -Q...
            ! Ria,P = sum_uv M(u,a) conjg(CholT(u,v,P) M(v,i)) 
            Ria(:,:) = (0.d0,0.d0)
            ! T1(a,v,Q) = sum_u conj(M(u,a)) CholT(u,v,Q)
            call zgemm('C','N',nvirK(ka,ispin),norbK(ki)*nQ,norbK(ka), &
                       (1.d0,0.d0),MOs(1,noccK(ka,ispin)+1,2),maxnorb,CholT(1,1),norbK(ki), &
                       (0.d0,0.d0),T1(1),noccK(ki,ispin))
            do n=1,nP
              ! Lia(i,a,P) = sum_v T1(i,v,P) * M(v,a)
              call zgemm('N','N',noccK(ki,ispin),nvirK(ka,ispin),norbK(ka), &
                       (1.d0,0.d0),T1((n-1)*noccK(ki,ispin)*norbK(ka)+1),noccK(ki,ispin), &
                       MOs(1,noccK(ka,ispin)+1,2),norbK(ka), &
                       (0.d0,0.d0),Lia(1,n),noccK(ki,ispin))            
            enddo
            
            !  3.d scale by eigenvalue factor  

            !  3.e accumulate Qpq 

            ! if get_sosex, accumulate Ypq

          enddo ! ispin
          !
          deallocate( CholT )
          !
        enddo

        deallocate( Chol )

        call mp_sum(Qpq,intra_image_comm)
        ! 4. keep copy if it is your turn, write Qpq to file if necessary 
        ! if(ionode) call esh5_posthf_write_Qmatrix(h5id_hamil%id,"Qpq",3,Q-1,iw-1,Qpq,error)
        if( MOD(iwQ,nproc_image) == me_image ) Qpq(:,:,2) = Qpq(:,:,1)

        iwQ = iwQ + 1
        if( iwQ == nproc_image ) then
          ! 5. calculate trace of Qpq and determinant of (1-Qpq) with LU factorization
          ! accumulate drpa energy
        endif    
        !
      enddo ! iq
    enddo ! iw
    !
    ! deallocate
    
    !
  END SUBROUTINE rpa
  !  
END MODULE rpa_module
