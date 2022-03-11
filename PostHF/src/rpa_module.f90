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
  USE mp,           ONLY: mp_root_sum, mp_sum, mp_max, mp_bcast, mp_barrier
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
  ! RPA calculation using cholesky decomposition of ERI tensor given in
  ! hamil_esh5 and rotation matrix to canonical orbitals/dft orbitals in
  ! orbs_esh5. Will read CanOrbMat from orbs_esh5
  !
  SUBROUTINE rpa_cholesky(edrpa,dfft,chol_type,hamil_esh5,orbmat_esh5,esosex) 
    !
    USE control_flags,        ONLY : use_gpu
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_esh5, orbmat_esh5, chol_type
    COMPLEX(DP), INTENT(OUT) :: edrpa
    COMPLEX(DP), INTENT(IN), OPTIONAL :: esosex
    !
#if defined(__CUDA)
!      if(use_gpu) then
!        call rpa_cholesky_gpu(edrpa,dfft,chol_type,hamil_esh5,orbmat_esh5,esosex)
!      else
        call rpa_cholesky_cpu(edrpa,dfft,chol_type,hamil_esh5,orbmat_esh5,esosex)
!      endif
#else
      call rpa_cholesky_cpu(edrpa,dfft,chol_type,hamil_esh5,orbmat_esh5,esosex)
#endif
  END SUBROUTINE rpa_cholesky
  
  SUBROUTINE rpa_cholesky_cpu(edrpa,dfft,chol_type,hamil_esh5,orbmat_esh5,esosex) 
    !
    USE qeh5_base_module
    USE constants, ONLY: tpi
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE gvect, ONLY : ecutrho
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_esh5, orbmat_esh5, chol_type
    COMPLEX(DP), INTENT(OUT) :: edrpa
    COMPLEX(DP), INTENT(IN), OPTIONAL :: esosex
    !
    TYPE(h5file_type) :: h5id_hamil, h5id_orbmat
    !
    LOGICAL :: get_sosex
    INTEGER :: i, j, n, iq, ip, Q, ka, ki, ia, ib, ii, jj
    INTEGER :: ispin, ik, ikk, iw, iwn, iwQ 
    INTEGER :: nP, maxnP, nPi, kk0
    INTEGER :: k_beg, k_end, nkloc, nia_max, ierror
    INTEGER :: nfreq, maxnorb, nchol_max, nabtot, noccmax
    REAL(DP) :: sx   ! spin degeneracy factor
    Complex(DP) :: ione
    Complex(DP) :: edmp2
    !
    INTEGER :: nel(2), maxocc, maxvir
    INTEGER, ALLOCATABLE :: noccK(:,:), nvirK(:,:)
    INTEGER, ALLOCATABLE :: norbK(:), ncholQ(:), Pbounds(:)
    REAL(DP), ALLOCATABLE :: weight(:,:),eigval(:,:)
    COMPLEX(DP), ALLOCATABLE :: eQ(:,:)    ! energy decomposition 
    COMPLEX(DP), ALLOCATABLE :: Tpq(:,:,:) ! local sector of Q matrix 
    COMPLEX(DP), ALLOCATABLE :: Qpq(:,:)   ! Q matrix 
    COMPLEX(DP), ALLOCATABLE :: Ypq(:,:,:) ! Y matrix 
    COMPLEX(DP), ALLOCATABLE :: MOs(:,:,:) ! Matrix of MO coefficients  
    COMPLEX(DP), ALLOCATABLE :: Chol(:,:) ! Cholesky matrix in SPO basis 
    COMPLEX(DP), ALLOCATABLE :: CholT(:,:) ! Cholesky matrix in SPO basis 
    COMPLEX(DP), ALLOCATABLE :: T1(:) ! Work space
    COMPLEX(DP), ALLOCATABLE :: T2(:,:) ! Work space 
    COMPLEX(DP), ALLOCATABLE :: Lia(:,:) ! left hand side of MO Cholesky matrix 
    COMPLEX(DP), ALLOCATABLE :: Ria(:,:) ! right hand side of MO Cholesky matrix 
    COMPLEX(DP), ALLOCATABLE :: Bia(:)   ! eigenvalue factors 
    REAL(DP), ALLOCATABLE :: xfreq(:), wfreq(:)
    TYPE(qeh5_file) :: qeh5_hamil 

    write(*,*)
    write(*,*) 'Starting RPA calculation. '
    CALL start_clock ( 'rpa' )

    if(chol_type.ne.'full' .and. chol_type.ne.'mo') &
      call errore('rpa_cholesky_cpu','Unknown chol_type',1)

    ione = (0.d0,1.d0)
    sx = 1.d0
    if(nspin == 1) sx = 2.d0

    get_sosex = present(esosex)
    get_sosex = .false.  ! turn off for now
    allocate( ncholQ(nksym), norbK(nksym) )
    if(chol_type == 'full') then
      CALL esh5_posthf_open_file_read(h5id_hamil%id,TRIM(hamil_esh5),  &
                                                  LEN_TRIM(hamil_esh5),ierror)
      if(ierror .ne. 0 ) &
        call errore('rpa','error opening hamil file',1)
      call esh5_posthf_read_nchol(h5id_hamil%id,ncholQ, ierror)
      if(ierror .ne. 0) call errore('rpa','Error: esh5_posthf_read_nchol',ierror)
      call esh5_posthf_read_norbk(h5id_hamil%id,norbK, ierror)
      if(ierror .ne. 0) call errore('rpa','Error: esh5_posthf_read_norbk',ierror)
    elseif(chol_type == 'mo') then
      call qeh5_openfile(qeh5_hamil,hamil_esh5,'read')
      call read_meta()
! read ncholQ, norbK, noccmax, nabtot        
    endif
    nchol_max = maxval(ncholQ(:))
    maxnorb = maxval(norbK(:))

    CALL esh5_posthf_open_file_read(h5id_orbmat%id,TRIM(orbmat_esh5),  &
                                                  LEN_TRIM(orbmat_esh5),ierror)
    if(ierror .ne. 0 ) &
      call errore('rpa','error opening orb mat file',1)
    call esh5_posthf_read_orbmat_info(h5id_orbmat%id,'CanOrbMat',9,i,n,ik,ispin,ierror)
    if(ierror .ne. 0 ) &
      call errore('rpa','error reading orbmat info',1)
    if( i .ne. maxnorb .or. n .ne. maxnorb .or. ik .ne. nksym .or. ispin .ne. nspin ) &
      call errore('rpa','inconsistent parameters in orbmat file',1)

    !  Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    ! Partition nchol_max 
    allocate(Pbounds(nproc_pool+1))
    Pbounds(1) = 1 
    maxnP=0
    do n=2,nproc_pool
      call fair_divide(Pbounds(n),i,n,nproc_pool,nchol_max)
    enddo
    Pbounds(nproc_pool+1)=nchol_max+1
    nP = Pbounds(me_pool+2) - Pbounds(me_pool+1) 
    maxnP = Pbounds(2) - Pbounds(1) 

    ! 1. Read eigenvalues, fermi weights and determinte occupied/virtual partitioning
    allocate( noccK(nksym,numspin), nvirK(nksym,numspin) )
    allocate(eigval(maxnorb,nkstot),weight(maxnorb,nkstot))
    weight(:,:) = 0.d0
    eigval(:,:) = 0.d0
    noccK(:,:)=0
    nvirK(:,:)=0
    do ikk=1,nkstot
      call esh5_posthf_read_et(h5id_orbmat%id,ikk-1,eigval(1,ikk), &
                             weight(1,ikk),ierror)
      if(ierror .ne. 0 ) &
        call errore('rpa','error reading weights',1)
    enddo
    ! find number of electrons
    call get_noccK(noccK,nel,maxnorb,nksym,numspin,weight,maxnorb)
    write(*,*) ' Number of electrons per spin channel: ',(nel(i),i=1,numspin)
    do ispin=1,numspin
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

    ! Generate frequency grid and weights
    ! for simplicity now, use transform Gauss-Legendre grid
    nfreq = 40 
    allocate( wfreq(nfreq), xfreq(nfreq) )
    call gaussleg_quad(nfreq,xfreq,wfreq)
    call transform_gaussleg(nfreq,0.5d0,xfreq,wfreq)
    call test_gaussleg_grid(nfreq,10.0d0,xfreq,wfreq)

    ! reading all cholesky vectors for now, until you have a routine that reads
    ! a section
    allocate( Qpq(nchol_max,nchol_max), Tpq(nP,nchol_max,nfreq) )
    allocate( T2(nchol_max,nchol_max) )
    allocate( Bia(nia_max), eQ(nksym,10) )
    eQ(:,:) = (0.d0,0.d0)
    if(chol_type == 'full') then
      allocate( MOs(maxnorb,maxnorb,2), T1(maxocc*maxnorb*nP))
      allocate( Lia(nia_max,maxnP), Ria(nia_max,maxnP) ) 
    elseif(chol_type == 'mo') then
      allocate( Lia(maxnP,nia_max), Ria(maxnP,nia_max) ) 
    endif

! MAM: put timers everywhere and improve performance where possible

    ! 3. Loop over Q and frequency 
    do iq=1,nQuniq

      Q = xQ(iq)
      if(chol_type == 'full') & 
        allocate( Chol(ncholQ(Q),maxnorb*maxnorb) )

      Tpq(:,:,:) = (0.d0,0.d0)
      do ik=1,nkloc

        ki = k_beg+ik-1
        ka = QKtoK2(Q,ki)  ! ka = ki-Q+G

        CALL start_clock ( 'rpa_io' )
        if(chol_type == 'full') then
          !   3.a Read Cholesky matrix 
          call esh5_posthf_read_cholesky(h5id_hamil%id,Q-1,ki-1,ncholQ(Q),&
                maxnorb*maxnorb,Chol(1,1),ierror)
          if(ierror .ne. 0 ) &
            call errore('rpa','error reading Chol',1)
          ! should I call blas?

          allocate( CholT(norbK(ki)*norbK(ka),nP) )
          ! improve this
          do n=1,nP
            do ii=1,norbK(ki)
            do jj=1,norbK(ka)
              CholT((jj-1)*norbK(ki)+ii,n) = &
                  Chol(Pbounds(me_pool+1)+n-1,(ii-1)*norbK(ka)+jj)
            enddo
            enddo
          enddo
        endif
        CALL stop_clock ( 'rpa_io' )

        do ispin=1,nspin

          kk0 = (ispin-1)*nksym

          if(chol_type == 'full') then
            !  3.b read MO matrix
            CALL start_clock ( 'rpa_io' )
            call esh5_posthf_read_orbmat(h5id_orbmat%id,"CanOrbMat",9,  &
                        ki-1,ispin-1,MOs(1,1,1),ierror)
            if(ierror .ne. 0 ) &
              call errore('rpa','error reading MOs',1)
            call esh5_posthf_read_orbmat(h5id_orbmat%id,"CanOrbMat",9,  &
                        ka-1,ispin-1,MOs(1,1,2),ierror)
            if(ierror .ne. 0 ) &
              call errore('rpa','error reading MOs',1)
            CALL stop_clock ( 'rpa_io' )

            CALL start_clock ( 'rpa_Lia' )
            ! process the occ-vir sector!
            ! Lia,P = sum_uv conj(M(u,i)) CholT(u,v,P) M(v,a)  
            Lia(:,:) = (0.d0,0.d0)
            ! T1(i,v,P) = sum_u conj(M(u,i)) CholT(u,v,P)
            call zgemm('C','N',noccK(ki,ispin),norbK(ka)*nP,norbK(ki), &
                       (1.d0,0.d0),MOs(1,1,1),size(MOs,1),  &
                                 CholT(1,1),norbK(ki), &
                       (0.d0,0.d0),T1(1),noccK(ki,ispin)) 
            do n=1,nP
              ! Lia(i,a,P) = sum_v T1(i,v,P) * M(v,a)
              call zgemm('N','N',noccK(ki,ispin),nvirK(ka,ispin),norbK(ka), &
                       (1.d0,0.d0),T1((n-1)*noccK(ki,ispin)*norbK(ka)+1),noccK(ki,ispin), &
                       MOs(1,noccK(ka,ispin)+1,2),size(MOs,1), &
                       (0.d0,0.d0),Lia(1,n),noccK(ki,ispin)) 
            enddo
            CALL stop_clock ( 'rpa_Lia' )
          elseif(chol_type == 'mo') then
            ! read Lp,ia and conjugate
            CALL start_clock ( 'rpa_io' )
            call read_Lia_conjugated(Q,ki,ka,ispin,Pbounds(me_pool+1),nP)
            CALL stop_clock ( 'rpa_io' )
          endif

          CALL start_clock ( 'rpa_Qpq' )
          do ip=1,nproc_pool

            if( Pbounds(ip) > ncholQ(Q) ) exit
            nPi = min( Pbounds(ip+1)-Pbounds(ip), ncholQ(Q)-Pbounds(ip)+1) 
            if(nPi == 0) exit ! nothing to do, finish loop
            ! broadcast Lia 
            if( me_pool+1 == ip ) Ria(:,:) = Lia(:,:) 
            call mp_bcast( Ria, ip-1, intra_pool_comm )

            do iw = 1, nfreq

              ! setup Bia
              ! MAM: This is wrong with partial occupations! Need to keep
              ! partially occupied states in both occ and virtual sets!
              n = 0
              do ia = noccK(ka,ispin)+1,noccK(ka,ispin)+nvirK(ka,ispin) 
                do ii = 1, noccK(ki,ispin)
                  ! fi / [ (e(a) - e(i)) - iw ]
                  n = n + 1
                  Bia(n) = weight(ii,ki+kk0) / ( (eigval(ia,ka+kk0) - eigval(ii,ki+kk0)) &
                                          - ione*xfreq(iw) )
                enddo
              enddo

              if(chol_type=='full') then
                !  scale by eigenvalue factor 
                do i=1,nPi
                  do ia=1,n
                    Ria(ia,i) = Ria(ia,i) * Bia(ia)
                  enddo
                enddo
                ! accumulate Tpq = sum conj(Liap) * Riap
                call zgemm('C','N',nP,nPi,n, &
                         (1.d0,0.d0)*sx,Lia,size(Lia,1),Ria,size(Ria,1), &
                         (1.d0,0.d0),Tpq(1,Pbounds(ip),iw),size(Tpq,1))
                !  remove eigenvalue factor 
                do i=1,nPi
                  do ia=1,n
                    Ria(ia,i) = Ria(ia,i) / Bia(ia)
                  enddo
                enddo
              elseif(chol_type == 'mo') then
                !  scale by eigenvalue factor 
                do ia=1,n
                  do i=1,nPi
                    Ria(i,ia) = Ria(i,ia) * conjg(Bia(ia))
                  enddo
                enddo
                ! accumulate Tpq = sum conj(Liap) * Riap
                call zgemm('N','C',nP,nPi,n, &
                         (1.d0,0.d0)*sx,Lia,size(Lia,1),Ria,size(Ria,1), &
                         (1.d0,0.d0),Tpq(1,Pbounds(ip),iw),size(Tpq,1))
                !  remove eigenvalue factor 
                do ia=1,n
                  do i=1,nPi
                    Ria(i,ia) = Ria(i,ia) / conjg(Bia(ia))
                  enddo
                enddo
              endif

            enddo ! iw
            !
          enddo ! ip - loop over procs in npool
          CALL stop_clock ( 'rpa_Qpq' )

          ! now repeat the process for the vir-occ sector!
          if(chol_type=='full') then
            CALL start_clock ( 'rpa_Lia' )
            ! Lia,P = sum_uv conj(M(u,a)) CholT(u,v,P) M(v,i)  
            Lia(:,:) = (0.d0,0.d0)
            ! T1(u,i,P) = sum_v CholT(u,v,P) M(v,i) 
            do n=1,nP  
              call zgemm('N','N',norbK(ki),noccK(ka,ispin),norbK(ka), &
                     (1.d0,0.d0),CholT(1,n),norbK(ki), &
                     MOs(1,1,2),size(MOs,1), &
                     (0.d0,0.d0),T1((n-1)*norbK(ki)*noccK(ka,ispin)+1),norbK(ki))
            enddo
            ! Lia(a,i,P) = sum_u conj(M(u,a)) T1(u,i,P)
            do n=1,nP  
              call zgemm('C','N',nvirK(ki,ispin),noccK(ka,ispin),norbK(ki), &
                     (1.d0,0.d0),MOs(1,noccK(ki,ispin)+1,1),size(MOs,1), &
                     T1((n-1)*norbK(ki)*noccK(ka,ispin)+1),norbK(ki), &
                     (0.d0,0.d0),Lia(1,n),nvirK(ki,ispin)) 
            enddo
            CALL stop_clock ( 'rpa_Lia' )
          elseif(chol_type == 'mo') then
            ! read Lp,ai and conjugate
            CALL start_clock ( 'rpa_io' )
            call read_Lai_conjugated(Q,ki,ka,ispin,Pbounds(me_pool+1),nP)
            CALL stop_clock ( 'rpa_io' )
          endif

          do ip=1,nproc_pool

            if( Pbounds(ip) > ncholQ(Q) ) exit
            nPi = min( Pbounds(ip+1)-Pbounds(ip), ncholQ(Q)-Pbounds(ip)+1) 
            if(nPi == 0) exit ! nothing to do, finish loop
            ! broadcast Lia 
            if( me_pool+1 == ip ) Ria(:,:) = Lia(:,:) 
            call mp_bcast( Ria, ip-1, intra_pool_comm )

            do iw = 1, nfreq

              ! setup Bia
              ! MAM: This is wrong with partial occupations! Need to keep
              ! partially occupied states in both occ and virtual sets!
              n = 0
              do ii = 1, noccK(ka,ispin)
                do ia = noccK(ki,ispin)+1,noccK(ki,ispin)+nvirK(ki,ispin) 
                  ! fi / [ (e(a) - e(i)) - iw ]
                  n = n + 1
                  Bia(n) = weight(ii,ka+kk0) / ( (eigval(ia,ki+kk0) - eigval(ii,ka+kk0)) &
                                          + ione*xfreq(iw) )
                enddo
              enddo

              if(chol_type=='full') then
                !  scale by eigenvalue factor 
                do i=1,nPi
                  do ia=1,n
                    Ria(ia,i) = Ria(ia,i) * Bia(ia)
                  enddo
                enddo
                ! accumulate Tpq = sum conj(Liap) * Riap
                call zgemm('C','N',nP,nPi,n, &
                         (1.d0,0.d0)*sx,Lia,size(Lia,1),Ria,size(Ria,1), &
                         (1.d0,0.d0),Tpq(1,Pbounds(ip),iw),size(Tpq,1))
                !  remove eigenvalue factor 
                do i=1,nPi
                  do ia=1,n
                    Ria(ia,i) = Ria(ia,i) / Bia(ia)
                  enddo
                enddo
              elseif(chol_type == 'mo') then
                !  scale by eigenvalue factor 
                do ia=1,n
                  do i=1,nPi
                    Ria(i,ia) = Ria(i,ia) * conjg(Bia(ia))
                  enddo
                enddo
                ! accumulate Tpq = sum Lpia * Rqia.T.conj()
                call zgemm('N','C',nP,nPi,n, &
                         (1.d0,0.d0)*sx,Lia,size(Lia,1),Ria,size(Ria,1), &
                         (1.d0,0.d0),Tpq(1,Pbounds(ip),iw),size(Tpq,1))
                !  remove eigenvalue factor 
                do ia=1,n
                  do i=1,nPi
                    Ria(i,ia) = Ria(i,ia) / conjg(Bia(ia))
                  enddo
                enddo
              endif

            enddo ! iw
            !
          enddo ! ip - loop over procs in npool

        enddo ! ispin
        !
        if(allocated(CholT)) deallocate( CholT )
        !
      enddo ! ik

      if(allocated(Chol)) deallocate( Chol )

      CALL start_clock ( 'rpa_process_Q' )
      ! accumulate Qpq and calculate energy contributions
      ! round-robin distribution over cores
      do iw = 1, nfreq
        !
        T2(:,:)=(0.d0,0.d0)
        T2(Pbounds(me_pool+1):Pbounds(me_pool+1)+nP-1,:) = Tpq(1:nP,:,iw)
        call mp_root_sum(T2,Qpq,mod(iw-1,nproc_image),intra_image_comm)
        ! record frequency 'assigned' to me_image 
        if(me_image == mod(iw-1,nproc_image)) iwn = iw 
        ! calculate contributions 
        if( mod(iw,nproc_image) == 0 .or. iw == nfreq) then
          ! nothing to do on last iteration for this rank   
          if(me_image > mod(iw-1,nproc_image)) exit  
          call rpa_process_Q(ncholQ(Q),wfreq(iwn),Qpq,T2,eQ(iq,1),eQ(iq,2))  
          !
        endif
        !
      enddo  
      call mp_barrier( intra_image_comm )
      CALL stop_clock ( 'rpa_process_Q' )
      !
    enddo ! iq
    ! done
    if(nproc_image > 1) CALL mp_sum ( eQ, intra_image_comm )
    edrpa = (0.d0,0.d0)
    edmp2 = (0.d0,0.d0)
    do iq=1,nQuniq
      edrpa = edrpa + wQ(iq)*eQ(iq,1)
      edmp2 = edmp2 + wQ(iq)*eQ(iq,2)
    enddo

    write(*,*)
    write(*,*) ' ************************************** '
    write(*,*) '  EdMP2 (Ha): ',edmp2
    write(*,*) '  EdRPA (Ha): ',edrpa

    IF ( ionode .and. verbose ) THEN
      !
      write(*,*)
      write(*,*) '  Q  EdRPA(Q) EdMP2(Q) '
      do iq=1,nQuniq
        write(*,'(i5,g14.4,"("g14.6,g14.6") "," ("g14.6,g14.6")")') iq,wQ(iq),&
            eQ(iq,1),eQ(iq,2)
      enddo
      !
    ENDIF
    write(*,*) ' ************************************** '

    if(chol_type=='full') then
      call esh5_posthf_close_file(h5id_hamil%id) 
    elseif(chol_type=='mo') then
      call qeh5_close(qeh5_hamil)
    endif

    CALL stop_clock ( 'rpa' )
    if(ionode) then
      write(*,*) 'Timers: '
      CALL print_clock ( 'rpa' )
      CALL print_clock ( 'rpa_io' )
      CALL print_clock ( 'rpa_Qpq' )
      if(chol_type == 'full' )CALL print_clock ( 'rpa_Lia' )
      CALL print_clock ( 'rpa_process_Q' )
    endif
    !
    ! deallocate
    if(allocated(noccK)) deallocate(noccK)
    if(allocated(nvirK)) deallocate(nvirK)
    if(allocated(norbK)) deallocate(norbK)
    if(allocated(ncholQ)) deallocate(ncholQ)
    if(allocated(weight)) deallocate(weight)
    if(allocated(eigval)) deallocate(eigval)
    if(allocated(eQ)) deallocate(eQ)
    if(allocated(Qpq)) deallocate(Qpq)
    if(allocated(Tpq)) deallocate(Tpq)
    if(allocated(Ypq)) deallocate(Ypq)
    if(allocated(MOs)) deallocate(MOs)
    if(allocated(Chol)) deallocate(Chol)
    if(allocated(CholT)) deallocate(CholT)
    if(allocated(T1)) deallocate(T1)
    if(allocated(T2)) deallocate(T2)
    if(allocated(Lia)) deallocate(Lia)
    if(allocated(Ria)) deallocate(Ria)
    if(allocated(Bia)) deallocate(Bia)
    if(allocated(xfreq)) deallocate(xfreq)
    if(allocated(wfreq)) deallocate(wfreq)
    if(allocated(Pbounds)) deallocate(Pbounds)
    !if(allocated()) deallocate()
 
    CONTAINS

    SUBROUTINE read_Lia_conjugated(Q,ki,ka,ispin,P_beg,nP)
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Q,ki,ka,ispin,P_beg,nP  
      INTEGER :: iaC,ia,a,i,n, nPi
      CHARACTER(LEN=4) :: Qstr,kstr,sstr
      TYPE(qeh5_file) :: ham,moL
      TYPE(qeh5_dataset) :: dset
      COMPLEX(DP), ALLOCATABLE :: buff(:,:)  
      COMPLEX(DP) :: eJ,eX

      nPi = min( ncholQ(Q)-P_beg+1, nP )  
      if(noccK(ki,ispin) > noccmax) &
        call errore('read_Lia','noccK(ki,ispin) > noccmax',1)  
      allocate(buff(nPi,noccK(ki,ispin)*norbK(ka)))  
      call qeh5_open_group(qeh5_hamil, "Hamiltonian", ham)
      call qeh5_open_group(ham, "MOCholesky", moL)
      write ( Qstr, '(I4)') Q-1
      write ( kstr, '(I4)') ki-1
      write ( sstr, '(I4)') ispin-1
      dset%name = "L"//trim(adjustl(Qstr))//"_k"//  &
                          trim(adjustl(kstr))//"_s"//trim(adjustl(sstr))
      call qeh5_set_space(dset, buff(1,1), 2, [nPi,size(buff,2)], 'm')
      call qeh5_open_dataset(moL, dset, ACTION = 'read')
      call qeh5_set_file_hyperslab(dset, OFFSET = [2*P_beg-2,0], &
                         COUNT = [2*nPi,size(buff,2)])
      call qeh5_read_dataset(buff, dset)

      ia=0
      ! careful here, ia is c-ordered in buff and fortran ordered in Lia  
      Lia(:,:)=(0.d0,0.d0)  
      do a=noccK(ka,ispin)+1,noccK(ka,ispin)+nvirK(ka,ispin)
      do i=1,noccK(ki,ispin)
        ia=ia+1
        iaC = (i-1)*norbK(ka) + a
        do n=1,nPi
          Lia(n,ia) = conjg(buff(n,iaC))  
        enddo   
      enddo   
      enddo   
      call qeh5_close(dset)
      call qeh5_close(moL)
      call qeh5_close(ham)
      deallocate(buff)
      !
    END SUBROUTINE read_Lia_conjugated

    SUBROUTINE read_Lai_conjugated(Q,ki,ka,ispin,P_beg,nP)
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Q,ki,ka,ispin,P_beg,nP
      INTEGER :: aiC,ai,a,i,n,nPi
      CHARACTER(LEN=4) :: Qstr,kstr,sstr
      TYPE(qeh5_file) :: ham,moL
      TYPE(qeh5_dataset) :: dset
      COMPLEX(DP), ALLOCATABLE :: buff(:,:)

      nPi = min( ncholQ(Q)-P_beg+1, nP )  
      if(noccK(ka,ispin) > noccmax) &
        call errore('read_Lai','noccK(ka,ispin) > noccmax',1)
      allocate(buff(nPi,nvirK(ki,ispin)*noccmax))
      call qeh5_open_group(qeh5_hamil, "Hamiltonian", ham)
      call qeh5_open_group(ham, "MOCholesky", moL)
      write ( Qstr, '(I4)') Q-1
      write ( kstr, '(I4)') ki-1
      write ( sstr, '(I4)') ispin-1
      dset%name = "L"//trim(adjustl(Qstr))//"_k"//  &
                          trim(adjustl(kstr))//"_s"//trim(adjustl(sstr))
      call qeh5_set_space(dset, buff(1,1), 2, [nPi,size(buff,2)], 'm')
      call qeh5_open_dataset(moL, dset, ACTION = 'read')
      call qeh5_set_file_hyperslab(dset, &
                         OFFSET = [2*P_beg-2,norbK(ka)*noccmax], &
                         COUNT = [2*nPi,size(buff,2)])
      call qeh5_read_dataset(buff, dset)
      Lia(:,:)=(0.d0,0.d0)  
      ai=0
      do i=1,noccK(ka,ispin)
      do a=1,nvirK(ki,ispin)
        ai=ai+1
        aiC = (a-1)*noccmax + i
        do n=1,nPi
          Lia(n,ai) = conjg(buff(n,aiC))
        enddo
      enddo
      enddo
      call qeh5_close(dset)
      call qeh5_close(moL)
      call qeh5_close(ham)
      deallocate(buff)
      !  
    END SUBROUTINE read_Lai_conjugated

    SUBROUTINE read_meta()
      USE h5lt
      !
      IMPLICIT NONE
      integer(HSIZE_T) :: dims(1)
      integer(HID_T) :: int_type_id, dp_type_id
      TYPE(qeh5_file) :: ham,moL
      INTEGER :: err, nmotot, d(1)

      call qeh5_open_group(qeh5_hamil, "Hamiltonian", ham)
      call qeh5_open_group(ham, "MOCholesky", moL)

      dims = [1]
      call h5ltread_dataset_f(moL%id, "noccmax", H5T_NATIVE_INTEGER, d, dims, err)
      call errore('read_meta','H5LTmake_dataset_f',abs(err))
      noccmax = d(1) 

      dims = [1]
      call h5ltread_dataset_f(moL%id, "nabtot", H5T_NATIVE_INTEGER, d, dims, err)
      call errore('read_meta','H5LTmake_dataset_f',abs(err))
      nabtot = d(1) 

      dims = [nksym]
      call h5ltread_dataset_f(moL%id, "NMOPerKP", H5T_NATIVE_INTEGER, norbK, & 
                        dims, err)
      call errore('read_meta','H5LTmake_dataset_f',abs(err))

      dims = [nksym]
      call h5ltread_dataset_f(moL%id, "NCholPerKP", H5T_NATIVE_INTEGER, ncholQ, &
                        dims, err)
      call errore('read_meta','H5LTmake_dataset_f',abs(err))

      call qeh5_close(moL)
      call qeh5_close(ham)

! read ncholQ, norbK, noccmax, nabtot 
    END SUBROUTINE read_meta

    !
  END SUBROUTINE rpa_cholesky_cpu
  !  
  ! Qpq = -<P| v^{1/2}(r) x(r,r') v^{1/2}(r') |Q>
  ! Edmp2 = -Tr[Qpq^2]
  ! Edrpa = -Tr[ ln(1+Qpq) - Qpq ] = - sum_n [ ln(en) + 1 - en ], 
  !        where en are the eigenvalues of 1+Qpq. We actually use an LU
  !        factorization of 1+Qpq to avoid the more expensive calculation of
  !        eigenvalues.
  SUBROUTINE rpa_process_Q(naux,iw,Qpq,Tpq,edrpa,edmp2) 
    !
    USE kinds, ONLY : DP
    USE constants, ONLY: tpi
    USE lsda_mod, ONLY: nspin
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: naux
    REAL(DP), INTENT(IN) :: iw
    COMPLEX(DP), INTENT(INOUT) :: edrpa, edmp2
    COMPLEX(DP), INTENT(INOUT) :: Qpq(:,:)
    COMPLEX(DP), INTENT(IN) :: Tpq(:,:)  
    INTEGER, ALLOCATABLE :: ipiv (:)
    !
    INTEGER :: i,j,k,n,info
    REAL(DP) :: dx, sg
    !
    ! weight of frequency and factor of 2*pi
    dx = iw / tpi
    ! E = -1/(2*pi) int_w Tr[ ln(1-Q) + Q ]
    ! evaluating emp2 energy for now to debug 
    do i=1,naux
      edrpa = edrpa - dx * Qpq(i,i) ! remember that Q has factor of -1.0 
      do j=1,naux
        edmp2 = edmp2 - 0.5d0 * dx * Qpq(i,j) * Qpq(j,i)  
      enddo 
    enddo 
    !
    allocate(ipiv(naux))
    !
    do i=1,naux
      Qpq(i,i) = Qpq(i,i)+1.d0
    enddo
    CALL zgetrf (naux, naux, Qpq, size(Qpq,1), ipiv, info)
    CALL errore ('rpa_process_Q', 'error in ZGETRF', abs (info) )
    sg = 1.d0
    do i=1,naux
      if(real(Qpq(i,i)) < 0.d0) then
        edrpa = edrpa + dx * log(-1.d0*Qpq(i,i)) 
        sg = -sg 
      else
        edrpa = edrpa + dx * log(Qpq(i,i)) 
      endif
      if(ipiv(i) .ne. i) sg = -sg 
    enddo
    if(sg < 0.d0) &
      write(*,*) ' *** WARNING: Negative determinant found in RPA Q matrix. *** '
    !
    deallocate(ipiv) 
    !
  END SUBROUTINE rpa_process_Q
  !  
END MODULE rpa_module
