!--------------------------------------------------------------------
! Written by Miguel A. Morales, LLNL, 2020 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE mp2_module
  !----------------------------------------------------------------------
  ! 
  ! Implements MP2 and some other utilities, e.g. pseudo-canonicalization
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY: omega, alat, tpiba, tpiba2, at, bg
  USE posthf_mod, ONLY: nksym,norb,numspin,e0,efc,e2Ha,&
                        nelec_tot,nup,ndown,ke_factorization,&
                        xksym,igksym,ngksym,QKtoK2,kminus,Qpts,DM,nmax_DM, &
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
  USE becmod,   ONLY : bec_type, becp, allocate_bec_type, deallocate_bec_type
  USE mp,           ONLY: mp_sum, mp_max, mp_bcast, mp_barrier
  USE mp_images, ONLY: intra_image_comm, me_image, root_image, nproc_image
  USE mp_pools,     ONLY: inter_pool_comm, intra_pool_comm, npool, nproc_pool, &
                          me_pool,root_pool,my_pool_id
  USE noncollin_module,     ONLY : noncolin, npol
  USE lsda_mod, ONLY: lsda, nspin
  USE read_orbitals_from_file, ONLY: get_orbitals_set
  use fft_interfaces,       ONLY : invfft, fwfft
  USE fft_types, ONLY: fft_type_descriptor
  USE posthf_mod, ONLY : ke_factorization
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
  SUBROUTINE mp2_g(emp2,dfft,esh5_file,reg_pow,reg_expo)
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE gvect, ONLY : ecutrho
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE paw_variables,           ONLY : okpaw
    USE becmod,  ONLY : bec_type, ALLOCATE_bec_type, DEALLOCATE_bec_type
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: esh5_file
    INTEGER, INTENT(IN), OPTIONAL :: reg_pow
    REAL(DP), INTENT(IN), OPTIONAL :: reg_expo
    COMPLEX(DP), INTENT(OUT) :: emp2
    !
    INTEGER :: Q, ka, kb, iab, ia, ib, ii, ij, ki, kj, ni, nj, iuv
    INTEGER :: ik, ibnd, nxxs, i, ikk, kk, jj, ia0, iu0, ispin, is1, is2, noa, nob, n
    INTEGER :: a_beg, a_end, b_beg, b_end, ab_beg, ab_end, nabpair  
    INTEGER :: naorb, nborb, iabmax
    INTEGER :: maxv_rank, error
    INTEGER :: nel(2), maxocc, maxvir, nke_dim
    INTEGER, ALLOCATABLE :: noccK(:,:), nvirK(:,:) 
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: regkappa
    INTEGER :: regp
    REAL(DP) :: dQ(3),dQ1(3),dk(3), dG(3, 27), scl, qcut(5)
    REAL(DP) :: residual, maxv, fXX, fac
    REAL(DP) :: etemp, etemp2, q2x
    LOGICAL :: first, regularize
    ! QPOINT stuff
    INTEGER, ALLOCATABLE :: norbK(:)
    REAL(DP), ALLOCATABLE :: weight(:,:),eigval(:,:)
    COMPLEX(DP), ALLOCATABLE :: Kia(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Vr(:),Vr2(:)    
    COMPLEX(DP), ALLOCATABLE :: Psiocc(:,:,:)   
    COMPLEX(DP), ALLOCATABLE :: Psivira(:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: Psivirb(:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Rgb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vab(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vba(:,:)  ! 
    TYPE(ke_factorization), ALLOCATABLE :: ke(:) ! for paw one-center terms
    COMPLEX(DP),ALLOCATABLE :: paw_one_center(:,:,:)
    TYPE(bec_type),ALLOCATABLE :: becpsia(:)
    TYPE(bec_type),ALLOCATABLE :: becpsib(:)
    TYPE(bec_type),ALLOCATABLE :: becocc(:)
    CHARACTER(len=10) read_type
    CHARACTER(len=256) h5name
    type(h5file_type) :: h5id_input_orbs
    !
    regp=1
    regkappa=0.d0
    if(present(reg_pow)) regp = max(1,reg_pow)
    if(present(reg_expo)) regkappa = max(0.d0,reg_expo) 
    regularize = (regkappa > 0.d0) 

    nxxs = dfft%nr1x*dfft%nr2x*dfft%nr3x

    if(nspin == 1 ) then
      fac = 2.d0 
    else
      fac = 1.d0 
    endif
    if(noncolin) call errore('mp2_g','No nonconlin yet.',1)

    IF(tqr) THEN
      tabxx => tabp
    ENDIF

    allocate( noccK(nksym,numspin), nvirK(nksym,numspin), norbK(nksym) )
    if(present(esh5_file)) then
      h5name = TRIM(esh5_file)
      call open_esh5_read(h5id_input_orbs,h5name)
      norbK(:) = h5id_input_orbs%norbK(:)
      allocate(eigval(h5id_input_orbs%maxnorb,nkstot), &
                      weight(h5id_input_orbs%maxnorb,nkstot))  
      weight(:,:) = 0.d0
      eigval(:,:) = 0.d0
      do ik=1,nksym
        call esh5_posthf_read_et(h5id_input_orbs%id,ik-1,eigval(1,ik), &
                                 weight(1,ik),error)  
        if(error .ne. 0 ) &
          call errore('mp2_g','error reading weights',1)
      enddo  
      ! find number of electrons
      call get_noccK(noccK,nel,h5id_input_orbs%maxnorb,nksym,numspin, &
                     weight,h5id_input_orbs%maxnorb)
    else
      if(npool > 1) call errore('mp2_g', &
             'Error: npool > 1 not allowed in mp2_g without esh5 orbitals',1)
      allocate(eigval(nbnd,nkstot))
      eigval(:,:) = et(:,:)*e2Ha
      ! find number of electrons
      call get_noccK(noccK,nel,minval(norbK(:)),nksym,numspin,wg,nbnd,wk)
      norbK(:) = norb
    endif
    write(*,*) ' Number of electrons per spin channel: ',(nel(i),i=1,numspin)
    nvirK(:,1) = norbK(:) - noccK(:,1)
    if(numspin==2) &
      nvirK(:,2) = norbK(:) - noccK(:,2)
    maxocc = maxval(noccK(1:nksym,1))
    maxvir = maxval(nvirK(1:nksym,1))

    ! limiting to insulators for now, need changes in case of metals
    first=.true.
    do ik=1,nksym
      do ii=1,norbK(ik)
        if( (weight(ii,ik) > 0.001d0) .and. &
            (weight(ii,ik) < 0.999d0) .and. &
              first) then 
          first = .false.
          if(ionode) then
            write(*,*) ''
            write(*,*) '**********************************************************************'
            write(*,*) '*** Warning: Partial occupations found. MP2 energy is not correct. ***'
            write(*,*) '**********************************************************************'
            write(*,*) ''
          endif
        endif
      enddo
    enddo

    ! MAM: parallelization over bands only, regardless of value of npools
    call find_2d_partition(maxvir,me_image+1,nproc_image,a_beg,a_end,b_beg,b_end)
    naorb   = a_end-a_beg+1
    nborb   = b_end-b_beg+1
    if( (naorb < 1) .or. (nborb < 1) ) &
      call errore('mp2_g','Error: Too many processors.',1)

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over nxxs and gather over proc grid
! MAM: NOTE: can store 3 components separately and multiply later on, 
! reduces memory cost from nx*ny*nz to 3*max(nx,ny,nz)
! application might be slow though
    allocate( phasefac(nxxs,27) )
    CALL calculate_phase_factor(dfft, phasefac, dG)

    if(ionode) write(*,*) 'Reading orbitals from file'
    !
    emp2 = (0.d0,0.d0)
    ispin = 1  ! alpha for now only

    ! eventually distribute nxxs and reduce matrices below within groups
    allocate( Kia(ngm,naorb,nksym), Kib(ngm,nborb,nksym) )
    allocate( Psivira(nxxs,naorb,nksym), Psiocc(nxxs,maxocc,nksym) )
    allocate( Psivirb(nxxs,nborb,nksym) )
    allocate( Vab(naorb,nborb), Vba(nborb,naorb), Rgb(ngm,max(naorb,nborb)) )
    allocate( vCoulG(ngm), Vr(nxxs), Vr2(nxxs) )

    call start_clock( 'mp2_io' )
    nke_dim = 0

    if(present(esh5_file)) then
      read_type = 'esh5'
    else
      read_type = 'davcio'
    endif
    if(okvan.or.okpaw) then
      ! store bec for Psia 
      ALLOCATE(becpsia(nksym))
      ALLOCATE(becpsib(nksym))
      ALLOCATE(becocc(nksym))
      do ik=1,nksym
        CALL allocate_bec_type( nkb, nborb, becpsib(ik))
        CALL allocate_bec_type( nkb, maxocc, becocc(ik))
      enddo
      if(okpaw) then
        allocate(ke(ntyp))
        call calculate_factorized_paw_one_center(ke,nke_dim)
        !allocate( paw_one_center(nke_dim,nabpair,nksym) )
      endif
      call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          1,Psiocc,1,maxocc,1,nksym,becpsi=becocc)
      do ki=1,nksym
        call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          1,Psivira(:,:,ki:ki),noccK(ki,1)+a_beg,naorb,ki,1, &
                          becpsi=becpsia(ki:ki))
        call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft, &
                          1,Psivirb(:,:,ki:ki),noccK(ki,1)+b_beg,nborb,ki,1, &
                          becpsi=becpsib(ki:ki))
      enddo
    else
      call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          1,Psiocc,1,maxocc,1,nksym)
      do ki=1,nksym
        call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          1,Psivira(:,:,ki:ki),noccK(ki,1)+a_beg,naorb,ki,1)
        call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          1,Psivirb(:,:,ki:ki),noccK(ki,1)+b_beg,nborb,ki,1)
      enddo
    endif
    if(present(esh5_file)) then
      call close_esh5_read(h5id_input_orbs)
    endif

    call stop_clock( 'mp2_io' )

    ! need eigenvalues, right now assuming they are given by pwscf files
    if(ionode) write(*,*) 'Starting loop over bands' 

    do ki=1,nksym

      if(ionode) write(*,*) 'K:',ki
      CALL start_clock ( 'mp2' )
      do ii=1,noccK(ki,ispin)
    
        call start_clock( 'mp2_Kia' )
        ! calculate Kia/Kib
        Kia(:,:,:)=(0.d0,0.d0)
        Kib(:,:,:)=(0.d0,0.d0)
        do ka=1,nksym

          ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(ka) - k(ki)  |^{-2}
          CALL g2_convolution(dfft%ngm, g, xksym (1:3, ka), xksym (1:3, ki), vCoulG)
        
          IF ( okvan .and..not.tqr ) &
            CALL qvan_init (dfft%ngm, xksym (1:3, ki), xksym (1:3, ka))  

          do ibnd=1,naorb

            ! Orbital pairs in R
            Vr(1:dfft%nnr) = CONJG(Psiocc(1:dfft%nnr,ii,ki)) * &
                                       Psivira(1:dfft%nnr,ibnd,ka) / omega 

            if(okvan .and. tqr) then
              ! Orbital pairs in R
              call addusxx_r(Vr(:),becocc(ki)%k(:,ii),&
                                         becpsia(ka)%k(:,ibnd))
            endif

            ! fwfft orbital pairs to G
            CALL fwfft ('Rho', Vr(:), dfft)

            if(okvan .and. .not.tqr) &
              CALL addusxx_g(dfft, Vr(:), xksym(1:3,ki), &
                  xksym(1:3,ka),'c',becphi_c=becocc(ki)%k(:,ii),&
                  becpsi_c=becpsia(ka)%k(:,ibnd))
  
            ! multiply by FFT[ 1/r ]
            Kia(1:ngm, ibnd, ka) = Vr(dfft%nl(1:ngm)) * vCoulG(1:ngm) * e2Ha / nksym

            if(okpaw) then
!              call contract_paw_one_center(ke,paw_one_center(:,ibnd,ka), &
!                  becpsia(ik)%k(:,ia-a_beg+1),becpsib(ik)%k(:,ib-b_beg+1))
!              paw_one_center(:,ibnd,ik) = paw_one_center(:,ibnd,ik) / (1.d0*nksym)
            endif

          enddo

          do ibnd=1,nborb

            ! Orbital pairs in R
            Vr(1:dfft%nnr) = CONJG(Psiocc(1:dfft%nnr,ii,ki)) * &
                                       Psivirb(1:dfft%nnr,ibnd,ka) / omega 

            if(okvan .and. tqr) then
              ! Orbital pairs in R
              call addusxx_r(Vr(:),becocc(ki)%k(:,ii),&
                                         becpsib(ka)%k(:,ibnd))
            endif

            ! fwfft orbital pairs to G
            CALL fwfft ('Rho', Vr(:), dfft)

            if(okvan .and. .not.tqr) &
              CALL addusxx_g(dfft, Vr(:), xksym(1:3,ki), &
                  xksym(1:3,ka),'c',becphi_c=becocc(ki)%k(:,ii),&
                  becpsi_c=becpsib(ka)%k(:,ibnd))

            ! multiply by FFT[ 1/r ]
            Kib(1:ngm, ibnd, ka) = Vr(dfft%nl(1:ngm)) * vCoulG(1:ngm) * e2Ha / nksym

            if(okpaw) then
!              call contract_paw_one_center(ke,paw_one_center(:,ibnd,ka), &
!                  becpsia(ik)%k(:,ia-a_beg+1),becpsib(ik)%k(:,ib-b_beg+1))
!              paw_one_center(:,ibnd,ik) = paw_one_center(:,ibnd,ik) / (1.d0*nksym)
            endif

          enddo

          IF ( okvan .and..not.tqr ) CALL qvan_clean ()

        end do
        call stop_clock( 'mp2_Kia' )

        do Q=1,nksym

          do kb=1,nksym

            ka = QKtoK2(Q, ki)
            kj = QKtoK2(Q, kb)

            ! (ia|jb) 
            dQ(1:3) = xksym(1:3, ka) - xksym(1:3, ki) + &
                      xksym(1:3, kb) - xksym(1:3, kj)
            kk=0
            do jj=1,27
              if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
                kk=jj
                exit
              endif
            enddo
            if(kk.lt.1) call errore('mp2','Can not find dQ in G list.',1)

            do ij=1,noccK(kj,ispin)

              IF ( okvan .and..not.tqr ) &
                CALL qvan_init (dfft%ngm, xksym (1:3, kj), xksym (1:3, kb))  

              ! Vab = (ia | jb) = sum_r Kia(r) * conjg(psi_j(r)) * psi_b(r) * exp(iQ) 
              do ib=1,nborb
                if(okvan) then
                  Vr(:) = (0.d0,0.d0)
                  if( tqr ) then
                    call addusxx_r(Vr(1:nxxs),becocc(kj)%k(:,ij),&
                                  becpsib(kb)%k(:,ib))
                  else
                    CALL addusxx_g(dfft, Vr(1:nxxs), xksym(1:3,kj), &
                      xksym(1:3,kb),'c',becphi_c=becocc(kj)%k(:,ij), &
                      becpsi_c=becpsib(kb)%k(:,ib))
                    CALL invfft ('Rho', Vr(1:nxxs), dfft)
                  endif
                  Vr2(1:nxxs) = Vr(1:nxxs) * omega &
                          + CONJG(Psiocc(1:nxxs,ij,kj)) * &
                            Psivirb(1:nxxs,ib,kb) * phasefac(1:nxxs, kk) 
                else
                  Vr2(1:nxxs) = CONJG(Psiocc(1:nxxs,ij,kj)) * &
                            Psivirb(1:nxxs,ib,kb) * phasefac(1:nxxs, kk) 
                endif

                Vr2(1:nxxs) = CONJG(Vr2(1:nxxs))

                ! fwfft orbital pairs to G
                CALL fwfft ('Rho', Vr2(:), dfft)

                Rgb(1:ngm,ib) = CONJG(Vr2(dfft%nl(1:ngm)))

              enddo

              IF ( okvan .and..not.tqr ) CALL qvan_clean()

              call start_clock( 'mp2_mm' )
              call zgemm('T','N',naorb,nborb,ngm,(1.d0,0.d0),Kia(1,1,ka),ngm,  &
                   Rgb(1,1),ngm,(0.d0,0.d0),Vab(1,1),naorb)
              call stop_clock( 'mp2_mm' )

              IF ( okvan .and..not.tqr ) &
                CALL qvan_init (dfft%ngm, xksym (1:3, kj), xksym (1:3, ka))  

              ! Vba = (ib | ja) = sum_r Kib(r) * conjg(psi_j(r)) * psi_a(r) * exp(iQ) 
              do ia=1,naorb
                if(okvan) then
                  Vr(:) = (0.d0,0.d0)
                  if( tqr ) then
                    call addusxx_r(Vr(1:nxxs),becocc(kj)%k(:,ij),&
                                  becpsia(ka)%k(:,ia))
                  else
                    CALL addusxx_g(dfft, Vr(1:nxxs), xksym(1:3,kj), &
                      xksym(1:3,kb),'c',becphi_c=becocc(kj)%k(:,ij), &
                      becpsi_c=becpsia(ka)%k(:,ia))
                    CALL invfft ('Rho', Vr(1:nxxs), dfft)
                  endif
                  Vr2(1:nxxs) = Vr(1:nxxs) * omega  &
                        + CONJG(Psiocc(1:nxxs,ij,kj)) * &
                          Psivira(1:nxxs,ia,ka) * phasefac(1:nxxs, kk) 
                else
                  Vr2(1:nxxs) = CONJG(Psiocc(1:nxxs,ij,kj)) * &
                          Psivira(1:nxxs,ia,ka) * phasefac(1:nxxs, kk) 
                endif

                Vr2(1:nxxs) = CONJG(Vr2(1:nxxs))

                ! fwfft orbital pairs to G
                CALL fwfft ('Rho', Vr2(:), dfft)

                Rgb(1:ngm,ia) = CONJG(Vr2(dfft%nl(1:ngm)))

              enddo

              IF ( okvan .and..not.tqr ) CALL qvan_clean()

              call start_clock( 'mp2_mm' )
              call zgemm('T','N',nborb,naorb,ngm,(1.d0,0.d0),Kib(1,1,kb),ngm,  &
                     Rgb(1,1),ngm,(0.d0,0.d0),Vba(1,1),nborb)
              call stop_clock( 'mp2_mm' )

              call start_clock( 'mp2_abij' )
              do ia=1,naorb
                if( noccK(ka,1) + a_beg - 1 + ia > norbK(ka) ) exit
                if(regularize) then
                  do ib=1,nborb
                    if( noccK(kb,1) + b_beg - 1 + ib > norbK(kb) ) cycle
! add one center terms here!
                    etemp = ( eigval(noccK(ka,1)+a_beg-1+ia,ka) + &
                                            eigval(noccK(kb,1)+b_beg-1+ib,kb) - &
                                            eigval(ii,ki) - eigval(ij,kj) )

                    etemp2 = ((1.d0 - exp(-regkappa*etemp))**(1.d0*regp)) / etemp
                    emp2 = emp2 - Vab(ia,ib) * etemp2 *  &
                                        conjg( 2.d0*Vab(ia,ib) - Vba(ib,ia) )
                  enddo
                else
                  do ib=1,nborb
                    if( noccK(kb,1) + b_beg - 1 + ib > norbK(kb) ) cycle
! add one center terms here!
                    etemp = 1.d0/( eigval(noccK(ka,1)+a_beg-1+ia,ka) + &
                                            eigval(noccK(kb,1)+b_beg-1+ib,kb) - &
                                            eigval(ii,ki) - eigval(ij,kj) )
                    emp2 = emp2 - etemp * Vab(ia,ib) * & 
                                        conjg( 2.d0*Vab(ia,ib) - Vba(ib,ia) ) 
                  enddo
                endif
              enddo
              call stop_clock( 'mp2_abij' )

            enddo

          enddo

        enddo

      enddo
      CALL stop_clock ( 'mp2' )
      CALL print_clock( 'mp2' )

    enddo

    if(nproc_image > 1) CALL mp_sum ( emp2, intra_image_comm ) 
    emp2 = emp2/(1.d0*nksym)
    write(*,*) '  EMP2 (Ha): ',emp2
    IF ( ionode .and. verbose ) THEN
      !
      WRITE( 6, * )
      !
      CALL print_clock ( 'mp2_Kia' )
      CALL print_clock ( 'mp2_io' )
      CALL print_clock ( 'mp2_mm' )
      CALL print_clock ( 'mp2_abij' )
      !
    ENDIF

    if(allocated(norbK)) deallocate(norbK)
    IF( ALLOCATED(Kia) ) DEALLOCATE (Kia)
    IF( ALLOCATED(Kib) ) DEALLOCATE (Kib)
    IF( ALLOCATED(Psiocc) ) DEALLOCATE (Psiocc)
    IF( ALLOCATED(Psivira) ) DEALLOCATE (Psivira)
    IF( ALLOCATED(Psivirb) ) DEALLOCATE (Psivirb)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(noccK) ) DEALLOCATE(noccK)
    IF( ALLOCATED(nvirK) ) DEALLOCATE(nvirK)
    IF( ALLOCATED(Vab) ) DEALLOCATE(Vab)
    IF( ALLOCATED(Vba) ) DEALLOCATE(Vba)
    IF( ALLOCATED(Vr) ) DEALLOCATE(Vr)
    IF( ALLOCATED(Vr2) ) DEALLOCATE(Vr2)
    IF( ALLOCATED(Rgb) ) DEALLOCATE(Rgb)
    if(allocated(eigval)) deallocate(eigval)
    if(allocated(weight)) deallocate(weight)
    if(allocated(paw_one_center)) deallocate(paw_one_center)
    if(okvan) then
      do ik=1,nksym
        CALL deallocate_bec_type(becpsia(ik))
        CALL deallocate_bec_type(becpsib(ik))
        CALL deallocate_bec_type(becocc(ik))
      enddo
      DEALLOCATE(becpsia)
      DEALLOCATE(becpsib)
      DEALLOCATE(becocc)
    endif
    if(okvan) then
      do i=1,ntyp
        deallocate(ke(i)%L)
      enddo
      deallocate(ke)
    endif

  END SUBROUTINE mp2_g

#if defined(__CUDA)
  SUBROUTINE mp2_gpu(emp2,dfft,esh5_file,reg_pow,reg_expo)
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE gvect, ONLY : ecutrho
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE paw_variables,           ONLY : okpaw
    USE becmod,  ONLY : bec_type, ALLOCATE_bec_type, DEALLOCATE_bec_type
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    !
    USE cublas
    USE cudafor
    USE device_util_m,           ONLY: dev_memcpy, dev_memset
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: esh5_file
    INTEGER, INTENT(IN), OPTIONAL :: reg_pow
    REAL(DP), INTENT(IN), OPTIONAL :: reg_expo
    COMPLEX(DP), INTENT(OUT) :: emp2
    !
    INTEGER :: iq, Q, ka, kb, iab, ia, ib, ii, ij, ki, kj, ni, nj, iuv
    INTEGER :: ik, ibnd, nxxs, i, j, kk0, ikk, kk, jj, ia0, iu0
    INTEGER :: ispin, is1, is2, noa, nob, n
    INTEGER :: a_beg, a_end, b_beg, b_end, ab_beg, ab_end, nabpair  
    INTEGER :: naorb, nborb, iabmax
    INTEGER :: maxv_rank, error
    INTEGER :: nel(2), maxocc, maxvir, nke_dim, maxocc2
    INTEGER, ALLOCATABLE :: noccK(:,:), nvirK(:,:) 
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: regkappa
    INTEGER :: regp
    REAL(DP) :: dQ(3),dQ1(3),dk(3), dG(3, 27), scl, qcut(5)
    REAL(DP) :: residual, maxv, fXX, fac
    REAL(DP) :: etemp, etemp2, q2x
    COMPLEX(DP) :: one, zero
    LOGICAL :: first,regularize
    ! QPOINT stuff
    INTEGER, ALLOCATABLE :: norbK(:)
    REAL(DP), ALLOCATABLE :: weight(:,:),eigval(:,:)
    COMPLEX(DP), ALLOCATABLE :: Kia(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Vr(:),Vr2(:)    
    COMPLEX(DP), ALLOCATABLE :: Psiocc(:,:,:)   
    COMPLEX(DP), ALLOCATABLE :: Psivira(:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: Psivirb(:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Rgb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vab(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vba(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Buff(:,:,:)  
    TYPE(ke_factorization), ALLOCATABLE :: ke(:) ! for paw one-center terms
    COMPLEX(DP),ALLOCATABLE :: paw_one_center(:,:,:)
    TYPE(bec_type),ALLOCATABLE :: becpsia(:)
    TYPE(bec_type),ALLOCATABLE :: becpsib(:)
    TYPE(bec_type),ALLOCATABLE :: becocc(:)
    CHARACTER(len=10) read_type
    CHARACTER(len=256) h5name
    type(h5file_type) :: h5id_input_orbs
    !
    INTEGER :: nfft, many_fft
    INTEGER, ALLOCATABLE :: nl_d(:)
    COMPLEX(DP), ALLOCATABLE :: eQ(:,:)
    COMPLEX(DP), ALLOCATABLE :: Kia_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: Kib_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: Vr_d(:)
    COMPLEX(DP), ALLOCATABLE :: Psiocc_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: Psivira_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: Psivirb_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: vCoulG_d(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: phasefac_d(:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Rgb_d(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vab_d(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vba_d(:,:)  
    !
    attributes(device) :: nl_d,Kia_d,Kib_d,Vr_d,Psiocc_d,Psivira_d,Psivirb_d
    attributes(device) :: vCoulG_d,phasefac_d,Rgb_d,Vab_d,Vba_d 
    !
    type(cublasHandle) :: handle_cublas
    integer :: istat
    istat = cublasCreate(handle_cublas)
    if (istat .ne. CUBLAS_STATUS_SUCCESS) print *,istat
    ! fix this!!!
    call start_clock( 'mp2_setup' )
    many_fft = 16
    one = (1.0,0.0)
    zero = (0.0,0.0)
    regp=1
    regkappa=0.d0
    if(present(reg_pow)) regp = max(1,reg_pow)
    if(present(reg_expo)) regkappa = max(0.d0,reg_expo) 
    regularize = (regkappa > 0.d0) 

    nxxs = dfft%nr1x*dfft%nr2x*dfft%nr3x

    if(nspin == 1 ) then
      fac = 2.d0 
    else
      fac = 1.d0 
    endif
    if(noncolin) call errore('mp2_gpu','No nonconlin yet.',1)

    IF(tqr) THEN
      tabxx => tabp
    ENDIF

    allocate( noccK(nksym,numspin), nvirK(nksym,numspin), norbK(nksym) )
    if(present(esh5_file)) then
      h5name = TRIM(esh5_file)
      call open_esh5_read(h5id_input_orbs,h5name)
      norbK(:) = h5id_input_orbs%norbK(:)
      allocate(eigval(h5id_input_orbs%maxnorb,nkstot), &
                      weight(h5id_input_orbs%maxnorb,nkstot))  
      weight(:,:) = 0.d0
      eigval(:,:) = 0.d0
      do ik=1,nkstot
        call esh5_posthf_read_et(h5id_input_orbs%id,ik-1,eigval(1,ik), &
                                 weight(1,ik),error)  
        if(error .ne. 0 ) &
          call errore('mp2_g','error reading weights',1)
      enddo  
      ! find number of electrons
      call get_noccK(noccK,nel,h5id_input_orbs%maxnorb,nksym,numspin, &
                     weight,h5id_input_orbs%maxnorb)
    else
      if(npool > 1) call errore('mp2_g', &
             'Error: npool > 1 not allowed in mp2_g without esh5 orbitals',1)
      allocate(eigval(nbnd,nkstot))
      eigval(:,:) = et(:,:)*e2Ha
      ! find number of electrons
      call get_noccK(noccK,nel,minval(norbK(:)),nksym,numspin,wg,nbnd,wk)
      norbK(:) = norb
    endif
    write(*,*) ' Number of electrons per spin channel: ',(nel(i),i=1,numspin)
    nvirK(:,1) = norbK(:) - noccK(:,1)
    if(numspin==2) &
      nvirK(:,2) = norbK(:) - noccK(:,2)
    maxocc = maxval(noccK(1:nksym,:))
    maxvir = maxval(nvirK(1:nksym,:))

    ! limiting to insulators for now, need changes in case of metals
    do ispin=1,numspin
      do ik=2,nksym
        if(noccK(ik,ispin) .ne. noccK(1,ispin)) &
          call errore('mp2_g','Error: Only insulators for now!!!',1) 
      enddo
    enddo

    ! MAM: parallelization over bands only, regardless of value of npools
    call find_2d_partition(maxvir,me_image+1,nproc_image,a_beg,a_end,b_beg,b_end)
    naorb   = a_end-a_beg+1
    nborb   = b_end-b_beg+1
    if( (naorb < 1) .or. (nborb < 1) ) &
      call errore('mp2_g','Error: Too many processors.',1)

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over nxxs and gather over proc grid
! MAM: NOTE: can store 3 components separately and multiply later on, 
! reduces memory cost from nx*ny*nz to 3*max(nx,ny,nz)
! application might be slow though
    allocate( phasefac(nxxs,27) )
    CALL calculate_phase_factor(dfft, phasefac, dG)

    if(verbose .and. ionode) then
      write(*,*) ' Summary - naorb, nborb, npwx, ngm, nxxs, nks, nspin:',  &
                        naorb,nborb,npwx,ngm,nxxs,nksym,nspin
      write(*,*) ' Host Allocations (in MB): '
      write(*,*) '     - Kia: ',ngm*naorb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Kib: ',ngm*nborb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Psivira:',nxxs*naorb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Psivirb:',nxxs*nborb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Psiocc: ',nxxs*maxocc*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Vab: ',2*naorb*nborb*16.0/1024.0/1024.0
      write(*,*) '     - Rgb: ',ngm*max(naorb,nborb)*16.0/1024.0/1024.0
      write(*,*) '     - Buff: ',nxxs*max(maxocc,nborb)*16.0/1024.0/1024.0
      write(*,*) ' Device Allocations (in MB): '
      write(*,*) '     - Kia: ',ngm*naorb*16.0/1024.0/1024.0
      write(*,*) '     - Kib: ',ngm*nborb*16.0/1024.0/1024.0
      write(*,*) '     - Psivira:',nxxs*naorb*16.0/1024.0/1024.0
      write(*,*) '     - Psivirb:',nxxs*nborb*16.0/1024.0/1024.0
      write(*,*) '     - Psiocc: ',nxxs*maxocc*16.0/1024.0/1024.0
      write(*,*) '     - Vab: ',2*naorb*nborb*16.0/1024.0/1024.0
      write(*,*) '     - Rgb: ',ngm*max(naorb,nborb)*16.0/1024.0/1024.0
      FLUSH( stdout )
    endif
    if(nproc_image > 1) call mp_barrier( intra_image_comm )

    if(ionode) write(*,*) 'Allocating arrays'
    !
    emp2 = (0.d0,0.d0)

    ! MAM: since naorb orbs are only updated inside the Q loop, 
    !      you can read them from file. This means that you can minimize
    !      memory usage by making nborb as small as possible 
    !      (instead of equal to naorb like you do now) 
    allocate( Kia(ngm,naorb,nksym), Kib(ngm,nborb,nksym) )
    allocate( Psivira(nxxs,naorb,nksym), Psiocc(nxxs,maxocc,nksym) )
    allocate( Psivirb(nxxs,nborb,nksym) )
    allocate( Vab(naorb,nborb), Vba(nborb,naorb), Rgb(ngm,max(naorb,nborb)) )
    allocate( vCoulG(ngm), Vr(nxxs), Vr2(nxxs), eQ(nQuniq,3) )
    allocate( Buff(nxxs,max(maxocc,nborb),1) )
    eQ(:,:) = (0.d0,0.d0)

    ! allocate gpu arrays
    allocate( Kia_d(ngm,naorb), Kib_d(ngm,nborb) )
    allocate( Psivira_d(nxxs,naorb), Psiocc_d(nxxs,maxocc) )
    allocate( Psivirb_d(nxxs,nborb) )
    allocate( Vab_d(naorb,nborb), Vba_d(nborb,naorb), Rgb_d(ngm,max(naorb,nborb)) )
    allocate( vCoulG_d(ngm), Vr_d(nxxs*many_fft), nl_d(ngm) )
     
    nl_d = dfft%nl
    !
    call stop_clock( 'mp2_setup' )

    if(ionode) write(*,*) 'Reading orbitals from file'
    nke_dim = 0
    if(present(esh5_file)) then
      read_type = 'esh5'
    else
      read_type = 'davcio'
    endif
    if(okvan.or.okpaw) then
      ! store bec for Psia 
      ALLOCATE(becpsia(nksym))
      ALLOCATE(becpsib(nksym))
      ALLOCATE(becocc(nksym))
      do ik=1,nksym
        CALL allocate_bec_type( nkb, nborb, becpsib(ik))
        CALL allocate_bec_type( nkb, maxocc, becocc(ik))
      enddo
      if(okpaw) then
        allocate(ke(ntyp))
        call calculate_factorized_paw_one_center(ke,nke_dim)
        !allocate( paw_one_center(nke_dim,nabpair,nksym) )
      endif
    endif

    ! need eigenvalues, right now assuming they are given by pwscf files
    if(ionode) write(*,*) 'Starting loop over bands' 

    do ispin=1,min(nspin,2)

      ! only for insulators, so all kpoints must have same occupation numbers
      maxocc = noccK(1,ispin)
      kk0 = nksym*(ispin-1)
      call start_clock( 'mp2_io' )
      if(okvan.or.okpaw) then
       call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          ispin,Psiocc,1,maxocc,1,nksym,becpsi=becocc)
       call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          ispin,Psivira,maxocc+a_beg,naorb,1,nksym,becpsi=becpsia)
       call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft, &
                          ispin,Psivirb,maxocc+b_beg,nborb,1,nksym,becpsi=becpsib)
      else
       call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          ispin,Psiocc,1,maxocc,1,nksym)
       call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          ispin,Psivira,maxocc+a_beg,naorb,1,nksym)
       call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          ispin,Psivirb,maxocc+b_beg,nborb,1,nksym)
      endif
      call stop_clock( 'mp2_io' )

      do ki=1,nksym

      if(ionode) then
        write(*,*) 'K,spin:',ki,ispin
        FLUSH(6)
      endif  
      CALL start_clock ( 'mp2' )
      do ii=1,noccK(ki,ispin)
    
        call start_clock( 'mp2_Kia' )
        ! calculate Kia/Kib
        Kia(:,:,:)=(0.d0,0.d0)
        Kib(:,:,:)=(0.d0,0.d0)
        Psiocc_d = Psiocc(:,:,ki)
        do ka=1,nksym

          ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(ka) - k(ki)  |^{-2}
          CALL g2_convolution(dfft%ngm, g, xksym (1:3, ka), xksym (1:3, ki), vCoulG)
          vCoulG_d = vCoulG
        
!          IF ( okvan .and..not.tqr ) &
!            CALL qvan_init (dfft%ngm, xksym (1:3, ki), xksym (1:3, ka))  

          Psivira_d = Psivira(:,:,ka)
          Psivirb_d = Psivirb(:,:,ka)
          Kia_d(:,:)=(0.d0,0.d0)
          Kib_d(:,:)=(0.d0,0.d0)
          do ibnd=1,naorb,many_fft

            nfft = min( many_fft, naorb-ibnd+1 )
            ! Orbital pairs in R
!$cuf kernel do(2)
            do j=0,nfft-1
              do i=1,nxxs
                Vr_d(i+nxxs*j) = CONJG(Psiocc_d(i,ii)) * &
                                     Psivira_d(i,ibnd+j) / omega 
              enddo
            enddo    

            ! fwfft orbital pairs to G
            CALL fwfft ('Rho', Vr_d, dfft, howmany=nfft)

            ! multiply by FFT[ 1/r ]
!$cuf kernel do(2)
            do j=0,nfft-1
              do i=1,ngm
                Kia_d(i, ibnd+j) = Vr_d(nl_d(i)+j*nxxs) * vCoulG_d(i) * e2Ha / nksym
              enddo
            enddo

          enddo
          Kia(:,:,ka) = Kia_d(:,:)  

          do ibnd=1,nborb,many_fft

            nfft = min( many_fft, nborb-ibnd+1 )
            ! Orbital pairs in R
!$cuf kernel do(2)
            do j=0,nfft-1
              do i=1,nxxs
                Vr_d(i+j*nxxs) = CONJG(Psiocc_d(i,ii)) * &
                                       Psivirb_d(i,ibnd+j) / omega 
              enddo
            enddo

            ! fwfft orbital pairs to G
            CALL fwfft ('Rho', Vr_d(:), dfft, howmany=nfft)

            ! multiply by FFT[ 1/r ]
!$cuf kernel do(2)
            do j=0,nfft-1
              do i=1,ngm
                Kib_d(i, ibnd+j) = Vr_d(nl_d(i)+j*nxxs) * vCoulG_d(i) * e2Ha / nksym
              enddo
            enddo

          enddo
          Kib(:,:,ka) = Kib_d(:,:)  

!          IF ( okvan .and..not.tqr ) CALL qvan_clean ()

        end do
        call stop_clock( 'mp2_Kia' )

        do iq=1,nQuniq

          Q = xQ(iq)
          ka = QKtoK2(Q, ki)
          ! move calculation of Kia(:,:,ka) here, do not store Kia for all ka!!!
          Kia_d = Kia(:,:,ka)
          Psivira_d = Psivira(:,:,ka)

          do kb=1,nksym

            kj = QKtoK2(Q, kb)

            ! (ia|jb) 
            dQ(1:3) = xksym(1:3, ka) - xksym(1:3, ki) + &
                      xksym(1:3, kb) - xksym(1:3, kj)
            kk=0
            do jj=1,27
              if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
                kk=jj
                exit
              endif
            enddo
            if(kk.lt.1) call errore('mp2','Can not find dQ in G list.',1)

            ! copy to gpu
            call start_clock( 'mp2_memcpy' )
            Kib_d = Kib(:,:,kb)
            Psiocc_d = Psiocc(:,:,kj)
            Psivirb_d = Psivirb(:,:,kb)
            phasefac_d = phasefac(:,kk)
            call stop_clock( 'mp2_memcpy' )

            do ij=1,noccK(kj,ispin)

              call start_clock( 'mp2_Rgb' )

              ! Vab = (ia | jb) = sum_r Kia(r) * conjg(psi_j(r)) * psi_b(r) * exp(iQ) 
              do ib=1,nborb,many_fft
                nfft = min( many_fft, nborb-ib+1 )
!$cuf kernel do(2)
                   do j=0,nfft-1
                     do i=1,nxxs
                       Vr_d(i+j*nxxs) = Psiocc_d(i,ij) * &
                            CONJG(Psivirb_d(i,ib+j) * phasefac_d(i)) 
                     enddo  
                   enddo  

                ! fwfft orbital pairs to G
                call start_clock( 'mp2_fwfft' )
                CALL fwfft ('Rho', Vr_d(:), dfft, howmany=nfft)
                call stop_clock( 'mp2_fwfft' )

!$cuf kernel do(2)
                do j=0,nfft-1
                  do i=1,ngm
                    Rgb_d(i,ib+j) = CONJG(Vr_d(nl_d(i)+j*nxxs))
                  enddo
                enddo

              enddo

!              IF ( okvan .and..not.tqr ) CALL qvan_clean()
              call stop_clock( 'mp2_Rgb' )

              call start_clock( 'mp2_mm' )
              istat = cublasZgemm_v2(handle_cublas,CUBLAS_OP_T,CUBLAS_OP_N,&
                    naorb,nborb,ngm,one,Kia_d,ngm,  &
                    Rgb_d,ngm,zero,Vab_d,naorb)
              if (istat .ne.CUBLAS_STATUS_SUCCESS) print *,istat
              call stop_clock( 'mp2_mm' )

              call start_clock( 'mp2_memcpy' )
              CALL dev_memcpy(Vab,Vab_d)
              call stop_clock( 'mp2_memcpy' )

              call start_clock( 'mp2_Rgb' )

              ! Vba = (ib | ja) = sum_r Kib(r) * conjg(psi_j(r)) * psi_a(r) * exp(iQ) 
              do ia=1,naorb,many_fft
                nfft = min( many_fft, naorb-ia+1 )
!$cuf kernel do(2)
                   do j=0,nfft-1
                     do i=1,nxxs
                       Vr_d(i+j*nxxs) = Psiocc_d(i,ij) * &
                            CONJG(Psivira_d(i,ia+j) * phasefac_d(i))
                     enddo
                   enddo

                ! fwfft orbital pairs to G
                call start_clock( 'mp2_fwfft' )
                CALL fwfft ('Rho', Vr_d, dfft, howmany=nfft)
                call stop_clock( 'mp2_fwfft' )

!$cuf kernel do(2)
                do j=0,nfft-1
                  do i=1,ngm
                    Rgb_d(i,ia+j) = CONJG(Vr_d(nl_d(i)+j*nxxs))
                  enddo
                enddo

              enddo

!              IF ( okvan .and..not.tqr ) CALL qvan_clean()
              call stop_clock( 'mp2_Rgb' )

              call start_clock( 'mp2_mm' )
              istat = cublasZgemm_v2(handle_cublas,CUBLAS_OP_T,CUBLAS_OP_N,&
                    nborb,naorb,ngm,one,Kib_d,ngm,  &
                    Rgb_d,ngm,zero,Vba_d,nborb)
              if (istat .ne.CUBLAS_STATUS_SUCCESS) print *,istat
              call stop_clock( 'mp2_mm' )

              call start_clock( 'mp2_memcpy' )
              CALL dev_memcpy(Vba,Vba_d)
              call stop_clock( 'mp2_memcpy' )

              call start_clock( 'mp2_abij' )
              do ia=1,naorb
                if( maxocc + a_beg - 1 + ia > norbK(ka) ) exit
                if(regularize) then
                  do ib=1,nborb
                    if( maxocc + b_beg - 1 + ib > norbK(kb) ) cycle
! add one center terms here!
                    etemp = ( eigval(maxocc+a_beg-1+ia,ka+kk0) + &
                                            eigval(maxocc+b_beg-1+ib,kb+kk0) - &
                                            eigval(ii,ki+kk0) - eigval(ij,kj+kk0) )

                    etemp2 = wQ(iq) * ((1.d0 - exp(-regkappa*etemp))**(1.d0*regp)) / etemp
                    emp2 = emp2 - Vab(ia,ib) * etemp2 *  &
                                        conjg( fac*Vab(ia,ib) - Vba(ib,ia) )
                    eQ(iq,1) = eQ(iq,1) - Vab(ia,ib) * etemp2 *  &
                                        conjg( fac*Vab(ia,ib) ) 
                    eQ(iq,2) = eQ(iq,2) + Vab(ia,ib) * etemp2 *  &
                                        conjg( Vba(ib,ia) )
                  enddo
                else
                  do ib=1,nborb
                    if( maxocc + b_beg - 1 + ib > norbK(kb) ) cycle
! add one center terms here!
                    etemp = wQ(iq)/( eigval(maxocc+a_beg-1+ia,ka+kk0) + &
                                            eigval(maxocc+b_beg-1+ib,kb+kk0) - &
                                            eigval(ii,ki+kk0) - eigval(ij,kj+kk0) )
                    emp2 = emp2 - etemp * Vab(ia,ib) * & 
                                        conjg( fac*Vab(ia,ib) - Vba(ib,ia) ) 
                    eQ(iq,1) = eQ(iq,1) - etemp * Vab(ia,ib) * & 
                                        conjg( fac*Vab(ia,ib) ) 
                    eQ(iq,2) = eQ(iq,2) + etemp * Vab(ia,ib) * & 
                                        conjg( Vba(ib,ia) ) 
                  enddo
                endif
              enddo
              call stop_clock( 'mp2_abij' )

            enddo  ! ij

            ! opposite spin contribution
            if( ispin==1 .and. nspin>=2 ) then

              maxocc2 = noccK(1,2)
              ! copy to gpu
              call start_clock( 'mp2_io' )
              ! no paw/uspp yet!!!
              call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          2,Buff,1,maxocc2,kj,1)
              call stop_clock( 'mp2_io' )
              call start_clock( 'mp2_memcpy' )
              Psiocc_d = Buff(:,1:maxocc2,1)
              call stop_clock( 'mp2_memcpy' )
              call start_clock( 'mp2_io' )
              call get_orbitals_set(h5id_input_orbs,read_type,'psir',dfft,&
                          2,Buff,maxocc2+b_beg,nborb,kb,1)
              call stop_clock( 'mp2_io' )
              call start_clock( 'mp2_memcpy' )
              Psivirb_d = Buff(:,1:nborb,1) 
              call stop_clock( 'mp2_memcpy' )

              do ij=1,noccK(kj,2)

                call start_clock( 'mp2_Rgb' )

                ! Vab = (ia | jb) = sum_r Kia(r) * conjg(psi_j(r)) * psi_b(r) * exp(iQ) 
                do ib=1,nborb,many_fft
                  nfft = min( many_fft, nborb-ib+1 )
!$cuf kernel do(2)
                  do j=0,nfft-1
                    do i=1,nxxs
                      Vr_d(i+j*nxxs) = Psiocc_d(i,ij) * &
                           CONJG(Psivirb_d(i,ib+j) * phasefac_d(i)) 
                    enddo  
                  enddo  

                  ! fwfft orbital pairs to G
                  call start_clock( 'mp2_fwfft' )
                  CALL fwfft ('Rho', Vr_d(:), dfft, howmany=nfft)
                  call stop_clock( 'mp2_fwfft' )

!$cuf kernel do(2)
                  do j=0,nfft-1
                    do i=1,ngm
                      Rgb_d(i,ib+j) = CONJG(Vr_d(nl_d(i)+j*nxxs))
                    enddo
                  enddo

                enddo

                call stop_clock( 'mp2_Rgb' )

                call start_clock( 'mp2_mm' )
                istat = cublasZgemm_v2(handle_cublas,CUBLAS_OP_T,CUBLAS_OP_N,&
                    naorb,nborb,ngm,one,Kia_d,ngm,  &
                    Rgb_d,ngm,zero,Vab_d,naorb)
                if (istat .ne.CUBLAS_STATUS_SUCCESS) print *,istat
                call stop_clock( 'mp2_mm' )

                call start_clock( 'mp2_memcpy' )
                CALL dev_memcpy(Vab,Vab_d)
                call stop_clock( 'mp2_memcpy' )

                call start_clock( 'mp2_abij' )
                do ia=1,naorb
                  if( maxocc + a_beg - 1 + ia > norbK(ka) ) exit
                  if(regularize) then
                    do ib=1,nborb
                      if( maxocc2 + b_beg - 1 + ib > norbK(kb) ) cycle
! add one center terms here!
                      etemp = ( eigval(maxocc+a_beg-1+ia,ka) + &
                                              eigval(maxocc2+b_beg-1+ib,kb+kk0) - &
                                              eigval(ii,ki) - eigval(ij,kj+kk0) )

                      etemp2 = wQ(iq) * ((1.d0 - exp(-regkappa*etemp))**(1.d0*regp)) / etemp
                      emp2 = emp2 - 2.d0 * DBLE(Vab(ia,ib) * etemp2 *  &
                                          conjg( Vab(ia,ib) ))
                      eQ(iq,1) = eQ(iq,1) - 2.d0 * DBLE(Vab(ia,ib) * etemp2 *  &
                                          conjg( Vab(ia,ib) )) 
                    enddo
                  else
                    do ib=1,nborb
                      if( maxocc2 + b_beg - 1 + ib > norbK(kb) ) cycle
! add one center terms here!
                      etemp = wQ(iq)/( eigval(maxocc+a_beg-1+ia,ka) + &
                                            eigval(maxocc2+b_beg-1+ib,kb+kk0) - &
                                            eigval(ii,ki) - eigval(ij,kj+kk0) )
                      emp2 = emp2 - 2.d0 * DBLE(etemp * Vab(ia,ib) * & 
                                          conjg( Vab(ia,ib) )) 
                      eQ(iq,1) = eQ(iq,1) - 2.d0 * DBLE(etemp * Vab(ia,ib) * & 
                                          conjg( Vab(ia,ib) )) 
                    enddo
                  endif
                enddo
                call stop_clock( 'mp2_abij' )

              enddo  ! ij
            
            endif ! opposite spin            

          enddo  ! kb

        enddo  !iq

      enddo ! ii
      CALL stop_clock ( 'mp2' )
      CALL print_clock( 'mp2' )

    enddo ! ki
      write(*,*) 'eQ: ',eQ(1,1),eQ(1,2)  
    enddo ! ispin

    if(nspin==2) then
      emp2 = emp2*0.5d0
      eQ(:,:) = eQ(:,:)*0.50
    endif
    if(nproc_image > 1) CALL mp_sum ( emp2, intra_image_comm ) 
    if(nproc_image > 1) CALL mp_sum ( eQ, intra_image_comm ) 
    write(*,*) '  EMP2 (Ha): ',emp2
    write(*,*) '  EJ (Ha): ',sum(eQ(:,1))
    write(*,*) '  EX (Ha): ',sum(eQ(:,2))
    IF ( ionode .and. verbose ) THEN
      write(*,*) '  Q  EMP2(Q) '
      do iq=1,nQuniq    
        write(*,'(i5,g14.8,"("g14.6,g14.6")","("g14.6,g14.6")")') iq,wQ(iq),&
            eQ(iq,1)/wQ(iq),eQ(iq,2)/wQ(iq)
      enddo
      !
      WRITE( 6, * )
      !
      CALL print_clock ( 'mp2_setup' )
      CALL print_clock ( 'mp2_Kia' )
      CALL print_clock ( 'mp2_io' )
      CALL print_clock ( 'mp2_memcpy' )
      CALL print_clock ( 'mp2_fwfft' )
      CALL print_clock ( 'mp2_Rgb' )
      CALL print_clock ( 'mp2_mm' )
      CALL print_clock ( 'mp2_abij' )
      !
    ENDIF

    if(allocated(norbK)) deallocate(norbK)
    IF( ALLOCATED(Kia) ) DEALLOCATE (Kia)
    IF( ALLOCATED(Kib) ) DEALLOCATE (Kib)
    IF( ALLOCATED(Psiocc) ) DEALLOCATE (Psiocc)
    IF( ALLOCATED(Psivira) ) DEALLOCATE (Psivira)
    IF( ALLOCATED(Psivirb) ) DEALLOCATE (Psivirb)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(noccK) ) DEALLOCATE(noccK)
    IF( ALLOCATED(nvirK) ) DEALLOCATE(nvirK)
    IF( ALLOCATED(Vab) ) DEALLOCATE(Vab)
    IF( ALLOCATED(Vba) ) DEALLOCATE(Vba)
    IF( ALLOCATED(Vr) ) DEALLOCATE(Vr)
    IF( ALLOCATED(Vr2) ) DEALLOCATE(Vr2)
    IF( ALLOCATED(Rgb) ) DEALLOCATE(Rgb)
    if(allocated(eigval)) deallocate(eigval)
    if(allocated(weight)) deallocate(weight)
    if(allocated(paw_one_center)) deallocate(paw_one_center)
    IF( ALLOCATED(eQ) ) DEALLOCATE(eQ)
    IF( ALLOCATED(Buff) ) DEALLOCATE(Buff)
    if(okvan) then
      do ik=1,nksym
        CALL deallocate_bec_type(becpsia(ik))
        CALL deallocate_bec_type(becpsib(ik))
        CALL deallocate_bec_type(becocc(ik))
      enddo
      DEALLOCATE(becpsia)
      DEALLOCATE(becpsib)
      DEALLOCATE(becocc)
    endif
    if(okvan) then
      do i=1,ntyp
        deallocate(ke(i)%L)
      enddo
      deallocate(ke)
    endif

    IF( ALLOCATED(Kia_d) ) DEALLOCATE (Kia_d)
    IF( ALLOCATED(Kib_d) ) DEALLOCATE (Kib_d)
    IF( ALLOCATED(Psiocc_d) ) DEALLOCATE (Psiocc_d)
    IF( ALLOCATED(Psivira_d) ) DEALLOCATE (Psivira_d)
    IF( ALLOCATED(Psivirb_d) ) DEALLOCATE (Psivirb_d)
    IF( ALLOCATED(vCoulG_d) ) DEALLOCATE (vCoulG_d)
    IF( ALLOCATED(phasefac_d) ) DEALLOCATE (phasefac_d)
    IF( ALLOCATED(Vab_d) ) DEALLOCATE(Vab_d)
    IF( ALLOCATED(Vba_d) ) DEALLOCATE(Vba_d)
    IF( ALLOCATED(Vr_d) ) DEALLOCATE(Vr_d)
    IF( ALLOCATED(Rgb_d) ) DEALLOCATE(Rgb_d)
    IF( ALLOCATED(nl_d) ) DEALLOCATE(nl_d)
    istat = cublasDestroy(handle_cublas)
    if (istat .ne.CUBLAS_STATUS_SUCCESS) print *,istat
    if(present(esh5_file)) then
      call close_esh5_read(h5id_input_orbs)
    endif
    !
  END SUBROUTINE mp2_gpu
#endif

  SUBROUTINE mp2no_g(dfft,mp2noFile,nskipvir,eigcut,esh5_file)
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE wavefunctions, ONLY : evc
    USE io_global, ONLY : stdout
    USE paw_variables,        ONLY : okpaw
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,  ONLY : bec_type, ALLOCATE_bec_type, DEALLOCATE_bec_type
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r, qgm
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
!    USE exx, ONLY: ecutfock
!    USE gvect,     ONLY : ecutrho
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    USE read_orbitals_from_file, ONLY: get_orbitals, close_esh5_write, open_esh5_write
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(IN) :: mp2noFile 
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: nskipvir
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: esh5_file
    REAL(DP), INTENT(IN) :: eigcut
    !
    CHARACTER(len=256) :: h5name,read_type 
    INTEGER :: Q, iab, ni, nj, iuv, ia0, ib0, ic0, h5len
    INTEGER :: ia, ib, ic, ii, ij, ka, kc, ki, kj, error
    INTEGER :: ik, ibnd, nxxs, i, ikk, kk, jj, iu0, ispin, is1, is2, noa, nob, n
    INTEGER :: c_beg, c_end, b_beg, b_end 
    INTEGER :: ncorb, nborb
    INTEGER :: maxv_rank, nv, nvmax, minvir, maxnumvir
    INTEGER :: nel(2), maxocc, maxauxvir
    INTEGER, ALLOCATABLE :: noccK(:,:), nauxvirK(:,:), numvir(:)
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: recv(3,3), at0(3,3)
    REAL(DP) :: dQ(3),dQ1(3),dk(3), dG(3, 27), scl
    REAL(DP) :: fac,pnorm, residual, maxv, fXX, ecut
    COMPLEX(DP) :: emp2
    COMPLEX(DP) :: ctemp, etemp, fij
    REAL(DP) :: dca, dcb, eaij, ebij
    ! QPOINT stuff
    INTEGER, ALLOCATABLE :: norbK(:)
    REAL(DP), ALLOCATABLE :: weight(:,:),eigval(:,:)
    COMPLEX(DP), ALLOCATABLE :: Dab(:,:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kic(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Psiocc(:,:,:,:)   
    COMPLEX(DP), ALLOCATABLE :: Psivir(:,:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Rgb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vca(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vcb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:)
    REAL(DP), ALLOCATABLE :: eig(:)   
    REAL(DP), ALLOCATABLE :: xkcart_(:,:)
    TYPE(h5file_type) :: h5id_input_orbs, h5id_output_orbs
    !

    nxxs = dfft%nr1x*dfft%nr2x*dfft%nr3x

    fac=1.d0
    if(nspin==1) fac=2.d0
    if(noncolin) call errore('mp2no_g','No nonconlin yet.',1)

    nvmax = 0
    if(eigcut < 0.0) nvmax = nint(-eigcut)

    if(nspin > 2 ) &
      call errore('mp2no','Error: nspin>2 in mp2no not implemented.',1)

    ! find number of electrons
    allocate( noccK(nksym,numspin), nauxvirK(nksym,numspin), norbK(nksym), &
              numvir(nksym) )
    if(present(esh5_file)) then
      h5name = TRIM(esh5_file)
      call open_esh5_read(h5id_input_orbs,h5name)
      norbK(:) = h5id_input_orbs%norbK(:)
      allocate(eigval(h5id_input_orbs%maxnorb,nkstot), &
               weight(h5id_input_orbs%maxnorb,nkstot))
      do ik=1,nksym
        call esh5_posthf_read_et(h5id_input_orbs%id,ik-1,eigval(1,ik),weight(1,ik),error)
        if(error .ne. 0 ) &
          call errore('mp2no_g','error reading weights',1)
      enddo
      ! find number of electrons
      call get_noccK(noccK,nel,h5id_input_orbs%maxnorb,nksym,numspin,&
                     weight,h5id_input_orbs%maxnorb)
    else
      if(npool > 1) call errore('mp2no', &
         'Error: npool > 1 not allowed in mp2no_g without esh5 orbitals.',1)
      allocate(eigval(nbnd,nkstot))
      eigval(:,:) = et(:,:)*e2Ha
      ! find number of electrons
      call get_noccK(noccK,nel,minval(norbK(:)),nksym,numspin,wg,nbnd,wk)
      norbK(:) = norb
    endif
    write(*,*) ' Number of electrons per spin channel: ',(nel(i),i=1,numspin)
    maxocc = maxval(noccK(1:nksym,1))
    if( maxval(noccK(1:nksym,numspin)) .gt. maxocc ) &
      call errore('mp2no','Error: Maximum occupation in down electron.',1) 

    ! lowset virtual orbital to be rotated
    ! range of a/b is [minvir,norb] and range of c is [noccK+1,norb] 
    nauxvirK(:,:)=0
    minvir = maxocc + max(1,nskipvir) 
    ! check that occupation of minvir is zero for all kpoint/spin
    do ispin=1,numspin
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        if(abs(wk(ikk))>1.d-10) then
          scl = 1.d0/wk(ikk)
        else
          scl = 1.d0
        endif
        if( abs(wg(minvir,ikk)*scl) > 1.d-6 ) & 
          call errore('mp2no','Error: Non zero occupation at minimum virtual state in mp2no.',1)
!        ! put some sort of warning if eigenvalue gap from maxocc to minvir is too small
!        if( abs(et(minvir,ikk) - et(noccK(ikk,ispin),ikk)) < 0.07d0  ) &   ! < 1eV
!          call errore('mp2_g','Error: Small eigenvalue gap, increase minvir.',1)
      enddo
    enddo
    nauxvirK(:,1) = norbK(:) - noccK(:,1)
    if(numspin>1) nauxvirK(:,2) = norbK(:) - noccK(:,2)
    maxauxvir = maxval(nauxvirK(1:nksym,:))
    numvir(:) = norbK(:) - minvir + 1
    maxnumvir = maxval(numvir(:))
    write(*,*) 'Starting virtual: ',minvir
    write(*,*) 'Maxocc, maxauxvir, maxnumvir:',maxocc,maxauxvir,maxnumvir

    ! should do this over rectangular grid maxauxvir x numvir
    call find_2d_partition(maxauxvir,me_image+1,nproc_image,c_beg,c_end,b_beg,b_end)
    ncorb   = c_end-c_beg+1
    nborb   = b_end-b_beg+1    ! careful here, since maxauxvir >= numvir 
    if( (ncorb < 1) .or. (nborb < 1) ) &
      call errore('mp2_g','Error: Too many processors.',1)

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over nxxs and gather over proc grid
    allocate( phasefac(nxxs,27) )
    CALL calculate_phase_factor(dfft, phasefac, dG)

    ! 
    ! Loop over Q points
    !
    ! final normalization of 1.0/nksym applied later to keep thresh consistent 
    ! with single k-point case 
    pnorm = 1.d0 / (omega*nxxs*nksym) 
    emp2 = (0.d0,0.d0)
    ispin = 1  ! alpha for now only

    if(verbose .and. ionode) then
      write(*,*) ' Summary - nborb, ncorb, nxxs, nks, nspin:',nborb,ncorb,nxxs,nksym,nspin
      write(*,*) ' Allocations (in MB): '
      write(*,*) '     - Dab: ',maxnumvir*nborb*nksym*nspin*16.0/1024.0/1024.0
      write(*,*) '     - Kic: ',nxxs*ncorb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Kib: ',nxxs*nborb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Psivir: ',nxxs*maxauxvir*nksym*numspin*16.0/1024.0/1024.0
      write(*,*) '     - Psiocc: ',nxxs*maxocc*nksym*numspin*16.0/1024.0/1024.0
      write(*,*) '     - Vca: ',ncorb*maxnumvir*16.0/1024.0/1024.0
      write(*,*) '     - Vcb: ',ncorb*nborb*16.0/1024.0/1024.0
      write(*,*) '     - Rgb: ',nxxs*maxauxvir*16.0/1024.0/1024.0
      FLUSH( stdout )
    endif
    if(nproc_image > 1) call mp_barrier( intra_image_comm ) 

    ! eventually distribute nxxs and reduce matrices below within groups
    allocate( Dab(maxnumvir,nborb,nksym,nspin) )
    allocate( Kic(nxxs,ncorb,nksym), Kib(nxxs,nborb,nksym) )
    allocate( Psivir(nxxs,maxauxvir,nksym,numspin) )  ! [noccK+1,norb]
    allocate( Psiocc(nxxs,maxocc,nksym,numspin) )
    allocate( Vca(ncorb,maxnumvir), Vcb(ncorb,nborb), Rgb(nxxs,maxauxvir) )
    allocate( vCoulG(nxxs) )

    call start_clock( 'mp2no_io' )
    ! need eigenvalues, right now assuming they are given by pwscf files
    if(present(esh5_file)) then
      read_type = 'esh5'
    else
      read_type = 'davcio'
    endif
    do ispin=1,numspin  
      do ik=1,nksym
        call get_orbitals(h5id_input_orbs,read_type,'psir',dfft,&
                    Psiocc(:,:,ik,ispin),1,maxocc,ik,ispin)
        call get_orbitals(h5id_input_orbs,read_type,'psir',dfft, &
                    Psivir(:,:,ik,ispin),noccK(ik,ispin)+1,numvir(ik),ik,ispin) 
      enddo
    enddo
    call stop_clock( 'mp2no_io' )

    if(ionode) write(*,*) 'Starting loop over K vectors.'

    Dab(:,:,:,:) = (0.d0,0.d0)
    do ki=1,nksym

      if(ionode) write(*,*) 'K:',ki
      CALL start_clock ( 'mp2no' )
    
      do ispin=1,numspin  

        do ii=1,noccK(ki,ispin)
    
          call start_clock( 'mp2no_Kia' )
          ! calculate Kic/Kib
          Kic(:,:,:)=(0.d0,0.d0)
          Kib(:,:,:)=(0.d0,0.d0)
          do kc=1,nksym

            ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(ka) - k(ki)  |^{-2}
            CALL g2_convolution(dfft%ngm, g, xksym (1:3, kc), xksym (1:3, ki),psic)
            vCoulG(:) = (0.d0,0.d0)
            vCoulG( dfft%nl(1:dfft%ngm) ) = e2Ha * psic(1:dfft%ngm)

            do ibnd=1,ncorb

              ! Orbital pairs in R
              Kic(1:dfft%nnr,ibnd,kc) = CONJG(Psiocc(1:dfft%nnr,ii,ki,ispin)) * &
                                       Psivir(1:dfft%nnr,c_beg-1+ibnd,kc,ispin) * pnorm

              ! fwfft orbital pairs to G
              CALL fwfft ('Rho', Kic(:,ibnd,kc), dfft)
  
              ! multiply by FFT[ 1/r ]
              Kic(1:nxxs,ibnd,kc) = Kic(1:nxxs,ibnd,kc) * vCoulG(1:nxxs)

              ! invfft to R
              CALL invfft ('Rho', Kic(:,ibnd,kc), dfft)

            enddo

            do ibnd=1,nborb
 
              ! Orbital pairs in R
              Kib(1:dfft%nnr,ibnd,kc) = CONJG(Psiocc(1:dfft%nnr,ii,ki,ispin)) * &
                                       Psivir(1:dfft%nnr,b_beg-1+ibnd,kc,ispin) * pnorm


              ! fwfft orbital pairs to G
              CALL fwfft ('Rho', Kib(:,ibnd,kc), dfft)

              ! multiply by FFT[ 1/r ]
              Kib(1:nxxs,ibnd,kc) = Kib(1:nxxs,ibnd,kc) * vCoulG(1:nxxs)

              ! invfft to R
              CALL invfft ('Rho', Kib(:,ibnd,kc), dfft)

            enddo

          end do
          call stop_clock( 'mp2no_Kia' )

          do Q=1,nksym

            do ka=1,nksym

              kc = QKtoK2(Q, ki)
              kj = QKtoK2(Q, ka)

              ! (ic|ja) 
              dQ(1:3) = xksym(1:3, kc) - xksym(1:3, ki) + &
                        xksym(1:3, ka) - xksym(1:3, kj)
              kk=0
              do jj=1,27
                if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
                  kk=jj
                  exit
                endif
              enddo
              if(kk.lt.1) call errore('mp2_g','Can not find dQ in G list.',1)

              do ij=1,noccK(kj,ispin)

                fij = wg(ii,ki + nksym*(ispin-1)) * wg(ij,kj + nksym*(ispin-1)) 
                if(abs(wk(ki + nksym*(ispin-1)))>1.d-10) fij = fij/wk(ki + nksym*(ispin-1))
                if(abs(wk(kj + nksym*(ispin-1)))>1.d-10) fij = fij/wk(kj + nksym*(ispin-1))

                ! Vca = (ic | ja) = sum_r Kic(r) * conjg(psi_j(r)) * psi_a(r) * exp(iQ) 
                call start_clock( 'mp2no_rgb' )
                do ia=1,numvir(ka)
                  ia0 = nauxvirK(ka,ispin) - numvir(ka)  ! index of ia on Psivir  
                  Rgb(1:dfft%nnr,ia) = CONJG(Psiocc(1:dfft%nnr,ij,kj,ispin)) * &
                                             Psivir(1:dfft%nnr,ia0+ia,ka,ispin) * &
                                             phasefac(1:dfft%nnr, kk) 
                enddo
                call stop_clock( 'mp2no_rgb' )

                call start_clock( 'mp2no_mm' )
                call zgemm('T','N',ncorb,numvir(ka),nxxs,(1.d0,0.d0),Kic(1,1,kc),nxxs,  &
                       Rgb(1,1),nxxs,(0.d0,0.d0),Vca(1,1),ncorb)
                call stop_clock( 'mp2no_mm' )

                ! Vcb = (jc | ib) = sum_r Kib(r) * conjg(psi_j(r)) * psi_c(r) * exp(iQ) 
                call start_clock( 'mp2no_rgb' )
                do ic=1,ncorb
                  Rgb(1:dfft%nnr,ic) = CONJG(Psiocc(1:dfft%nnr,ij,kj,ispin)) * &
                                             Psivir(1:dfft%nnr,c_beg-1+ic,kc,ispin) * &
                                             phasefac(1:dfft%nnr, kk)
                enddo
                call stop_clock( 'mp2no_rgb' )

                call start_clock( 'mp2no_mm' )
                call zgemm('T','N',ncorb,nborb,nxxs,(1.d0,0.d0),Rgb(1,1),nxxs, &
                      Kib(1,1,ka),nxxs,(0.d0,0.d0),Vcb(1,1),ncorb)
                call stop_clock( 'mp2no_mm' )

                call start_clock( 'mp2no_abij' )
                ! emp2
                do ib=1,nborb
                  ib0 = noccK(ka,ispin) + nauxvirK(ka,ispin) - numvir(ka)
                  if( ib0 + b_beg - 1 + ib > norbK(ka) ) cycle
                  ebij = eigval(ib0+b_beg-1+ib,ka+ nksym*(ispin-1)) - &
                         eigval(ii,ki + nksym*(ispin-1)) - eigval(ij,kj + nksym*(ispin-1))
                  do ic=1,ncorb
                    if( c_beg - 1 + ic > nauxvirK(kc,ispin) ) cycle
                    dcb = 1.d0 / (ebij + eigval(noccK(kc,ispin)+c_beg-1+ic,kc+nksym*(ispin-1)))
                    !dcb = min( dcb, 13.605d0 )
                    emp2 = emp2 - fij * dcb * Vca(ic,b_beg-1+ib) * &
                      conjg( fac*Vca(ic,b_beg-1+ib) - Vcb(ic,ib) )
                  enddo
                enddo

                do ia=1,numvir(ka)
                  ia0 = noccK(ka,ispin) + nauxvirK(ka,ispin) - numvir(ka)
                  eaij = eigval(ia0+ia,ka + nksym*(ispin-1)) - &
                         eigval(ii,ki + nksym*(ispin-1)) - eigval(ij,kj + nksym*(ispin-1)) 
                  do ib=1,nborb
                    if( ia0 + b_beg - 1 + ib > norbK(ka) ) cycle
                    ebij = eigval(ia0+b_beg-1+ib,ka + nksym*(ispin-1)) - &
                        eigval(ii,ki + nksym*(ispin-1)) - eigval(ij,kj + nksym*(ispin-1)) 
                    do ic=1,ncorb
                      if( c_beg - 1 + ic > nauxvirK(kc,ispin) ) cycle
                      dca = 1.d0 / (eaij + eigval(noccK(kc,ispin)+c_beg-1+ic,kc + nksym*(ispin-1)))
                      dcb = 1.d0 / (ebij + eigval(noccK(kc,ispin)+c_beg-1+ic,kc + nksym*(ispin-1)))
!                      dca = min( dca, 13.605d0 )
!                      dcb = min( dcb, 13.605d0 )
                      Dab(ia,ib,ka,ispin) = Dab(ia,ib,ka,ispin) + fij * dca * dcb * Vca(ic,ia) * &
                        conjg( fac*Vca(ic,b_beg-1+ib) - Vcb(ic,ib) ) / nksym / nksym 
                    enddo
                  enddo
                enddo  
                call stop_clock( 'mp2no_abij' )

                if(nspin == 2) then
                  ! opposite spin contribution
  
                endif

              enddo

            enddo  ! ka

          enddo  ! Q

        enddo  ! ii

      enddo  ! ispin
      CALL stop_clock ( 'mp2no' )
      CALL print_clock( 'mp2no' )

    enddo

    if(nproc_image > 1) CALL mp_sum ( emp2, intra_image_comm ) 

    IF( ALLOCATED(Kic) ) DEALLOCATE (Kic)
    IF( ALLOCATED(Kib) ) DEALLOCATE (Kib)
    IF( ALLOCATED(Psiocc) ) DEALLOCATE (Psiocc)
    IF( ALLOCATED(Psivir) ) DEALLOCATE (Psivir)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(noccK) ) DEALLOCATE(noccK)
    IF( ALLOCATED(Vca) ) DEALLOCATE(Vca)
    IF( ALLOCATED(Vcb) ) DEALLOCATE(Vcb)
    IF( ALLOCATED(Rgb) ) DEALLOCATE(Rgb)

    ! now generate NO
    if(ionode) write(*,*) 'Calculating natural orbitals'
    if(nproc_image > 1) call mp_barrier( intra_image_comm ) 
    allocate( Vca(maxnumvir,maxnumvir) )

    if( me_image == root_image) then
      !
      call open_esh5_write(h5id_output_orbs,dfft,mp2noFile,.false.)  
      allocate( Orbitals(npwx,1), eig(maxnumvir), Vcb(npwx,maxnumvir) )
      !
    endif

    do ka=1,nksym
    
      Vca(:,:) = (0.d0,0.d0)
      Vca(1:numvir(ka),b_beg:b_end) = Dab(1:numvir(ka),1:nborb,ka,1)  
      call mp_sum( Vca, intra_image_comm )

      if(me_image == root_image) then

        call start_clock( 'mp2no_io' )
        call get_orbitals(h5id_input_orbs,read_type,'psig',dfft, &
                    Vcb(:,:),minvir,numvir(ka),ka,1)
        call stop_clock( 'mp2no_io' )

        call eigsys('V', 'U', .true., numvir(ka), maxnumvir, Vca, eig) 

        if(nvmax > 0) then
          nv=nvmax
          if(verbose) then
            do ia=1,nv
              write(*,*) ka,ia,eig(ia)
            enddo
          endif
        else
          do ia=1,numvir(ka)
            if(verbose) write(*,*) ka,ia,eig(ia)
            if( dble(eig(ia)) < eigcut ) exit
            nv = ia
          enddo
        endif

        if( size( Orbitals, 2) < nv ) then
          deallocate(Orbitals)
          allocate( Orbitals(npwx,nv) )
        endif
        Orbitals(:,:) = (0.d0,0.d0)
        call zgemm('N','N',ngksym(ka),nv,numvir(ka),(1.d0,0.d0),Vcb(1,1),npwx,  &
                     Vca,maxnumvir,(0.d0,0.d0),Orbitals,npwx)
        ! Write to file
        numvir(ka) = nv
        call esh5_posthf_write(h5id_output_orbs%id,orbsG,5,ka-1,npwx, &
                               numvir(ka),Orbitals,npwx,error)
        if(error .ne. 0 ) &
          call errore('mp2no_g','error writing orbital',1)
        !
        if(numspin .eq. 2) then
          call errore(' Error: lsda not yet implemented in mp2no.',1)
        endif
      endif

    enddo

    if( me_image == root_image) then
      call esh5_posthf_write_norb(h5id_output_orbs%id,orbsG,5,nksym,numvir(1),error)
      if(error .ne. 0 ) &
        call errore('mp2no_g','error writing nvirK OrbsG',1)
      call close_esh5_write(h5id_output_orbs)
    endif
    if(present(esh5_file)) then
      call close_esh5_read(h5id_input_orbs)
    endif
    
    IF( ALLOCATED(numvir) ) DEALLOCATE(numvir)
    IF( ALLOCATED(norbK) ) DEALLOCATE(norbK)
    if(allocated(eigval)) deallocate(eigval)
    if(allocated(weight)) deallocate(weight)
    if(allocated(eig)) deallocate(eig)
    if(allocated(Orbitals)) deallocate(Orbitals)
    IF( ALLOCATED(nauxvirK) ) DEALLOCATE(nauxvirK)
    IF( ALLOCATED(Vca) ) DEALLOCATE(Vca)
    IF( ALLOCATED(Vcb) ) DEALLOCATE(Vcb)
    IF( ALLOCATED(Dab) ) DEALLOCATE(Dab)

    if(ionode) write(*,*) 'EMP2 (Ha): ',emp2/(1.d0*nksym)
    IF ( ionode .and. verbose ) THEN
      !
      WRITE( 6, * )
      !
      CALL print_clock ( 'mp2no_Kia' )
      CALL print_clock ( 'mp2no_io' )
      CALL print_clock ( 'mp2no_mm' )
      CALL print_clock ( 'mp2no_rgb' )
      CALL print_clock ( 'mp2no_abij' )
      !
    ENDIF

  END SUBROUTINE mp2no_g

  SUBROUTINE mp2no(h5file,dfft,nskipvir,eigcut,esh5_file)
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE wavefunctions, ONLY : evc
    USE io_global, ONLY : stdout
    USE paw_variables,        ONLY : okpaw
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,  ONLY : bec_type, ALLOCATE_bec_type, DEALLOCATE_bec_type
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r, qgm
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
!    USE exx, ONLY: ecutfock
!    USE gvect,     ONLY : ecutrho
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*) :: h5file 
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: nskipvir
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: esh5_file
    REAL(DP), INTENT(IN) :: eigcut
    !
    CHARACTER(len=256) :: tmp
    INTEGER :: Q, iab, ni, nj, iuv, ia0, ib0, ic0, h5len
    INTEGER :: ia, ib, ic, ii, ij, ka, kc, ki, kj, error
    INTEGER :: ik, ibnd, nxxs, i, ikk, kk, jj, iu0, ispin, is1, is2, noa, nob, n
    INTEGER :: c_beg, c_end, b_beg, b_end 
    INTEGER :: ncorb, nborb
    INTEGER :: maxv_rank, nv, nvmax, minvir, numvir 
    INTEGER :: nel(2), maxocc, maxauxvir
    INTEGER, ALLOCATABLE :: noccK(:,:), nauxvirK(:,:)
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: recv(3,3), at0(3,3)
    REAL(DP) :: dQ(3),dQ1(3),dk(3), dG(3, 27), scl
    REAL(DP) :: fac,pnorm, residual, maxv, fXX, ecut
    COMPLEX(DP) :: emp2
    COMPLEX(DP) :: ctemp, etemp, fij
    REAL(DP) :: dca, dcb, eaij, ebij
    ! QPOINT stuff
    COMPLEX(DP), ALLOCATABLE :: Dab(:,:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kic(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Psiocc(:,:,:,:)   
    COMPLEX(DP), ALLOCATABLE :: Psivir(:,:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Rgb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vca(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vcb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:)
    REAL(DP), ALLOCATABLE :: eig(:)   
    REAL(DP), ALLOCATABLE :: xkcart_(:,:)
    TYPE(h5file_type) :: h5id_input_orbs, h5id_output_orbs
    !

    nxxs = dfft%nr1x*dfft%nr2x*dfft%nr3x

    if(npool > 1) &
      call errore('mp2no','Error: npool > 1 not allowed in mp2no',1)
    if(noncolin) call errore('mp2no','No nonconlin yet.',1)

    fac=1.d0
    if(nspin==1) fac=2.d0

    nvmax = 0
    if(eigcut < 0.0) nvmax = nint(-eigcut)

    if(nspin > 2 ) &
      call errore('mp2no','Error: nspin>2 in mp2no not implemented.',1)

    ! find number of electrons
    nel(:) = 0
    allocate( noccK(nksym,numspin), nauxvirK(nksym,numspin) )
    noccK(:,:)=0
    nauxvirK(:,:)=0
    do ispin=1,numspin
      pnorm=0.d0
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        if(abs(wk(ikk))>1.d-10) then
            scl = 1.d0/wk(ikk)
        else
            scl = 1.d0
        endif
        do ia=1,norb
          pnorm = pnorm + wg(ia,ikk)*scl
          if( abs(wg(ia,ikk)*scl) > 0.01d0 ) then
            noccK(ik,ispin) = noccK(ik,ispin) + 1
          endif
        enddo
      enddo
      nel(ispin) = nint(pnorm)  
    enddo
    maxocc = maxval(noccK(1:nksym,1))
    if( maxval(noccK(1:nksym,numspin)) .gt. maxocc ) &
      call errore('mp2no','Error: Maximum occupation in down electron.',1) 

    ! lowset virtual orbital to be rotated
    ! range of a/b is [minvir,norb] and range of c is [noccK+1,norb] 
    minvir = maxocc + max(1,nskipvir) 
    ! check that occupation of minvir is zero for all kpoint/spin
    do ispin=1,numspin
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        if(abs(wk(ikk))>1.d-10) then
          scl = 1.d0/wk(ikk)
        else
          scl = 1.d0
        endif
        if( abs(wg(minvir,ikk)*scl) > 1.d-6 ) & 
          call errore('mp2no','Error: Non zero occupation at minimum virtual state in mp2no.',1)
!        ! put some sort of warning if eigenvalue gap from maxocc to minvir is too small
!        if( abs(et(minvir,ikk) - et(noccK(ikk,ispin),ikk)) < 0.07d0  ) &   ! < 1eV
!          call errore('pw2posthf','Error: Small eigenvalue gap, increase minvir.',1)
      enddo
    enddo
    nauxvirK(:,:) = norb - noccK(:,:)
    maxauxvir = maxval(nauxvirK(1:nksym,:))
    numvir = norb - minvir + 1

    ! should do this over rectangular grid maxauxvir x numvir
    call find_2d_partition(maxauxvir,me_pool+1,nproc_pool,c_beg,c_end,b_beg,b_end)
    ncorb   = c_end-c_beg+1
    nborb   = b_end-b_beg+1    ! careful here, since maxauxvir >= numvir 
    if( (ncorb < 1) .or. (nborb < 1) ) &
      call errore('mp2_g','Error: Too many processors.',1)

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over nxxs and gather over proc grid
    allocate( phasefac(nxxs,27) )
    CALL calculate_phase_factor(dfft, phasefac, dG)

    if(ionode) write(*,*) 'Starting loop over Q vectors.'
    ! 
    ! Loop over Q points
    !
    ! final normalization of 1.0/nksym applied later to keep thresh consistent 
    ! with single k-point case 
    pnorm = 1.d0 / (omega*nxxs*nksym) 
    emp2 = (0.d0,0.d0)
    ispin = 1  ! alpha for now only

    if(verbose .and. ionode) then
      write(*,*) ' Summary - nborb, ncorb, nxxs, nks, nspin:',nborb,ncorb,nxxs,nksym,nspin
      write(*,*) ' Allocations (in MB): '
      write(*,*) '     - Dab: ',numvir*nborb*nksym*nspin*16.0/1024.0/1024.0
      write(*,*) '     - Kic: ',nxxs*ncorb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Kib: ',nxxs*nborb*nksym*16.0/1024.0/1024.0
      write(*,*) '     - Psivir: ',nxxs*maxauxvir*nksym*numspin*16.0/1024.0/1024.0
      write(*,*) '     - Psiocc: ',nxxs*maxocc*nksym*numspin*16.0/1024.0/1024.0
      write(*,*) '     - Vca: ',ncorb*numvir*16.0/1024.0/1024.0
      write(*,*) '     - Vcb: ',ncorb*nborb*16.0/1024.0/1024.0
      write(*,*) '     - Rgb: ',nxxs*maxauxvir*16.0/1024.0/1024.0
      FLUSH( stdout )
    endif
    if(nproc_image > 1) call mp_barrier( intra_image_comm ) 

    ! eventually distribute nxxs and reduce matrices below within groups
    allocate( Dab(numvir,nborb,nksym,nspin) )
    allocate( Kic(nxxs,ncorb,nksym), Kib(nxxs,nborb,nksym) )
    allocate( Psivir(nxxs,maxauxvir,nksym,numspin) )  ! [noccK+1,norb]
    allocate( Psiocc(nxxs,maxocc,nksym,numspin) )
    allocate( Vca(ncorb,numvir), Vcb(ncorb,nborb), Rgb(nxxs,maxauxvir) )
    allocate( vCoulG(nxxs) )

    ! need eigenvalues, right now assuming they are given by pwscf files

    do ispin=1,numspin  
      call get_orbitals_set(h5id_input_orbs,'davcio','psir',dfft,ispin,&
                    Psiocc(:,:,:,ispin),1,maxocc,1,nksym)
      call get_orbitals_set(h5id_input_orbs,'davcio','psir',dfft,ispin,&
                    Psivir(:,:,:,ispin),1,norb,1,nksym,nminK=noccK(:,ispin))
    enddo

    Dab(:,:,:,:) = (0.d0,0.d0)
    do ki=1,nksym

      if(ionode) write(*,*) 'K:',ki
      CALL start_clock ( 'mp2no' )
    
      do ispin=1,numspin  

        do ii=1,noccK(ki,ispin)
    
          ! calculate Kic/Kib
          Kic(:,:,:)=(0.d0,0.d0)
          Kib(:,:,:)=(0.d0,0.d0)
          do kc=1,nksym

            ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(ka) - k(ki)  |^{-2}
            CALL g2_convolution(dfft%ngm, g, xksym (1:3, kc), xksym (1:3, ki),psic)
            vCoulG(:) = (0.d0,0.d0)
            vCoulG( dfft%nl(1:dfft%ngm) ) = e2Ha * psic(1:dfft%ngm)

            do ibnd=1,ncorb

              ! Orbital pairs in R
              Kic(1:dfft%nnr,ibnd,kc) = CONJG(Psiocc(1:dfft%nnr,ii,ki,ispin)) * &
                                       Psivir(1:dfft%nnr,c_beg-1+ibnd,kc,ispin) * pnorm

              ! fwfft orbital pairs to G
              CALL fwfft ('Rho', Kic(:,ibnd,kc), dfft)
  
              ! multiply by FFT[ 1/r ]
              Kic(1:nxxs,ibnd,kc) = Kic(1:nxxs,ibnd,kc) * vCoulG(1:nxxs)

              ! invfft to R
              CALL invfft ('Rho', Kic(:,ibnd,kc), dfft)

            enddo

            do ibnd=1,nborb
 
              ! Orbital pairs in R
              Kib(1:dfft%nnr,ibnd,kc) = CONJG(Psiocc(1:dfft%nnr,ii,ki,ispin)) * &
                                       Psivir(1:dfft%nnr,b_beg-1+ibnd,kc,ispin) * pnorm


              ! fwfft orbital pairs to G
              CALL fwfft ('Rho', Kib(:,ibnd,kc), dfft)

              ! multiply by FFT[ 1/r ]
              Kib(1:nxxs,ibnd,kc) = Kib(1:nxxs,ibnd,kc) * vCoulG(1:nxxs)

              ! invfft to R
              CALL invfft ('Rho', Kib(:,ibnd,kc), dfft)

            enddo

          end do

          do Q=1,nksym

            do ka=1,nksym

              kc = QKtoK2(Q, ki)
              kj = QKtoK2(Q, ka)

              ! (ic|ja) 
              dQ(1:3) = xksym(1:3, kc) - xksym(1:3, ki) + &
                        xksym(1:3, ka) - xksym(1:3, kj)
              kk=0
              do jj=1,27
                if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
                  kk=jj
                  exit
                endif
              enddo
              if(kk.lt.1) call errore('mp2no','Can not find dQ in G list.',1)

              do ij=1,noccK(kj,ispin)

                fij = wg(ii,ki + nksym*(ispin-1)) * wg(ij,kj + nksym*(ispin-1)) 
                if(abs(wk(ki + nksym*(ispin-1)))>1.d-10) fij = fij/wk(ki + nksym*(ispin-1))
                if(abs(wk(kj + nksym*(ispin-1)))>1.d-10) fij = fij/wk(kj + nksym*(ispin-1))

                ! Vca = (ic | ja) = sum_r Kic(r) * conjg(psi_j(r)) * psi_a(r) * exp(iQ) 
                do ia=1,numvir
                  ia0 = nauxvirK(ka,ispin) - numvir  ! index of ia on Psivir  
                  Rgb(1:dfft%nnr,ia) = CONJG(Psiocc(1:dfft%nnr,ij,kj,ispin)) * &
                                             Psivir(1:dfft%nnr,ia0+ia,ka,ispin) * &
                                             phasefac(1:dfft%nnr, kk) 
                enddo

                call zgemm('T','N',ncorb,numvir,nxxs,(1.d0,0.d0),Kic(1,1,kc),nxxs,  &
                       Rgb(1,1),nxxs,(0.d0,0.d0),Vca(1,1),ncorb)

                ! Vcb = (jc | ib) = sum_r Kib(r) * conjg(psi_j(r)) * psi_c(r) * exp(iQ) 
                do ic=1,ncorb
                  Rgb(1:dfft%nnr,ic) = CONJG(Psiocc(1:dfft%nnr,ij,kj,ispin)) * &
                                             Psivir(1:dfft%nnr,c_beg-1+ic,kc,ispin) * &
                                             phasefac(1:dfft%nnr, kk)
                enddo

                call zgemm('T','N',ncorb,nborb,nxxs,(1.d0,0.d0),Rgb(1,1),nxxs, &
                      Kib(1,1,ka),nxxs,(0.d0,0.d0),Vcb(1,1),ncorb)

                ! emp2
                do ib=1,nborb
                  ib0 = noccK(ka,ispin) + nauxvirK(ka,ispin) - numvir
                  if( ib0 + b_beg - 1 + ib > norb ) cycle
                  ebij = et(ib0+b_beg-1+ib,ka+ nksym*(ispin-1)) - &
                         et(ii,ki + nksym*(ispin-1)) - et(ij,kj + nksym*(ispin-1))
                  do ic=1,ncorb
                    if( c_beg - 1 + ic > nauxvirK(kc,ispin) ) cycle
                    dcb = 1.d0 / e2Ha / (ebij + et(noccK(kc,ispin)+c_beg-1+ic,kc+nksym*(ispin-1)))
                    !dcb = min( dcb, 13.605d0 )
                    emp2 = emp2 - fij * dcb * Vca(ic,b_beg-1+ib) * &
                      conjg( fac*Vca(ic,b_beg-1+ib) - Vcb(ic,ib) )
                  enddo
                enddo

                do ia=1,numvir
                  ia0 = noccK(ka,ispin) + nauxvirK(ka,ispin) - numvir
                  eaij = et(ia0+ia,ka + nksym*(ispin-1)) - &
                         et(ii,ki + nksym*(ispin-1)) - et(ij,kj + nksym*(ispin-1)) 
                  do ib=1,nborb
                    if( ia0 + b_beg - 1 + ib > norb ) cycle
                    ebij = et(ia0+b_beg-1+ib,ka + nksym*(ispin-1)) - &
                        et(ii,ki + nksym*(ispin-1)) - et(ij,kj + nksym*(ispin-1)) 
                    do ic=1,ncorb
                      if( c_beg - 1 + ic > nauxvirK(kc,ispin) ) cycle
                      dca = 1.d0 / e2Ha / (eaij + et(noccK(kc,ispin)+c_beg-1+ic,kc + nksym*(ispin-1)))
                      dcb = 1.d0 / e2Ha / (ebij + et(noccK(kc,ispin)+c_beg-1+ic,kc + nksym*(ispin-1)))
                      dca = min( dca, 13.605d0 )
                      dcb = min( dcb, 13.605d0 )
                      Dab(ia,ib,ka,ispin) = Dab(ia,ib,ka,ispin) + fij * dca * dcb * Vca(ic,ia) * &
                        conjg( fac*Vca(ic,b_beg-1+ib) - Vcb(ic,ib) ) / nksym / nksym 
                    enddo
                  enddo
                enddo  

                if(nspin == 2) then
                  ! opposite spin contribution
  
                endif

              enddo

            enddo  ! ka

          enddo  ! Q

        enddo  ! ii

      enddo  ! ispin
      CALL stop_clock ( 'mp2no' )
      CALL print_clock( 'mp2no' )

    enddo

    if(nproc_image > 1) CALL mp_sum ( emp2, intra_image_comm ) 
    if(ionode) write(*,*) 'EMP2 (Ha): ',emp2/(1.d0*nksym)

    IF( ALLOCATED(Kic) ) DEALLOCATE (Kic)
    IF( ALLOCATED(Kib) ) DEALLOCATE (Kib)
    IF( ALLOCATED(Psiocc) ) DEALLOCATE (Psiocc)
    IF( ALLOCATED(Psivir) ) DEALLOCATE (Psivir)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(noccK) ) DEALLOCATE(noccK)
    IF( ALLOCATED(Vca) ) DEALLOCATE(Vca)
    IF( ALLOCATED(Vcb) ) DEALLOCATE(Vcb)
    IF( ALLOCATED(Rgb) ) DEALLOCATE(Rgb)

    ! now generate NO
    if(ionode) write(*,*) 'Calculating natural orbitals'
    if(nproc_image > 1) call mp_barrier( intra_image_comm ) 
    allocate( Vca(numvir,numvir) )

    if( me_image == root_image) then
      !
      h5len = LEN_TRIM(h5file) 
      allocate( xkcart_(3,nksym) )
      do ik=1,nksym
        xkcart_(1:3,ik) = xksym(1:3,ik)*tpiba
      enddo
      recv(1:3,1:3) = bg(1:3,1:3) * tpiba
      at0(1:3,1:3) = at(1:3,1:3) * alat
#if defined(__HDF5) || defined(__HDF5_C)
      CALL esh5_posthf_open_write(h5id_output_orbs%id,h5file,h5len, error)
      if(error .ne. 0 ) &
          call errore('mp2no','error opening orbital file for write',1)
      CALL esh5_posthf_write_meta(h5id_output_orbs%id,orbsG,5,nksym,1,npwx,xkcart_, &
                          1,dfft%nr1,dfft%nr2,dfft%nr3,at0,recv,alat,error)
      if(error .ne. 0 ) &
          call errore('mp2no','error writing meta data OrbG',1)
#endif
      deallocate(xkcart_)
      !
      allocate( Orbitals(npwx,1), eig(numvir) )
      !
    endif

    do ka=1,nksym
    
      Vca(:,:) = (0.d0,0.d0)
      Vca(1:numvir,b_beg:b_end) = Dab(1:numvir,1:nborb,ka,1)  
      call mp_sum( Vca, intra_image_comm )

      if(me_image == root_image) then

        CALL davcio (evc, 2*nwordwfc, iunwfc, ka, - 1)

        call eigsys('V', 'U', .true., nauxvirK(ka,1), numvir, Vca, eig) 

        if(nvmax > 0) then
          nv=nvmax
          if(verbose) then
            do ia=1,nv
              write(*,*) ka,ia,eig(ia)
            enddo
          endif
        else
          do ia=1,nauxvirK(ka,1)
            if(verbose) write(*,*) ka,ia,eig(ia)
            if( dble(eig(ia)) < eigcut ) exit
            nv = ia
          enddo
        endif

        if( size( Orbitals, 2) < nv ) then
          deallocate(Orbitals)
          allocate( Orbitals(npwx,nv) )
        endif
        Orbitals(:,:) = (0.d0,0.d0)
        call zgemm('N','N',ngksym(ka),nv,nauxvirK(ka,1),(1.d0,0.d0),evc(1,maxocc+1),npwx,  &
                     Vca,numvir,(0.d0,0.d0),Orbitals,npwx)
        ! Write to file
        nauxvirK(ka,1) = nv
        call esh5_posthf_write(h5id_output_orbs%id,orbsG,5,ka-1,npwx,nauxvirK(ka,1),Orbitals,npwx,error)
        if(error .ne. 0 ) &
          call errore('mp2no','error writing orbital',1)
        !
        if(numspin .eq. 2) then
          call errore(' Error: lsda not yet implemented in mp2no.',1)
        endif
      endif

    enddo

    if( me_image == root_image) then
#if defined(__HDF5) || defined(__HDF5_C)
      call esh5_posthf_write_norb(h5id_output_orbs%id,orbsG,5,nksym,nauxvirK(1,1),error)
      if(error .ne. 0 ) &
        call errore('mp2no','error writing nvirK OrbsG',1)
      call esh5_posthf_close_write(h5id_output_orbs%id)
#endif
    endif
    
    if(allocated(eig)) deallocate(eig)
    if(allocated(Orbitals)) deallocate(Orbitals)
    IF( ALLOCATED(nauxvirK) ) DEALLOCATE(nauxvirK)
    IF( ALLOCATED(Vca) ) DEALLOCATE(Vca)
    IF( ALLOCATED(Dab) ) DEALLOCATE(Dab)

  END SUBROUTINE mp2no

  SUBROUTINE approx_mp2no(h5file,dfft,nskipvir,eigcut,low_memory,esh5_file)
    USE parallel_include
    USE wvfct, ONLY: wg, et
    USE klist, ONLY: wk
    USE paw_variables,        ONLY : okpaw
    USE uspp,       ONLY : okvan
    USE uspp,     ONLY : vkb, nkb
    USE wavefunctions, ONLY : evc
    USE io_global, ONLY : stdout
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE paw_variables,           ONLY : okpaw
    USE becmod,  ONLY : bec_type, becp, ALLOCATE_bec_type, DEALLOCATE_bec_type
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r, qgm
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
!    USE exx, ONLY: ecutfock
!    USE gvect,     ONLY : ecutrho
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    USE read_orbitals_from_file, ONLY: get_psi_esh5
    USE orbital_generators, ONLY: get_spanning_basis
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(IN) :: h5file 
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: nskipvir
    REAL(DP), INTENT(IN) :: eigcut
    LOGICAL, INTENT(IN) :: low_memory
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: esh5_file
    !
    CHARACTER(len=256) :: tmp
    INTEGER :: Q, iab, ni, nj, iuv, ia0, ib0, ic0, grid_type_, h5len
    INTEGER :: ia, ib, ic, ii, ka, kc, ki
    INTEGER :: ik, ibnd, nxxs, i, ikk, kk, jj, iu0, ispin, is1, is2, noa, nob, n
    INTEGER :: c_beg, c_end, b_beg, b_end, npw 
    INTEGER :: ncorb, nborb, error
    INTEGER :: maxv_rank, nv, nvb, nvmax, minvir, numvir 
    INTEGER :: nel(2), maxocc, maxauxvir
    INTEGER, ALLOCATABLE :: noccK(:,:), nauxvirK(:,:)
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: recv(3,3), at0(3,3)
    REAL(DP) :: dQ(3),dQ1(3),dk(3), dG(3, 27), scl
    REAL(DP) :: fac,pnorm, residual, maxv, fXX, ecut
    COMPLEX(DP) :: ctemp, etemp, fii
    COMPLEX(DP) :: emp2
    REAL(DP) :: dca, dcb, eaii, ebii
    ! QPOINT stuff
    COMPLEX(DP), ALLOCATABLE :: Dab(:,:,:,:)    
    COMPLEX(DP), ALLOCATABLE :: Kic(:,:)    
    COMPLEX(DP), ALLOCATABLE :: Psiocc(:,:)   
    COMPLEX(DP), ALLOCATABLE :: Psivir(:,:,:,:)  
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Rgb(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Vca(:,:)  ! 
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:)
    COMPLEX(DP), ALLOCATABLE :: AlphaOrbs(:,:)
    COMPLEX(DP), ALLOCATABLE :: BetaOrbs(:,:)
    COMPLEX(DP), ALLOCATABLE :: S(:,:) 
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:) 
    COMPLEX(DP), ALLOCATABLE :: buff(:) 
    REAL(DP), ALLOCATABLE :: eig(:)   
    REAL(DP), ALLOCATABLE :: xkcart_(:,:)
    TYPE(h5file_type) :: h5id_input_orbs, h5id_output_orbs
    !

    nxxs = dfft%nr1x*dfft%nr2x*dfft%nr3x

    if(npool > 1) &
      call errore('approx_mp2no','Error: npool > 1 not allowed in appmp2no',1)
    if(noncolin) call errore('app_mp2no','No nonconlin yet.',1)

    fac=1.d0
    if(nspin==1) fac=2.d0

    nvmax = 0
    if(eigcut < 0.0) nvmax = nint(-eigcut)

    if(nspin > 2 ) &
      call errore('approx_mp2no','Error: nspin>2 in mp2no not implemented.',1)

    ! find number of electrons
    nel(:) = 0
    allocate( noccK(nksym,nspin), nauxvirK(nksym,nspin) )
    noccK(:,:)=0
    nauxvirK(:,:)=0
    do ispin=1,numspin
      pnorm=0.d0
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        if(abs(wk(ikk))>1.d-10) then
            scl = 1.d0/wk(ikk)
        else
            scl = 1.d0
        endif
        do ia=1,norb
          pnorm = pnorm + wg(ia,ikk)*scl
          if( abs(wg(ia,ikk)*scl) > 0.01d0 ) then
            noccK(ik,ispin) = noccK(ik,ispin) + 1
          endif
        enddo
      enddo
      nel(ispin) = nint(pnorm)  
    enddo
    maxocc = maxval(noccK(1:nksym,1))
    if( maxval(noccK(1:nksym,numspin)) .gt. maxocc ) &
      call errore('approx_mp2no','Error: Maximum occupation in down electron.',1) 

    ! lowset virtual orbital to be rotated
    ! range of a/b is [minvir,norb] and range of c is [noccK+1,norb] 
    minvir = maxocc + max(1,nskipvir)
    ! check that occupation of minvir is zero for all kpoint/spin
    do ispin=1,numspin
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        if(abs(wk(ikk))>1.d-10) then
          scl = 1.d0/wk(ikk)
        else
          scl = 1.d0
        endif
        if( abs(wg(minvir,ikk)*scl) > 1.d-6 ) & 
          call errore('approx_mp2no','Error: Non zero occupation at minimum virtual state in mp2no.',1)
!        ! put some sort of warning if eigenvalue gap from maxocc to minvir is too small
!        if( abs(et(minvir,ikk) - et(noccK(ikk,ispin),ikk)) < 0.07d0  ) &   ! < 1eV
!          call errore('pw2posthf','Error: Small eigenvalue gap, increase minvir.',1)
      enddo
    enddo
    nauxvirK(:,:) = norb - noccK(:,:)
    maxauxvir = maxval(nauxvirK(1:nksym,:))
    numvir = norb - minvir + 1
    if(numvir < nvmax) &
      call errore('approx_mp2no','Error: Requesting too many states, increase norb.',1)

    write(*,*) 'Approximate MP2 NO: '
    write(*,*) 'Maximum # of auxiliary virtuals: ',maxauxvir
    write(*,*) 'Range of virtuals:',minvir,minvir+numvir
    write(*,*) 'Number of virtuals:',numvir

    ! should do this over rectangular grid maxauxvir x numvir
    call fair_divide(c_beg,c_end,me_pool+1,nproc_pool,maxauxvir)
    ncorb   = c_end-c_beg+1

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over nxxs and gather over proc grid
    allocate( phasefac(nxxs,27) )
    CALL calculate_phase_factor(dfft, phasefac, dG)

    if(ionode) write(*,*) 'Starting loop over Q vectors.'
    ! 
    ! Loop over Q points
    !
    ! final normalization of 1.0/nksym applied later to keep thresh consistent 
    ! with single k-point case 
    pnorm = 1.d0 / (omega*nxxs*nksym) 
    emp2 = (0.d0,0.d0)
    ispin = 1  ! alpha for now only

    if(verbose .and. ionode) then
      write(*,*) ' Summary - nborb, ncorb, nxxs, nks, nspin:',nborb,ncorb,nxxs,nksym,nspin
      write(*,*) ' Allocations (in MB): '
      write(*,*) '     - Dab: ',numvir*numvir*nksym*nspin*16.0/1024.0/1024.0
      write(*,*) '     - Kic: ',nxxs*ncorb*16.0/1024.0/1024.0
      if(.not.low_memory) & 
      write(*,*) '     - Psivir:',nxxs*maxauxvir*nksym*numspin*16.0/1024.0/1024.0
      write(*,*) '     - Psiocc: ',2*nxxs*16.0/1024.0/1024.0
      write(*,*) '     - Vca: ',ncorb*numvir*16.0/1024.0/1024.0
      write(*,*) '     - Rgb: ',nxxs*numvir*16.0/1024.0/1024.0
      FLUSH( stdout )
    endif
    if(nproc_image > 1) call mp_barrier( intra_image_comm )

    ! eventually distribute nxxs and reduce matrices below within groups
    allocate( Dab(numvir,numvir,nksym,nspin) )
    allocate( Kic(nxxs,ncorb) )
    if(low_memory) then
      allocate( buff(nxxs) )  ! [noccK+1,norb]
    else
      allocate( Psivir(nxxs,maxauxvir,nksym,numspin) )  ! [noccK+1,norb]
    endif
    allocate( Psiocc(nxxs,2) )
    allocate( Vca(ncorb,numvir), Rgb(nxxs,numvir) )
    allocate( vCoulG(nxxs) )

    ! need eigenvalues, right now assuming they are given by pwscf files

    if(low_memory) then
      ! write Psivir to file and read on demand
      tmp = TRIM( h5file ) // '.psivir.h5'
      h5len = LEN_TRIM(tmp)   
      allocate( xkcart_(3,2*nksym)  )
      if(ionode) then
        recv(1:3,1:3) = bg(1:3,1:3) * tpiba
        at0(1:3,1:3) = at(1:3,1:3) * alat
        CALL esh5_posthf_open_write(h5id_output_orbs%id,tmp,h5len, error)
        CALL esh5_posthf_write_meta(h5id_output_orbs%id,orbsR,5,2*nksym,1,0,xkcart_, &
                        0,dfft%nr1,dfft%nr2,dfft%nr3,at0,recv,alat,error)
        if(error .ne. 0 ) &
          call errore('approx_mp2no','error writing meta data OrbR',1)

        do ispin=1,numspin
          do ik=1,nksym
            ikk = ik + nksym*(ispin-1)

            CALL gk_sort (xksym (1:3, ik), ngm, g, ecutwfc / tpiba2, &
                   npw, igksym(1), g2kin)

            CALL davcio (evc, 2*nwordwfc, iunwfc, ikk, - 1)

            ni = norb - noccK(ik,ispin)
            do ibnd=1,ni
               !
               psic(:) = (0.d0,0.d0) 
               psic(dfft%nl(igksym(1:npw)))=evc(1:npw,noccK(ik,ispin)+ibnd)
               if(gamma_only) &
                 psic(dfft%nlm(igksym(1:npw)))=CONJG(evc(1:npw,noccK(ik,ispin)+ibnd))
               ! 
              CALL invfft ('Wave', psic(:), dfft)
              !
              call esh5_posthf_write_band(h5id_output_orbs%id, orbsR, 5, ikk-1, &
                    nxxs, ibnd-1, psic, error)
              if(error .ne. 0 ) &
                call errore('approx_mp2no','error writing orbital in real space',1)  
              !
            enddo
          enddo
        enddo
        call esh5_posthf_write_norb(h5id_output_orbs%id,orbsR,5,2*nksym,nauxvirK, error)
        if(error .ne. 0 ) &
          call errore('approx_mp2no','error writing norb OrbsR',1)
        call esh5_posthf_close_write(h5id_output_orbs%id)
      endif
      if(nproc_image > 1) call mp_barrier( intra_image_comm ) 
      ! now everybody opens read only  
      call esh5_posthf_open_read_r(h5id_output_orbs%id,tmp,h5len,error)
      deallocate( xkcart_ )
    else
      do ispin=1,numspin
        call get_orbitals_set(h5id_input_orbs,'davcio','psir',dfft,ispin,&
                    Psivir(:,:,:,ispin),1,norb,1,nksym,nminK=noccK(:,ispin))
      enddo  
    endif

    Dab(:,:,:,:) = (0.d0,0.d0)
    do ki=1,nksym

      if(ionode) write(*,*) 'K:',ki
      CALL start_clock ( 'mp2no' )
    
      do ispin=1,numspin  
        !
        CALL gk_sort (xksym (1:3, ki + nksym*(ispin-1)), ngm, g, ecutwfc / tpiba2, &
               npw, igksym(1), g2kin)
        CALL davcio (evc, 2*nwordwfc, iunwfc, ki + nksym*(ispin-1), - 1)
        !
        do ii=1,noccK(ki,ispin)
          !
          Psiocc(:,:) = (0.d0,0.d0)
          Psiocc(dfft%nl(igksym(1:npw)),1)=evc(1:npw,ii)
          if(gamma_only) &
            Psiocc(dfft%nlm(igksym(1:npw)),1)=CONJG(evc(1:npw,ii))
          ! 
          CALL invfft ('Wave', Psiocc(:,1), dfft)
          !
          Psiocc(:,1) = CONJG(Psiocc(:,1))
          !
          do Q=1,nksym

            kc = QKtoK2(Q, ki)
            ka = -1
            do ikk = 1,nksym
              if(QKtoK2(Q, ikk) == ki) then
                ka = ikk
                exit
              endif  
            enddo
            if(ka.lt.1) call errore('approx_mp2no','Can not find ka.',1)
            
            ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(ka) - k(ki)  |^{-2}
            CALL g2_convolution(dfft%ngm, g, xksym (1:3, kc), xksym (1:3,ki),psic)
            vCoulG(:) = (0.d0,0.d0)
            vCoulG( dfft%nl(1:dfft%ngm) ) = e2Ha * psic(1:dfft%ngm)

            Kic(:,:)=(0.d0,0.d0)
            do ibnd=1,ncorb
        
              if(c_beg-1+ibnd > norb-noccK(kc,ispin)) cycle
              if(low_memory) then
                CALL start_clock ( 'band_io' )
                call get_psi_esh5(h5id_input_orbs,c_beg-1+ibnd,kc,ispin,buff)  
                CALL stop_clock ( 'band_io' )
                ! Orbital pairs in R
                Kic(1:dfft%nnr,ibnd) = Psiocc(1:dfft%nnr,1) * &
                                      buff(1:dfft%nnr) * pnorm
              else
                ! Orbital pairs in R
                Kic(1:dfft%nnr,ibnd) = Psiocc(1:dfft%nnr,1) * &
                                      Psivir(1:dfft%nnr,c_beg-1+ibnd,kc,ispin) * pnorm
              endif  

              ! fwfft orbital pairs to G
              CALL start_clock ( 'fwfft' )
              CALL fwfft ('Rho', Kic(:,ibnd), dfft)
              CALL stop_clock ( 'fwfft' )

              ! multiply by FFT[ 1/r ]
              Kic(1:nxxs,ibnd) = Kic(1:nxxs,ibnd) * vCoulG(1:nxxs)

              ! invfft to R
              CALL start_clock ( 'invfft' )
              CALL invfft ('Rho', Kic(:,ibnd), dfft)
              CALL stop_clock ( 'invfft' )

            enddo

            ! (ic|ia) 
            dQ(1:3) = xksym(1:3, kc) - xksym(1:3, ki) + &
                      xksym(1:3, ka) - xksym(1:3, ki)
            kk=0
            do jj=1,27
              if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
                kk=jj
                exit
              endif
            enddo
            if(kk.lt.1) call errore('approx_mp2no','Can not find dQ in G list.',1)

            fii = wg(ii,ki + nksym*(ispin-1)) * wg(ii,ki + nksym*(ispin-1)) 
            if(abs(wk(ki + nksym*(ispin-1)))>1.d-10) &
                fii = fii/wk(ki + nksym*(ispin-1))/wk(ki + nksym*(ispin-1))

            Psiocc(1:dfft%nnr,2) = Psiocc(1:dfft%nnr,1) * phasefac(1:dfft%nnr, kk)

            CALL start_clock ( 'Rgb' )
            if(low_memory) then
              ! Vca = (ic | ia) = sum_r Kic(r) * conjg(psi_i(r)) * psi_a(r) * exp(iQ) 
              do ia=1,numvir
                ia0 = nauxvirK(ka,ispin) - numvir  ! index of ia on Psivir  
                if(ia0+ia > norb-noccK(ka,ispin)) then
                   Rgb(1:dfft%nnr,ia) = (0.d0,0.d0)
                   cycle 
                endif 
                CALL start_clock ( 'band_io' )
                call get_psi_esh5(h5id_input_orbs,ia0+ia,ka,ispin,buff)  
                CALL stop_clock ( 'band_io' )
                Rgb(1:dfft%nnr,ia) = Psiocc(1:dfft%nnr,2) * &
                                           buff(1:dfft%nnr) 
              enddo
            else
              ! Vca = (ic | ia) = sum_r Kic(r) * conjg(psi_i(r)) * psi_a(r) * exp(iQ) 
              do ia=1,numvir
                ia0 = nauxvirK(ka,ispin) - numvir  ! index of ia on Psivir  
                Rgb(1:dfft%nnr,ia) = Psiocc(1:dfft%nnr,2) * &
                                           Psivir(1:dfft%nnr,ia0+ia,ka,ispin) 
              enddo
            endif
            CALL stop_clock ( 'Rgb' )
            CALL start_clock ( 'gemm' )
            call zgemm('T','N',ncorb,numvir,nxxs,(1.d0,0.d0),Kic(1,1),nxxs,  &
                   Rgb(1,1),nxxs,(0.d0,0.d0),Vca(1,1),ncorb)
            CALL stop_clock ( 'gemm' )

            CALL start_clock ( 'Dab' )
            do ia=1,numvir
              ia0 = noccK(ka,ispin) + nauxvirK(ka,ispin) - numvir
              eaii = et(ia0+ia,ka + nksym*(ispin-1)) - &
                     2.d0*et(ii,ki + nksym*(ispin-1)) 
              do ib=1,numvir
                if( ia0+ib > norb ) cycle
                ebii = et(ia0+ib,ka + nksym*(ispin-1)) - &
                    2.d0*et(ii,ki + nksym*(ispin-1)) 
                do ic=1,ncorb
                  if( c_beg - 1 + ic > nauxvirK(kc,ispin) ) cycle
                  dca = 1.d0 / e2Ha / (eaii + et(noccK(kc,ispin)+c_beg-1+ic,kc + nksym*(ispin-1)))
                  dcb = 1.d0 / e2Ha / (ebii + et(noccK(kc,ispin)+c_beg-1+ic,kc + nksym*(ispin-1)))
                  !dca = min( dca, 13.605d0 )
                  !dcb = min( dcb, 13.605d0 )
                  Dab(ia,ib,ka,ispin) = Dab(ia,ib,ka,ispin) + fii * dca * dcb * Vca(ic,ia) * &
                    conjg( Vca(ic,ib) ) / nksym  
                enddo
              enddo
            enddo  
            CALL stop_clock ( 'Dab' )

          enddo  ! Q

        enddo  ! ii

      enddo  ! ispin
      CALL stop_clock ( 'mp2no' )
      CALL print_clock( 'mp2no' )

    enddo
    CALL print_clock( 'fwfft' )
    CALL print_clock( 'invfft' )
    if(low_memory) CALL print_clock( 'band_io' )
    CALL print_clock( 'Rgb' )
    CALL print_clock( 'gemm' )
    CALL print_clock( 'Dab' )

    if(nproc_image > 1) CALL mp_sum ( emp2, intra_image_comm ) 
    if(ionode) write(*,*) 'EMP2: ',emp2/(1.d0*nksym)

    IF( ALLOCATED(Kic) ) DEALLOCATE (Kic)
    IF( ALLOCATED(Psiocc) ) DEALLOCATE (Psiocc)
    IF( ALLOCATED(Psivir) ) DEALLOCATE (Psivir)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(noccK) ) DEALLOCATE(noccK)
    IF( ALLOCATED(Vca) ) DEALLOCATE(Vca)
    IF( ALLOCATED(Rgb) ) DEALLOCATE(Rgb)
    IF( ALLOCATED(buff) ) DEALLOCATE(buff)

    if(low_memory) call esh5_posthf_close_read(h5id_output_orbs%id) 

    ! now generate NO
    if(ionode) write(*,*) 'Calculating natural orbitals'
    if(nproc_image > 1) call mp_barrier( intra_image_comm ) 
    allocate( Vca(numvir,numvir) )

    if( me_image == root_image) then
      !
      h5len = LEN_TRIM(h5file) 
      allocate( xkcart_(3,nksym) )
      do ik=1,nksym
        xkcart_(1:3,ik) = xksym(1:3,ik)*tpiba
      enddo
      recv(1:3,1:3) = bg(1:3,1:3) * tpiba
      at0(1:3,1:3) = at(1:3,1:3) * alat
#if defined(__HDF5) || defined(__HDF5_C)
      CALL esh5_posthf_open_write(h5id_output_orbs%id,h5file,h5len, error)
      if(error .ne. 0 ) &
        call errore('approx_mp2no','error opening orbital file for write',1)
      CALL esh5_posthf_write_meta(h5id_output_orbs%id,orbsG,5,nksym,1,npwx,xkcart_, &
                          1,dfft%nr1,dfft%nr2,dfft%nr3,at0,recv,alat,error)
      if(error .ne. 0 ) &
          call errore('approx_mp2no','error writing meta data OrbG',1)
#endif
      deallocate(xkcart_)
      !
      allocate( Orbitals(npwx,1), eig(2*numvir), S(2*numvir,2*numvir) )
      allocate( AlphaOrbs(npwx,1) )
      allocate( BetaOrbs(npwx,1) )
      if(okpaw .or. okvan) then
        allocate( spsi(npwx,1) )
        CALL allocate_bec_type ( nkb, numvir, becp )
      endif
      !
    endif
    
    do ka=1,nksym

      Vca(:,:) = (0.d0,0.d0)
      Vca(1:numvir,1:numvir) = Dab(1:numvir,1:numvir,ka,1)  
      call mp_sum( Vca, intra_image_comm )

      if(me_image == root_image) then

        CALL davcio (evc, 2*nwordwfc, iunwfc, ka, - 1)

        call eigsys('V', 'U', .true., numvir, numvir, Vca, eig) 

        if(nvmax > 0) then
          nv=nvmax
          if(verbose) then
            do ia=1,nv
              write(*,*) ka,ia,eig(ia)
            enddo
          endif
        else
          do ia=1,numvir
            if(verbose) write(*,*) ka,ia,eig(ia)
            if( dble(eig(ia)) < eigcut ) exit
            nv = ia
          enddo
        endif

        if( size( AlphaOrbs, 2) < nv ) then
          deallocate(AlphaOrbs)
          allocate( AlphaOrbs(npwx,nv) )
        endif
        AlphaOrbs(:,:) = (0.d0,0.d0)
        call zgemm('N','N',ngksym(ka),nv,numvir,(1.d0,0.d0),evc(1,minvir),npwx,  &
                   Vca,numvir,(0.d0,0.d0),AlphaOrbs,npwx)

      endif

  
      if( nspin .eq. 2 ) then

        Vca(:,:) = (0.d0,0.d0)
        Vca(1:numvir,1:numvir) = Dab(1:numvir,1:numvir,ka,2)
        call mp_sum( Vca, intra_image_comm )

        if(me_image == root_image) then

          CALL davcio (evc, 2*nwordwfc, iunwfc, ka+nksym, - 1)

          call eigsys('V', 'U', .true., numvir, numvir, Vca, eig)

          ! for simplicity, forcing nv states  

          if(verbose) then
            do ia=1,nv
              write(*,*) 'Beta:',ka,ia,eig(ia)
            enddo
          endif

          if( size( BetaOrbs, 2) < nv ) then
            deallocate(BetaOrbs)
            allocate( BetaOrbs(npwx,nv) )
            if(okpaw .or. okvan) then
              call errore('approx_mp2nok: finish paw in mp2no.',1)  
              ! call allocate_becp, init_us_2  
              deallocate(spsi)
              allocate( spsi(npwx,nv) )
            endif
          endif
          BetaOrbs(:,:) = (0.d0,0.d0)
          call zgemm('N','N',ngksym(ka),nv,numvir,(1.d0,0.d0),evc(1,minvir),npwx,  &
                     Vca,numvir,(0.d0,0.d0),BetaOrbs,npwx)

          ! safe
          if( size( Orbitals, 2) < 2*nv ) then
            deallocate(Orbitals)
            allocate( Orbitals(npwx,2*nv) )
          endif

          if(okvan .or. okpaw) then
            CALL init_us_2 (ngk(ka), igk_k(1,ka), xk (1,ka), vkb)
          endif

          nvb = 0
          call get_spanning_basis(ka, nv, nvb, AlphaOrbs,npwx,BetaOrbs,npwx,Orbitals, &
                      npwx, S, 2*numvir, spsi, eig, 1.d-3)
          if(verbose) write(*,*) '# orbitals:',nv,nvb

          ! Write to file
          nauxvirK(ka,1) = min(nv,nvb)
          call esh5_posthf_write(h5id_output_orbs%id,orbsG,5,ka-1,npwx,nauxvirK(ka,1), &
                Orbitals,npwx,error)
          if(error .ne. 0 ) &
            call errore('approx_mp2no','error writing orbital',1)
          !
        endif

      else
        if(me_image == root_image) then
          ! Write to file
          nauxvirK(ka,1) = nv
          call esh5_posthf_write(h5id_output_orbs%id,h5id_output_orbs,orbsG,5,ka-1,npwx,&
                nauxvirK(ka,1),AlphaOrbs,npwx,error)
          if(error .ne. 0 ) &
            call errore('approx_mp2no','error writing orbital',1)
          !
        endif  
      endif  

    enddo

    if( me_image == root_image) then
#if defined(__HDF5) || defined(__HDF5_C)
      call esh5_posthf_write_norb(h5id_output_orbs%id,orbsG,5,nksym,nauxvirK(1,1),error)
      if(error .ne. 0 ) &
        call errore('approx_mp2no','error writing nvirK OrbsG',1)
      call esh5_posthf_close_write(h5id_output_orbs%id)
#endif
    endif
    
    if(allocated(eig)) deallocate(eig)
    if(allocated(S)) deallocate(S)
    if(allocated(Orbitals)) deallocate(Orbitals)
    if(allocated(AlphaOrbs)) deallocate(AlphaOrbs)
    if(allocated(BetaOrbs)) deallocate(BetaOrbs)
    IF( ALLOCATED(nauxvirK) ) DEALLOCATE(nauxvirK)
    IF( ALLOCATED(Vca) ) DEALLOCATE(Vca)
    IF( ALLOCATED(Dab) ) DEALLOCATE(Dab)
    if(allocated(spsi)) deallocate(spsi)
    if(okvan .or. okpaw) CALL deallocate_bec_type (becp)

  END SUBROUTINE approx_mp2no


END MODULE mp2_module
