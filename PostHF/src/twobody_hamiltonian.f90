!-------------------------------------------------------------------i
! Written by Miguel A. Morales, LLNL, 2020 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE twobody_hamiltonian 
  !----------------------------------------------------------------------
  ! 
  ! Implements MP2 and some other utilities, e.g. pseudo-canonicalization
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY: omega, alat, tpiba, tpiba2, at, bg
  USE posthf_mod, ONLY: nksym,norb,numspin,e0,efc,e2Ha,&
                        nelec_tot,nup,ndown,ke_factorization,&
                        xksym,igksym,ngksym,QKtoK2,kminus,Qpts,DM,DM_mf,nmax_DM
  USE read_orbitals_from_file, ONLY: get_psi_esh5,&
                    open_esh5_read,close_esh5_read,h5file_type
  USE wavefunctions, ONLY : evc, psic
  USE control_flags, ONLY : gamma_only
  USE gvect, ONLY: ngm, ngm_g, g, gstart, gg 
  USE io_files, ONLY: nwordwfc, iunwfc 
  USE io_global, ONLY: stdout, ionode,  ionode_id
  USE wvfct, ONLY: nbnd, npwx, g2kin
  USE klist,  ONLY : nkstot, wk
  USE becmod,   ONLY : bec_type, becp, allocate_bec_type, deallocate_bec_type
  USE mp,           ONLY: mp_sum, mp_max, mp_bcast, mp_barrier
  USE mp_images, ONLY: intra_image_comm, me_image, root_image, nproc_image
  USE mp_pools,     ONLY: inter_pool_comm, intra_pool_comm, npool, nproc_pool, &
                          me_pool,root_pool,my_pool_id
  USE mp_bands,     ONLY: nproc_bgrp
  USE noncollin_module,     ONLY : noncolin, npol
  USE lsda_mod, ONLY: lsda, nspin
  USE read_orbitals_from_file, ONLY: get_orbitals, get_orbitals_set
  use fft_interfaces,       ONLY : invfft, fwfft
  USE fft_types, ONLY: fft_type_descriptor
  USE posthf_mod, ONLY : ke_factorization
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: verbose
  ! original code, could be a bit faster but uses more memory
  !
  CONTAINS
  !  
  SUBROUTINE cholesky_r(ncmax,thresh,dfft,hamil_file,orb_file)
    USE parallel_include
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE paw_variables,           ONLY : okpaw
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist 
    USE exx, ONLY: ecutfock
    USE gvect,     ONLY : ecutrho
    USE uspp_param,         ONLY : nh 
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_file,orb_file
    REAL(DP), INTENT(IN) :: ncmax, thresh
    !
    CHARACTER(len=256) :: h5name
    INTEGER :: nchol_max, nchol, h5len, oldh5
    INTEGER :: Q, ka, kb, iab, ia, ib, iu, iv, ku, kv, ni, nj, iuv
    INTEGER :: i, ii,jj,kk, ia0, iu0, ispin, is1, is2, noa, nob
    INTEGER :: a_beg, a_end, b_beg, b_end, ab_beg, ab_end, nabpair   ! assigned orbitals/pairs
    INTEGER :: naorb, nborb, iabmax,nkloc,k_beg,k_end 
    INTEGER :: maxv_rank,nke_dim,error,ibnd,ik,j, maxnorb
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: dQ(3),dQ1(3),dk(3), dG(3, 27)
    REAL(DP) :: residual, maxv, fXX, fac, max_scl, rtemp, sqDuv
    REAL(DP), ALLOCATABLE :: maxres(:)
    REAL(DP), ALLOCATABLE :: comm(:)
    COMPLEX(DP) :: ctemp, ctemp2, etemp, E1, max_scl_intg, err_(3)
    LOGICAL :: more
    ! QPOINT stuff
    COMPLEX(DP), ALLOCATABLE :: Chol(:,:,:)    ! Cholesky Matrix 
    COMPLEX(DP), ALLOCATABLE :: GChol(:,:,:)   ! Collect Chol for energy eval  
    COMPLEX(DP), ALLOCATABLE :: Psia(:,:,:)    ! orbitals in real space
    COMPLEX(DP), ALLOCATABLE :: Psib(:,:,:)    ! orbitals in real space
    COMPLEX(DP), ALLOCATABLE :: Diag(:,:)      ! Diagonal of Residual Matrix 
    COMPLEX(DP), ALLOCATABLE :: Kqab(:,:,:)    ! Exchange potential in real space
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)
    COMPLEX(DP), ALLOCATABLE :: Vuv(:)         ! 
    COMPLEX(DP), ALLOCATABLE :: Vuv2(:)        ! 
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: v1(:,:)  ! 
    INTEGER, ALLOCATABLE :: ncholQ(:)
    TYPE(ke_factorization), ALLOCATABLE :: ke(:) ! for paw one-center terms
    COMPLEX(DP),ALLOCATABLE :: paw_one_center(:,:,:)
    TYPE(bec_type),ALLOCATABLE :: becpsia(:)
    TYPE(bec_type),ALLOCATABLE :: becpsib(:)
    COMPLEX(DP)::eX,eJ,eX_mf,eJ_mf,eX2,eJ2
    COMPLEX(DP) :: PAW_J_energy
    !
    CHARACTER(len=4) :: str_me_image
    TYPE(h5file_type) :: h5id_orbs, h5id_hamil 
    !

    ! open hamiltonian file
    if(ionode) then
      oldh5=1
      h5name = TRIM(hamil_file) 
    else
      oldh5=0
      WRITE ( str_me_image, '(I4)') me_image
      h5name = TRIM(hamil_file) //"_part"//adjustl(trim(str_me_image))
    endif
    h5len = LEN_TRIM(h5name)
    CALL esh5_posthf_open_file(h5id_hamil%id,h5name,h5len,oldh5)
    if(oldh5 .ne. 0 ) &
      call errore('cholesky','error opening hamil file',1)

    ! open orbital file 
    h5name = TRIM( orb_file ) 
    call open_esh5_read(h5id_orbs,h5name)
    if( h5id_orbs%grid_type .ne. 1 ) &
      call errore('cholesky','grid_type ne 1',1)
    if( h5id_orbs%nspin .ne. 1 ) &
      call errore('cholesky','nspin ne 1',1)
    maxnorb = h5id_orbs%maxnorb

    nmax_DM = h5id_orbs%nmax_DM
    if( nmax_DM <= 0 ) &
      call errore('cholesky','error nmax_DM <= 0',1)
    if( allocated(DM) ) then
      if( (size(DM,1).ne.nmax_DM) .or. &  
          (size(DM,2).ne.nmax_DM) .or. &
          (size(DM,3).ne.nksym) .or. &
          (size(DM,4).ne.nspin) ) then  
        deallocate(DM)
        allocate( DM(nmax_DM, nmax_DM, nksym, nspin) )
      endif 
    else
      allocate( DM(nmax_DM, nmax_DM, nksym, nspin) )
    endif
    if( allocated(DM_mf) ) then
      if( (size(DM_mf,1).ne.nmax_DM) .or. &  
          (size(DM_mf,2).ne.nmax_DM) .or. &
          (size(DM_mf,3).ne.nksym) .or. &
          (size(DM_mf,4).ne.nspin) ) then  
        deallocate(DM_mf)
        allocate( DM_mf(nmax_DM, nmax_DM, nksym, nspin) )
      endif  
    else
      allocate( DM_mf(nmax_DM, nmax_DM, nksym, nspin) )
    endif
    do i=1,nspin
      do ik=1,nksym
        call esh5_posthf_read_dm(h5id_orbs%id,"DM",2,ik-1,i-1,DM(1,1,ik,i),error)
        if(error .ne. 0 ) &
          call errore('cholesky','error reading DM',1)
        call esh5_posthf_read_dm(h5id_orbs%id,"DM_mf",5,ik-1,i-1,DM_mf(1,1,ik,i),error)
        if(error .ne. 0 ) &
          call errore('cholesky','error reading DM',1)
      enddo
    enddo

    if(nspin == 1 ) then
      fac = 2.d0 
    else
      fac = 1.d0 
    endif

    IF(tqr) THEN
!       IF(ecutfock==ecutrho) THEN
!          WRITE(stdout,'(5x,"Real-space augmentation: EXX grid -> DENSE grid")')
          tabxx => tabp
!       ELSE
!          WRITE(stdout,'(5x,"Real-space augmentation: initializing EXX grid")')
!          CALL qpointlist(dfft, tabxx)
!       ENDIF
    ENDIF

    ! parallelization has restrictions
    if(nproc_pool > 1 .and. npool .ne. nksym) &
      call errore('pw2posthf','Error: nb > 1 only allowed if npools == nkstot',1) 

    !   Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    !
    ! Data Distribution setup
    !   Partition orbital pairs:  Simple for now, find more memory friendly version
    call fair_divide(ab_beg,ab_end,me_pool+1,nproc_pool,maxnorb*maxnorb)
    nabpair = ab_end - ab_beg + 1
    !   Determine assigned orbitals
! fix this!!!!
!    a_beg   = (ab_beg-1)/maxnorb + 1
!    a_end   = (ab_end-1)/maxnorb + 1
    a_beg   = 1 
    a_end   = maxnorb 
    naorb   = a_end-a_beg+1
    b_beg   = 1
    b_end   = maxnorb   ! keeping all for simplicity now 
    nborb   = b_end-b_beg+1

    ! maximum number of cholesky vectors allowed, 
    ! must be sufficient to reach thresh
    nchol_max = int(ncmax*maxnorb)

    !  Allocate space for assigned orbitals in real space and corresponding
    !  orbital pairs
    allocate( Kqab(dfft%nnr,nabpair,nkloc), Chol(nchol_max,nabpair,nkloc) )
    allocate( Psia(dfft%nnr,naorb,nkloc), Psib(dfft%nnr,nborb,nkloc) )
    allocate( vCoulG(dfft%nnr) )
    allocate( Diag(nabpair,nkloc) )
    allocate( phasefac(dfft%nnr,27) )
    allocate( maxres(nchol_max), comm(nproc_image*4) )
    allocate( ncholQ(nksym),  xkcart(3,nksym) )
    allocate( GChol(nchol_max,nmax_DM,nmax_DM) )

    ! fix with symmetry 
    do ik=1,nksym
      xkcart(1:3,ik) = xksym(1:3,ik)*tpiba
    enddo

    if(ionode) write(*,*) 'Generating orbitals in real space'

    nke_dim = 0
    ! store bec for Psia 
    if(okvan.or.okpaw) then 
      ALLOCATE(becpsia(nkloc))
      ALLOCATE(becpsib(nkloc))
      do ik=1,nkloc
        CALL allocate_bec_type( nkb, naorb, becpsia(ik))
        CALL allocate_bec_type( nkb, nborb, becpsib(ik))
      enddo
      if(okpaw) then
        allocate(ke(ntyp)) 
        call calculate_factorized_paw_one_center(ke,nke_dim)
        !allocate( paw_one_center(nke_dim,nabpair,2*nkloc) )
        allocate( paw_one_center(nke_dim,nabpair,nkloc) )
      endif  
      call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,&
                    1,Psia,a_beg,naorb,k_beg,nkloc,becpsi=becpsia)
    else
      call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,&
                    1,Psia,a_beg,naorb,k_beg,nkloc)
    endif

    ! allocate Vuv/Vuv2 now
    allocate( Vuv(dfft%nnr+nchol_max+nke_dim+1), Vuv2(dfft%nnr) )

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over dfft%nnr and gather over proc grid
    CALL calculate_phase_factor(dfft, phasefac, dG)

    CALL start_clock ( 'orb_cholesky' )
    if(ionode) write(*,*) 'Starting loop over Q vectors.'
    !
    ! 
    ! Loop over Q points
    !
    ! final normalization of 1.0/nksym applied later to keep thresh consistent 
    ! with single k-point case 
    eX = (0.d0,0.d0)
    eJ = (0.d0,0.d0)
    eX_mf = (0.d0,0.d0)
    eJ_mf = (0.d0,0.d0)
    do Q=1,nksym

      if( Q .gt. kminus(Q) ) cycle 
      if( Q .eq. kminus(Q) ) then
        fXX=1.d0
      else
        fXX=2.d0
      endif
      max_scl_intg = (0.d0,0.d0) 
      max_scl = 0.d0

      if(ionode) write(*,*) ' Q: ',Q
      if(ionode) write(*,*) '  --  Generating orbitals in real space'
      ! 
      ! 0. Generate 'b' orbitals in real space
      !
      if(okvan.or.okpaw) then 
        call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,&
                      1,Psib,b_beg,nborb,k_beg,nkloc,Q,QKtoK2,becpsib)
      else
        call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,&
                      1,Psib,b_beg,nborb,k_beg,nkloc,Q,QKtoK2)
      endif  

      ! 
      ! Generate exchange matrix Kab(r) 
      !
      if(ionode) write(*,*) '  --  Generating exchange matrix Kab(r)'

      ! 1. Construct orbital pairs in R 
      ! 2. fwfft to G 
      ! 3. Multiply orbital pairs by FFT[ 1/r ]  
      ! 4. invfft back to R
      Kqab(:,:,:)=(0.d0,0.d0)
      do ik=1,nkloc

        ka = k_beg + ik - 1  
        kb = QKtoK2(Q, ka)
        iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)

        ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
        CALL g2_convolution(dfft%ngm, g, xksym (1:3, kb), xksym (1:3, ka), psic)  
        vCoulG(:) = (0.d0,0.d0)
        vCoulG( dfft%nl(1:dfft%ngm) ) = e2Ha * psic(1:dfft%ngm) 
        if(gamma_only) vCoulG( dfft%nlm(1:dfft%ngm) ) = e2Ha * conjg(psic(1:dfft%ngm))

        IF ( okvan .and..not.tqr ) &
          CALL qvan_init (dfft%ngm, xksym (1:3, ka), xksym (1:3, kb))

        do ibnd=1,nabpair

          ! iab = (ia-1)*norbK(kb) + ib 
          iab = ab_beg + ibnd - 1 
          if(iab > iabmax) exit 
          ia = (iab-1)/h5id_orbs%norbK(kb) + 1 
          ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1 

          do j = 1,dfft%nnr
            Kqab(j,ibnd,ik) = CONJG(Psia(j,ia-a_beg+1,ik)) * Psib(j,ib-b_beg+1,ik) / omega
          enddo
          if(okvan .and. tqr) then  
            ! Orbital pairs in R
            call addusxx_r(Kqab(:,ibnd,ik),becpsia(ik)%k(:,ia-a_beg+1),&
                                         becpsib(ik)%k(:,ib-b_beg+1)) 
          endif  

          ! fwfft orbital pairs to G
          CALL fwfft ('Rho', Kqab(:,ibnd,ik), dfft)

          if(okvan .and. .not.tqr) &
            CALL addusxx_g(dfft, Kqab(:,ibnd,ik), xksym(1:3,ka), &
                xksym(1:3,kb),'c',becphi_c=becpsia(ik)%k(:,ia-a_beg+1),&
                becpsi_c=becpsib(ik)%k(:,ib-b_beg+1))

          ! multiply by FFT[ 1/r ]
          do j = 1,dfft%nnr   
            Kqab(j,ibnd,ik) = Kqab(j,ibnd,ik) * vCoulG(j) 
          enddo

          ! invfft kab to R
          CALL invfft ('Rho', Kqab(:,ibnd,ik), dfft)

          ! add normalization here, why not!
          Kqab(1:dfft%nnr,ibnd,ik) = Kqab(1:dfft%nnr,ibnd,ik)/  & 
                (dfft%nr1*dfft%nr2*dfft%nr3*nksym) 

          if(okpaw) then
            call contract_paw_one_center(ke,paw_one_center(:,ibnd,ik), &
                becpsia(ik)%k(:,ia-a_beg+1),becpsib(ik)%k(:,ib-b_beg+1))
            paw_one_center(:,ibnd,ik) = paw_one_center(:,ibnd,ik) / (1.d0*nksym)
          endif
        
        enddo

        IF ( okvan .and..not.tqr ) CALL qvan_clean ()

      end do

    if(.false.) then
    do ik=1,nkloc

      ka = k_beg + ik - 1
      kb = QKtoK2(Q, ka)

      IF ( okvan .and..not.tqr ) &
        CALL qvan_init (dfft%ngm, xksym (1:3, kb), xksym (1:3, ka))

      do ia=1,nup
        do ib=1,nup

          ibnd = (ia-1) * h5id_orbs%norbK(kb) + ib

          Vuv(:) = (0.d0,0.d0)
          if( tqr ) then  
            call addusxx_r(Vuv(1:dfft%nnr),becpsib(ik)%k(:,ib),&
                                  becpsia(ik)%k(:,ia))
          else
            ! augmentation charge in 'G'
            CALL addusxx_g(dfft, Vuv(1:dfft%nnr), xksym(1:3,kb), &
              xksym(1:3,ka),'c',becphi_c=becpsib(ik)%k(:,ib-b_beg+1),&
              becpsi_c=becpsia(ik)%k(:,ia-a_beg+1))
            ! fft to real space 
            CALL invfft ('Rho', Vuv(1:dfft%nnr), dfft)
          endif
          do j = 1,dfft%nnr
            Vuv(j) = Vuv(j) * omega + CONJG(Psib(j,ib,ik)) * Psia(j,ia,ik) 
          enddo
          do j=1,dfft%nnr
            eX = eX - fXX * fac * Kqab(j, ibnd, ik) * Vuv(j) 
          enddo
          ctemp = PAW_J_energy(becpsia(ik)%k(:,ia),becpsib(ik)%k(:,ib), &
                             becpsib(ik)%k(:,ib),becpsia(ik)%k(:,ia)) / (1.d0*nksym) 
          eX = eX - fXX * fac * ctemp
          eX2 = eX2 - fXX * fac * ctemp

        enddo
      enddo

      IF ( okvan .and..not.tqr ) CALL qvan_clean ()

    enddo

    if(Q==1) then
    do ka=1,nkloc
    do ku=1,nkloc

      IF ( okvan .and..not.tqr ) &
        CALL qvan_init (dfft%ngm, xksym (1:3, ku), xksym (1:3, ku))

      do ia=1,nup
      do iu=1,nup

          iab = (ia-1) * h5id_orbs%norbK(ka) + ia

          Vuv(:) = (0.d0,0.d0)  
          if(tqr) then
            call addusxx_r(Vuv(1:dfft%nnr),becpsia(ku)%k(:,iu),&
                                  becpsia(ku)%k(:,iu))
          else
            CALL addusxx_g(dfft, Vuv(1:dfft%nnr), xksym(1:3,ku), &
              xksym(1:3,ku),'c',becphi_c=becpsia(ku)%k(:,iu),&
              becpsi_c=becpsia(ku)%k(:,iu))
            ! fft to real space 
            CALL invfft ('Rho', Vuv(1:dfft%nnr), dfft)
          endif
          do j = 1,dfft%nnr
            Vuv(j) = Vuv(j) * omega + CONJG(Psia(j,iu,ku)) * Psia(j,iu,ku) 
          enddo
          do j=1,dfft%nnr
            eJ = eJ + fXX * fac * fac * Kqab(j, iab, ka) * Vuv(j) 
          enddo
          ctemp = PAW_J_energy(becpsia(ka)%k(:,ia),becpsia(ka)%k(:,ia), &
                             becpsia(ku)%k(:,iu),becpsia(ku)%k(:,iu)) / (1.d0*nksym)
          eJ = eJ + fXX * fac * fac * ctemp  
          eJ2 = eJ2 + fXX * fac * fac * ctemp  

      enddo
      enddo

      IF ( okvan .and..not.tqr ) CALL qvan_clean ()

    enddo
    enddo
    endif
    goto 401
    endif

      if(ionode) write(*,*) '  --  Constructing Diagonal of Residual Matrix '
      Diag(:,:)=(0.d0,0.d0)
      iu=-1 
      iv=-1 
      ku=-1 
      maxv=0.d0
      do ik=1,nkloc

        ka = k_beg + ik - 1  
        kb = QKtoK2(Q, ka)
        iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)

        IF ( okvan .and..not.tqr ) &
          CALL qvan_init (dfft%ngm, xksym (1:3, kb), xksym (1:3, ka))

        do ibnd=1,nabpair

          iab = ab_beg + ibnd - 1
          if(iab > iabmax) exit 
          ia = (iab-1)/h5id_orbs%norbK(kb) + 1
          ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1

          ! (ab|cd) = kab(r) * conjg(Psi(b,r)) * Psi(a,r)
          ! Diag(a,b) = (ab|ba)  
          if(okvan) then
            Vuv(1:dfft%nnr) = (0.d0,0.d0)
            ! add charge augmentation  
            if( tqr ) then
              call addusxx_r(Vuv(1:dfft%nnr),becpsib(ik)%k(:,ib-b_beg+1),&
                                  becpsia(ik)%k(:,ia-a_beg+1))
            else
              ! augmentation charge in 'G'
              CALL addusxx_g(dfft, Vuv(1:dfft%nnr), xksym(1:3,kb), &
                xksym(1:3,ka),'c',becphi_c=becpsib(ik)%k(:,ib-b_beg+1),&
                becpsi_c=becpsia(ik)%k(:,ia-a_beg+1))                
              ! fft to real space 
              CALL invfft ('Rho', Vuv(1:dfft%nnr), dfft)
            endif
            do j = 1,dfft%nnr
              Diag(ibnd,ik) = Diag(ibnd,ik) + Kqab(j, ibnd, ik) * & 
                ( omega * Vuv(j) + CONJG(Psib(j,ib-b_beg+1,ik)) * Psia(j,ia-a_beg+1,ik)) 
            enddo        
            if(okpaw) then
              call contract_paw_one_center(ke,Vuv(1:nke_dim), &
                becpsib(ik)%k(:,ib-b_beg+1),becpsia(ik)%k(:,ia-a_beg+1)) 
              do j=1,nke_dim
                Diag(ibnd,ik) = Diag(ibnd,ik) + paw_one_center(j,ibnd,ik) * Vuv(j)
              enddo
            endif
          else
            do j = 1,dfft%nnr
              Diag(ibnd,ik) = Diag(ibnd,ik) + Kqab(j, ibnd, ik) * & 
                            CONJG(Psib(j,ib-b_beg+1,ik)) * Psia(j,ia-a_beg+1,ik) 
            enddo        
          endif
          if( Aimag(Diag(ibnd,ik)) .gt. 1.d-9 ) &
            write(*,*) 'WARNING: Complex diagonal term: ',Q,ik, &
                                ia,ib,Diag(ibnd,ik)
          if( Dble(Diag(ibnd,ik)) .lt. -1.d-12 ) &
            write(*,*) 'WARNING: Negative diagonal term: ',Q,ik, &
                                ia,ib,Diag(ibnd,ik)
          ! enforce positive, real diagonal matrix  
          Diag(ibnd,ik) = CMPLX( abs(Dble(Diag(ibnd,ik))), 0.d0, kind=DP)   
          if( abs(Diag(ibnd,ik)) .gt. maxv ) then
            maxv = abs(Diag(ibnd,ik))
            iu = ia
            iv = ib
            ku = ka
          endif

        enddo

        IF ( okvan .and..not.tqr ) CALL qvan_clean ()

      enddo

      maxv_rank=0
      if(nproc_image > 1) then
        comm( 4*me_image+1 ) = maxv
        comm( 4*me_image+2 ) = iu
        comm( 4*me_image+3 ) = iv
        comm( 4*me_image+4 ) = ku
#if defined (__MPI)
        call MPI_ALLGATHER( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                            comm, 4, MPI_DOUBLE_PRECISION, intra_image_comm, error) 
        if(error .ne. 0) &
          call errore('pw2posthf','Error: mpi_allgather',error)
#else
        call errore('pw2posthf','Error: undefined __MPI',error)
#endif
        maxv=0.d0
        do i=0,nproc_image-1
          if( comm( 4*i+1 ) > maxv ) then
            maxv_rank=i
            maxv = comm( 4*i+1 )
          endif  
        enddo
        iu = comm( 4*maxv_rank+2 )
        iv = comm( 4*maxv_rank+3 )
        ku = comm( 4*maxv_rank+4 )
        if( maxv_rank < 0 .or. maxv_rank >= nproc_image ) then
          write(*,*) 'maxv_rank out of bounds: ',maxv_rank
          call errore('pw2posthf',' maxv_rank out of bounds ',1)
        endif
        if(maxv_rank == me_image) then
          if( ku < k_beg .or. ku > k_end) then
            write(*,*) 'ku out of bounds: ',ku,k_beg,k_end,me_image
            call errore('pw2posthf',' ku out of bounds ',1)
          endif
          iuv = (iu-1)*h5id_orbs%norbK( QKtoK2(Q, ku) ) + iv   
          if( iuv < ab_beg .or. iuv > ab_end) then
            write(*,*) 'iuv out of bounds: ',iuv,ab_beg,ab_end,me_image
            call errore('pw2posthf',' iuv out of bounds ',1)
          endif
        endif
      endif

      if(ionode) write(*,*) '  --  Constructing Cholesky Vectors' 
      Chol(:,:,:)=(0.d0,0.d0)
      more=.true.
      nchol = 0
      do while(more) 

        nchol = nchol + 1
        if( nchol .gt. nchol_max ) & 
          call errore('pw2posthf.f90','too many cholesky vectors, increase nchol_max',1)
        if(verbose .and. ionode) write(*,*) 'it: ',nchol,maxv

        CALL start_clock ( 'orb_comm_ovr' )
        if(maxv_rank == me_image) then
          ! construct Vuv(r) = conjg(Psi(r,iv,kv)) * Psi(r,iu,ku)
          ik = ku-k_beg+1
          ibnd = (iu-1)*h5id_orbs%norbK( QKtoK2(Q, ku) ) + iv - ab_beg + 1
          ! add charge augmentation  
          Vuv(1:dfft%nnr) = (0.d0,0.d0)
          if(okvan) then
            if(tqr) then
              call addusxx_r(Vuv(1:dfft%nnr),becpsib(ik)%k(:,iv-b_beg+1),&
                                  becpsia(ik)%k(:,iu-a_beg+1))
            else
              ! can write your own versions of qvan that don't
              ! allocate/deallocate
              CALL qvan_init (dfft%ngm, xksym (1:3, QKtoK2(Q,ku)), xksym (1:3, ku))  
              ! augmentation charge in 'G'
              CALL addusxx_g(dfft, Vuv(1:dfft%nnr), xksym(1:3,QKtoK2(Q,ku)), &
                xksym(1:3,ku),'c',becphi_c=becpsib(ik)%k(:,iv-b_beg+1),&
                becpsi_c=becpsia(ik)%k(:,iu-a_beg+1)) 
              ! augmentation charge in real space  
              CALL invfft ('Rho', Vuv(1:dfft%nnr), dfft)
              ! add orbital contribution
              call qvan_clean()
            endif
            if(okpaw) then 
              call contract_paw_one_center(ke,Vuv(dfft%nnr+nchol:dfft%nnr+nchol-1+nke_dim), &
                becpsib(ik)%k(:,iv-b_beg+1),becpsia(ik)%k(:,iu-a_beg+1))
              Vuv(dfft%nnr+nchol-1+nke_dim+1) = sqrt(abs(dble(Diag(ibnd,ik)))) 
            endif
          endif  
          do j = 1,dfft%nnr
            Vuv(j) = Vuv(j) * omega + CONJG(Psib(j,iv-b_beg+1,ik)) * &
                       Psia(j,iu-a_beg+1,ik)
          enddo  
          do j = 1,nchol-1
            Vuv(dfft%nnr+j) = CONJG(Chol(j,ibnd,ik)) 
          enddo  
        endif
        if(nproc_image > 1) call mp_bcast(Vuv(1:dfft%nnr+nchol-1+nke_dim+1),&
                                    maxv_rank,intra_image_comm) 
        if(okpaw) sqDuv = dble(Vuv(dfft%nnr+nchol-1+nke_dim+1))
        CALL stop_clock ( 'orb_comm_ovr' )

        CALL start_clock ( 'orb_2el' )
        ! contruct new cholesky vector L(n,iab,ik)
        ! L(ab,n) = (ab|vu) - sum_{1,n-1} L(ab,p) * conjg(L(uv,p))
        do ik=1,nkloc

          ka = k_beg + ik - 1  
          kb = QKtoK2(Q, ka)
          iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)

          ! dQ = Qba + Qdc = kb - ka + kd - kc = -Qpts(Q) + ku - kv
          dQ(1:3) = xksym(1:3, kb) - xksym(1:3, ka) + &
                  xksym(1:3, ku) - xksym(1:3, QKtoK2(Q,ku))
          ii=0
          do jj=1,27
            if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
              ii=jj
              exit 
            endif
          enddo
          if(ii.lt.1) call errore('pw2posthf','Can not find dQ in G list.',1)
          Vuv2(1:dfft%nnr) = Vuv(1:dfft%nnr) * phasefac(1:dfft%nnr, ii)

          ! possibly wasting effort here!
          call zgemv('T',dfft%nnr,nabpair,(1.d0,0.d0),Kqab(1,1,ik),dfft%nnr, &
                       Vuv2(1),1,(0.d0,0.d0),Chol(nchol,1,ik),nchol_max) 

          if(okpaw) then
            ! add onsite term: Chol(nchol,nabpair,ik) += (ab|uv)_onsite
            call zgemv('T',nke_dim,nabpair,(1.d0,0.d0),  &
                   paw_one_center(1,1,ik),nke_dim,   &
                   Vuv(dfft%nnr+nchol:dfft%nnr+nchol-1+nke_dim),1,  &
                   (1.d0,0.d0),Chol(nchol,1,ik),nchol_max) 
          endif

        enddo
        CALL stop_clock ( 'orb_2el' )

        CALL start_clock ( 'orb_cholvgen' )
        if(nchol > 1) &
          call zgemv('T',nchol-1,nabpair*nkloc,(-1.d0,0.d0),Chol(1,1,1),nchol_max, &
                       Vuv(dfft%nnr+1),1,(1.d0,0.d0),Chol(nchol,1,1),nchol_max)

        ! enforce positive definite in residual matrix
        if(okpaw) then
          do ik=1,nkloc
            ka = k_beg + ik - 1
            kb = QKtoK2(Q, ka)
            iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)
            do ibnd=1,nabpair
              iab = ab_beg + ibnd - 1
              if(iab > iabmax) exit
              rtemp = sqrt(abs(dble(Diag(ibnd,ik))))*sqDuv/abs(Chol(nchol,ibnd,ik))
              if( rtemp < 1.d0 ) then 
                if( rtemp > max_scl ) max_scl = rtemp
                if( (1.d0-rtemp)*abs(Chol(nchol,ibnd,ik)) > abs(max_scl_intg) ) &
                  max_scl_intg = (1.d0-rtemp)*Chol(nchol,ibnd,ik)    
                Chol(nchol,ibnd,ik) = (1.d0 - 1.e-8)*rtemp*Chol(nchol,ibnd,ik)
              endif
            enddo
          enddo
        endif
        ctemp = 1.d0/sqrt(maxv)
        call zscal(nabpair*nkloc,ctemp,Chol(nchol,1,1),nchol_max)
        CALL stop_clock ( 'orb_cholvgen' )

        ! update diagonal: Diag(ibnd,ik) -= L(nchol,ibnd,ik) * conjg(L(nchol,ibnd,ik))
        CALL start_clock ( 'orb_diagupd' )
        maxres(nchol) = maxv
        iu=-1
        iv=-1
        ku=-1
        maxv=0.d0
        do ik=1,nkloc

          ka = k_beg + ik - 1
          kb = QKtoK2(Q, ka)
          iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)

          do ibnd=1,nabpair

            iab = ab_beg + ibnd - 1
            if(iab > iabmax) exit

            Diag(ibnd,ik) = Diag(ibnd,ik) - Chol(nchol,ibnd,ik)*CONJG(Chol(nchol,ibnd,ik))
            if( Dble(Diag(ibnd,ik)) .lt. -thresh*0.01 ) then 
              iab = ab_beg + ibnd - 1
              ia = (iab-1)/h5id_orbs%norbK(kb) + 1
              ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1
              write(*,*) 'WARNING: Negative diagonal term: ',Q,ik,ibnd,ia,ib,Diag(ibnd,ik)
            endif
            if( abs(Diag(ibnd,ik)) .gt. maxv ) then
              maxv = abs(Diag(ibnd,ik))
              iab = ab_beg + ibnd - 1
              iu = (iab-1)/h5id_orbs%norbK(kb) + 1
              iv = MOD((iab-1), h5id_orbs%norbK(kb)) + 1
              ku = k_beg + ik - 1  
            endif

          enddo

        enddo
        CALL stop_clock ( 'orb_diagupd' )

        CALL start_clock ( 'orb_comm_ovr' )
        if(nproc_image > 1) then
          comm( 4*me_image+1 ) = maxv
          comm( 4*me_image+2 ) = iu
          comm( 4*me_image+3 ) = iv
          comm( 4*me_image+4 ) = ku
#if defined (__MPI)
          call MPI_ALLGATHER( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                            comm, 4, MPI_DOUBLE_PRECISION, intra_image_comm, error)
          if(error .ne. 0) &
            call errore('pw2posthf','Error: mpi_allgather',error)
#endif
          maxv=0.d0
          maxv_rank=0
          do i=0,nproc_image-1
            if( comm( 4*i+1 ) > maxv ) then
              maxv_rank=i
              maxv = comm( 4*i+1 )
            endif
          enddo
          iu = comm( 4*maxv_rank+2 )
          iv = comm( 4*maxv_rank+3 )
          ku = comm( 4*maxv_rank+4 )
          if( maxv_rank < 0 .or. maxv_rank >= nproc_image ) then
            write(*,*) 'maxv_rank out of bounds: ',maxv_rank
            call errore('pw2posthf',' maxv_rank out of bounds ',1)
          endif
          if(maxv_rank == me_image) then
            if( ku < k_beg .or. ku > k_end) then
              write(*,*) 'ku out of bounds: ',ku,k_beg,k_end,me_image
              call errore('pw2posthf',' ku out of bounds ',1)
            endif
            iuv = (iu-1)*h5id_orbs%norbK( QKtoK2(Q, ku) ) + iv
            if( iuv < ab_beg .or. iuv > ab_end) then
              write(*,*) 'iuv out of bounds: ',iuv,ab_beg,ab_end,me_image
              call errore('pw2posthf',' iuv out of bounds ',1)
            endif
          endif
        endif
        CALL stop_clock ( 'orb_comm_ovr' )

        if( maxv .gt. maxres(nchol) ) then
          call errore('pw2posthf','Cholesky residual increase.',1)
        endif

        if(maxv .lt. thresh) then
          more=.false.  
        endif  

      end do

      ncholQ(Q) = nchol
      CALL start_clock ( 'orb_cholwrt' )
      if(root_image == me_image) then
        call esh5_posthf_cholesky_root(h5id_hamil%id,Q,k_beg-1,nkloc,ab_beg-1,nabpair,nchol,  &
                  nksym,maxnorb*maxnorb,nchol_max,Chol,error)
      else
        call esh5_posthf_cholesky(h5id_hamil%id,Q,k_beg-1,nkloc,ab_beg-1,nabpair,nchol,  &
                  nksym,maxnorb*maxnorb,nchol_max,Chol,error)
      endif
      CALL stop_clock ( 'orb_cholwrt' )
      if(error.ne.0) &
        call errore('pw2posthf','Error in esh5_posthf_cholesky',error)
      if(ionode) write(*, *) ' Done with Q:',Q,'  nchol:',nchol
      if(ionode .and. okpaw) &
        write(*,*) '    -- Largest PAW approx: ',max_scl,max_scl_intg

      if( Q==1 ) then
        allocate( v1(nchol_max,2) )
        v1(:,:) = (0.d0,0.d0)
      endif

      ! very inefficient but simple for now
      ! terrible for large supercells!!!!
      ! eX = sum_ka sum_a,b,u,v (ab|vu) G[ka](a,u) G[kb](v,b)
      ! eJ = sum_s1_s2 sum_nc v(nc, s1) * conjg(v(nc, s2))
      ! v(nc, is) = sum_ka sum_ab Chol(n,ab,ka) * DM(a,b,ka,is) 
      do ik=1,nkloc

        ka = k_beg + ik - 1
        kb = QKtoK2(Q, ka)
        iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)
        noa = min(nmax_DM,h5id_orbs%norbK(ka))
        nob = min(nmax_DM,h5id_orbs%norbK(kb))

        GChol(:,:,:) = (0.d0,0.d0)
        do ibnd=1,nabpair

          iab = ab_beg + ibnd - 1
          if(iab > iabmax) exit 
          ia = (iab-1)/h5id_orbs%norbK(kb) + 1
          ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1
          if( (ia <= noa) .and. (ib <= nob) ) GChol(:,ia,ib) = Chol(:,ibnd,ik)

        enddo
        if( nproc_pool > 1 ) call mp_sum(GChol,intra_pool_comm)

        ! round robin for simplicity
        iab = -1 
        do ia=1,noa
        do ib=1,nob

          iab = iab + 1
          if(MOD(iab,nproc_pool) .ne. me_pool ) cycle

          ! EJ  
          if( Q == 1 ) then
            do ispin=1,min(2,nspin)
              !  
              etemp = DM(ia,ib,ka,ispin)
              iuv = (ia-1) * h5id_orbs%norbK(kb) + ib  
              v1(1:nchol,1) = v1(1:nchol,1) + GChol(1:nchol,ia,ib)*etemp
              etemp = DM_mf(ia,ib,ka,ispin)
              v1(1:nchol,2) = v1(1:nchol,2) + GChol(1:nchol,ia,ib)*etemp
              !
            enddo
          endif  

          ! EX  
          ! can call gemm to calculate (ab|uv) for all uv
          do iu=1,noa
          do iv=1,nob

            if(noncolin) then
              ctemp = DM(ia,iu,ka,1) * DM(iv,ib,kb,1) + &
                      DM(ia,iu,ka,2) * DM(iv,ib,kb,2) + &
                      DM(ia,iu,ka,3) * DM(iv,ib,kb,4) + &
                      DM(ia,iu,ka,4) * DM(iv,ib,kb,3) 
              ctemp2 = DM_mf(ia,iu,ka,1) * DM_mf(iv,ib,kb,1) + &
                      DM_mf(ia,iu,ka,2) * DM_mf(iv,ib,kb,2) + &
                      DM_mf(ia,iu,ka,3) * DM_mf(iv,ib,kb,4) + &
                      DM_mf(ia,iu,ka,4) * DM_mf(iv,ib,kb,3)
            elseif(lsda) then
              ctemp = DM(ia,iu,ka,1) * DM(iv,ib,kb,1) + &
                      DM(ia,iu,ka,2) * DM(iv,ib,kb,2)
              ctemp2 = DM_mf(ia,iu,ka,1) * DM_mf(iv,ib,kb,1) + &
                      DM_mf(ia,iu,ka,2) * DM_mf(iv,ib,kb,2)
            else
              ctemp = 2.d0 * DM(ia,iu,ka,1) * DM(iv,ib,kb,1)
              ctemp2 = 2.d0 * DM_mf(ia,iu,ka,1) * DM_mf(iv,ib,kb,1)
            endif
            if( (abs(ctemp) < 1e-8) .and. (abs(ctemp2) < 1e-8) ) CYCLE 

            ! (ab|vu) = sum_nc L(nc, ia, ib) * conjg(L(nc, iu, iv))
            ! can get all ib and iv together with DGEMM
            etemp = (0.d0,0.d0) 
            do j=1,nchol
              etemp = etemp + GChol(j,ia,ib)*CONJG(GChol(j,iu,iv))
            enddo
            eX = eX - fXX * etemp * ctemp 
            eX_mf = eX_mf - fXX * etemp * ctemp2

          enddo
          enddo

        enddo
        enddo

      enddo

      if( Q==1 ) then
        call mp_sum(v1,intra_image_comm)
        etemp = (0.d0,0.d0)
        do j=1,nchol
          etemp = etemp + v1(j,1) * conjg(v1(j,1))
        enddo
        eJ = fXX * fac * fac * etemp
        etemp = (0.d0,0.d0)
        do j=1,nchol
          etemp = etemp + v1(j,2) * conjg(v1(j,2))
        enddo
        eJ_mf = fXX * fac * fac * etemp
        deallocate( v1 )
      endif

401   CONTINUE

    enddo
402 CONTINUE
    CALL stop_clock ( 'orb_cholesky' )
    call mp_sum(eX,intra_image_comm)
    call mp_sum(eX_mf,intra_image_comm)
    if(ionode) then
      write(*,*) 'EJ(1Det),EJ(MF) (Ha):',0.5d0*eJ/(nksym*1.0), &
                0.5d0*eJ_mf/(nksym*1.0)
      write(*,*) 'EXX(1Det),EXX(MF) (Ha):',0.5d0*eX/(nksym*1.0), &
                0.5d0*eX_mf/(nksym*1.0)
    endif

    if(noncolin) then
      CALL esh5_posthf_kpoint_info(h5id_hamil%id,nup+ndown,0,e0*nksym,efc*nksym,nksym,xkcart, &
                                kminus,QKtoK2,h5id_orbs%norbK,ncholQ)
    else
      CALL esh5_posthf_kpoint_info(h5id_hamil%id,nup,ndown,e0*nksym,efc*nksym,nksym,xkcart, &
                                kminus,QKtoK2,h5id_orbs%norbK,ncholQ)
    endif

    if(ionode) then
      write(*,*) 'Timers: '
      CALL print_clock ( 'orb_cholesky' )
      CALL print_clock ( 'orb_comm_ovr' )
      CALL print_clock ( 'orb_2el' )
      CALL print_clock ( 'orb_cholvgen' )
      CALL print_clock ( 'orb_diagupd' )
      CALL print_clock ( 'orb_cholwrt' )
    ENDIF

    IF( ALLOCATED(Kqab) ) DEALLOCATE (Kqab)
    IF( ALLOCATED(Vuv) ) DEALLOCATE (Vuv)
    IF( ALLOCATED(Vuv2) ) DEALLOCATE (Vuv2)
    IF( ALLOCATED(Diag) ) DEALLOCATE (Diag)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(maxres) ) DEALLOCATE (maxres)
    IF( ALLOCATED(xkcart) ) DEALLOCATE (xkcart)
    IF( ALLOCATED(comm) ) DEALLOCATE (comm)
    IF( ALLOCATED(ncholQ) ) DEALLOCATE (ncholQ)
    if(allocated(paw_one_center)) deallocate(paw_one_center)
    if(okvan) then
      do ik=1,nkloc
        CALL deallocate_bec_type(becpsia(ik))
        CALL deallocate_bec_type(becpsib(ik))
      enddo
      DEALLOCATE(becpsia)
      DEALLOCATE(becpsib)
    endif
    if(okvan) then
      do i=1,ntyp
        deallocate(ke(i)%L)  
      enddo
      deallocate(ke)
    endif

    if(me_image .ne. root_image) &
      CALL esh5_posthf_close_file(h5id_hamil%id)
    if(nproc_image > 1) call mp_barrier( intra_image_comm )

    if(ionode) then
      h5name = TRIM( hamil_file ) 
      h5len = LEN_TRIM(h5name)
      call esh5_posthf_join_all(h5id_hamil%id,h5name,h5len,nksym,nproc_image,error)
      if(error .ne. 0) then
        write(*,*) 'Error: ',error
        call errore('pw2posthf','Error in esh5_posthf_join_all',1)
      endif
      CALL esh5_posthf_close_file(h5id_hamil%id)
    endif
    call close_esh5_read(h5id_orbs)
    if(nproc_image > 1) call mp_barrier( intra_image_comm )

  END SUBROUTINE cholesky_r
  !
  SUBROUTINE eri_df(thresh,dfft,hamil_file,orb_file,df_file)
    USE parallel_include
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp,                    ONLY : okvan
    USE paw_variables,           ONLY : okpaw
    USE exx, ONLY: ecutfock
    USE gvect,     ONLY : ecutrho
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_file,orb_file,df_file
    REAL(DP), INTENT(IN) :: thresh
    !
    CHARACTER(len=256) :: h5name
    INTEGER :: h5len, oldh5
    INTEGER :: numdf, naux, nqdf
    INTEGER :: n1, n2, n3
    INTEGER :: Q, Qp, ka, kb, iab, ia, ib, iu, iv, ku, kv, ni, nj, iuv
    INTEGER :: i,n,m,ii,jj,kk, ia0, iu0, ispin, is1, is2, noa, nob
    INTEGER :: a_beg, a_end, b_beg, b_end, ab_beg, ab_end, nabpair   ! assigned orbitals/pairs
    INTEGER :: naorb, nborb, iabmax,nkloc,k_beg,k_end 
    INTEGER :: error,ibnd,ik,j, maxnorb, many_ab
    INTEGER, ALLOCATABLE :: Q2xk(:)
    INTEGER, ALLOCATABLE :: g2g(:)
    REAL(DP), ALLOCATABLE :: xkcart(:,:)
    REAL(DP) :: fXX, fac, rtemp
    REAL(DP) :: xk0(3), dQ(3), dG(3, 27)
    COMPLEX(DP) :: ctemp, ctemp2, etemp, E1
    LOGICAL :: more
    ! Allocatables 
    COMPLEX(DP), ALLOCATABLE :: orb(:) 
    COMPLEX(DP), ALLOCATABLE :: Vpq(:,:)    
    COMPLEX(DP), ALLOCATABLE :: Vinvpq(:,:)    
    COMPLEX(DP), ALLOCATABLE :: Pab(:,:)    ! 3-center integrals in df basis 
    COMPLEX(DP), ALLOCATABLE :: LPab(:,:)   ! factorized ERI  
    COMPLEX(DP), ALLOCATABLE :: CPab(:,:)   ! factorized ERI  
    ! Work Arrays
    COMPLEX(DP), ALLOCATABLE :: T1(:,:)    
    COMPLEX(DP), ALLOCATABLE :: T2(:,:)    
    COMPLEX(DP), ALLOCATABLE :: T3(:,:)    
    COMPLEX(DP), ALLOCATABLE :: PsiP(:,:)    ! DF basis in G space  
    COMPLEX(DP), ALLOCATABLE :: Psia(:,:,:)    ! orbitals in real space
    COMPLEX(DP), ALLOCATABLE :: Psib(:,:,:)    ! orbitals in real space
    COMPLEX(DP), ALLOCATABLE :: Kp(:,:)        ! invfft[Psip(G) * v(G+Q)] * exp(iQr) 
    COMPLEX(DP), ALLOCATABLE :: Tab(:,:)       ! Tab(r,pa) = conj(Psia(r)) * Psib(r) 
    COMPLEX(DP), ALLOCATABLE :: phasefac(:,:)
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:)      ! 
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
    COMPLEX(DP), ALLOCATABLE :: v1(:,:)  
    INTEGER, ALLOCATABLE :: nauxQ(:)
    REAL(DP), ALLOCATABLE :: sval(:)
    COMPLEX(DP)::eX,eJ,eX_mf,eJ_mf,eX2,eJ2
    !
    CHARACTER(len=4) :: str_me_image
    TYPE(h5file_type) :: h5id_orbs, h5id_hamil, h5id_df 
    !
    if(okvan .or. okpaw) &
      call errore('eri_df','No PAW or US allowed yet.',1)

    if(thresh > 1.d-8) then
      write(*,*) '\n*************************************************'
      write(*,*) 'WARNING: SVD threshold in eri_df seems large.'
      write(*,*) 'Consider using a smaller value, e.g. 1e-8.'
      write(*,*) '*************************************************\n'
    endif

    ! open hamiltonian file
    if(ionode) then
      oldh5=1
      h5name = TRIM(hamil_file) 
    else
      oldh5=0
      WRITE ( str_me_image, '(I4)') me_image
      h5name = TRIM(hamil_file) //"_part"//adjustl(trim(str_me_image))
    endif
    h5len = LEN_TRIM(h5name)
    CALL esh5_posthf_open_file(h5id_hamil%id,h5name,h5len,oldh5)
    if(oldh5 .ne. 0 ) &
      call errore('eri_df','error opening hamil file',1)

    ! open orbital file 
    h5name = TRIM( orb_file ) 
    call open_esh5_read(h5id_orbs,h5name)
    if( h5id_orbs%grid_type .ne. 1 ) &
      call errore('eri_df','grid_type ne 1',1)
    if( h5id_orbs%nspin .ne. 1 ) &
      call errore('eri_df','nspin ne 1',1)
    maxnorb = h5id_orbs%maxnorb

    ! open df file 
    h5name = TRIM( df_file ) 
    call open_esh5_read(h5id_df,h5name,.false.)
    if( h5id_df%nspin .ne. 1 ) &
      call errore('eri_df','nspin ne 1 in df file',1)
    if( h5id_df%nkpts .ne. nksym ) &
      call errore('eri_df','nkpts ne nksym in df file',1)
    if( h5id_df%nr1*h5id_df%nr2*h5id_df%nr3 .ne. dfft%nnr ) &
      call errore('eri_df','FFT grids do not match.',1)

    allocate( Q2xk(nksym) )
    do i=1,nksym
      ! find h5id_df%xk in Qpts list, make sure it exists
      ik=-1
      do jj=1,nksym
        if(sum( (h5id_df%xk(1:3,jj)-Qpts(1:3,i))**2 ) < 1.d-6) then
          ik = jj 
          exit
        endif
      enddo
      if(ik > 0) then
        Q2xk(i) = ik
      else
        write(*,*) '\n********************************************'
        write(*,*) '                   ERROR'
        write(*,*) '********************************************'
        write(*,*) 'missing Qpt:',i,Qpts(1:3,i)*tpiba
        write(*,*) 'Qpt grid in file: '
          write(*,*) '   k,xk(:):'
        do n=1,h5id_df%nkpts
          write(*,*) '   ',n,h5id_df%xk(1:3,n)*tpiba
        enddo
        call errore('eri_df',' error: Q-point grids do not agree',1)
      endif
    enddo
    numdf = maxval(h5id_df%norbK(:))
    if(h5id_df%grid_type == 0) then
      ! setup g2g mapping
      allocate(g2g(dfft%nnr), orb(dfft%nnr))
      n=1
      do kk=0,h5id_df%nr3-1
        n3 = kk
        if( n3 > (h5id_df%nr3-1)/2 ) n3 = n3 - h5id_df%nr3 + dfft%nr3
        do jj=0,h5id_df%nr2-1
          n2 = jj
          if( n2 > (h5id_df%nr2-1)/2 ) n2 = n2 - h5id_df%nr2 + dfft%nr2
          do ii=0,h5id_df%nr1-1
            n1 = ii
            if( n1 > (h5id_df%nr1-1)/2 ) n1 = n1 - h5id_df%nr1 + dfft%nr1
            ik = 1 + n1 + n2 * dfft%nr1x + n3 * dfft%nr1x * dfft%nr2x
            g2g(n) = ik
            n = n + 1
          enddo
        enddo
      enddo
    elseif(h5id_df%grid_type == 2) then
      ! nothing to do
    else
      call errore('eri_df','unknown grid_type in df file',1)
    endif

    nmax_DM = h5id_orbs%nmax_DM
    if( nmax_DM <= 0 ) &
      call errore('cholesky','error nmax_DM <= 0',1)
    if( allocated(DM) ) then
      if( (size(DM,1).ne.nmax_DM) .or. &  
          (size(DM,2).ne.nmax_DM) .or. &
          (size(DM,3).ne.nksym) .or. &
          (size(DM,4).ne.nspin) ) then  
        deallocate(DM)
        allocate( DM(nmax_DM, nmax_DM, nksym, nspin) )
      endif 
    else
      allocate( DM(nmax_DM, nmax_DM, nksym, nspin) )
    endif
    if( allocated(DM_mf) ) then
      if( (size(DM_mf,1).ne.nmax_DM) .or. &  
          (size(DM_mf,2).ne.nmax_DM) .or. &
          (size(DM_mf,3).ne.nksym) .or. &
          (size(DM_mf,4).ne.nspin) ) then  
        deallocate(DM_mf)
        allocate( DM_mf(nmax_DM, nmax_DM, nksym, nspin) )
      endif  
    else
      allocate( DM_mf(nmax_DM, nmax_DM, nksym, nspin) )
    endif
    do i=1,nspin
      do ik=1,nksym
        call esh5_posthf_read_dm(h5id_orbs%id,"DM",2,ik-1,i-1,DM(1,1,ik,i),error)
        if(error .ne. 0 ) &
          call errore('cholesky','error reading DM',1)
        call esh5_posthf_read_dm(h5id_orbs%id,"DM_mf",5,ik-1,i-1,DM_mf(1,1,ik,i),error)
        if(error .ne. 0 ) &
          call errore('cholesky','error reading DM',1)
      enddo
    enddo

    if(nspin == 1 ) then
      fac = 2.d0 
    else
      fac = 1.d0 
    endif

    ! parallelization has restrictions
    if(nproc_pool > 1 .and. npool .ne. nksym) &
      call errore('pw2posthf','Error: nb > 1 only allowed if npools == nkstot',1) 

    !   Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    !
    ! Data Distribution setup
    !   Partition orbital pairs:  Simple for now, find more memory friendly version
    call fair_divide(ab_beg,ab_end,me_pool+1,nproc_pool,maxnorb*maxnorb)
    nabpair = ab_end - ab_beg + 1

    ! somewhat arbitrary
    !many_ab = min(128, nabpair)  
    many_ab = nabpair

    !   Determine assigned orbitals
! fix this
!    a_beg   = (ab_beg-1)/maxnorb + 1
!    a_end   = (ab_end-1)/maxnorb + 1
    a_beg   = 1 
    a_end   = maxnorb 
    naorb   = a_end-a_beg+1
    b_beg   = 1
    b_end   = maxnorb   ! keeping all for simplicity now 
    nborb   = b_end-b_beg+1

    !  Allocate space for assigned orbitals in real space and corresponding
    !  orbital pairs
    allocate( Kp(dfft%nnr,numdf), Tab(dfft%nnr,many_ab), PsiP(dfft%nnr,numdf) ) 
    allocate( Psia(dfft%nnr,naorb,1), Psib(dfft%nnr,nborb,1) )
    allocate( T1(numdf,numdf), T2(numdf,numdf) ) 
    allocate( Vinvpq(numdf,numdf), CPab(numdf, nabpair) )
    allocate( Vpq(numdf,numdf), Pab(numdf,nabpair), LPab(numdf,nabpair) ) 
    allocate( vCoulG(dfft%nnr) )
    allocate( phasefac(dfft%nnr,27) )
    allocate( nauxQ(nksym),  xkcart(3,nksym) )
    allocate( spsi(1, 1), sval(numdf) )

    ! fix with symmetry 
    do ik=1,nksym
      xkcart(1:3,ik) = xksym(1:3,ik)*tpiba
    enddo

    ! generate phase factors in real space, phasefac(ir,G) = exp(i*G*ir), 
    ! where G are the Fourier frequencies with indexes {-1,0,1}, 27 in total 
    ! if memory cost is too much, distribute over dfft%nnr and gather over proc
    ! grid
    CALL calculate_phase_factor(dfft, phasefac, dG)

    if(ionode) write(*,*) 'Starting loop over Q vectors.'
    !
    ! 
    ! Loop over Q points
    !
    eX = (0.d0,0.d0)
    eJ = (0.d0,0.d0)
    eX_mf = (0.d0,0.d0)
    eJ_mf = (0.d0,0.d0)
    do Q=1,1
    !do Q=1,nksym

      if( Q .gt. kminus(Q) ) cycle 
      if( Q .eq. kminus(Q) ) then
        fXX=1.d0
      else
        fXX=2.d0
      endif
      ! the Q-point of the basis function is -Q
      Qp = kminus(Q)    

      if(ionode) write(*,*) ' Q: ',Q
      nqdf = h5id_df%norbK(Q2xk(Qp))

      !!! MAM: right now there is no parallelization on the generation of
      !!!      of V^{-1/2}. Distribute work among all procs 
  
      ! read df basis for kpt=Q  
      ! need to use get_psi_esh5 since the orbital might need modifications
      PsiP(:,:) = (0.d0,0.d0)
      do ia=1,nqdf
        if(h5id_df%grid_type == 0) then
          call get_psi_esh5(h5id_df,ia,Q2xk(Qp),1,orb(:)) 
          PsiP(g2g(1:dfft%nnr),ia) = orb(1:dfft%nnr)
        elseif(h5id_df%grid_type == 2) then
          call get_psi_esh5(h5id_df,ia,Q2xk(Qp),1,PsiP(:,ia)) 
        endif
        ! normalized in G-space on dfft grid.
        call Overlap(1,1,dfft%nnr,(1.d0,0.d0),PsiP(1,ia),size(PsiP,1), &
                       PsiP(1,ia),size(PsiP,1),(0.d0,0.d0),T1(1,1),1,  &
                       .false.,spsi(:,1))
        rtemp = 1.d0 / sqrt(dble(T1(1,1)))
        PsiP(:,ia) = PsiP(:,ia) * rtemp
      enddo

      ! Setting up density fitting basis 
      ! 1. calculate overlap matrix
      T1(:,:) = (0.d0,0.d0)
      call Overlap(nqdf,nqdf,dfft%nnr,(1.d0,0.d0),PsiP,size(PsiP,1), &
                       PsiP,size(PsiP,1),(0.d0,0.d0),T1,size(T1,1),  &
                       .false.,spsi)

      ! 2. orthonormalize basis, assume linear dependencies 
      call svdh(nqdf, t1, size(T1, 1), sval) 

      ! build orthonormal basis 
      naux=0
      do ia = 1, nqdf 
        if(sval(ia) .gt. thresh) then  
          T1(:,ia)=T1(:,ia)/sqrt(sval(ia))
          naux = naux + 1 
        else
          exit
        endif
      enddo
      nauxQ(Q) = naux

      ! rotate PsiP(:,q) = T1(p,q) PsiP(:,q)
      allocate( T3(dfft%nnr, naux) )
      call zgemm('N','N',dfft%nnr,naux,nqdf,(1.d0,0.d0),PsiP,dfft%nnr,  &
                 T1,size(T1,1),(0.d0,0.d0),T3,dfft%nnr)
      PsiP(:,1:naux) = T3
      deallocate( T3 )

      ! sanity check
      call Overlap(naux,naux,dfft%nnr,(1.d0,0.d0),PsiP,size(PsiP,1), &
                       PsiP,size(PsiP,1),(0.d0,0.d0),T1,size(T1,1),  &
                       .false.,spsi)
      do ii=1,naux
        if( abs(T1(ii,ii)-1.d0) > 1.0d-8 )  then 
          write(*,*) ii,T1(ii,ii)
          call errore('eri_df','DF basis not orthonormal after rotation.',1)
        endif
        do jj=ii+1,naux
          if( abs(T1(ii,jj)) > 1.0d-8 )  then 
            write(*,*) ii,jj,T1(ii,jj)
            call errore('eri_df','DF basis not orthonormal after rotation.',1)
          endif
        enddo
      enddo

      ! 3. Calculate coulomb matrix elements 
      ! T1(p,q) = <PsiP(r,p) | v(r,r') | PsiP(r',q)> 
      !         = sum_G conj(PsiP(G,p)) * v(G+Qp) * PsiP(G,q)
      ! vCoulG( iG ) = | G + Qp |^{-2} = | G + Qp - 0  |^{-2}
      xk0(1:3) = 0.d0
      CALL g2_convolution(dfft%ngm, g, Qpts (1:3, Qp), xk0, psic) 

      vCoulG(:) = (0.d0,0.d0)
      vCoulG( dfft%nl(1:dfft%ngm) ) = e2Ha * psic(1:dfft%ngm)
      if(gamma_only) vCoulG( dfft%nlm(1:dfft%ngm) ) = e2Ha * conjg(psic(1:dfft%ngm))

      T1(:,:) = (0.d0,0.d0)
      do ii=1,naux
        do jj=ii,naux
          ctemp=(0.d0,0.d0)
          do n=1,dfft%nnr
            ctemp=ctemp+conjg(PsiP(n,ii))*vCoulG(n)*PsiP(n,jj)
          enddo
          T1(ii,jj) = ctemp
          if(jj > ii) T1(jj,ii) = conjg(ctemp)
        enddo
      enddo

      ! 4. Calculate Vpq = T1^{-1/2} = T1^{-1} * T1^{1/2} = U s^{-1/2} U^+ 
      call eigsys('V', 'U', .false., naux, size(T1,1), T1, sval)
      T2(1:naux,1:naux) = T1(1:naux,1:naux)
      do ia = 1, naux
        T1(:,ia)=T1(:,ia)/sqrt(sval(ia))
      enddo
      call zgemm('N','C',naux,naux,naux,(1.d0,0.d0),T1,size(T1,1),  &
                 T2,size(T2,1),(0.d0,0.d0),Vpq,size(Vpq,1))

      ! Vinv
      do ia = 1, naux
        T1(:,ia)=T2(:,ia)/sval(ia)
      enddo
      call zgemm('N','C',naux,naux,naux,(1.d0,0.d0),T1,size(T1,1),  &
                 T2,size(T2,1),(0.d0,0.d0),Vinvpq,size(Vinvpq,1))

      if( Q==1 ) then
        allocate( v1(naux,2) )
        v1(:,:) = (0.d0,0.d0)
      endif

      ! calculate Kp(r) = invfft[ v(G+Qp) * PsiP(G) ] 
      ! vCoulG( iG ) already calculated, reuse 
      Kp(:,:)=(0.d0,0.d0)
      do n=1,naux
        !
        ! multiply by FFT[ 1/r ]
        do j = 1,dfft%nnr
          Kp(j,n) = PsiP(j,n) * vCoulG(j) 
        enddo

        ! invfft kab to R
        CALL invfft ('Rho', Kp(:,n), dfft)

        ! add normalization 
        Kp(1:dfft%nnr,n) = Kp(1:dfft%nnr,n) /  &
              (dfft%nr1*dfft%nr2*dfft%nr3*sqrt(nksym*omega))
        !
      enddo

      ! remove later on  
      do n=1,naux  
        CALL fwfft ('Rho', PsiP(:,n), dfft)
      enddo  

      do ik=1,nkloc

        ! Q = ka - kb  !!! careful not to confuse Q with Qp 
        ka = k_beg + ik - 1
        kb = QKtoK2(Q, ka)
        iabmax = h5id_orbs%norbK(ka)*h5id_orbs%norbK(kb)

        ! dQ = Qba - Qp = kb - ka - Qpts(Qp) 
        dQ(1:3) = xksym(1:3, kb) - xksym(1:3, ka) - Qpts(1:3,Qp) 
        ii=0
        do jj=1,27
          if(sum( (dQ(1:3)-dG(1:3,jj))**2 ) .lt. 1.d-8) then
            ii=jj
            exit
          endif
        enddo
        if(ii.lt.1) call errore('eri_df','Can not find dQ in G list.',1)

        ! read orbitals, why read them all?
        call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,&
                    1,Psia,a_beg,naorb,ka,1)
        call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,&
                    1,Psib,b_beg,nborb,kb,1)
 
        ! 5. Calculate Pab = (P|ab) 
        !    (P|ab) = sum_r conjg(Kp(r)) * conjg(Psia(r)) * Psib(r) * exp(i*(kb-ka-Qp)*r),

        ! add phase factor
        do n=1,naux 
          Kp(1:dfft%nnr,n) = Kp(1:dfft%nnr,n) * &
                             phasefac(1:dfft%nnr, ii)
        enddo

        Pab(:,:)=(0.d0,0.d0)
        do ibnd=1,nabpair,many_ab 

          n = 0
          do m=1,min(many_ab, nabpair-ibnd+1)

            ! iab = (ia-1)*norbK(kb) + ib 
            iab = ab_beg + ibnd - 1 + m - 1 
            if(iab > iabmax) exit
            ia = (iab-1)/h5id_orbs%norbK(kb) + 1
            ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1

            n=n+1
            do j = 1,dfft%nnr
              Tab(j,n) = conjg(Psia(j,ia-a_beg+1,1)) *  &
                               Psib(j,ib-b_beg+1,1)        
            enddo

          enddo 

          ! (P|ab) = Pab = sum_r' conj(Kp(r)) Tab(r) 
          if(n > 0) then 
            call zgemm('C','N',naux,n,dfft%nnr,               &
                   (1.d0,0.d0),Kp,size(Kp,1),Tab,size(Tab,1), &
                   (0.d0,0.d0),Pab(1,ibnd),size(Pab,1))
          endif

        enddo

!        Pab(:,:)=(0.d0,0.d0)
!        do ibnd=1,nabpair
!          iab = ab_beg + ibnd - 1 
!          if(iab > iabmax) exit
!          ia = (iab-1)/h5id_orbs%norbK(kb) + 1
!          ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1
!          do n=1,naux
!          do j=1,dfft%nnr
!            Pab(n,ibnd) = Pab(n,ibnd) + conjg(Kp(j,n) * Psia(j,ia-a_beg+1,1) ) * Psib(j,ib-b_beg+1,1) 
!          enddo
!          enddo
!        enddo

        ! remove phase factor
        do n=1,naux 
          Kp(1:dfft%nnr,n) = Kp(1:dfft%nnr,n) / &
                               phasefac(1:dfft%nnr, ii)
        enddo

        ! 6. Get ERI factorization: LPab = V^{-1/2} * (P|ab) = Vpq * Pab  
        call zgemm('N','N',naux,nabpair,naux,                  &
                 (1.d0,0.d0),Vpq,size(Vpq,1),Pab,size(Pab,1),  &
                 (0.d0,0.d0),LPab,size(LPab,1))

        call zgemm('N','N',naux,nabpair,naux,                  &
                 (1.d0,0.d0),Vinvpq,size(Vinvpq,1),Pab,size(Pab,1),  &
                 (0.d0,0.d0),CPab,size(CPab,1))

!now output only the segment associated with nkloc

!        write(*,*) 'Q,ki,kl,i,k,j,l,V:'
!        do ia=1,2
!        do ib=1,2
!        do iu=1,2
!        do iv=1,2
            
!          iab = (ia-1)*h5id_orbs%norbK(kb)+ib
!          iuv = (iu-1)*h5id_orbs%norbK(kb)+iv
!          ctemp = (0.d0,0.d0) 
!          do n=1,naux
!            ctemp = ctemp + LPab(n,iab) * conjg(LPab(n,iuv))
!          enddo
!          write(100+me_image,'(7i4,2g20.8)') Q,ka,ka,ia,ib,iv,iu,real(ctemp),imag(ctemp) 
      
!        enddo
!        enddo
!        enddo
!        enddo
        if(ka==1) then
          ia=1
          ib=1
          iu=1
          iv=1
          iab = (ia-1)*h5id_orbs%norbK(kb)+ib
          iuv = (iu-1)*h5id_orbs%norbK(kb)+iv
          ctemp = (0.d0,0.d0) 
          do n=1,naux
            ctemp = ctemp + LPab(n,iab) * conjg(LPab(n,iuv))
          enddo
          write(*,'(7i4,2g20.8)') &
          Q,ka,ka,ia,ib,iv,iu,real(ctemp),imag(ctemp)

          T1(:,:) = (0.d0,0.d0)
          do n=1,naux
            PsiP(:,n) = PsiP(:,n)*sqrt(dfft%nnr*1.d0)
            do j=1,dfft%nnr
              T1(n,1) = T1(n,1) + PsiP(j,n)
              T1(n,2) = T1(n,2) + conjg(PsiP(j,n))*PsiP(j,n)
            enddo          
            T1(n,1) = T1(n,1)
            T1(n,2) = T1(n,2)
            write(*,*) 'PsiP: ',n,T1(n,2)
          enddo          

          CPab(:,:) = (0.d0,0.d0)
          do ibnd=1,nabpair
            iab = ab_beg + ibnd - 1 
            if(iab > iabmax) exit
            ia = (iab-1)/h5id_orbs%norbK(kb) + 1
            ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1
            do j = 1,dfft%nnr
              Tab(j,ibnd) = conjg(Psia(j,ia-a_beg+1,1)) *  &
                                  Psib(j,ib-b_beg+1,1)
            enddo
          enddo
          call zgemm('C','N',naux,nabpair,dfft%nnr,                  &
                     (1.d0,0.d0),PsiP,size(PsiP,1),Tab,size(Tab,1),  &
                     (0.d0,0.d0),CPab,size(CPab,1))

          do ia=1,h5id_orbs%norbK(ka)
          do ib=ia,ia
            iab = (ia-1)*h5id_orbs%norbK(kb)+ib
            ctemp=(0.d0,0.d0)
            do n=1,naux
              ctemp = ctemp + CPab(n,iab)*T1(n,1)
            enddo
            ctemp = ctemp/dfft%nnr
            write(*,*) '   ',ia,ib,real(ctemp),imag(ctemp)    
          enddo
          enddo

          do ia=1,h5id_orbs%norbK(ka)
          do ib=ia,ia
            iab = (ia-1)*h5id_orbs%norbK(kb)+ib
            rtemp=(0.d0,0.d0)      
            do j=1,dfft%nnr
              ctemp2 = conjg( Psia(j,ia-a_beg+1,1) ) * Psib(j,ib-b_beg+1,1) 
              etemp = (0.d0,0.d0)
              do n=1,naux
                etemp = etemp + CPab(n,iab) * PsiP(j,n) 
              enddo
write(200+ia,*) j,real(etemp),imag(etemp),real(ctemp2),imag(ctemp2)
              ctemp2 = ctemp2 - etemp 
              rtemp = rtemp + ctemp2 * conjg(ctemp2) 
            enddo
            rtemp = sqrt(rtemp/dfft%nnr)
            write(*,*) '   ',ia,ib,rtemp
          enddo
          enddo

        elseif(ka==2) then
          ia=2
          ib=2
          iu=2
          iv=2
          iab = (ia-1)*h5id_orbs%norbK(kb)+ib
          iuv = (iu-1)*h5id_orbs%norbK(kb)+iv
          ctemp = (0.d0,0.d0) 
          do n=1,naux
            ctemp = ctemp + LPab(n,iab) * conjg(LPab(n,iuv))
          enddo

          rtemp=0.d0        
          do j=1,dfft%nnr
            ctemp2 = conjg( Psia(j,ia-a_beg+1,1) ) * Psib(j,ib-b_beg+1,1)         
            do n=1,naux
              ctemp2 = ctemp2 - CPab(n,iab) * PsiP(j,n)
            enddo
            rtemp = rtemp + abs(ctemp2)
          enddo
          rtemp = rtemp / dfft%nnr

          write(*,'(7i4,3g20.8)') &
          Q,ka,ka,ia,ib,iv,iu,real(ctemp),imag(ctemp),rtemp
        elseif(ka>2) then
          goto 991
        endif

        ! EJ  
        if( Q == 1 ) then
          ! accumulate EJ
          noa = min(nmax_DM,h5id_orbs%norbK(ka))
          nob = min(nmax_DM,h5id_orbs%norbK(kb))

          do ibnd=1,nabpair
            iab = ab_beg + ibnd - 1
            if(iab > iabmax) exit
            ia = (iab-1)/h5id_orbs%norbK(kb) + 1
            ib = MOD((iab-1), h5id_orbs%norbK(kb)) + 1
            if(ia > noa .or. ib > nob) cycle

            do ispin=1,min(2,nspin)
              !  
              etemp = DM(ia,ib,ka,ispin)
              v1(1:naux,1) = v1(1:naux,1) + LPab(1:naux,ibnd)*etemp
              etemp = DM_mf(ia,ib,ka,ispin)
              v1(1:naux,2) = v1(1:naux,2) + LPab(1:naux,ibnd)*etemp
              !
            enddo
            !
          enddo ! ibnd
          !
        endif ! Q==1

      enddo ! ik

      if(Q==1) then
        !
        call mp_sum(v1,intra_image_comm)
        etemp = (0.d0,0.d0)
        do j=1,naux
          etemp = etemp + v1(j,1) * conjg(v1(j,1))
        enddo
        eJ = fXX * fac * fac * etemp
        etemp = (0.d0,0.d0)
        do j=1,naux
          etemp = etemp + v1(j,2) * conjg(v1(j,2))
        enddo
        eJ_mf = fXX * fac * fac * etemp
        !
        deallocate(v1)
      endif ! Q==1

      if(ionode) write(*, *) ' Done with Q, naux(Q):',Q,nauxQ(Q)

    enddo !Q

    call mp_sum(eX,intra_image_comm)
    call mp_sum(eX_mf,intra_image_comm)
    if(ionode) then
      write(*,*) 'EJ(1Det),EJ(MF) (Ha):',0.5d0*eJ/(nksym*1.0), &
                0.5d0*eJ_mf/(nksym*1.0)
      write(*,*) 'EXX(1Det),EXX(MF) (Ha):',0.5d0*eX/(nksym*1.0), &
                0.5d0*eX_mf/(nksym*1.0)
    endif

    if(ionode) then
      write(*,*) 'Timers: '
    ENDIF
991 continue


    IF( ALLOCATED(Pab) ) DEALLOCATE (Pab)
    IF( ALLOCATED(Vpq) ) DEALLOCATE (Vpq)
    IF( ALLOCATED(Vinvpq) ) DEALLOCATE (Vinvpq)
    IF( ALLOCATED(LPab) ) DEALLOCATE (LPab)
    IF( ALLOCATED(CPab) ) DEALLOCATE (CPab)
    IF( ALLOCATED(Tab) ) DEALLOCATE (Tab)
    IF( ALLOCATED(T1) ) DEALLOCATE (T1)
    IF( ALLOCATED(T2) ) DEALLOCATE (T2)
    IF( ALLOCATED(T3) ) DEALLOCATE (T3)
    IF( ALLOCATED(Psia) ) DEALLOCATE (Psia)
    IF( ALLOCATED(Psib) ) DEALLOCATE (Psib)
    IF( ALLOCATED(PsiP) ) DEALLOCATE (PsiP)
    IF( ALLOCATED(Kp) ) DEALLOCATE (Kp)
    IF( ALLOCATED(Q2xk) ) DEALLOCATE (Q2xk)
    IF( ALLOCATED(vCoulG) ) DEALLOCATE (vCoulG)
    IF( ALLOCATED(phasefac) ) DEALLOCATE (phasefac)
    IF( ALLOCATED(xkcart) ) DEALLOCATE (xkcart)
    IF( ALLOCATED(nauxQ) ) DEALLOCATE (nauxQ)
    IF( ALLOCATED(g2g) ) DEALLOCATE (g2g)
    IF( ALLOCATED(orb) ) DEALLOCATE (orb)
    IF( ALLOCATED(spsi) ) DEALLOCATE (spsi)
    IF( ALLOCATED(spsi) ) DEALLOCATE (spsi)
    IF( ALLOCATED(v1) ) DEALLOCATE (v1)

    if(me_image .ne. root_image) &
      CALL esh5_posthf_close_file(h5id_hamil%id)
    if(nproc_image > 1) call mp_barrier( intra_image_comm )

    if(ionode) then
      h5name = TRIM( hamil_file ) 
      h5len = LEN_TRIM(h5name)
      call esh5_posthf_join_all(h5id_hamil%id,h5name,h5len,nksym,nproc_image,error)
      if(error .ne. 0) then
        write(*,*) 'Error: ',error
        call errore('pw2posthf','Error in esh5_posthf_join_all',1)
      endif
      CALL esh5_posthf_close_file(h5id_hamil%id)
    endif
    call close_esh5_read(h5id_orbs)
    if(nproc_image > 1) call mp_barrier( intra_image_comm )

  END SUBROUTINE eri_df

  ! calculate position dependent range parameter used in basis set correction scheme  
  ! assumes that Psia, DM and Chol have been calculated
  ! not yet distributed  
  SUBROUTINE calculate_KS_bscorr(dfft,hamil_file,orb_file)
    USE scf, ONLY: rho
    USE constants, ONLY: sqrtpi, sqrt2
    USE wvfct, ONLY: btype
    ! 
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_file,orb_file
    !
    CHARACTER(len=256) :: h5name
    INTEGER :: h5len, oldh5
    INTEGER :: ik,ki,kj,kk,kl,ka,kb,i,j,k,l,a,b,n,p,iik,ilj,Q
    INTEGER :: nc0, nc1, nchol_max, nchol, error
    REAL(DP) :: nel1,nel2, tmp, ecloc
    REAL(DP) :: s, bscorr, rhox, arhox, ec, vc(2), rs, beta, zeta
    COMPLEX(DP) :: n2_
    COMPLEX(DP), ALLOCATABLE :: Gp(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: n1(:,:),buff(:,:)
    COMPLEX(DP), ALLOCATABLE :: v(:,:)
    COMPLEX(DP), ALLOCATABLE :: Chol(:,:,:)    ! Cholesky Matrix 
    COMPLEX(DP), ALLOCATABLE :: Psi(:,:,:)    ! orbitals in real space
    COMPLEX(DP), ALLOCATABLE :: KS_n2(:)   
    COMPLEX(DP), ALLOCATABLE :: KS_ur(:)   
    INTEGER, ALLOCATABLE :: ncholQ(:)
    REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                           vanishing_mag    = 1.D-20, &
                           small = 1.E-10_DP, & 
                           third = 1.0_DP / 3.0_DP, &
                           pi34 = 0.6203504908994_DP, &
                           fact = -2.04311121652593_DP 
    TYPE(h5file_type) :: h5id_orbs, h5id_hamil
    !
    
    h5name = TRIM(hamil_file)
    h5len = LEN_TRIM(h5name)
    CALL esh5_posthf_open_file_read(h5id_hamil%id,h5name,h5len,error)
    if(error .ne. 0 ) &
      call errore('cholesky','error opening hamil file',1)

    ! open orbital file 
    h5name = TRIM( orb_file )
    call open_esh5_read(h5id_orbs,h5name)
    if( h5id_orbs%grid_type .ne. 1 ) &
      call errore('cholesky','grid_type ne 1',1)

    nmax_DM = h5id_orbs%nmax_DM
    if( nmax_DM <= 0 ) &
      call errore('cholesky','error nmax_DM <= 0',1)
    allocate( DM(nmax_DM, nmax_DM, nksym, nspin) )
    do i=1,nspin
      do ik=1,nksym
        call esh5_posthf_read_dm(h5id_orbs%id,"DM",2,ik-1,i-1,DM(1,1,ik,i),error)
        if(error .ne. 0 ) &
          call errore('cholesky','error reading DM',1)
      enddo
    enddo

    if(noncolin) &
      call errore('pw2posthf','noncolin not yet in calculate_KS_ur',1)
    
    ! eventually setup an independent grid here and transform Psi to this grid
    ! not efficient for large cells

    allocate( ncholQ(nksym) )
    call esh5_posthf_read_nchol(h5id_hamil%id,ncholQ, error) 
    if(error .ne. 0) call errore('pw2posthf','Error: esh5_posthf_read_nchol',error)
    nchol_max = maxval(ncholQ(:))

    allocate( Chol(nchol_max,h5id_orbs%maxnorb*h5id_orbs%maxnorb,2), &
              Psi(dfft%nnr,h5id_orbs%maxnorb,nksym) )
    allocate( KS_ur(dfft%nnr), KS_n2(dfft%nnr), Gp(dfft%nnr,nmax_DM,2)  )   
    allocate( v(dfft%nnr,2), n1(dfft%nnr,2) )

    call get_orbitals_set(h5id_orbs, 'esh5','psir',dfft,1,Psi,1,h5id_orbs%maxnorb,1,nksym)
                
    n1(:,:) = (0.d0,0.d0)
    KS_ur(:) = (0.d0,0.d0)
    KS_n2(:) = (0.d0,0.d0)

    ! calculate n2
    do ki=1,nksym
      call zgemm('N','T',dfft%nnr,nmax_DM,nmax_DM,(1.d0,0.d0),Psi(1,1,ki),dfft%nnr,  &
                     DM(1,1,ki,1),nmax_DM,(0.d0,0.d0),Gp(1,1,1),dfft%nnr)
      do i=1,nmax_DM
        n1(1:dfft%nnr,1) = n1(1:dfft%nnr,1) + Gp(1:dfft%nnr,i,1) * &
                                    conjg(Psi(1:dfft%nnr,i,ki))
      enddo
    enddo
    if( nspin > 1 ) then
      do ki=1,nksym
        call zgemm('N','T',dfft%nnr,nmax_DM,nmax_DM,(1.d0,0.d0),Psi(1,1,ki),dfft%nnr,  &
                     DM(1,1,ki,2),nmax_DM,(0.d0,0.d0),Gp(1,1,1),dfft%nnr)
        do i=1,nmax_DM
          n1(1:dfft%nnr,2) = n1(1:dfft%nnr,2) + Gp(1:dfft%nnr,i,2) * &
                                      conjg(Psi(1:dfft%nnr,i,ki))
        enddo
      enddo
    else
      n1(1:dfft%nnr,2) = n1(1:dfft%nnr,1)
    endif

    ! conjugating Psi
    do k=1,size(Psi,3)
      do j=1,size(Psi,2)
        do i=1,size(Psi,1)
          Psi(i,j,k) = conjg(Psi(i,j,k))
        enddo
      enddo
    enddo

    do Q=1,nksym
     write(*,*) ' KS u(r), Q:',Q
     if(nksym > 1) &
       call errore('pw2posthf',' Error: need to read Chol for every Q point.',1)

     ! temporary parallelization, terrible for large supercells !!!
     call fair_divide(nc0,nc1,me_pool+1,nproc_pool,ncholQ(Q))

     do ki=1,nksym

      kk = QKtoK2(Q,ki) 
      ka = kk

      ! get Gp(up)
      call zgemm('N','N',dfft%nnr,nmax_DM,nmax_DM,(1.d0,0.d0),Psi(1,1,ka),dfft%nnr,  &
                         DM(1,1,ka,1),nmax_DM,(0.d0,0.d0),Gp(1,1,1),dfft%nnr)      

      call esh5_posthf_read_cholesky(h5id_hamil%id,Q-1,ki-1,ncholQ(Q),&
                h5id_orbs%maxnorb*h5id_orbs%maxnorb,Chol(1,1,1),error)
      if(error .ne. 0) &
        call errore('pw2posthf',' Error: esh5_posthf_read_cholesky.',1)
    
      do kl=1,nksym
      
       kj = QKtoK2(Q,kl) 
       kb = kl  

       ! Gp(down) = Gpl 
       call zgemm('N','N',dfft%nnr,nmax_DM,nmax_DM,(1.d0,0.d0),Psi(1,1,kb),dfft%nnr,  &
                          DM(1,1,kb,nspin),nmax_DM,(0.d0,0.d0),Gp(1,1,2),dfft%nnr)      

       call esh5_posthf_read_cholesky(h5id_hamil%id,Q-1,kl-1,ncholQ(Q),&
                h5id_orbs%maxnorb*h5id_orbs%maxnorb,Chol(1,1,2),error)
       if(error .ne. 0) &
         call errore('pw2posthf',' Error: esh5_posthf_read_cholesky.',1)

       do n=nc0,nc1 

        ! v(p,1) = sum_i,k Lik_n * Gp(up) * conjg(Psi(p,i))  
        ! v(p,2) = sum_l,j conjg(Llj_n) * Gp(down) * conjg(Psi(p,j))  
        ! KS_ur(p) = v(p,1)*v(p,2)
        ! use gemm later on
        v(:,:)=(0.d0,0.d0)
        do i=1,h5id_orbs%norbK(ki) 
         iik = (i-1)*h5id_orbs%norbK(kk)
         do k=1,nmax_DM 
          iik=iik+1
          do p=1,dfft%nnr
            v(p,1) = v(p,1) + Chol(n,iik,1)*Gp(p,k,1)*conjg(Psi(p,i,ki))  
          enddo
         enddo
        enddo 
        do l=1,nmax_DM
         ilj = (l-1)*h5id_orbs%norbK(kj)
         do j=1,h5id_orbs%norbK(kj)
          ilj=ilj+1
          do p=1,dfft%nnr
            v(p,2) = v(p,2) + conjg(Chol(n,ilj,2))*Gp(p,l,2)*conjg(Psi(p,j,kj))
          enddo
         enddo
        enddo
        KS_ur(1:dfft%nnr) = KS_ur(1:dfft%nnr) + v(1:dfft%nnr,1) * v(1:dfft%nnr,2)
 
       enddo !n

      enddo !kl
     enddo !ki
    enddo !Q

    KS_n2(1:dfft%nnr) = 2.d0 * n1(1:dfft%nnr,1) * n1(1:dfft%nnr,2)
    if( nproc_image > 1 ) call mp_sum ( KS_ur, intra_image_comm ) 
    KS_ur(1:dfft%nnr) = sqrtpi * 0.5d0 * 2.d0 * KS_ur(1:dfft%nnr) / KS_n2(1:dfft%nnr)
    n1(:,:) = n1(:,:) / omega
    KS_n2(:) = KS_n2(:) / omega / omega

    nel1=0.d0
    nel2=0.d0  
    do i=1,dfft%nnr
      nel1 = nel1 + n1(i,1) + n1(i,2) 
    enddo
    if(ionode) write(*,*) 'nel:',nel1*omega/(dfft%nr1 * dfft%nr2 * dfft%nr3)/(nksym*1.d0) 

    bscorr=0.d0
    ecloc=0.d0
    do i=1,dfft%nnr
      !
      rhox = dble(n1(i,1) + n1(i,2))
      !
      arhox = ABS( rhox )        
      !
      if ( arhox > vanishing_charge ) then 
        !
        s=0.d0
        zeta = 0.d0
        if(lsda) zeta = dble(n1(i,1) - n1(i,2)) / arhox 
        !
        rs = pi34 / arhox**third 
!        call pw (rs, 1, ec, vc(1))
        !
        n2_ = KS_n2(i) / (1.d0 + 2.d0 / (sqrtpi * KS_ur(i))) 
        beta = fact*ec*arhox/n2_
        bscorr = bscorr + arhox * ec / (1.d0 + beta * KS_ur(i)**3.d0 ) 
        ecloc = ecloc + arhox * ec
        !
      endif  
      !  
    enddo 

    if(ionode) then
      write(*,*) ' Basis set correction: ', &
        bscorr*omega/(dfft%nr1 * dfft%nr2 * dfft%nr3)
      write(*,*) ' Ec (lda): ', &
        ecloc*omega/(dfft%nr1 * dfft%nr2 * dfft%nr3)

      write(*,*) 'i,n2,ur:'
      do i=0,dfft%nr3-1
        j = 1 + i + dfft%nr1x * ( i + dfft%nr2x*i )

        !rhox = rho%of_r(j,1) ! + rho_core(i)
        rhox = dble(n1(j,1) + n1(j,2))
        !
        arhox = ABS( rhox )
        !
          s=0.d0
          zeta = 0.d0
          if(lsda) zeta = dble(n1(j,1) - n1(j,2)) / arhox 
          !
          rs = pi34 / arhox**third
!          call pw (rs, 1, ec, vc(1))
          !
          n2_ = KS_n2(j) / (1.d0 + 2.d0 / (sqrtpi * KS_ur(j)))
          beta = fact*ec*arhox/n2_
          bscorr = bscorr + arhox * ec / (1.d0 + beta * KS_ur(j)**3.d0 )
          !
          write(*,*) i,abs(KS_n2(j)),abs(KS_ur(j)),arhox,ec,beta, &
                    dble(rhox * ec / (1.d0 + beta * KS_ur(j)**3.d0 )),rs
          !
      enddo
    endif

    call close_esh5_read(h5id_orbs)

    CALL esh5_posthf_close_file(h5id_hamil%id)

    if( allocated(Chol) ) deallocate(Chol)
    if( allocated(Psi) ) deallocate(Psi)
    if( allocated(KS_n2) ) deallocate(KS_n2)
    if( allocated(KS_ur) ) deallocate(KS_ur)
    if( allocated(Gp) ) deallocate(Gp)
    if( allocated(v) ) deallocate(v)
    if( allocated(n1) ) deallocate(n1)
    if( allocated(buff) ) deallocate(buff)
    IF( ALLOCATED(ncholQ) ) DEALLOCATE (ncholQ)

  END SUBROUTINE calculate_KS_bscorr

  ! specialized routine for the case where SCFOrbMat is a delta function
  ! only makes sense for nspin==1 case
  SUBROUTINE add_FockM_rhf(h5id,h5id_occ,ki,ispin,dfft,noccK,FockM) 
    USE fft_types, ONLY: fft_type_descriptor
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE gvecw, ONLY : ecutwfc
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE paw_variables,           ONLY : okpaw
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
    USE exx, ONLY: ecutfock
    USE gvect,     ONLY : ecutrho
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id, h5id_occ
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: ki,ispin    
    INTEGER, INTENT(IN) :: noccK(:,:)    
    COMPLEX(DP), INTENT(INOUT) :: FockM(:,:)   
    !
    INTEGER :: b,kb,i,j,M,M2,i_beg,i_end,j_beg,j_end,ni,nj
    INTEGER :: k_beg,k_end,nkloc,npw,n
    COMPLEX(DP) :: ctemp
    REAL(DP) :: fac
    !
    INTEGER, ALLOCATABLE :: igki(:)
    COMPLEX(DP), ALLOCATABLE :: psi_i(:)   
    COMPLEX(DP), ALLOCATABLE :: psi_b(:)   
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:,:)     
    COMPLEX(DP), ALLOCATABLE :: Pbb(:)     
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:) 
    TYPE(ke_factorization), ALLOCATABLE :: ke(:) 
    COMPLEX(DP),ALLOCATABLE :: paw_one_center(:,:,:)
    TYPE(bec_type),ALLOCATABLE :: becpsii
    TYPE(bec_type),ALLOCATABLE :: becpsib
    TYPE(bec_type),ALLOCATABLE :: becpsij(:)
    !

    call start_clock('add_FockM')

    if(nspin == 1 ) then
      fac = 2.d0
    else
      call errore('add_FockM_rhf','nspin>1',1)
    endif

    if(okvan.or.okpaw) &
      call errore('add_FockM','No PAW/USPP yet',1)

    !  Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    ! Partition M*M 
    M = h5id%norbK(ki)
    call find_2d_partition(M,me_pool+1,nproc_pool,i_beg,i_end,j_beg,j_end)
    ni = i_end-i_beg+1
    nj = j_end-j_beg+1

    allocate( Kib(dfft%ngm,ni+nj,1) ) 
    allocate( psi_b(dfft%nnr), psi_i(dfft%nnr), Pbb(dfft%nnr)) 
    allocate( vCoulG(dfft%ngm), igki(npwx) )

    Pbb(:) = (0.d0,0.d0)  
    CALL gk_sort (xksym (1:3, ki), ngm, g, ecutwfc / tpiba2, &
             npw, igki(1), g2kin)

    do kb=k_beg,k_end

      ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
      CALL g2_convolution(dfft%ngm, g, xksym (1:3, kb), xksym (1:3, ki), vCoulG)
      ! 
      CALL gk_sort (xksym (1:3, kb), ngm, g, ecutwfc / tpiba2, &
              npw, igksym(1), g2kin)
      !
      do b=1,noccK(kb,1)
        ! read psi_b 
        call start_clock('diag_io')
        psi_b(:) = (0.d0,0.d0)
        call get_psi_esh5(h5id_occ, b,kb,1,psic)
        psi_b(dfft%nl(igksym(1:npw)))=psic(1:npw)
        if(gamma_only) &
          psi_b(dfft%nlm(igksym(1:npw)))=CONJG(psic(1:npw))
        call stop_clock('diag_io')
        call start_clock('diag_invfft')
        CALL invfft ('Wave', psi_b(:), dfft)
        call start_clock('diag_invfft')

        Pbb(1:dfft%nnr) = Pbb(1:dfft%nnr) + CONJG(psi_b(1:dfft%nnr)) * psi_b(1:dfft%nnr)  

        do i=i_beg,i_end
          !
          ! read psi_i 
          ! can eliminate this read/invfft if you are willing to save the bands
          !  
          call start_clock('diag_io')
          psi_i(:) = (0.d0,0.d0)
          call get_psi_esh5(h5id, i,ki,1,psic)
          psi_i(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
          if(gamma_only) &
            psi_i(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
          call stop_clock('diag_io')
          call start_clock('diag_invfft')
          CALL invfft ('Wave', psi_i(:), dfft)
          call stop_clock('diag_invfft')
          !
          psic(1:dfft%nnr) = CONJG(psi_i(1:dfft%nnr)) * psi_b(1:dfft%nnr) / omega
          !
          ! fwfft orbital pairs to G
          call start_clock('diag_fwfft')
          CALL fwfft ('Rho', psic(:), dfft)
          call stop_clock('diag_fwfft')  
!          ! multiply by FFT[ 1/r ]
          Kib(1:dfft%ngm,i-i_beg+1,1) = psic(dfft%nl(1:dfft%ngm)) * vCoulG(1:dfft%ngm)
! incomplete in gamma_only!!!
        enddo

        do j=j_beg,j_end
          !
          ! read psi_i 
          ! can eliminate this read/invfft if you are willing to save the bands
          !  
          call start_clock('diag_io')
          psi_i(:) = (0.d0,0.d0)
          call get_psi_esh5(h5id, j,ki,1,psic)
          psi_i(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
          if(gamma_only) &
            psi_i(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
          call stop_clock('diag_io')
          call start_clock('diag_invfft')
          CALL invfft ('Wave', psi_i(:), dfft)
          call stop_clock('diag_invfft')
          !
          psic(1:dfft%nnr) = CONJG(psi_i(1:dfft%nnr)) * psi_b(1:dfft%nnr) 
          !
          ! fwfft orbital pairs to G
          call start_clock('diag_fwfft')
          CALL fwfft ('Rho', psic(:), dfft)
          call stop_clock('diag_fwfft')
!          ! multiply by FFT[ 1/r ]
          Kib(1:dfft%ngm,j-j_beg+1+ni,1) = CONJG(psic(dfft%nl(1:dfft%ngm))) 
! incomplete in gamma_only!!!
        enddo

        do i=i_beg,i_end
          do j=j_beg,j_end

            ! exchage contribution
            ctemp=(0.d0,0.d0)
            do n=1,dfft%ngm
              ctemp = ctemp + Kib(n,i-i_beg+1,1) * Kib(n,j-j_beg+1+ni,1)
! incomplete in gamma_only!!!
            enddo 
            FockM(i,j) = FockM(i,j) - e2Ha*ctemp/(nksym)

          enddo ! j
        enddo ! j

      enddo  ! b

    enddo  ! kb

    ! add coulomb contribution here
    ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
    CALL g2_convolution(dfft%ngm, g, xksym (1:3, ki), xksym (1:3, ki), vCoulG)

    ! fwfft orbital pairs to G
    psic(1:dfft%nnr) = Pbb(1:dfft%nnr)
    CALL fwfft ('Rho', psic(:), dfft)
    Pbb(1:dfft%ngm) = vCoulG(1:dfft%ngm) * CONJG(psic(dfft%nl(1:dfft%ngm)))
    !
    do i=i_beg,i_end
      !
      call start_clock('diag_io')
      psi_i(:) = (0.d0,0.d0)
      call get_psi_esh5(h5id, i,ki,1,psic)
      psi_i(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
      if(gamma_only) &
        psi_i(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
      call stop_clock('diag_io')
      call start_clock('diag_invfft')
      CALL invfft ('Wave', psi_i(:), dfft)
      call stop_clock('diag_invfft')
      !
      do j=j_beg,j_end
        !
        call start_clock('diag_io')
        psi_b(:) = (0.d0,0.d0)
        call get_psi_esh5(h5id, j,ki,1,psic)
        psi_b(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
        if(gamma_only) &
          psi_b(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
        call stop_clock('diag_io')
        call start_clock('diag_invfft')
        CALL invfft ('Wave', psi_b(:), dfft)
        call stop_clock('diag_invfft')
        !
        psic(1:dfft%nnr) = CONJG(psi_i(1:dfft%nnr)) * psi_b(1:dfft%nnr) / omega

        call start_clock('diag_fwfft')
        ! fwfft orbital pairs to G 
        CALL fwfft ('Rho', psic(:), dfft)
        call stop_clock('diag_fwfft')

        ! Coulomb contribution
        ctemp=(0.d0,0.d0)
        do n=1,dfft%ngm
          ctemp = ctemp + psic(dfft%nl(n)) * Pbb(n)
! incomplete in gamma_only!!!
        enddo
        FockM(i,j) = FockM(i,j) + fac * e2Ha * ctemp / nksym
        !
      enddo ! j
      !
    enddo ! i 

    if(nproc_image>1) call mp_sum( FockM, intra_image_comm )

    if(allocated(igki)) deallocate(igki)    
    if(allocated(psi_i)) deallocate(psi_i)    
    if(allocated(psi_b)) deallocate(psi_b)    
    if(allocated(Kib)) deallocate(Kib)    
    if(allocated(Pbb)) deallocate(Pbb)    
    if(allocated(vCoulG)) deallocate(vCoulG)    
    call stop_clock('add_FockM')

    if(ionode) then
      write(*,*) 'Timers: '
      CALL print_clock ( 'diag_io' )
      CALL print_clock ( 'diag_invfft' )
      CALL print_clock ( 'diag_fwfft' )
      CALL print_clock ( 'add_FockM' )
    endif

  END SUBROUTINE add_FockM_rhf

  SUBROUTINE add_FockM(h5id,h5id_occ,ki,ispin,dfft,noccK,FockM) 
    USE fft_types, ONLY: fft_type_descriptor
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE gvecw, ONLY : ecutwfc
    USE uspp,                    ONLY : okvan
    USE paw_variables,           ONLY : okpaw
    USE exx, ONLY: ecutfock
    USE gvect,     ONLY : ecutrho
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id, h5id_occ
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: ki,ispin    
    INTEGER, INTENT(IN) :: noccK(:,:)    
    COMPLEX(DP), INTENT(INOUT) :: FockM(:,:)   
    !
    INTEGER :: b,kb,i,j,M,M2,i_beg,i_end,j_beg,j_end,ni,nj
    INTEGER :: k_beg,k_end,nkloc,npw,n,p,maxb,b1,b2
    INTEGER :: nOrbM1,nOrbM2,nk,ns,is 
    COMPLEX(DP) :: ctemp
    REAL(DP) :: fac
    !
    INTEGER, ALLOCATABLE :: igki(:)
    COMPLEX(DP), ALLOCATABLE :: Psi(:,:)   
    COMPLEX(DP), ALLOCATABLE :: psi_b(:)   
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:)     
    COMPLEX(DP), ALLOCATABLE :: Pbb(:)     
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:) 
    COMPLEX(DP), ALLOCATABLE :: Fm(:,:)
    COMPLEX(DP), ALLOCATABLE :: Psi_d(:,:)
    !
    INTEGER :: nfft, ierror
    COMPLEX(DP), ALLOCATABLE :: SCFOrbMat(:,:)
    COMPLEX(DP), ALLOCATABLE :: psi_scf(:,:)
    !
    call start_clock('add_FockM')
    
    if(nspin == 1 ) then
      fac = 2.d0
    else
      fac = 1.d0
    endif

    ! for mixed_basis=.false., you should use add_Fock__rhf routine 
    call esh5_posthf_read_orbmat_info(h5id_occ%id,'SCFOrbMat',9, &
             nOrbM1,nOrbM2,nk,ns,ierror)
    if(ierror.ne.0) &
      call errore('add-FockM','Error in esh5_posthf_read_orbmat_info',1)
    if((nk.ne.nksym) .or. (ns.ne.nspin) .or. &
        (nOrbM2 < maxval(noccK(:,:))) ) &
      call errore('add-FockM','Inconsistent SCFOrbMat',1)
    maxb = (maxval(noccK(:,:))+nproc_pool-1)/nproc_pool

    if(okvan.or.okpaw) &
      call errore('add_FockM','No PAW/USPP yet',1)

    !  Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    ! Partition M*M 
    M = h5id%norbK(ki)
    call find_2d_partition(M,me_pool+1,nproc_pool,i_beg,i_end,j_beg,j_end)
    ni = i_end-i_beg+1
    nj = j_end-j_beg+1

    allocate( Psi(npwx,min(16,maxval(h5id%norbK(:)))) ) 
    allocate( psi_b(dfft%nnr), Pbb(dfft%nnr)) 
    allocate( vCoulG(dfft%ngm), igki(npwx) )
    allocate( SCFOrbMat(nOrbM1,nOrbM2), psi_scf(npwx,maxb) )
    allocate( Kib(dfft%ngm,ni+nj) )
    allocate( Psi_d(dfft%nnr,ni+nj) )
    allocate( Fm(ni,nj) )
    !

    Fm(:,:) = FockM(i_beg:i_end,j_beg:j_end) ! set Fm to incoming value 
    !

    Pbb(:) = (0.d0,0.d0)  
    CALL gk_sort (xksym (1:3, ki), dfft%ngm, g, ecutwfc / tpiba2, &
             npw, igki(1), g2kin)

    call start_clock('diag_io')
    ! basis orbitals are spin-independent
    do i=1,ni
      call get_psi_esh5(h5id,i+i_beg-1,ki,1,psic)
      psi_b(:)=(0.d0,0.d0)
      psi_b(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
      if(gamma_only) &
        psi_b(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
      Psi_d(:,i) = psi_b(:)  
    enddo
    ! basis orbitals are spin-independent
    do j=1,nj
      call get_psi_esh5(h5id,j+j_beg-1,ki,1,psic)
      psi_b(:)=(0.d0,0.d0)
      psi_b(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
      if(gamma_only) &
        psi_b(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
      Psi_d(:,ni+j) = psi_b(:)  
    enddo
    ! copy wfns to device
    do i=1,ni+nj
      CALL invfft ('Wave', Psi_d(:,i), dfft)
    enddo
    call stop_clock('diag_io')

    do kb=k_beg,k_end

      ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
      CALL g2_convolution(dfft%ngm, g, xksym (1:3, kb), xksym (1:3, ki), vCoulG)
      ! 
      CALL gk_sort (xksym (1:3, kb), dfft%ngm, g, ecutwfc / tpiba2, &
              npw, igksym(1), g2kin)
      !
      ! read MF overlap matrix and generate (hopefully!!!) HF orbitals  
      ! spin dependence of HF state in SCFOrbMat  
      !
      call esh5_posthf_read_orbmat(h5id_occ%id,'SCFOrbMat',9,kb-1,ispin-1, &
               SCFOrbMat,ierror) 
      if(ierror.ne.0) &
        call errore('add-FockM','Error in esh5_posthf_read_orbmat',1)

      !
      ! distributing psi_scf round-robin for simpler implementation  
      ! 
      call start_clock('diag_io')
      ! modify SCFOrbMat matrix
      maxb=0
      do b=1,noccK(kb,ispin)
        b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
        b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
        if(me_pool == b1) then
          maxb=b2
          SCFOrbMat(:,b2) = SCFOrbMat(:,b)
        endif
      enddo

      ! using Psi as a working array 
      ! 
      psi_scf(:,:) = (0.d0,0.d0)
      if(maxb > 0) then
        b1 = min(nOrbM1,h5id%norbK(kb))
        do i=1,b1,size(Psi,2)
          n = min( size(Psi,2), b1-i+1 )
          do j=1,n
            ! basis set is spin independent
            call get_psi_esh5(h5id, j+i-1,kb,1,Psi(:,j))
          enddo  
          call zgemm('N','N',npw,maxb,n,(1.d0,0.d0),Psi(1,1),npwx,  &
                     SCFOrbMat(i,1),nOrbM1,(1.d0,0.d0),psi_scf,npwx)
        enddo
      endif
      call stop_clock('diag_io')
      !
      do b=1,noccK(kb,ispin)
        ! read psi_b 
        psi_b(:) = (0.d0,0.d0)
        b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
        b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
        if(me_pool == b1) psic(1:npw) = psi_scf(1:npw,b2)
        call mp_bcast( psic, b1, intra_pool_comm ) 
        psi_b(dfft%nl(igksym(1:npw)))=psic(1:npw)
        if(gamma_only) &
          psi_b(dfft%nlm(igksym(1:npw)))=CONJG(psic(1:npw))
        call start_clock('diag_invfft')
        CALL invfft ('Wave', psi_b(:), dfft)
        call stop_clock('diag_invfft')

        if(me_pool == b1) &
          Pbb(1:dfft%nnr) = Pbb(1:dfft%nnr) + CONJG(psi_b(1:dfft%nnr)) * psi_b(1:dfft%nnr)  

        do i=1,ni
          !
          do n = 1,dfft%nnr  
            psic(n) = CONJG(Psi_d(n,i)) * &
                                  psi_b(n) / omega
          enddo
          !
          ! fwfft orbital pairs to G
          call start_clock('diag_fwfft')
          CALL fwfft ('Rho', psic, dfft)
          call stop_clock('diag_fwfft')  

!          ! multiply by FFT[ 1/r ]
          do n = 1,dfft%ngm  
            Kib(n,i) = psic(dfft%nl(n)) * vCoulG(n)
          enddo
! incomplete in gamma_only!!!
        enddo

        do j=1,nj
          !
          do n = 1,dfft%nnr
            psic(n) = CONJG(Psi_d(n,j+ni)) * psi_b(n) 
          enddo
          !
          ! fwfft orbital pairs to G
          call start_clock('diag_fwfft')
          CALL fwfft ('Rho', psic, dfft)
          call stop_clock('diag_fwfft')

!          ! multiply by FFT[ 1/r ]
          do n = 1,dfft%ngm
            Kib(n,j+ni) = CONJG(psic(dfft%nl(n))) 
          enddo
! incomplete in gamma_only!!!
        enddo

        ctemp = -e2Ha/(1.d0*nksym)
        call zgemm('T','N',& 
            ni,nj,dfft%ngm,ctemp,Kib,dfft%ngm,Kib(1,ni+1),dfft%ngm,&
            (1.d0,0.d0),Fm(1,1),ni) 

      enddo  ! b

    enddo  ! kb

    ! copy back to host
    FockM(i_beg:i_end,j_beg:j_end) = Fm(:,:)

    ! add other spin contribution to Pbb
    if(nspin > 1) then
      is = 3-ispin
      do kb=k_beg,k_end
        ! 
        CALL gk_sort (xksym (1:3, kb), dfft%ngm, g, ecutwfc / tpiba2, &
                npw, igksym(1), g2kin)
        !
        ! read MF overlap matrix and generate (hopefully!!!) HF orbitals  
        !
        call esh5_posthf_read_orbmat(h5id_occ%id,'SCFOrbMat',9,kb-1,is-1, &
               SCFOrbMat,ierror) 
        if(ierror.ne.0) &
          call errore('add-FockM','Error in esh5_posthf_read_orbmat',1)

        !
        ! distributing psi_scf round-robin for simpler implementation  
        ! 
        call start_clock('diag_io')
        ! modify SCFOrbMat matrix
        maxb=0
        do b=1,noccK(kb,is)
          b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
          b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
          if(me_pool == b1) then
            maxb=b2
            SCFOrbMat(:,b2) = SCFOrbMat(:,b)
          endif
        enddo

        ! in this implementation, h5id_occ has no reason to exist!!!
        ! using Psi as a working array 
        ! 
        psi_scf(:,:) = (0.d0,0.d0)
        if(maxb > 0) then
          b1 = min(nOrbM1,h5id%norbK(kb))
          do i=1,b1,size(Psi,2)
            n = min( size(Psi,2), b1-i+1 )
            do j=1,n
              call get_psi_esh5(h5id, j+i-1,kb,1,Psi(:,j))
            enddo  
            call zgemm('N','N',npw,maxb,n,(1.d0,0.d0),Psi(1,1),npwx,  &
                       SCFOrbMat(i,1),nOrbM1,(1.d0,0.d0),psi_scf,npwx)
          enddo
        endif
        call stop_clock('diag_io')
        !
        do b=1,noccK(kb,is)
          psi_b(:) = (0.d0,0.d0)
          b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
          b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
          if(me_pool == b1) then 
            psi_b(dfft%nl(igksym(1:npw)))=psi_scf(1:npw,b2)
            if(gamma_only) &
              psi_b(dfft%nlm(igksym(1:npw)))=CONJG(psi_scf(1:npw,b2))
            call start_clock('diag_invfft')
            CALL invfft ('Wave', psi_b(:), dfft)
            call stop_clock('diag_invfft')
            Pbb(1:dfft%nnr) = Pbb(1:dfft%nnr) + CONJG(psi_b(1:dfft%nnr)) * psi_b(1:dfft%nnr)  
          endif
        enddo
        !
      enddo
      !
    endif
    ! reduce Pbb within pool
    call mp_sum( Pbb, intra_pool_comm ) 

    ! add coulomb contribution here
    ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
    CALL g2_convolution(dfft%ngm, g, xksym (1:3, ki), xksym (1:3, ki), vCoulG)
    ! fwfft orbital pairs to G
    psic(1:dfft%nnr) = Pbb
    CALL fwfft ('Rho', psic(1:dfft%nnr), dfft)
    do n=1,dfft%ngm    
      psi_b(n) = vCoulG(n) * CONJG(psic(dfft%nl(n)))
    enddo
    !
    do i=1,ni
      !
      do j=1,nj
        !
        do n = 1,dfft%nnr
          psic(n) = CONJG(Psi_d(n,i)) * &
                           Psi_d(n,j+ni) / omega
        enddo
        !
        call start_clock('diag_fwfft')
        ! fwfft orbital pairs to G 
        CALL fwfft ('Rho', psic, dfft)
        call stop_clock('diag_fwfft')

        ! Coulomb contribution
        ctemp=(0.d0,0.d0)
        do n=1,dfft%ngm
          ctemp = ctemp + psic(dfft%nl(n)) * psi_b(n)
! incomplete in gamma_only!!!
        enddo
        FockM(i+i_beg-1,j+j_beg-1) = FockM(i+i_beg-1,j+j_beg-1) + fac * e2Ha * ctemp / nksym
        !
      enddo ! j
      !
    enddo ! i 

    if(nproc_image>1) call mp_sum( FockM, intra_image_comm )

    if(allocated(igki)) deallocate(igki)    
    if(allocated(Psi)) deallocate(Psi)    
    if(allocated(psi_b)) deallocate(psi_b)    
    if(allocated(Pbb)) deallocate(Pbb)    
    if(allocated(vCoulG)) deallocate(vCoulG)    
    !
    if(allocated(Psi_d)) deallocate(Psi_d)
    if(allocated(Kib)) deallocate(Kib)
    if(allocated(Fm)) deallocate(Fm)
    IF( ALLOCATED(SCFOrbMat) ) DEALLOCATE (SCFOrbMat)
    IF( ALLOCATED(psi_scf) ) DEALLOCATE (psi_scf)
    call stop_clock('add_FockM')

    if(ionode) then
      write(*,*) 'Timers: '
      CALL print_clock ( 'diag_io' )
      CALL print_clock ( 'diag_invfft' )
      CALL print_clock ( 'diag_fwfft' )
      CALL print_clock ( 'add_FockM' )
    endif

    !
  END SUBROUTINE add_FockM

#if defined(__CUDA)
  ! you should only calculate upper/lower triangular part and hermitize 
  ! adds the 2-electron interaction part of the Fock matric to FockM
  ! As usual, this assumes a spin-independent basis.
  ! The spin dependence is encapsulated on the HF matrix found in
  ! h5id_occ::SCFOrbMat
  SUBROUTINE add_FockM_gpu(h5id,h5id_occ,ki,ispin,dfft,noccK,FockM) 
    USE fft_types, ONLY: fft_type_descriptor
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE gvecw, ONLY : ecutwfc
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE paw_variables,           ONLY : okpaw
    USE us_exx,         ONLY : qvan_init, qvan_clean, addusxx_r, addusxx_g, &
                                newdxx_g, newdxx_r
    USE realus, ONLY: tabxx,tabp,generate_qpointlist,qpointlist
    USE exx, ONLY: ecutfock
    USE gvect,     ONLY : ecutrho
    USE uspp_param,         ONLY : nh
    USE control_flags, ONLY : tqr
    USE paw_exx, ONLY : PAW_xx_energy,PAW_init_fock_kernel,PAW_clean_fock_kernel
    !
    USE cublas
    USE cudafor
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id, h5id_occ
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: ki,ispin    
    INTEGER, INTENT(IN) :: noccK(:,:)    
    COMPLEX(DP), INTENT(INOUT) :: FockM(:,:)   
    !
    INTEGER :: b,kb,i,j,M,M2,i_beg,i_end,j_beg,j_end,ni,nj
    INTEGER :: k_beg,k_end,nkloc,npw,n,p,maxb,b1,b2
    INTEGER :: nOrbM1,nOrbM2,nk,ns,is 
    COMPLEX(DP) :: ctemp
    REAL(DP) :: fac
    !
    INTEGER, ALLOCATABLE :: igki(:)
    COMPLEX(DP), ALLOCATABLE :: Psi(:,:)   
    COMPLEX(DP), ALLOCATABLE :: psi_b(:)   
    COMPLEX(DP), ALLOCATABLE :: Kib(:,:)     
    COMPLEX(DP), ALLOCATABLE :: Pbb(:)     
    COMPLEX(DP), ALLOCATABLE :: vCoulG(:) 
    TYPE(ke_factorization), ALLOCATABLE :: ke(:) 
    COMPLEX(DP),ALLOCATABLE :: paw_one_center(:,:,:)
    TYPE(bec_type),ALLOCATABLE :: becpsii
    TYPE(bec_type),ALLOCATABLE :: becpsib
    TYPE(bec_type),ALLOCATABLE :: becpsij(:)
    !
    INTEGER :: nfft, many_fft, ierror
    COMPLEX(DP), ALLOCATABLE :: Psi_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: psib_d(:)
    COMPLEX(DP), ALLOCATABLE :: psic_d(:)
    COMPLEX(DP), ALLOCATABLE :: Kib_d(:,:)
    COMPLEX(DP), ALLOCATABLE :: Fm_d(:,:)
    INTEGER, ALLOCATABLE :: nl_d(:)
    COMPLEX(DP), ALLOCATABLE :: vCoulG_d(:) 
    COMPLEX(DP), ALLOCATABLE :: SCFOrbMat(:,:)
    COMPLEX(DP), ALLOCATABLE :: psi_scf(:,:)
    attributes(DEVICE) :: Psi_d,psic_d,psib_d,Kib_d,Fm_d,nl_d,vCoulG_d
    !
    type(cublasHandle) :: handle_cublas
    integer :: istat
    istat = cublasCreate(handle_cublas)
    if (istat .ne. CUBLAS_STATUS_SUCCESS) print *,istat
    !
    call start_clock('add_FockM')
    
    many_fft=16
    if(nspin == 1 ) then
      fac = 2.d0
    else
      fac = 1.d0
    endif

    ! for mixed_basis=.false., this is just a delta function, so avoid
    ! all of this with a flag? or maybe write in the hdf5 the type of basis?
    ! read dimensions of SCFOrbMat
    call esh5_posthf_read_orbmat_info(h5id_occ%id,'SCFOrbMat',9, &
             nOrbM1,nOrbM2,nk,ns,ierror)
    if(ierror.ne.0) &
      call errore('add-FockM','Error in esh5_posthf_read_orbmat_info',1)
    if((nk.ne.nksym) .or. (ns.ne.nspin) .or. &
        (nOrbM2 < maxval(noccK(:,:))) ) &
      call errore('add-FockM','Inconsistent SCFOrbMat',1)
    maxb = (maxval(noccK(:,:))+nproc_pool-1)/nproc_pool

    if(okvan.or.okpaw) &
      call errore('add_FockM','No PAW/USPP yet',1)

    !  Partition k-points among MPI tasks
    call fair_divide(k_beg,k_end,my_pool_id+1,npool,nksym)
    nkloc   = k_end - k_beg + 1

    ! Partition M*M 
    M = h5id%norbK(ki)
    call find_2d_partition(M,me_pool+1,nproc_pool,i_beg,i_end,j_beg,j_end)
    ni = i_end-i_beg+1
    nj = j_end-j_beg+1

    allocate( Psi(npwx,min(16,maxval(h5id%norbK(:)))) ) 
    allocate( psi_b(dfft%nnr), Pbb(dfft%nnr)) 
    allocate( vCoulG(dfft%ngm), igki(npwx) )
    allocate( SCFOrbMat(nOrbM1,nOrbM2), psi_scf(npwx,maxb) )
    !
    allocate( Kib_d(dfft%ngm,ni+nj) )
    allocate( Psi_d(dfft%nnr,ni+nj) )
    allocate( Fm_d(ni,nj) )
    allocate( psib_d(dfft%nnr), psic_d(dfft%nnr*many_fft))
    allocate( vCoulG_d(dfft%ngm), nl_d(dfft%ngm) )

    nl_d = dfft%nl
    Fm_d(:,:) = FockM(i_beg:i_end,j_beg:j_end) ! set Fm to incoming value 
    !

    Pbb(:) = (0.d0,0.d0)  
    CALL gk_sort (xksym (1:3, ki), dfft%ngm, g, ecutwfc / tpiba2, &
             npw, igki(1), g2kin)

    call start_clock('diag_io')
    ! basis orbitals are spin-independent
    do i=1,ni
      call get_psi_esh5(h5id,i+i_beg-1,ki,1,psic)
      psi_b(:)=(0.d0,0.d0)
      psi_b(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
      if(gamma_only) &
        psi_b(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
      Psi_d(:,i) = psi_b(:)  
    enddo
    ! basis orbitals are spin-independent
    do j=1,nj
      call get_psi_esh5(h5id,j+j_beg-1,ki,1,psic)
      psi_b(:)=(0.d0,0.d0)
      psi_b(dfft%nl(igki(1:ngksym(ki))))=psic(1:ngksym(ki))
      if(gamma_only) &
        psi_b(dfft%nlm(igki(1:ngksym(ki))))=CONJG(psic(1:ngksym(ki)))
      Psi_d(:,ni+j) = psi_b(:)  
    enddo
    ! copy wfns to device
!    Psi_d = Psi
    do i=1,ni+nj,many_fft
      nfft = min( many_fft, ni+nj-i+1 )  
      CALL invfft ('Wave', Psi_d(:,i), dfft, howmany=nfft)
    enddo
    call stop_clock('diag_io')

    do kb=k_beg,k_end

      ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
      CALL g2_convolution(dfft%ngm, g, xksym (1:3, kb), xksym (1:3, ki), vCoulG)
      vCoulG_d = vCoulG
      ! 
      CALL gk_sort (xksym (1:3, kb), dfft%ngm, g, ecutwfc / tpiba2, &
              npw, igksym(1), g2kin)
      !
      ! read MF overlap matrix and generate (hopefully!!!) HF orbitals  
      ! spin dependence of HF state in SCFOrbMat  
      !
      call esh5_posthf_read_orbmat(h5id_occ%id,'SCFOrbMat',9,kb-1,ispin-1, &
               SCFOrbMat,ierror) 
      if(ierror.ne.0) &
        call errore('add-FockM','Error in esh5_posthf_read_orbmat',1)

      !
      ! distributing psi_scf round-robin for simpler implementation  
      ! 
      call start_clock('diag_io')
      ! modify SCFOrbMat matrix
      maxb=0
      do b=1,noccK(kb,ispin)
        b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
        b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
        if(me_pool == b1) then
          maxb=b2
          SCFOrbMat(:,b2) = SCFOrbMat(:,b)
        endif
      enddo

      ! in this implementation, h5id_occ has no reason to exist!!!
      ! using Psi as a working array 
      ! 
      psi_scf(:,:) = (0.d0,0.d0)
      if(maxb > 0) then
        b1 = min(nOrbM1,h5id%norbK(kb))
        do i=1,b1,size(Psi,2)
          n = min( size(Psi,2), b1-i+1 )
          do j=1,n
            ! basis set is spin independent
            call get_psi_esh5(h5id, j+i-1,kb,1,Psi(:,j))
          enddo  
          call zgemm('N','N',npw,maxb,n,(1.d0,0.d0),Psi(1,1),npwx,  &
                     SCFOrbMat(i,1),nOrbM1,(1.d0,0.d0),psi_scf,npwx)
        enddo
      endif
      call stop_clock('diag_io')
      !
      do b=1,noccK(kb,ispin)
        ! read psi_b 
!        call start_clock('diag_io')
        psi_b(:) = (0.d0,0.d0)
!        call get_psi_esh5(h5id_occ, b,kb,1,psic)
        b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
        b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
        if(me_pool == b1) psic(1:npw) = psi_scf(1:npw,b2)
        call mp_bcast( psic, b1, intra_pool_comm ) 
        psi_b(dfft%nl(igksym(1:npw)))=psic(1:npw)
        if(gamma_only) &
          psi_b(dfft%nlm(igksym(1:npw)))=CONJG(psic(1:npw))
!        call stop_clock('diag_io')
        call start_clock('diag_invfft')
        CALL invfft ('Wave', psi_b(:), dfft)
        call stop_clock('diag_invfft')

        if(me_pool == b1) &
          Pbb(1:dfft%nnr) = Pbb(1:dfft%nnr) + CONJG(psi_b(1:dfft%nnr)) * psi_b(1:dfft%nnr)  

        psib_d = psi_b

        do i=1,ni,many_fft
          !
          nfft = min( many_fft, ni-i+1 )
          !
!$cuf kernel do(2)
          do p = 0,nfft-1  
            do n = 1,dfft%nnr  
              psic_d(n+p*dfft%nnr) = CONJG(Psi_d(n,i+p)) * &
                                        psib_d(n) / omega
            enddo
          enddo
          !
          ! fwfft orbital pairs to G
          call start_clock('diag_fwfft')
          CALL fwfft ('Rho', psic_d, dfft, howmany=nfft)
          call stop_clock('diag_fwfft')  

!          ! multiply by FFT[ 1/r ]
!$cuf kernel do(2)
          do p = 0,nfft-1  
            do n = 1,dfft%ngm  
              Kib_d(n,i+p) = psic_d(nl_d(n)+p*dfft%nnr) * vCoulG_d(n)
            enddo
          enddo
! incomplete in gamma_only!!!
        enddo

        do j=1,nj,many_fft
          !
          nfft = min( many_fft, nj-j+1 )
          !
!$cuf kernel do(2)
          do p = 0,nfft-1
            do n = 1,dfft%nnr
              psic_d(n+p*dfft%nnr) = CONJG(Psi_d(n,j+ni+p)) * psib_d(n) 
            enddo
          enddo
          !
          ! fwfft orbital pairs to G
          call start_clock('diag_fwfft')
          CALL fwfft ('Rho', psic_d, dfft, howmany=nfft)
          call stop_clock('diag_fwfft')

!          ! multiply by FFT[ 1/r ]
!$cuf kernel do(2)
          do p = 0,nfft-1
            do n = 1,dfft%ngm
              Kib_d(n,j+ni+p) = CONJG(psic_d(nl_d(n)+p*dfft%nnr)) 
            enddo
          enddo
! incomplete in gamma_only!!!
        enddo

        ctemp = -e2Ha/(1.d0*nksym)
        istat = cublasZgemm_v2(handle_cublas,CUBLAS_OP_T,CUBLAS_OP_N,& 
            ni,nj,dfft%ngm,ctemp,Kib_d,dfft%ngm,Kib_d(1,ni+1),dfft%ngm,&
            (1.d0,0.d0),Fm_d(1,1),ni) 
        if (istat .ne.CUBLAS_STATUS_SUCCESS) print *,istat

      enddo  ! b

    enddo  ! kb

    ! copy back to host
    FockM(i_beg:i_end,j_beg:j_end) = Fm_d(:,:)

    ! add other spin contribution to Pbb
    if(nspin > 1) then
      is = 3-ispin
      do kb=k_beg,k_end
        ! 
        CALL gk_sort (xksym (1:3, kb), dfft%ngm, g, ecutwfc / tpiba2, &
                npw, igksym(1), g2kin)
        !
        ! read MF overlap matrix and generate (hopefully!!!) HF orbitals  
        !
        call esh5_posthf_read_orbmat(h5id_occ%id,'SCFOrbMat',9,kb-1,is-1, &
               SCFOrbMat,ierror) 
        if(ierror.ne.0) &
          call errore('add-FockM','Error in esh5_posthf_read_orbmat',1)

        !
        ! distributing psi_scf round-robin for simpler implementation  
        ! 
        call start_clock('diag_io')
        ! modify SCFOrbMat matrix
        maxb=0
        do b=1,noccK(kb,is)
          b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
          b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
          if(me_pool == b1) then
            maxb=b2
            SCFOrbMat(:,b2) = SCFOrbMat(:,b)
          endif
        enddo

        ! in this implementation, h5id_occ has no reason to exist!!!
        ! using Psi as a working array 
        ! 
        psi_scf(:,:) = (0.d0,0.d0)
        if(maxb > 0) then
          b1 = min(nOrbM1,h5id%norbK(kb))
          do i=1,b1,size(Psi,2)
            n = min( size(Psi,2), b1-i+1 )
            do j=1,n
              call get_psi_esh5(h5id, j+i-1,kb,1,Psi(:,j))
            enddo  
            call zgemm('N','N',npw,maxb,n,(1.d0,0.d0),Psi(1,1),npwx,  &
                       SCFOrbMat(i,1),nOrbM1,(1.d0,0.d0),psi_scf,npwx)
          enddo
        endif
        call stop_clock('diag_io')
        !
        do b=1,noccK(kb,is)
          psi_b(:) = (0.d0,0.d0)
          b1 = mod(b-1,nproc_pool)   ! index of mpi task with state
          b2 = (b-1)/nproc_pool+1      ! location of state in matrix 
          if(me_pool == b1) then 
            psi_b(dfft%nl(igksym(1:npw)))=psi_scf(1:npw,b2)
            if(gamma_only) &
              psi_b(dfft%nlm(igksym(1:npw)))=CONJG(psi_scf(1:npw,b2))
            call start_clock('diag_invfft')
            CALL invfft ('Wave', psi_b(:), dfft)
            call stop_clock('diag_invfft')
            Pbb(1:dfft%nnr) = Pbb(1:dfft%nnr) + CONJG(psi_b(1:dfft%nnr)) * psi_b(1:dfft%nnr)  
          endif
        enddo
        !
      enddo
      !
    endif
    ! reduce Pbb within pool
    call mp_sum( Pbb, intra_pool_comm ) 

    ! add coulomb contribution here
    ! vCoulG( iG ) = | G - Q |^{-2} = | G + k(kb) - k(ka)  |^{-2}
    CALL g2_convolution(dfft%ngm, g, xksym (1:3, ki), xksym (1:3, ki), vCoulG)
    vCoulG_d = vCoulG
    
    ! fwfft orbital pairs to G
    psic_d(1:dfft%nnr) = Pbb
    CALL fwfft ('Rho', psic_d(1:dfft%nnr), dfft)
!$cuf kernel do
    do n=1,dfft%ngm    
      psib_d(n) = vCoulG_d(n) * CONJG(psic_d(nl_d(n)))
    enddo
    !
    do i=1,ni
      !
      do j=1,nj,many_fft
        !
        nfft = min( many_fft, nj-j+1 )
        !
!$cuf kernel do(2)
        do p = 0,nfft-1
          do n = 1,dfft%nnr
            psic_d(n+p*dfft%nnr) = CONJG(Psi_d(n,i)) * &
                           Psi_d(n,j+ni+p) / omega
          enddo
        enddo
        !
        call start_clock('diag_fwfft')
        ! fwfft orbital pairs to G 
        CALL fwfft ('Rho', psic_d, dfft, howmany=nfft)
        call stop_clock('diag_fwfft')

        ! Coulomb contribution
        do p = 0,nfft-1
          ctemp=(0.d0,0.d0)
!$cuf kernel do
          do n=1,dfft%ngm
            ctemp = ctemp + psic_d(nl_d(n)+p*dfft%nnr) * psib_d(n)
! incomplete in gamma_only!!!
          enddo
          FockM(i+i_beg-1,j+j_beg-1+p) = FockM(i+i_beg-1,j+j_beg-1+p) + fac * e2Ha * ctemp / nksym
        enddo
        !
      enddo ! j
      !
    enddo ! i 

    if(nproc_image>1) call mp_sum( FockM, intra_image_comm )

    if(allocated(igki)) deallocate(igki)    
    if(allocated(Psi)) deallocate(Psi)    
    if(allocated(psi_b)) deallocate(psi_b)    
    if(allocated(Pbb)) deallocate(Pbb)    
    if(allocated(vCoulG)) deallocate(vCoulG)    
    !
    if(allocated(Psi_d)) deallocate(Psi_d)
    if(allocated(Kib_d)) deallocate(Kib_d)
    if(allocated(psib_d)) deallocate(psib_d)
    if(allocated(psic_d)) deallocate(psic_d)
    if(allocated(Fm_d)) deallocate(Fm_d)
    IF( ALLOCATED(nl_d) ) DEALLOCATE(nl_d)
    IF( ALLOCATED(vCoulG_d) ) DEALLOCATE (vCoulG_d)
    IF( ALLOCATED(SCFOrbMat) ) DEALLOCATE (SCFOrbMat)
    IF( ALLOCATED(psi_scf) ) DEALLOCATE (psi_scf)
    call stop_clock('add_FockM')

    if(ionode) then
      write(*,*) 'Timers: '
      CALL print_clock ( 'diag_io' )
      CALL print_clock ( 'diag_invfft' )
      CALL print_clock ( 'diag_fwfft' )
      CALL print_clock ( 'add_FockM' )
    endif

    istat = cublasDestroy(handle_cublas)
    if (istat .ne.CUBLAS_STATUS_SUCCESS) print *,istat
    !
  END SUBROUTINE add_FockM_gpu
#endif

END MODULE twobody_hamiltonian 
