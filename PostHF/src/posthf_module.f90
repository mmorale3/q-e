!
! Written by Miguel A. Morales, LLNL, 2020
!
! Store common variables here 
!
MODULE posthf_mod
  !  
  USE kinds, ONLY: DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info
  !  
  IMPLICIT NONE
  !  
  TYPE ke_factorization
    COMPLEX(DP), ALLOCATABLE :: L(:,:)
  END TYPE ke_factorization
  !
  INTEGER :: nksym, nkfull, numspin, norb
  INTEGER :: nelec_tot, nup, ndown, nQuniq
  REAL(DP), ALLOCATABLE :: xksym(:,:)
  INTEGER, ALLOCATABLE :: ngksym(:) 
  INTEGER, ALLOCATABLE :: igksym(:) 
  REAL(DP), ALLOCATABLE :: xkfull(:,:)
  INTEGER, ALLOCATABLE :: QKtoK2(:,:)
  REAL(DP), ALLOCATABLE :: Qpts(:,:)
  INTEGER, ALLOCATABLE :: kminus(:)
  REAL(DP), ALLOCATABLE :: wQ(:)
  INTEGER, ALLOCATABLE :: xQ(:)
  REAL(DP) :: e2Ha, e0, efc
  INTEGER :: nmax_DM
  LOGICAL :: use_symm 
  LOGICAL :: use_regularization, use_coulomb_vcut_ws, use_coulomb_vcut_spheric  
  REAL(DP) :: ecutvcut
  TYPE(vcut_type) :: vcut
  REAL(DP) :: exxdiv = 0._dp
  COMPLEX(DP), ALLOCATABLE :: DM(:,:,:,:)  ! density matrix consistent with
                                           ! single determinant trial wave function 
  COMPLEX(DP), ALLOCATABLE :: DM_mf(:,:,:,:)  ! mean field density matrix 
  !
  CONTAINS
  !
  ! right now only works without symmetry
  ! so xk and xksym are identical. 
  ! eventually will implement symmetry and these will be different
  !    
  SUBROUTINE init_posthf(norb_,symm,exxdiv_treat)
    !  
    USE constants, ONLY: e2
    USE control_flags, ONLY: gamma_only
    USE vlocal, ONLY : strf
    USE io_global, ONLY : stdout
    USE ions_base, ONLY : nat, nsp, ityp, tau, zv
    USE cell_base, ONLY: omega, alat, tpiba, tpiba2, at, bg
    USE wvfct, ONLY: npwx, nbnd, g2kin
    USE klist , ONLY: nks, nelec, nelup, neldw, xk, nkstot
    USE gvect, ONLY: ngm, ngm_g, g, ig_l2g, gstart, gg, gcutm
    USE gvecw, ONLY : ecutwfc
    USE io_global, ONLY: ionode,  ionode_id
    USE mp,           ONLY: mp_bcast, mp_barrier
    USE mp_images, ONLY: intra_image_comm, me_image, root_image
    USE symm_base, ONLY: set_sym_bl, find_sym
    USE lsda_mod, ONLY: lsda, nspin
    USE noncollin_module,     ONLY : noncolin, npol
    USE exx_base, ONLY: nq1,nq2,nq3, &  
            x_gamma_extrapolation, exx_divergence
    USE moire, ONLY: lmoire
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: norb_
    LOGICAL, INTENT(IN), OPTIONAL :: symm
    CHARACTER(len=*) , INTENT(IN), OPTIONAL :: exxdiv_treat
    !
    INTEGER :: ik,j
    ! ... external funcitons
    REAL(DP), EXTERNAL :: ewald
    REAL(DP), ALLOCATABLE :: m_loc(:,:)
    LOGICAL :: magnetic, gate
    REAL(DP) :: atws(3,3)
    !
    if(gamma_only .and. noncolin) &
      call errore('init_posthf','Error: gamma_only and noncolin not yet implemented. ',1)
    !
    use_symm = .true.
    if(present(symm)) use_symm=symm
    !
    use_regularization = .false.
    use_coulomb_vcut_spheric = .false.
    use_coulomb_vcut_ws = .false.
    x_gamma_extrapolation = .false.
    if(present(exxdiv_treat)) then
      SELECT CASE ( TRIM(exxdiv_treat) )
      CASE ( "gygi-baldereschi", "gygi-bald", "g-b", "gb" )
        !
        use_regularization = .TRUE.
        !
      CASE ( "vcut_ws" )
        !
        use_regularization = .TRUE.
        use_coulomb_vcut_ws = .TRUE.
        !
      CASE ( "vcut_spherical" )
        !
        use_regularization = .TRUE.
        use_coulomb_vcut_spheric = .TRUE.
        !
      CASE DEFAULT
        !
      END SELECT
      !
    endif
    !
    if(ionode) call print_freemem()
    e2Ha = 1.d0/e2
    norb = norb_
    if(norb .lt. 0) norb=nbnd
    if(norb .gt. nbnd) then
      write(*,*) 'WARNING: Requesting too many orbitals: ',norb
      write(*,*) 'Setting number_of_orbitals to number of orbitals in file:',nbnd
      norb=nbnd
    endif

    numspin = nspin
    if( lsda ) then
      nksym = nkstot/2
      nkfull = nksym   ! modify when symmetry is enabled
      if(noncolin) call errore('pw2qmcpack','lsda and nonclin (??)',1)
    else
      nksym = nkstot
      nkfull = nksym   ! modify when symmetry is enabled
      if(noncolin) numspin=1  ! treating as a single Slater Matrix with 2*norb orbitals
    endif
    !
    ! get the number of electrons
    nelec_tot= NINT(nelec)
    nup=NINT(nelup)
    ndown=NINT(neldw)
 
    if(nup .eq. 0) then
      ndown=nelec_tot/2
      nup=nelec_tot-ndown
    endif
    write(*,*) 'nelec,nup,ndown:',nelec_tot,nup,ndown
    !
    efc = 0.d0
    if (lmoire) then
    e0 = 0.d0
    else
    e0 = e2Ha * ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
    endif
    !
!    allocate( norbK(nksym) )
!    norbK(1:nksym) = norb   ! by default set to norb
!    maxnorb=nbnd
    !  
    allocate( igksym(npwx), ngksym(nkstot), xksym(3,nkstot) )
    if( ionode ) xksym(1:3,1:nkstot) = xk(1:3,1:nkstot)
    call mp_bcast(xksym, root_image, intra_image_comm)
    do ik=1,nkstot
        CALL gk_sort (xksym (1:3, ik), ngm, g, ecutwfc / tpiba2, &
                  ngksym(ik), igksym(1), g2kin)
    enddo
    write(*,*) 'Maximum/Minimum number of plane waves: ', &
            maxval(ngksym(:)),minval(ngksym(:))

    allocate( QKtoK2(nksym,nksym), kminus(nksym), Qpts(3,nksym), &
              xQ(nksym), wQ(nksym) )

    ! generate QMCPACK maps between {Ka, Kb} <--> {Q=(Ka-Kb)+G, K=Ka}
    CALL calculate_Qpt_map(nksym,xksym,Qpts,QKtoK2,kminus,nq1,nq2,nq3,use_regularization)
    !
    if(use_regularization) then
      ! find nqX, call 
      IF ( use_coulomb_vcut_ws .OR. use_coulomb_vcut_spheric ) THEN
         !
         ! build the superperiodicity direct lattice
         !
         atws = alat * at
         !
         atws(:,1) = atws(:,1) * nq1
         atws(:,2) = atws(:,2) * nq2
         atws(:,3) = atws(:,3) * nq3
         !
         CALL vcut_init( vcut, atws, ecutvcut )
         !
         IF (ionode) CALL vcut_info( stdout , vcut )
         !
      ENDIF
      !
      exxdiv = exx_divergence()  
      !
    endif
    ! from PW/src/setup.f90
    !
    !  ... generate transformation matrices for the crystal point group
    !  ... First we generate all the symmetry matrices of the Bravais lattice
    !
    if(use_symm) then
      !
      call set_sym_bl ( )
      !
      magnetic = .false.
      gate = .false.
      ALLOCATE( m_loc( 3, nat ) )
      CALL find_sym ( nat, tau, ityp, magnetic, m_loc, gate ) 
      DEALLOCATE(m_loc)
      !
      CALL find_Q_symm(nksym, xksym, Qpts, QKtoK2, nQuniq, xQ, wQ)
      !
      write(*,*) 'Found ',nQuniq,' Qpts in the irreducible subgroup.'
      do ik=1,nQuniq
        write(*,'(5x,2i5,f12.7)') ik,xQ(ik),wQ(ik)
        write(*,'(8x,3f12.7)') (Qpts(j,xQ(ik)),j=1,3)
      enddo
    else
      write(*,*) 'Not using symmetry.'
      nQuniq=nksym
      wQ(:) = 1.d0/dble(nksym)
      do ik=1,nksym
        xQ(ik)=ik
      enddo
    endif
    !
  END SUBROUTINE init_posthf

  SUBROUTINE clean_posthf()
    !  
    IMPLICIT NONE
    !
    ! add clean up routines in all posthf modules and call them here!
    ! right now there are arrays that are not deallocated!
    !
    if(allocated(igksym)) deallocate(igksym)  
    if(allocated(ngksym)) deallocate(ngksym)
    if(allocated(xksym)) deallocate(xksym)
    if(allocated(xkfull)) deallocate(xkfull)
    if(allocated(QKtoK2)) deallocate(QKtoK2)
    if(allocated(Qpts)) deallocate(Qpts)
    if(allocated(kminus)) deallocate(kminus)
    !
  END SUBROUTINE clean_posthf
  !
END MODULE posthf_mod
