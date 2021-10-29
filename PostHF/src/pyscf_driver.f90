!--------------------------------------------------------------------
! Written by Miguel A. Morales, LLNL, 2020 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!
!  !!! should rename pyscf_driver!!!
!
! called first to initialize QE
subroutine pyscf_driver_init(nim, npl, nb, norb_, pref, outd,verbose, &
                nktot_,natm_,nsp_,npwx_,ngm_,pwmesh_)
  !
  USE KINDS, ONLY: DP
  USE io_files,  ONLY : prefix, tmp_dir
  USE io_global, ONLY : ionode 
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start
  USE command_line_options, ONLY : set_command_line
  USE lsda_mod, ONLY: lsda
  USE klist, ONLY: nkstot,qnorm,nks,xk
  USE gvect, ONLY: ngm, ngm_g
  USE ions_base,    ONLY : nat, nsp
  USE wvfct, ONLY: npwx
  USE klist,  ONLY : nkstot
  use fft_base,             ONLY : dffts
  USE posthf_mod, ONLY: init_posthf
  USE mp_bands,     ONLY: nproc_bgrp
  USE realus,   ONLY : real_space
  USE paw_variables, ONLY : okpaw
  USE uspp,       ONLY : okvan
  USE control_flags, ONLY : tqr
  USE gvecs, ONLY : doublegrid
  USE orbital_generators, ONLY: orb_verbose => verbose
  USE twobody_hamiltonian, ONLY: chol_verbose => verbose
  USE mp2_module, ONLY: mp2_verbose => verbose
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nim, npl, nb, norb_
  LOGICAL, INTENT(IN) :: verbose
  CHARACTER(len=256), INTENT(IN) :: pref, outd
  INTEGER, INTENT(OUT) :: nktot_,natm_,nsp_,npwx_,ngm_,pwmesh_(3)
  !
!  INTEGER :: natm_,nsp_,npwx_,ngm_,pwmesh_(3)
  INTEGER :: ik,iq,nxxs
  REAL(DP) :: qnorm_
  LOGICAL :: start_images
  !
  ! needed for PAW/USPP and kpoints
  qnorm = 3.5d0
  tmp_dir = TRIM(outd)
  prefix = TRIM(pref)
  chol_verbose = verbose
  orb_verbose = verbose
  mp2_verbose = verbose
  !
  start_images = (nim>1)  
  !
  CALL set_command_line(nimage=nim, npool=npl, nband=nb)
  CALL mp_startup ( start_images=start_images )
  CALL environment_start ( 'pyscf_driver' )
  !
  CALL start_clock ( 'read_file' )
  CALL read_file
  CALL stop_clock ( 'read_file' )
  qnorm_ = 0.0_dp
  DO iq = 1,nks
     DO ik = iq+1,nks
        qnorm_ = max(qnorm_, sqrt( sum((xk(:,ik)-xk(:,iq))**2) ))
     ENDDO
  ENDDO
  if(qnorm_ > qnorm) call errore('pyscf_driver: qnorm guess is wrong.',1)
  !
  CALL openfil_pp
  !
  ! some checks
  ! 
  if(nproc_bgrp > 1) & 
    call errore('pyscf_driver', &
                'Error: Only pool and band parallelization allowed.',1)
  if(doublegrid) & 
    call errore('pyscf_driver','Error: doublegrid not implemented',1)
  if(real_space) &
    call errore('pyscf_driver','real_space not allowed',1)
  if(tqr .and. nkstot>1) &
    call errore('pyscf_driver','tqr with nkpots>1 not working',1)
  if( ngm .ne. ngm_g ) &
    call errore('pyscf_driver','Error: No PW parallelization allowed.',1) 
  !
  if(ionode) then
    if(okvan) write(*,*) ' Found USPP. '
    if(okpaw) write(*,*) ' Found PAW. '
  endif
  !
  nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x   ! global size of FFT grid
  write(*,*) 'nr1, nr2, nr3, nnr:',dffts%nr1,dffts%nr2,dffts%nr3,dffts%nnr
  write(*,*) 'ngm: ',dffts%ngm
  write(*,*) 'nkstot: ',nkstot

  if( nxxs .ne. dffts%nnr ) then
    CALL errore('pw2posthf', ' nnr .ne. nr1x*nr2x*nr3x ', 1)
  endif
  !
  ! initialize posthf module
  !
  call init_posthf(norb_)  
  !
  ! returned values
  !
  nktot_ = nkstot
  if(lsda) nktot_ = nktot_/2 
  natm_ = nat 
  nsp_ = nsp
  npwx_ = npwx
  ngm_ = ngm_g
  pwmesh_(1) = dffts%nr1  
  pwmesh_(2) = dffts%nr2  
  pwmesh_(3) = dffts%nr3 
  !
end subroutine pyscf_driver_init
!
! called at the end to finalize QE
!
subroutine pyscf_driver_end()
  !
  USE environment,ONLY : environment_end
  USE posthf_mod, ONLY: clean_posthf
  USE io_files, ONLY: iunwfc
  USE mp_global, ONLY: mp_global_end
  !
  IMPLICIT NONE
  !
  INTEGER :: info
  LOGICAL :: op
  ! 
  call clean_posthf()  
  !
  CALL environment_end ( 'pyscf_driver' )
  ! stop_pp without the call to STOP
  !
  INQUIRE ( iunwfc, opened = op )
  IF ( op ) CLOSE (unit = iunwfc, status = 'delete')
   
  CALL mp_global_end()
  ! 
  !
end subroutine pyscf_driver_end
!
! Returns basis info needed to setup cell object in pyscf
!
!subroutine pyscf_driver_get_info(na_,ns_,nk_,species_ids, atom_ids, atom_pos, kpts, lattice)
subroutine pyscf_driver_get_info(na_,ns_,nk_,atom_ids,atom_pos,kpts,lattice)
  !
  USE ions_base,    ONLY : nat, nsp, tau, ityp
  USE cell_base,    ONLY : at, tpiba, alat
  USE lsda_mod, ONLY: lsda
  USE klist, ONLY: nkstot,xk
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na_,nk_,ns_
!  CHARACTER(len=3), INTENT(INOUT) :: species_ids(ns_)  
  INTEGER, INTENT(OUT) :: atom_ids(na_)  
  REAL(kind=8), INTENT(OUT) :: atom_pos(3,na_)  
  REAL(kind=8), INTENT(OUT) :: kpts(3,nk_)  
  REAL(kind=8), INTENT(OUT) :: lattice(3,3)  
  !
  INTEGER :: nk
  !  
!  do i =1,nsp
!    species_ids(i) = TRIM(atm(i))
!  enddo  
  atom_ids(1:nat) = ityp(1:nat)
  atom_pos(1:3,1:nat) = tau(1:3,1:nat)*alat
  nk = nkstot
  if(lsda) nk = nk/2
  kpts(1:3,1:nk) = xk(1:3,1:nk)*tpiba  
  lattice(:,:) = at(:,:)*alat
  !
end subroutine pyscf_driver_get_info
!
! driver for mp2 calculation  
!  - orbitals are taken from QE and/or from esh5 files
!  - resulting orbital set can be pseudo-canonicalized (is this a word?) 
!  - mp2 calculation follows
! 
! This driver is meant to be called many times with similar settings,
!   e.g.: norb, etc...
!
subroutine pyscf_driver_mp2no(out_prefix_,diag,diag_type,appnos,nread_from_h5, &
                         h5_add_orbs_,nskip,eigcut,nextracut,mp2noecut,&
                         regkappa,regp)
  !
!  USE posthf_mod, ONLY: 
  USE KINDS, ONLY: DP
  USE io_files,  ONLY : tmp_dir 
  USE mp_global, ONLY : nimage, nproc_image, me_image, root_image, &
                        intra_image_comm 
  USE mp, ONLY: mp_barrier  
  USE mp_pools, ONLY: npool 
  USE orbital_generators, ONLY: generate_orbitals, &
                                diag_hf, davcio_to_esh5, generate_full_pw_basis
  USE onebody_hamiltonian, ONLY: getH1
  USE twobody_hamiltonian, ONLY: cholesky_r, calculate_KS_bscorr
  USE mp2_module, ONLY: mp2no_g,approx_mp2no
  use fft_base,             ONLY : dffts
  USE read_orbitals_from_file, ONLY: open_esh5_read,close_esh5_read,h5file_type 
  USE posthf_mod, ONLY: norb
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: h5_add_orbs_,diag_type,out_prefix_
  REAL(kind=8), INTENT(IN) :: eigcut, nextracut, regkappa, mp2noecut 
  LOGICAL, INTENT(IN) :: diag,appnos
  INTEGER, INTENT(IN) :: nread_from_h5, regp, nskip 
  !
  REAL(DP) :: occeigcut
  LOGICAL :: diagonalize  
  INTEGER :: h5len, oldh5  
  CHARACTER(len=256) :: out_prefix, h5_add_orbs, h5qeorbs
  CHARACTER(len=256) :: hamilFile,orbsFile,canOrbsFile,mp2noFile
  COMPLEX(DP) :: energy
  TYPE(h5file_type) :: h5id_hamil,h5id_orbs
   
  occeigcut = 1.d-8
  !
  out_prefix = TRIM(out_prefix_) 
  h5_add_orbs = TRIM(h5_add_orbs_)
  hamilFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.hamil.h5'
  orbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
  mp2noFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.mp2no.h5'
  canOrbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
  h5qeorbs = TRIM( tmp_dir )//TRIM( out_prefix ) // '.qe.orbs.h5'
  !
  !
  diagonalize = diag
  if( (TRIM(diag_type) == 'fullpw') .or. (nread_from_h5 > 0 ) ) & 
    diagonalize = .true.
  !    
!      all orbitals are accessible through davcio
  if( (TRIM(diag_type) == 'fullpw') .or. &
      (npool > 1) ) call davcio_to_esh5(h5qeorbs,norb,dffts)
  !  
  if(me_image == root_image) then
    !
    ! 1. generate orbitals
    !
    write(*,*) 'Generating orbitals' 
    if( TRIM(diag_type) == 'fullpw' ) then
      CALL generate_full_pw_basis(out_prefix,dffts)  
    else if(npool > 1) then
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.,esh5_orbs=h5qeorbs)
    else
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.)
    endif
    !
    ! 2. diagonalize HF is requested
    !
    if(diagonalize) then
      h5len = LEN_TRIM(hamilFile)
      oldh5=0
      CALL esh5_posthf_open_file(h5id_hamil%id,hamilFile,h5len,oldh5)
      ! all tasks open orbital file in read_only mode
      call open_esh5_read(h5id_orbs,orbsFile)
      if( h5id_orbs%grid_type .ne. 1) &
        call errore('pyscf_driver','grid_type ne 1',1)
      write(*,*) 'Generating 1-Body Hamiltonian'
      CALL getH1(h5id_orbs,h5id_hamil,dffts)
      CALL esh5_posthf_close_file(h5id_hamil%id)
      call close_esh5_read(h5id_orbs)

      call mp_barrier(intra_image_comm)  
      write(*,*) ' Diagonalizing HF hamiltonian' 
      ! orbitals are read from esh5 only right now
      if(TRIM(diag_type) == 'full') then
        call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile)
      else if(TRIM(diag_type) == 'keep_occ') then
        call diag_hf(dffts,'vir_keep_occ_et',hamilFile,orbsFile,orbsFile,canOrbsFile)
      else if(TRIM(diag_type) == 'fullpw') then
        call diag_hf(dffts,'full',hamilFile,h5qeorbs,orbsFile,canOrbsFile)
      endif
    endif
  else
    if(diagonalize) then  
      call mp_barrier(intra_image_comm)  
      ! orbitals are read from esh5 only right now
      if(TRIM(diag_type) == 'full') then
        call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile)
      elseif(TRIM(diag_type) == 'keep_occ') then
        call diag_hf(dffts,'vir_keep_occ_et',hamilFile,orbsFile,orbsFile,canOrbsFile)
      else if(TRIM(diag_type) == 'fullpw') then
        call diag_hf(dffts,'full',hamilFile,h5qeorbs,orbsFile,canOrbsFile)
      endif
    endif
  endif
  call mp_barrier(intra_image_comm)  
  !
  ! 3. calculate mp2nos
  !
  if(diagonalize) orbsFile = canOrbsFile
  write(*,*) ' Starting MP2NO'
  if(appnos) then
!    call approx_mp2no(energy,dffts,orbsFile,reg_expo=regkappa,reg_pow=regp)
    call errore('pyscf_driver_mp2no','Finish appnos',1)
  else  
    call mp2no_g(dffts,mp2noFile,nskip,mp2noecut,  &
                 esh5_file=orbsFile)  !,reg_expo=regkappa,reg_pow=regp)
  endif  
end subroutine pyscf_driver_mp2no

subroutine pyscf_driver_mp2(out_prefix_,diag,diag_type,nread_from_h5, &
                         h5_add_orbs_,eigcut,nextracut,regkappa,regp,emp2)
  !
!  USE posthf_mod, ONLY: 
  USE KINDS, ONLY: DP
  USE io_files,  ONLY : tmp_dir 
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : nimage, nproc_image, me_image, root_image, &
                        intra_image_comm 
  USE mp, ONLY: mp_barrier  
  USE mp_pools, ONLY: npool 
  USE orbital_generators, ONLY: generate_orbitals,write_trial_wavefunction, &
                                diag_hf, davcio_to_esh5, generate_full_pw_basis
  USE onebody_hamiltonian, ONLY: getH1
  USE twobody_hamiltonian, ONLY: cholesky_r, calculate_KS_bscorr
  USE mp2_module, ONLY: mp2_g,mp2no,approx_mp2no
#if defined(__CUDA)
  USE mp2_module, ONLY: mp2_gpu
#endif
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: nbnd, wg
  USE klist,  ONLY : wk
  use fft_base,             ONLY : dffts
  USE read_orbitals_from_file, ONLY: open_esh5_read,close_esh5_read,h5file_type 
  USE posthf_mod, ONLY: norb,nksym,nmax_DM,DM,DM_mf,e0, efc, numspin
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: h5_add_orbs_,diag_type,out_prefix_
  REAL(kind=8), INTENT(IN) :: eigcut, nextracut, regkappa 
  LOGICAL, INTENT(IN) :: diag
  INTEGER, INTENT(IN) :: nread_from_h5, regp 
  REAL(kind=8), INTENT(OUT) :: emp2
  !
  REAL(DP) :: occeigcut
  LOGICAL :: diagonalize  
  INTEGER :: h5len, oldh5  
  CHARACTER(len=256) :: out_prefix, h5_add_orbs, h5qeorbs
  CHARACTER(len=256) :: hamilFile,orbsFile,canOrbsFile
  COMPLEX(DP) :: energy
  TYPE(h5file_type) :: h5id_hamil,h5id_orbs
  INTEGER :: maxnorb, nelmax, error, ndet
  COMPLEX(DP), ALLOCATABLE :: M(:,:,:,:)  ! Overlap between basis states and
                                            ! occupied KS states.    
  ndet=1
  occeigcut=1.d-8
  !
  out_prefix = TRIM(out_prefix_) 
  h5_add_orbs = TRIM(h5_add_orbs_)
  hamilFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.hamil.h5'
  orbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
  canOrbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
  h5qeorbs = TRIM( tmp_dir )//TRIM( out_prefix ) // '.qe.orbs.h5'
  !
  !
  diagonalize = diag
  if( (TRIM(diag_type) == 'fullpw') .or. (nread_from_h5 > 0 ) ) & 
    diagonalize = .true.
!      all orbitals are accessible through davcio
  if( (TRIM(diag_type) == 'fullpw') .or. &
      (npool > 1) ) call davcio_to_esh5(h5qeorbs,norb,dffts)
  !  
  if(me_image == root_image) then
    !
    ! 1. generate orbitals
    !
    write(*,*) 'Generating orbitals' 
    if( TRIM(diag_type) == 'fullpw' ) then
      CALL generate_full_pw_basis(out_prefix,dffts)  
    else if(npool > 1) then
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.,esh5_orbs=h5qeorbs)
    else
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.)
    endif
    !
    ! 2. diagonalize HF is requested
    !
    if(diagonalize) then
      h5len = LEN_TRIM(hamilFile)
      oldh5=0
      CALL esh5_posthf_open_file(h5id_hamil%id,hamilFile,h5len,oldh5)
      ! all tasks open orbital file in read_only mode
      call open_esh5_read(h5id_orbs,orbsFile)
      if( h5id_orbs%grid_type .ne. 1) &
        call errore('pyscf_driver','grid_type ne 1',1)

      ! mp2 now requires Overlap matrix between basis and KS/HF states
      call get_nelmax(nbnd,nksym,numspin,wk,size(wg,1),wg,nelmax)
      maxnorb = maxval(h5id_orbs%norbK(:))
      write(*,*)'maxnorb, nelmax:',maxnorb, nelmax
      allocate( M(maxnorb,nelmax,nksym,min(2,nspin)) )
      if(npool > 1) then
        call write_trial_wavefunction(h5id_orbs,h5id_hamil,&
                    dffts,nelmax,M,ndet,0.05d0,0.95d0,h5wfn=h5qeorbs)
      else
        call write_trial_wavefunction(h5id_orbs,h5id_hamil,&
                    dffts,nelmax,M,ndet,0.05d0,0.95d0)
      endif

      write(*,*) 'Generating 1-Body Hamiltonian'
      CALL getH1(h5id_orbs,h5id_hamil,dffts)
      CALL esh5_posthf_close_file(h5id_hamil%id)
      call close_esh5_read(h5id_orbs)

      h5len = LEN_TRIM(orbsFile)
      oldh5 = 1
      CALL esh5_posthf_open_file(h5id_hamil%id,TRIM(orbsFile),h5len,oldh5)
      call esh5_posthf_write_orbmat(h5id_hamil%id,"SCFOrbMat",9, &
            maxnorb,nelmax,nksym,size(M,4),M,error)
      if(error .ne. 0 ) &
        call errore('pyscf_driver::mp2_driver','error writing orbmat',1)
      CALL esh5_posthf_close_file(h5id_hamil%id)
      if(allocated(M)) deallocate(M)

      call mp_barrier(intra_image_comm)  
      write(*,*) ' Diagonalizing HF hamiltonian' 
      ! orbitals are read from esh5 only right now
      if(TRIM(diag_type) == 'full') then
        call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile)
      else if(TRIM(diag_type) == 'keep_occ') then
        call diag_hf(dffts,'vir_keep_occ_et',hamilFile,orbsFile,orbsFile,canOrbsFile)
      else if(TRIM(diag_type) == 'fullpw') then
        call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile)
      endif
    endif
  else
    if(diagonalize) then  
      call mp_barrier(intra_image_comm)  
      ! orbitals are read from esh5 only right now
      if(TRIM(diag_type) == 'full') then
        call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile)
      elseif(TRIM(diag_type) == 'keep_occ') then
        call diag_hf(dffts,'vir_keep_occ_et',hamilFile,orbsFile,orbsFile,canOrbsFile)
      else if(TRIM(diag_type) == 'fullpw') then
        call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile)
      endif
    endif
  endif
  call mp_barrier(intra_image_comm)  
  !
  ! 3. calculate mp2
  !
  if(diagonalize) orbsFile = canOrbsFile 
  write(*,*) ' Starting MP2' 
  FLUSH( stdout )
#if defined(__CUDA)
  call mp2_gpu(energy,dffts,orbsFile,reg_expo=regkappa,reg_pow=regp)
#else
  call mp2_g(energy,dffts,orbsFile,reg_expo=regkappa,reg_pow=regp)
#endif
  emp2 = real(dble(energy),kind=8)
  !
  call print_clock('')    
  !
end subroutine pyscf_driver_mp2

! MAM: This routine right now forces diagonalization, to generate a file with
! the correct structure. If you want to allow calculations directly on the QE
! basis, generate a pass through that skips the first step and simply writes the
! qeorbs in the correct h5 format. Then call rpa routine.
subroutine pyscf_driver_rpa(out_prefix_,diag_type,nread_from_h5, &
                         h5_add_orbs_,eigcut,nextracut,thresh, ncholmax,erpa)
  !
!  USE posthf_mod, ONLY: 
  USE KINDS, ONLY: DP
  USE io_files,  ONLY : tmp_dir 
  USE io_global, ONLY : stdout
  USE mp_global, ONLY : nimage, nproc_image, me_image, root_image, &
                        intra_image_comm 
  USE mp, ONLY: mp_barrier, mp_bcast 
  USE mp_pools, ONLY: npool 
  USE orbital_generators, ONLY: generate_orbitals,write_trial_wavefunction, &
                                diag_hf, davcio_to_esh5, generate_full_pw_basis
  USE onebody_hamiltonian, ONLY: getH1
  USE twobody_hamiltonian, ONLY: cholesky_r, cholesky_MO
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: nbnd, wg
  USE klist,  ONLY : wk
  use fft_base,             ONLY : dffts
  USE read_orbitals_from_file, ONLY: open_esh5_read,close_esh5_read,h5file_type 
  USE posthf_mod, ONLY: norb,nksym,nmax_DM,DM,DM_mf,e0, efc, numspin
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: h5_add_orbs_,diag_type,out_prefix_
  REAL(kind=8), INTENT(IN) :: eigcut, nextracut, thresh, ncholmax
  INTEGER, INTENT(IN) :: nread_from_h5 
  REAL(kind=8), INTENT(OUT) :: erpa
  !
  REAL(DP) :: occeigcut
  LOGICAL :: diagonalize  
  INTEGER :: h5len, oldh5  
  CHARACTER(len=256) :: out_prefix, h5_add_orbs, h5qeorbs
  CHARACTER(len=256) :: hamilFile,orbsFile,canOrbsFile
  COMPLEX(DP) :: energy
  TYPE(h5file_type) :: h5id_hamil,h5id_orbs
  INTEGER :: maxnorb, nelmax, error, ndet, nel(2)
  INTEGER, ALLOCATABLE :: noccK(:,:)
  REAL(DP), ALLOCATABLE :: weights(:,:), eigval(:,:)
  COMPLEX(DP), ALLOCATABLE :: M(:,:,:,:)  ! Overlap between basis states and
                                            ! occupied KS states.    
  ndet=1
  occeigcut=1.d-8
  !
  out_prefix = TRIM(out_prefix_) 
  h5_add_orbs = TRIM(h5_add_orbs_)
  hamilFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.MOCholesky.h5'
  orbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
  canOrbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
  h5qeorbs = TRIM( tmp_dir )//TRIM( out_prefix ) // '.qe.orbs.h5'
  !
  !
  diagonalize = .true. 
  if( (TRIM(diag_type) == 'fullpw') .or. (nread_from_h5 > 0 ) ) & 
    diagonalize = .true.
!      all orbitals are accessible through davcio
  if( (TRIM(diag_type) == 'fullpw') .or. &
      (npool > 1) ) call davcio_to_esh5(h5qeorbs,norb,dffts)
  !  
  if(me_image == root_image) then
    !
    ! 1. generate orbitals
    !
    write(*,*) 'Generating orbitals' 
    if( TRIM(diag_type) == 'fullpw' ) then
      CALL generate_full_pw_basis(out_prefix,dffts)  
    else if(npool > 1) then
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.,esh5_orbs=h5qeorbs)
    else
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.)
    endif
    !
    ! 2. diagonalize HF is requested
    !
    if(diagonalize) then
      h5len = LEN_TRIM(hamilFile)
      oldh5=0
      CALL esh5_posthf_open_file(h5id_hamil%id,hamilFile,h5len,oldh5)
      ! all tasks open orbital file in read_only mode
      call open_esh5_read(h5id_orbs,orbsFile)
      if( h5id_orbs%grid_type .ne. 1) &
        call errore('pyscf_driver','grid_type ne 1',1)

      ! mp2 now requires Overlap matrix between basis and KS/HF states
      call get_nelmax(nbnd,nksym,numspin,wk,size(wg,1),wg,nelmax)
      maxnorb = maxval(h5id_orbs%norbK(:))
      call mp_bcast(nelmax, root_image, intra_image_comm)
      write(*,*)'maxnorb, nelmax:',maxnorb, nelmax
      allocate( M(maxnorb,nelmax,nksym,min(2,nspin)) )
      if(npool > 1) then
        call write_trial_wavefunction(h5id_orbs,h5id_hamil,&
                    dffts,nelmax,M,ndet,0.05d0,0.95d0,h5wfn=h5qeorbs)
      else
        call write_trial_wavefunction(h5id_orbs,h5id_hamil,&
                    dffts,nelmax,M,ndet,0.05d0,0.95d0)
      endif

      write(*,*) 'Generating 1-Body Hamiltonian'
      CALL getH1(h5id_orbs,h5id_hamil,dffts)
      CALL esh5_posthf_close_file(h5id_hamil%id)
      call close_esh5_read(h5id_orbs)

      h5len = LEN_TRIM(orbsFile)
      oldh5 = 1
      CALL esh5_posthf_open_file(h5id_hamil%id,TRIM(orbsFile),h5len,oldh5)
      call esh5_posthf_write_orbmat(h5id_hamil%id,"SCFOrbMat",9, &
            maxnorb,nelmax,nksym,size(M,4),M,error)
      if(error .ne. 0 ) &
        call errore('pyscf_driver::rpa_driver','error writing orbmat',1)
      CALL esh5_posthf_close_file(h5id_hamil%id)
      if(allocated(M)) deallocate(M)

      call mp_barrier(intra_image_comm)  
      write(*,*) ' Diagonalizing HF hamiltonian' 
      ! orbitals are read from esh5 only right now
      call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile,.true.)
    endif
  else
    if(diagonalize) then  
      call mp_bcast(nelmax, root_image, intra_image_comm)
      call mp_barrier(intra_image_comm)  
      ! orbitals are read from esh5 only right now
      call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile,.true.)
    endif
  endif
  call mp_barrier(intra_image_comm)  
  !
  if(diagonalize) orbsFile = canOrbsFile 
  !
  ! calculate cholesky matrix  
! can use the resulting cholesky decomposition for the diaghf step  
  call cholesky_MO(nelmax,ncholmax,thresh,dffts,hamilFile,orbsFile)
  !
  ! 3. calculate rpa 
  !
  write(*,*) ' Starting RPA' 
  FLUSH( stdout )
!  call rpa_cholesky_MO(energy,nelmax,dffts,hamilFile,orbsFile)
  erpa = real(dble(energy),kind=8)
  !
  call print_clock('')    
  !
end subroutine pyscf_driver_rpa

subroutine pyscf_driver_hamil(out_prefix_, nread_from_h5, h5_add_orbs_, &
       ndet, eigcut, nextracut, thresh, ncholmax, &
       get_hf, get_mp2, get_rpa, update_qe_bands, ehf, emp2, erpa)
  !
  USE mp_global, ONLY : nimage, nproc_image, me_image, root_image, &
                        intra_image_comm
  USE mp, ONLY: mp_barrier  
  USE mp_pools, ONLY: npool
  use fft_base,             ONLY : dffts
  USE KINDS, ONLY: DP
  USE io_files,  ONLY : tmp_dir
  use fft_base,             ONLY : dffts
  USE lsda_mod, ONLY: nspin
  USE wvfct, ONLY: nbnd, wg 
  USE klist,  ONLY : wk
  USE orbital_generators, ONLY: generate_orbitals,write_trial_wavefunction, &
                                diag_hf, verbose, davcio_to_esh5
  USE onebody_hamiltonian, ONLY: getH1
  USE twobody_hamiltonian, ONLY: cholesky_r, calculate_KS_bscorr
  USE mp2_module, ONLY: mp2_g,mp2no,approx_mp2no
  USE rpa_module, ONLY: rpa_cholesky
#if defined(__CUDA)
  USE mp2_module, ONLY: mp2_gpu
#endif
  USE read_orbitals_from_file, ONLY: open_esh5_read,close_esh5_read, &
                                     update_bands,h5file_type
  USE posthf_mod, ONLY: norb,nksym,nmax_DM,DM,DM_mf,e0, efc, numspin
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: h5_add_orbs_,out_prefix_
  REAL(kind=8), INTENT(IN) :: eigcut, nextracut, thresh, ncholmax
  LOGICAL, INTENT(IN) :: get_hf, get_mp2, get_rpa, update_qe_bands
  INTEGER, INTENT(IN) :: nread_from_h5, ndet
  REAL(kind=8), INTENT(OUT) :: ehf, emp2, erpa
  !
  INTEGER :: h5len, oldh5, error
  INTEGER :: maxnorb, nelmax
  COMPLEX(DP), ALLOCATABLE :: M(:,:,:,:)  ! Overlap between basis states and
                                            ! occupied KS states. 
  CHARACTER(len=256) :: out_prefix, h5_add_orbs, h5qeorbs
  CHARACTER(len=256) :: hamilFile,orbsFile,canOrbsFile
  COMPLEX(DP) :: e1, e1_so, e2_mp2, e2_rpa, e1_mf, e1_so_mf
  TYPE(h5file_type) :: h5id_hamil,h5id_orbs
  REAL(DP) :: occeigcut  

  occeigcut = 1.d-8
  !
  out_prefix = TRIM(out_prefix_)
  h5_add_orbs = TRIM(h5_add_orbs_)
  hamilFile = TRIM( tmp_dir )//TRIM( out_prefix ) // '.hamil.h5'
  orbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
  canOrbsFile = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
  h5qeorbs = TRIM( tmp_dir )//TRIM( out_prefix ) // '.qe.orbs.h5'
  !
  ehf=0.d0
  emp2=0.d0  
  erpa=0.d0  
  e1=(0.d0,0.d0) 
  e1_so=(0.d0,0.d0)  
  e1_mf=(0.d0,0.d0) 
  e1_so_mf=(0.d0,0.d0)  
  e2_mp2=(0.d0,0.d0) 
  !
  if(npool > 1) call davcio_to_esh5(h5qeorbs,norb,dffts)
  !
  if(me_image == root_image) then
    !
    ! 1. generate orbitals
    !
    if(npool > 1) then
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.,esh5_orbs=h5qeorbs)
    else
      CALL generate_orbitals(out_prefix,dffts,nread_from_h5,h5_add_orbs, &
            eigcut,occeigcut,nextracut,.false.)
    endif   
    !
    h5len = LEN_TRIM(hamilFile)
    oldh5=0
    CALL esh5_posthf_open_file(h5id_hamil%id,hamilFile,h5len,oldh5)

    ! all tasks open orbital file in read_only mode
    call open_esh5_read(h5id_orbs,orbsFile)
    if( h5id_orbs%grid_type .ne. 1) &
      call errore('pyscf_driver_hamil','grid_type ne 1',1)

    ! calculates expansion of occupied KS states in the generated  
    ! (spin-independent) basis and writes trial wfn to file
    call get_nelmax(nbnd,nksym,numspin,wk,size(wg,1),wg,nelmax)
    maxnorb = maxval(h5id_orbs%norbK(:))
    allocate( M(maxnorb,nelmax,nksym,min(2,nspin)) )
    if(npool > 1) then 
      call write_trial_wavefunction(h5id_orbs,h5id_hamil,dffts, &
                    nelmax,M,ndet,0.05d0,0.95d0,h5wfn=h5qeorbs)
    else
      call write_trial_wavefunction(h5id_orbs,h5id_hamil,dffts,&
                    nelmax,M,ndet,0.05d0,0.95d0)
    endif

    ! calculate and write 1 body hamiltonian
    ! make DM and optional argument to H1
    CALL getH1(h5id_orbs,h5id_hamil,dffts,e1,e1_so,e1_mf,e1_so_mf)

    CALL esh5_posthf_close_file(h5id_hamil%id)
    call close_esh5_read(h5id_orbs)

    h5len = LEN_TRIM(orbsFile)
    oldh5 = 1
    CALL esh5_posthf_open_file(h5id_hamil%id,TRIM(orbsFile),h5len,oldh5)

    call esh5_posthf_write_orbmat(h5id_hamil%id,"SCFOrbMat",9,maxnorb,nelmax,nksym,size(M,4),M,error)
    if(error .ne. 0 ) &
      call errore('pw2posthf','error writing orbmat',1)

    call esh5_posthf_write_dm(h5id_hamil%id,"DM",2,nmax_DM,nksym,nspin,DM,error)
    if(error .ne. 0 ) &
      call errore('pyscf_driver_hamil','error writing DM',1)
    call esh5_posthf_write_dm(h5id_hamil%id,"DM_mf",5,nmax_DM,nksym,nspin,DM_mf,error)
    if(error .ne. 0 ) &
      call errore('pyscf_driver_hamil','error writing DM_mf',1)

    write(*,*) 'E0, E1(1Det),E_SO(1Det) (Ha):',e0,e1/(nksym*1.0d0),e1_so/(nksym*1.0d0)
    write(*,*) '    E1,E_SO (Ha):',e1_mf/(nksym*1.0),e1_so_mf/(nksym*1.0)

    CALL esh5_posthf_close_file(h5id_hamil%id)
    if(allocated(M)) deallocate(M)
  endif
  call mp_barrier(intra_image_comm)
  if(update_qe_bands) then
    call update_bands(orbsFile,dffts)
  endif
  call mp_barrier(intra_image_comm)

  ! calculate cholesky matrix  
! note: you can construct FockM here if you want, it would be 
!       quite cheap. Add it as an optional
  call cholesky_r(ncholmax,thresh,dffts,hamilFile,orbsFile)

  ! mp2 requires hf rediag
  if(get_hf .or. get_mp2 .or. get_rpa) then
    call mp_barrier(intra_image_comm)
    ! orbitals are read from esh5 only right now
    call diag_hf(dffts,'full',hamilFile,orbsFile,orbsFile,canOrbsFile,.true.)
  endif

  if(get_mp2) then
    call mp_barrier(intra_image_comm)
#if defined(__CUDA)
    call mp2_gpu(e2_mp2,dffts,canOrbsFile)
#else
    call mp2_g(e2_mp2,dffts,canOrbsFile)
#endif
    emp2 = real(dble(e2_mp2),kind=8)
  endif
  call mp_barrier(intra_image_comm)
  !

  if(get_rpa) then
    call mp_barrier(intra_image_comm)
    call rpa_cholesky(e2_rpa,dffts,hamilFile,canOrbsFile)
    erpa = real(dble(e2_rpa),kind=8)
  endif
  call mp_barrier(intra_image_comm)
  !
  call print_clock('')    
  !
end subroutine pyscf_driver_hamil

! some of these drivers need to call pwscf from python. figure out how to do it!
! for some of these drivers
! to write: include option to calculate hf and mp2 energy on given basis
! 1. get_hamil  (ks + external basis, h5_add)
! 2. mp2no + rediag + get_hamil
! 3. mp2, write routine to diagonalize fock matrix on full npwx space
!    from a previously converged set of states (hopefully HF states) 
!    Use this to get mp2 energy with all pws on a given cutoff from a converged
!    hf calculation.
!
