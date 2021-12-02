!
! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM posthf
  !----------------------------------------------------------------------- 

  ! This subroutine writes the file "prefix".pwscf.xml and "prefix".pwscf.h5
  ! containing the  plane wave coefficients and other stuff needed by QMCPACK. 

  USE io_files,  ONLY : nd_nmbr, prefix, tmp_dir
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, npool, nimage, nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : gamma_only
  USE mp_world,   ONLY : world_comm, nproc
  USE environment,ONLY : environment_start, environment_end
  USE klist, ONLY: qnorm,nks,xk
  USE input_parameters,  ONLY : lmoire, amoire_in_ang, vmoire_in_mev, pmoire_in_deg, mstar, epsmoire
  USE KINDS, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: ios, number_of_orbitals, ndet, nskipvir, read_from_h5
  INTEGER :: iq, ik, regp
  REAL(DP) :: ncholmax, nextracut, thresh, eigcut_occ, eigcut, regkappa
  LOGICAL :: write_psir, expand_kp, debug, verbose 
  LOGICAL :: low_memory, update_qe_bands,get_hf,get_mp2,use_symm 
  REAL(DP) :: qnorm_
  ! directory for temporary files
  CHARACTER(len=256) :: outdir, out_prefix, h5_add_orbs, run_type, diag_type
  CHARACTER(len=256) :: exxdiv_treatment 
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck

  NAMELIST / inputpp / prefix, outdir, write_psir, expand_kp, debug, & 
            ndet, ncholmax,thresh, number_of_orbitals,eigcut_occ, eigcut, out_prefix, verbose, & 
            h5_add_orbs, nskipvir,low_memory, get_hf,get_mp2,use_symm, &
            regp,regkappa,read_from_h5,nextracut, update_qe_bands, run_type, diag_type, exxdiv_treatment, &
            lmoire, amoire_in_ang, vmoire_in_mev, pmoire_in_deg, mstar, epsmoire
#ifdef __MPI
  CALL mp_startup ( )
#endif

  CALL environment_start ( 'posthf' )
#if defined(__HDF5) || defined(__HDF5_C)
  IF ( nimage > 1) &
     CALL errore('posthf', ' image parallelization not (yet) implemented', 1)

  !   CALL start_postproc(nd_nmbr)
  ! 
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf'
  out_prefix = 'pwscf'
  write_psir = .false.
  expand_kp = .false.
  number_of_orbitals = 0
  ndet=1
  ncholmax=15.d0
  thresh=1.d-5
  eigcut=1.d-2
  eigcut_occ=1.d-8
  debug = .false.
  verbose = .false.
  run_type = ''
  diag_type = 'full'
  low_memory = .false.
  update_qe_bands = .false.
  h5_add_orbs = '' 
  nskipvir=0
  read_from_h5=-1
  nextracut=-1.0
  get_hf=.true.
  get_mp2=.true.
  regp=0
  regkappa=0.d0
  use_symm=.true.  
  exxdiv_treatment = 'none'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  ios = 0
  IF ( ionode )  THEN 
     !
     CALL input_from_file ( )
     !READ (5, inputpp, err=200, iostat=ios)
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  CALL mp_bcast( ios, ionode_id, world_comm ) 
  IF ( ios/=0 ) CALL errore('posthf', 'reading inputpp namelist', ABS(ios))
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast(prefix, ionode_id, world_comm ) 
  CALL mp_bcast(out_prefix, ionode_id, world_comm ) 
  CALL mp_bcast(tmp_dir, ionode_id, world_comm ) 
  CALL mp_bcast(write_psir, ionode_id, world_comm ) 
  CALL mp_bcast(expand_kp, ionode_id, world_comm ) 
  CALL mp_bcast(number_of_orbitals, ionode_id, world_comm ) 
  CALL mp_bcast(ndet, ionode_id, world_comm ) 
  CALL mp_bcast(ncholmax, ionode_id, world_comm ) 
  CALL mp_bcast(thresh, ionode_id, world_comm ) 
  CALL mp_bcast(eigcut, ionode_id, world_comm ) 
  CALL mp_bcast(eigcut_occ, ionode_id, world_comm ) 
  CALL mp_bcast(debug, ionode_id, world_comm ) 
  CALL mp_bcast(run_type, ionode_id, world_comm ) 
  CALL mp_bcast(diag_type, ionode_id, world_comm )
  CALL mp_bcast(h5_add_orbs, ionode_id, world_comm ) 
  CALL mp_bcast(nskipvir, ionode_id, world_comm ) 
  CALL mp_bcast(low_memory, ionode_id, world_comm ) 
  CALL mp_bcast(read_from_h5, ionode_id, world_comm ) 
  CALL mp_bcast(nextracut, ionode_id, world_comm ) 
  CALL mp_bcast(update_qe_bands, ionode_id, world_comm ) 
  CALL mp_bcast(get_hf, ionode_id, world_comm ) 
  CALL mp_bcast(get_mp2, ionode_id, world_comm ) 
  CALL mp_bcast(use_symm, ionode_id, world_comm ) 
  CALL mp_bcast(exxdiv_treatment, ionode_id, world_comm ) 
  CALL mp_bcast(regp, ionode_id, world_comm ) 
  CALL mp_bcast(regkappa, ionode_id, world_comm ) 
  CALL mp_bcast(lmoire, ionode_id, world_comm)
  CALL mp_bcast(amoire_in_ang, ionode_id, world_comm)
  CALL mp_bcast(vmoire_in_mev, ionode_id, world_comm)
  CALL mp_bcast(pmoire_in_deg, ionode_id, world_comm)
  CALL mp_bcast(mstar, ionode_id, world_comm)
  CALL mp_bcast(epsmoire, ionode_id, world_comm)
  !
  ! MAM: Problematic situation, I need to modify qnorm before init_us_1 is 
  !      called in read_file, but to modify qnorm accurately I need to know
  !      the list of kpoints
  ! from exx_base.f90
  !  qnorm = 0.0_dp
  !  DO iq = 1,nkqs
  !     DO ik = 1,nks
  !        qnorm = max(qnorm, sqrt( sum((xk(:,ik)-xkq_collect(:,iq))**2) ))
  !     ENDDO
  !  ENDDO
  ! guessing now
  qnorm = 3.5d0  
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
  if(qnorm_ > qnorm) call errore('posthf: qnorm guess is wrong.',1) 
  !
  CALL openfil_pp
  !
  ! passing arguments is becoming a real problem, make a module and pass an
  ! "options" object
  CALL pp_posthf(out_prefix,number_of_orbitals, expand_kp, thresh, &
       eigcut, eigcut_occ, nextracut, ncholmax, ndet, nskipvir, run_type, diag_type, write_psir, update_qe_bands, &
       h5_add_orbs, read_from_h5, get_hf, get_mp2, use_symm, exxdiv_treatment, regp, regkappa, &
       low_memory, verbose, debug)
  !
#else
!#error HDF5 flag neither enabled during configure nor added manually in make.inc
  CALL errore('posthf', ' HDF5 flag not enabled during configure',1)
#endif
  CALL environment_end ( 'posthf' )
  CALL stop_pp
  STOP
END PROGRAM posthf

! Note: Hartrees are used throughout this routine!!!
SUBROUTINE pp_posthf(out_prefix, norb_, expand_kp, thresh, eigcut, &
                          occeigcut, nextracut,  &
                          ncmax, ndet, nskipvir, run_type, diag_type, &
                          write_psir, update_qe_bands, h5_add_orbs, read_from_h5, &
                          get_hf, get_mp2, use_symm, exxdiv_treatment, regp, regkappa, low_memory, &
                          verbose, debug) 

  USE kinds, ONLY: DP
  USE ions_base, ONLY : nat, nsp, ityp, tau, zv, atm
  USE cell_base, ONLY: omega, alat, tpiba, tpiba2, at, bg
  USE constants, ONLY: tpi, fpi, e2
  USE run_info,  ONLY: title
  USE gvect, ONLY: ngm, ngm_g, g, gstart, gg
  USE vlocal,               ONLY : strf
  USE klist , ONLY: nks, nelec, nelup, neldw, xk, nkstot
  USE lsda_mod, ONLY: lsda, nspin, isk
  USE scf, ONLY: rho, rho_core, rhog_core, vnew
  USE wvfct, ONLY: npw, npwx, nbnd, g2kin, et, wg
  USE klist, ONLY: igk_k, ngk, wk
  USE gvecw, ONLY : ecutwfc
  USE gvecs, ONLY : doublegrid
  USE control_flags, ONLY: gamma_only
  USE io_global, ONLY: stdout, ionode,  ionode_id
  USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc, iun => iunsat, tmp_dir, prefix
  USE wavefunctions, ONLY : psic
  USE noncollin_module,     ONLY : noncolin, npol
  use scatter_mod,          ONLY : gather_grid, scatter_grid 
  use fft_base,             ONLY : dffts
  use fft_interfaces,       ONLY : invfft, fwfft
  USE symm_base,            ONLY : nsym, s
  USE ldaU,     ONLY : lda_plus_u

  USE mp,           ONLY: mp_sum, mp_max, mp_bcast, mp_barrier
  USE mp_world,     ONLY: world_comm, nproc, mpime, root
  USE mp_pools,     ONLY: inter_pool_comm, intra_pool_comm, npool, nproc_pool, &
                          me_pool,root_pool,my_pool_id
  USE mp_global, ONLY: intra_image_comm, me_image, root_image, nproc_image
  USE mp_bands,     ONLY: nproc_bgrp
  USE realus,   ONLY : real_space
  USE paw_variables, ONLY : okpaw
  USE uspp,       ONLY : okvan
  USE control_flags, ONLY : tqr
  USE posthf_mod, ONLY: nksym, nkfull, norb, xkfull, xksym, &
                        igksym, ngksym, kminus, QKtoK2, numspin, &
                        nelec_tot, nup, ndown, e0, efc, &
                        init_posthf, clean_posthf, ke_factorization, &
                        nmax_DM,DM,DM_mf 
  USE read_orbitals_from_file, ONLY: h5file_type,update_bands,  &
                        open_esh5_read, close_esh5_read, get_orbitals_set
  USE orbital_generators, ONLY: mixed_basis, generate_orbitals, &
                                get_spanning_basis, write_trial_wavefunction, &
                                orb_verbose => verbose, diag_hf,&
                                davcio_to_esh5
  USE onebody_hamiltonian, ONLY: getH1
  USE twobody_hamiltonian, ONLY: cholesky_r, calculate_KS_bscorr, &
                             chol_verbose => verbose
  USE mp2_module, ONLY: mp2_g,mp2no,approx_mp2no,mp2_verbose => verbose  
#if defined(__CUDA)
  USE mp2_module, ONLY: mp2_gpu
#endif

  IMPLICIT NONE
  !  
  INTEGER, INTENT(IN) :: norb_, regp
  INTEGER, INTENT(INOUT) :: ndet, nskipvir, read_from_h5
  REAL(DP), INTENT(IN) :: ncmax, thresh, regkappa
  REAL(DP), INTENT(INOUT) :: occeigcut, eigcut, nextracut
  LOGICAL, INTENT(IN) :: expand_kp, debug, verbose, write_psir, use_symm
  LOGICAL, INTENT(IN) :: low_memory, update_qe_bands, get_hf, get_mp2
  CHARACTER(len=256), INTENT(IN) :: h5_add_orbs, run_type, diag_type, exxdiv_treatment
  INTEGER :: ibnd, ik, j
  INTEGER :: ios, ierr, h5len,oldh5,ig_c,save_complex
  CHARACTER(256) :: tmp,tmp2,h5name,out_prefix,h5qeorbs
  INTEGER :: error
  REAL(DP) :: lowcut, highcut
  COMPLEX(DP) :: CONE,CZERO,CMINUSONE
  INTEGER :: i, nxxs, maxnorb, nelmax
  COMPLEX(DP), ALLOCATABLE :: M(:,:,:,:)  ! Overlap between basis states and
                                            ! occupied KS states. 
  COMPLEX(DP) :: e1,e1_so,e1_mf,e1_so_mf,emp2
  TYPE(h5file_type) :: h5id_input_orbs, h5id_output_orbs, h5id_hamil
! **********************************************************************

  if(real_space) &
    call errore('posthf','real_space not allowed',1)
  if(tqr .and. nkstot>1) &
    call errore('posthf','tqr with nkpots>1 not working',1)

  if(ionode) then
    if(okvan) write(*,*) ' Found USPP. '
    if(okpaw) write(*,*) ' Found PAW. '
  endif

  if(nextracut < 0) nextracut = 1.d-4

  if(nproc_bgrp > 1) then
    write(*,*) ' Error: Only pool and band parallelization allowed.'
    call errore('posthf','Error: Only pool and band parallelization allowed.',1)
  endif

  if(doublegrid) then
    write(*,*) ' Error: doublegrid not implemented yet ' 
    call errore('posthf','Error: doublegrid not implemented',1)
  endif

  CONE  = (1.d0,0.d0)
  CMINUSONE  = (-1.d0,0.d0)
  CZERO = (0.d0,0.d0)
  lowcut = 0.2d0
  highcut = 0.8d0
  chol_verbose = verbose
  orb_verbose = verbose
  mp2_verbose = verbose

  ! initialize variables in posthf_module
  call init_posthf(norb_,symm=use_symm,exxdiv_treat=exxdiv_treatment)

  ! this routine is designed to split over bands and kpoints, no PW
  ! parallelization is allowed here. A different routine will be written for PW 
  ! parallelization
  nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x   ! global size of FFT grid
  write(*,*) 'nr1, nr2, nr3, nnr:',dffts%nr1,dffts%nr2,dffts%nr3,dffts%nnr
  write(*,*) 'ngm: ',dffts%ngm
  write(*,*) 'nkstot, nksym, nkfull: ',nkstot,nksym,nkfull
  write(*,*) 'nspin (QE), nspin (posthf): ',nspin, numspin

  if( nxxs .ne. dffts%nnr ) then
    CALL errore('posthf', ' nxxs .ne. nnr ', 1) 
  endif

  if( TRIM(run_type) == 'one_body' ) then

    ! right now one_body is serial, so make esh5 file in parallel 
    if(npool > 1) then
      h5qeorbs = TRIM( tmp_dir )//TRIM( out_prefix ) // '.qe.orbs.h5'
      call davcio_to_esh5(h5qeorbs,norb,dffts) 
    endif

    if(me_image == root_image) then
      if(npool > 1) then
        CALL generate_orbitals(out_prefix,dffts,read_from_h5,h5_add_orbs, &
              eigcut,occeigcut,nextracut,write_psir,esh5_orbs=h5qeorbs)
      else
        CALL generate_orbitals(out_prefix,dffts,read_from_h5,h5_add_orbs, &
              eigcut,occeigcut,nextracut,write_psir)
      endif
      
      ! open hdf5 file 
      h5name = TRIM( tmp_dir )//TRIM( out_prefix ) // '.hamil.h5'
      h5len = LEN_TRIM(h5name)
      oldh5=0
      CALL esh5_posthf_open_file(h5id_hamil%id,h5name,h5len,oldh5)

      ! all tasks open orbital file in read_only mode
      h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
      call open_esh5_read(h5id_input_orbs,h5name)
      if( h5id_input_orbs%grid_type .ne. 1) &
        call errore('posthf','grid_type ne 1',1)

      ! calculates expansion of occupied KS states in the generated  
      ! (spin-independent) basis and writes trial wfn to file
      call get_nelmax(nbnd,nksym,numspin,wk,size(wg,1),wg,nelmax)  
      maxnorb = maxval(h5id_input_orbs%norbK(:))  
      write(*,*)'maxnorb, nelmax:',maxnorb, nelmax
      allocate( M(maxnorb,nelmax,nksym,min(2,nspin)) )
      if(npool > 1) then
        call write_trial_wavefunction(h5id_input_orbs,h5id_hamil,&
                    dffts,nelmax,M,ndet,lowcut,highcut,h5wfn=h5qeorbs)
      else
        call write_trial_wavefunction(h5id_input_orbs,h5id_hamil,&
                    dffts,nelmax,M,ndet,lowcut,highcut)
      endif
  
      ! calculate and write 1 body hamiltonian
      ! make DM and optional argument to H1
      CALL getH1(h5id_input_orbs,h5id_hamil,dffts,e1,e1_so,e1_mf,e1_so_mf)

      CALL esh5_posthf_close_file(h5id_hamil%id)
      call close_esh5_read(h5id_input_orbs)


      h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
      h5len = LEN_TRIM(h5name)
      oldh5 = 1
      CALL esh5_posthf_open_file(h5id_hamil%id,h5name,h5len,oldh5)

      call esh5_posthf_write_orbmat(h5id_hamil%id,"SCFOrbMat",9, &
            maxnorb,nelmax,nksym,size(M,4),M,error)
      if(error .ne. 0 ) &
        call errore('posthf','error writing orbmat',1)
      call esh5_posthf_write_dm(h5id_hamil%id,"DM",2,nmax_DM,nksym,nspin,DM,error)
      if(error .ne. 0 ) &
        call errore('posthf','error writing DM',1)
      call esh5_posthf_write_dm(h5id_hamil%id,"DM_mf",5,nmax_DM,nksym,nspin,DM_mf,error)
      if(error .ne. 0 ) &
        call errore('posthf','error writing DM_mf',1)
      write(*,*) 'E0, E1(1Det), E_SO(1Det) (Ha):',e0,e1/(nkfull*1.0),e1_so/(nkfull*1.0)
      write(*,*) '    E1, E_SO (Ha):',e1_mf/(nkfull*1.0),e1_so_mf/(nkfull*1.0)


      CALL esh5_posthf_close_file(h5id_hamil%id)
      if(allocated(M)) deallocate(M)

    endif

    if(update_qe_bands) then
      h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
      call update_bands(h5name,dffts) 
    endif

  elseif(TRIM(run_type) == 'hf_diag_vir') then

    ! orbitals are read from esh5 only right now
    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.hamil.h5'
    tmp = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
    tmp2 = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
    call diag_hf(dffts,'vir_update_occ_et',h5name,tmp,tmp,tmp2)
    if(update_qe_bands) call update_bands(tmp2,dffts)

  elseif(TRIM(run_type) == 'hf_diag_full') then

    ! orbitals are read from esh5 only right now
    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.hamil.h5'
    tmp = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
    tmp2 = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
    call diag_hf(dffts,'full',h5name,tmp,tmp,tmp2)
    if(update_qe_bands) call update_bands(tmp2,dffts)

  elseif(TRIM(run_type) == 'mp2') then

#if defined(__CUDA)
    call mp2_gpu(emp2,dffts,reg_expo=regkappa,reg_pow=regp)
#else
    call mp2_g(emp2,dffts,reg_expo=regkappa,reg_pow=regp)
#endif

  elseif(TRIM(run_type) == 'mp2_esh5') then

    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.canonical.orbitals.h5'
#if defined(__CUDA)
    call mp2_gpu(emp2,dffts,h5name,reg_expo=regkappa,reg_pow=regp)
#else
    call mp2_g(emp2,dffts,h5name,reg_expo=regkappa,reg_pow=regp)
#endif

  elseif(TRIM(run_type) == 'mp2no') then

    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.mp2no.h5'
    call mp2no(h5name,dffts,nskipvir,eigcut)
    if(update_qe_bands) call update_bands(h5name,dffts)

  elseif(TRIM(run_type) == 'appmp2no') then

    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.appmp2no.h5'
    call approx_mp2no(h5name,dffts,nskipvir,eigcut,low_memory)
    if(update_qe_bands) call update_bands(h5name,dffts)

  elseif(TRIM(run_type) == 'cholesky')  then !  

    ! calculate cholesky matrix  
    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.hamil.h5'
    tmp = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
    call cholesky_r(ncmax,thresh,dffts,h5name,tmp)

  elseif(TRIM(run_type) == 'bscorr') then

    h5name = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.hamil.h5'
    tmp = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
    call calculate_KS_bscorr(dffts,h5name,tmp)

  ! now can call pyscf drivers from here
  elseif(TRIM(run_type) == 'hamil') then  

    call pyscf_driver_hamil(out_prefix, read_from_h5, h5_add_orbs, &
       ndet, eigcut, nextracut, thresh, ncmax, &
       get_hf, get_mp2, update_qe_bands, e1, emp2)
    if(get_hf) &
      write(*,*) 'E0, E1, E1_SO (Ha):',e0,e1/(nkfull*1.0),e1_so/(nkfull*1.0)
    if(get_mp2) &
      write(*,*) 'EMP2 (Ha):',emp2

  elseif(TRIM(run_type) == 'mp2_fullpw') then

    call pyscf_driver_mp2(out_prefix,.true.,'fullpw',0,'',0.d0,0.d0,&
         regkappa,regp,emp2)
    write(*,*) 'EMP2 (Ha):',emp2

  elseif(TRIM(run_type) == 'mp2_driver') then

    call pyscf_driver_mp2(out_prefix,.true.,TRIM(diag_type),read_from_h5,h5_add_orbs,&
            eigcut,nextracut,regkappa,regp,emp2)
    write(*,*) 'EMP2 (Ha):',emp2

  else
    call errore('posthf','Error: Unknown run_type',1)
  endif

END SUBROUTINE pp_posthf

