!--------------------------------------------------------------------
! Written by Miguel A. Morales, LLNL, 2020 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE orbital_generators 
  !----------------------------------------------------------------------
  ! 
  ! Routines used to generate spanning spin-independent basis sets. 
  !
  USE kinds, ONLY : DP
  USE posthf_mod, ONLY: nksym,norb,numspin,xksym,ngksym,igksym,&
                nelec_tot,nup,ndown,nmax_DM,DM,DM_mf,e2Ha
  USE read_orbitals_from_file, ONLY: h5file_type, &
                        open_esh5_read, close_esh5_read,&
                        open_esh5_write,close_esh5_write, &
                        get_orbitals,get_orbitals_set
  USE wavefunctions, ONLY : psic
  USE control_flags, ONLY : gamma_only
  USE cell_base, ONLY: tpiba2
  USE gvect, ONLY: ngm, g, gstart
  USE gvecw, ONLY : ecutwfc
  USE io_files, ONLY: nwordwfc, iunwfc, tmp_dir, prefix
  USE wvfct, ONLY: nbnd, npwx, wg, et, g2kin
  USE klist,  ONLY : nkstot, wk, nks
  USE becmod,   ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE mp_images, ONLY: intra_image_comm, me_image, root_image
  USE noncollin_module,     ONLY : noncolin, npol
  USE lsda_mod, ONLY: lsda, nspin
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: mixed_basis, verbose
  CHARACTER(6), PARAMETER  :: orbsG = 'OrbsG'
  CHARACTER(6), PARAMETER  :: orbsR = 'OrbsR'
  !
  CONTAINS
  !
  ! Generates spin-independent basis.
  !  
  SUBROUTINE generate_orbitals(out_prefix,dfft,nread_from_h5_,h5_add_orbs,&
            eigcut,occeigcut,nextracut,write_psir,esh5_orbs)
    !
    USE cell_base, ONLY: alat, tpiba, at, bg
    USE uspp,     ONLY : vkb, nkb
    USE paw_variables, ONLY : okpaw
    USE uspp,       ONLY : okvan
    USE fft_types, ONLY: fft_type_descriptor
    USE io_global, ONLY: ionode
    ! 
    IMPLICIT NONE
    !
    CHARACTER(len=256), INTENT(IN) :: out_prefix, h5_add_orbs
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: nread_from_h5_
    REAL(DP), INTENT(IN) :: eigcut, nextracut, occeigcut
    LOGICAL, INTENT(IN) :: write_psir 
    CHARACTER(len=256), INTENT(IN), OPTIONAL :: esh5_orbs
    !
    CHARACTER(len=256) :: tmp,ftype
    INTEGER :: i, j, k, ik, ka, ikk, ie, ia, ib, n, is1, is2, npw2
    INTEGER :: minnorb, maxnextra, error, ir, nnr_
    INTEGER :: lda, nxxs, h5len, ibnd, npw
    INTEGER :: n1, n2, n3, nocc, nvir, nread_from_h5
    COMPLEX(DP) :: ctemp
    REAL(DP) :: recv(3,3), at0(3,3)
    REAL(DP) :: rtemp, scl
    REAL(DP), ALLOCATABLE :: eig(:)   ! eigenvalues of overlap matrix 
    COMPLEX(DP), ALLOCATABLE :: S(:,:)  ! Overlap between basis states and
                                        ! occupied KS states. 
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
    COMPLEX(DP), ALLOCATABLE :: Orbs(:,:,:) ! old orbitals 
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:) 
    INTEGER, ALLOCATABLE :: norb_extra(:)
    INTEGER, ALLOCATABLE :: g2g(:)
    COMPLEX(DP), ALLOCATABLE :: orb_extra(:)  
    COMPLEX(DP), ALLOCATABLE :: virtuals(:,:)  
    COMPLEX(DP), ALLOCATABLE :: Psir(:,:)  
    REAL(DP), ALLOCATABLE :: eigval(:),weights(:),gk(:,:)
    !
    TYPE(h5file_type) :: h5id_add_orbs, h5id_input_orbs, h5id_output_orbs
    !

    ! only roots of each pool generate orbitals 
    if(me_image .ne. root_image) return
    if(ionode) call print_freemem()

    nxxs = dfft%nr1x*dfft%nr2x*dfft%nr3x

    tmp = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
    call open_esh5_write(h5id_output_orbs,dfft,tmp,write_psir)

    ! setup h5id_input_orbs if reading from esh5 
    ftype = 'davcio'
    if( present(esh5_orbs) ) then
      ! assuming file is correct for now
      call open_esh5_read(h5id_input_orbs,TRIM(esh5_orbs)) 
      ftype = 'esh5'
    endif

    nread_from_h5 = nread_from_h5_
    maxnextra=0
#if defined(__HDF5) || defined(__HDF5_C)
    if( LEN(TRIM(h5_add_orbs)) > 0 ) then

      if(okpaw .or. okvan) call errore('generate_orbitals','h5_add_orbs with paw/uspp',1)  

      write(*,*) ' Reading additional orbitals '
      call open_esh5_read(h5id_add_orbs,TRIM(h5_add_orbs))
      allocate(norb_extra(nksym))
      norb_extra(1:nksym) = h5id_add_orbs%norbK(1:nksym)  
      do i=1,nksym
        maxnextra = max(maxnextra,norb_extra(i)) 
      enddo  
      write(*,*) 'Adding (maxnextra)',maxnextra,'  orbitals in files.'
      if(nread_from_h5 <= 0) nread_from_h5 = maxnextra 

      if(h5id_add_orbs%nr1 > dfft%nr1 .or. &
         h5id_add_orbs%nr2 > dfft%nr2 .or. &
         h5id_add_orbs%nr3 > dfft%nr3) &
        call errore('generate_orbitals',' error: grid dimensions in file too large',1) 

      if(h5id_add_orbs%grid_type == 0) then
        nnr_ = h5id_add_orbs%nr1*h5id_add_orbs%nr2*h5id_add_orbs%nr3
        allocate(orb_extra(nnr_), g2g(nnr_))
        n=1  
        do k=0,h5id_add_orbs%nr3-1  
          n3 = k
          if( n3 > (h5id_add_orbs%nr3-1)/2 ) n3 = n3 - h5id_add_orbs%nr3 + dfft%nr3
          do j=0,h5id_add_orbs%nr2-1
            n2 = j  
            if( n2 > (h5id_add_orbs%nr2-1)/2 ) n2 = n2 - h5id_add_orbs%nr2 + dfft%nr2
            do i=0,h5id_add_orbs%nr1-1
              n1 = i  
              if( n1 > (h5id_add_orbs%nr1-1)/2 ) n1 = n1 - h5id_add_orbs%nr1 + dfft%nr1
              ir = 1 + n1 + n2 * dfft%nr1x + n3 * dfft%nr1x * dfft%nr2x
              g2g(n) = ir 
              n = n + 1
            enddo  
          enddo  
        enddo  
      elseif(h5id_add_orbs%grid_type==1) then
        nnr_ = h5id_add_orbs%npwx
        if( nnr_ .gt. npwx .or. nnr_ < 1) then
          write(*,*) ' In file npwx > npwx:',nnr_,npwx 
          call errore('generate_orbitals','npwx(file) > npwx:',1)
        endif 
        allocate(orb_extra(nnr_), g2g(1))
      else
        write(*,*) ' Unknown grid_type (1):',h5id_add_orbs%grid_type
        call errore('generate_orbitals',' Unknown grid_type:',1)
      endif
        
    endif
#endif
    if(norb == 0) then
      ! just add orbitals from file

      if( maxnextra == 0 ) &
        call errore('generate_orbitals',' norb==0 and maxnextra==0',1) 

      mixed_basis = .true.
      lda = maxnextra
      if(okvan .or. okpaw) then
        allocate( spsi(npol*npwx, lda) )
        CALL allocate_bec_type ( nkb, lda, becp )
      else
        allocate( spsi(1, 1) )
      endif

      allocate( Orbitals(npwx,lda), S(lda,lda) )

      do ik=1,nksym
        !
        Orbitals(:,:) = (0.d0,0.d0)
        !
!        if(okvan .or. okpaw .or. write_psir) &
            CALL gk_sort (xksym (1:3, ik), ngm, g, ecutwfc / tpiba2, &
                  npw, igksym(1), g2kin)
        !
        if(okvan .or. okpaw) CALL init_us_2 (ngksym(ik), igksym(1), xksym (1, ik), vkb)
        !
        h5id_output_orbs%norbK(ik) = 0 
        !
        call add_extra_orbitals_from_file(h5id_add_orbs, ik, h5id_output_orbs%norbK(ik), & 
            min(norb_extra(ik),nread_from_h5), Orbitals, &
            orb_extra, S(:,1), spsi, dfft, g2g, igksym, nextracut )
        !
        ! check overlaps
        call Overlap(h5id_output_orbs%norbK(ik),h5id_output_orbs%norbK(ik),npwx,&
                            (1.d0,0.d0),Orbitals(1,1),npwx, &
                            Orbitals(1,1),npwx,(0.d0,0.d0),S(1,1),lda,.true.,spsi)
        do ia=1,h5id_output_orbs%norbK(ik)
          if( abs(S(ia,ia)-1.d0) > 1.d-7 ) then
            write(*,*) 'WARNING: Non-orthogonal orbitals (diag) ik,ia,S:',&
                            ik,ia,S(ia,ia)
          endif
          do ib=ia+1,h5id_output_orbs%norbK(ik)
            if( abs(S(ia,ib)) > 1.d-7 .or. abs(S(ib,ia)) > 1.d-7) then
              write(*,*) 'WARNING: Non-orthogonal orbitals (offdiag) ik,ia,ib,S:',&
                            ik,ia,ib,S(ia,ib),S(ib,ia)
            endif

          enddo
        enddo
        !
        if( h5id_output_orbs%norbK(ik) == 0 ) &
          call errore('generate_orbitals',' Did not find basis states. (norbK(ik)==0) ',1) 
        !
        if(gamma_only .AND. gstart == 2 ) &
          Orbitals(1,:) = CMPLX( DBLE( Orbitals(1,:) ), 0.D0 ,kind=DP)
        call esh5_posthf_write(h5id_output_orbs%id,orbsG, 5, ik-1, npwx, &
                        h5id_output_orbs%norbK(ik), Orbitals, npwx, error)
        if(error .ne. 0 ) &
          call errore('pw2postha','error writing orbital',1)
        !
        if(write_psir) then
          if(noncolin) &
            call errore('generate_orbitals','noncolin not yet allowed with write_psir',1)
          allocate( Psir(nxxs,h5id_output_orbs%norbK(ik)) )
          if(okvan .or. okpaw) call errore('missing aug charge!!!',1)
          call orbitals_g2r(ngksym(ik),igksym(1),dfft,nxxs,h5id_output_orbs%norbK(ik), &
                Orbitals,Psir)
          ! MAM: this is wrong if nrix != nri, since the meta data has nri as 
          !      the grid dims
          call esh5_posthf_write(h5id_output_orbs%id,orbsR, 5, ik-1, nxxs, &
                h5id_output_orbs%norbK(ik), Psir, nxxs, error)
          if(error .ne. 0 ) &
            call errore('generate_orbitals','error writing orbital in real space',1)
          deallocate(Psir)
        endif
        !
      enddo
      !
      if(allocated(Orbitals)) deallocate(Orbitals)
      if(allocated(S)) deallocate(S)

    elseif(nspin.gt.1) then

      if( noncolin .and. nspin.ne.4 ) &
        call errore('generate_orbitals','noncolin and nspin!=4 (??)',1)
      if( noncolin .and. npol.ne.2 ) &
        call errore('generate_orbitals','noncolin and npol!=2 (??)',1)
      if( noncolin .and. norb < nelec_tot .and. ionode ) &
        write(*,*) ' ******* Warning: norb < nelec in noncolinear calculation. ********'
      if( lsda .and. nspin.ne.2 ) &
        call errore('generate_orbitals','lsda and nspin!=2 (??)',1)
      if( lsda .and. npol.ne.1 ) &
        call errore('generate_orbitals','lsda and npol!=1 (??)',1)
      mixed_basis = .true.
      lda = 2*norb+maxnextra
      if(okvan .or. okpaw) then
        allocate( spsi(npol*npwx, lda) )
        CALL allocate_bec_type ( nkb, lda, becp )
      else
        allocate( spsi(1, 1) )
      endif

      ! Orbitals will contain the spin-dependent basis
      ! Orbs is used to read orbitals from file.
      allocate( Orbitals(npwx,lda), Orbs(npol*npwx,norb,2) )

      h5id_output_orbs%norbK(1:nksym) = 0 

      write(*,*) 'Generating mixed basis set.' 

      ! Working arrays for overlap and eigenvalues
      allocate( S(lda,lda), eig(lda) )

      do ik=1,nksym
        !
        Orbitals(:,:) = (0.d0,0.d0)
        !
        npw = ngksym(ik)
        !
!        if(okvan .or. okpaw .or. write_psir) &
            CALL gk_sort (xksym (1:3, ik), ngm, g, ecutwfc / tpiba2, &
                  npw, igksym(1), g2kin)
        !
        if(okvan .or. okpaw) CALL init_us_2 (npw, igksym(1), xksym (1, ik), vkb)
        !
        call get_orbitals(h5id_input_orbs,ftype,'psig',dfft,Orbs(:,:,1), &
                    1,norb,ik,1)
        Orbs(npw+1:npwx,:,1) = (0.d0,0.d0)
        if(noncolin) Orbs(npwx+npw+1:2*npwx,:,1) = (0.d0,0.d0)
        nocc=0
        scl = 1.d0
        if(abs(wk(ik))>1.d-10) scl = 1.d0/wk(ik)
        do ia=1,nbnd
          if(abs(wg(ia,ik)*scl) > 0.01d0) nocc = ia
        enddo
        !
        if(lsda) then
          ikk = ik + nksym
          npw2 = ngksym(ikk)
          call get_orbitals(h5id_input_orbs,ftype,'psig',dfft,Orbs(:,:,2), &
                    1,norb,ik,2)
          Orbs(npw2+1:npwx,:,2) = (0.d0,0.d0)
          scl = 1.d0
          if(abs(wk(ikk))>1.d-10) scl = 1.d0/wk(ikk)
          do ia=1,nbnd
            if(abs(wg(ia,ikk)*scl) > 0.01d0 .and. ia .gt. nocc) nocc = ia
          enddo
        endif
        !
        if(norb < nocc) &
            call errore('generate_orbitals',  &
                        'Error: # occupied states > number_of_orbitals. Increase orbitals.',1)
        ! 
        ! add occupied states with high small cutoff
        !
        h5id_output_orbs%norbK(ik) = 0
        call get_spanning_basis(ik, nocc, h5id_output_orbs%norbK(ik), Orbs(1,1,1), &
                npol*npwx, Orbs(1,1,2), npol*npwx, Orbitals, npwx, S, lda, spsi, eig, occeigcut)
        if(verbose) write(*,*) '# orbitals after adding occupied only:',nocc,&
                h5id_output_orbs%norbK(ik)
        ! 
        ! add virtual states with requested cutoff
        !
        if( norb > nocc ) then
          call get_spanning_basis(ik, norb-nocc, nvir, Orbs(1,nocc+1,1),    & 
            npol*npwx,Orbs(1,nocc+1,2),npol*npwx, &
            Orbitals(:,h5id_output_orbs%norbK(ik)+1),npwx,S,lda,spsi,eig,eigcut)
          call add_extra_orbitals(h5id_output_orbs%norbK(ik), nvir, Orbitals, &
                    S(:,1), spsi, nextracut ) 
          if(verbose) &
            write(*,*) '# orbitals after adding virtual:',nvir,h5id_output_orbs%norbK(ik)
        endif
        !
        ! Add second set of virtuals here if desired!!!
        !
        ! Add states from file if requested
        if(maxnextra > 0) &
            call add_extra_orbitals_from_file(h5id_add_orbs, ik,h5id_output_orbs%norbK(ik), & 
                min(norb_extra(ik),nread_from_h5), Orbitals, &
                orb_extra, S(:,1), spsi, dfft, g2g, igksym, nextracut ) 

        ! check overlaps
        call Overlap(h5id_output_orbs%norbK(ik),h5id_output_orbs%norbK(ik),npwx,&
                (1.d0,0.d0),Orbitals(1,1),npwx, &
                Orbitals(1,1),npwx,(0.d0,0.d0),S(1,1),lda,.true.,spsi)
        do ia=1,h5id_output_orbs%norbK(ik)
          if( abs(S(ia,ia)-1.d0) > 1.d-7 ) then
            write(*,*) 'WARNING: Non-orthogonal orbitals (diag) ik,ia,S:',&
                            ik,ia,S(ia,ia)
          endif
          do ib=ia+1,h5id_output_orbs%norbK(ik)
            if( abs(S(ia,ib)) > 1.d-7 .or. abs(S(ib,ia)) > 1.d-7) then
              write(*,*) 'WARNING: Non-orthogonal orbitals (offdiag) ik,ia,ib,S:',&
                            ik,ia,ib,S(ia,ib),S(ib,ia)
            endif
            
          enddo
        enddo
        !
        if(gamma_only .AND. gstart == 2 ) &
          Orbitals(1,:) = CMPLX( DBLE( Orbitals(1,:) ), 0.D0 ,kind=DP)
        call esh5_posthf_write(h5id_output_orbs%id,orbsG, 5, ik-1, npwx, &
            h5id_output_orbs%norbK(ik), Orbitals, npwx, error)
        if(error .ne. 0 ) & 
          call errore('generate_orbitals','error writing orbital',1)
        if(write_psir) then
          if(okvan .or. okpaw) call errore('missing aug charge!!!',1)
          allocate( Psir(nxxs,h5id_output_orbs%norbK(ik)) )
          call orbitals_g2r(npw,igksym(1),dfft,nxxs,h5id_output_orbs%norbK(ik), &
                Orbitals,Psir)
          ! MAM: this is wrong if nrix != nri, since the meta data has nri as 
          !      the grid dims
          call esh5_posthf_write(h5id_output_orbs%id,orbsR, 5, ik-1, nxxs, &
                    h5id_output_orbs%norbK(ik), Psir, nxxs, error)
          if(error .ne. 0 ) &
            call errore('generate_orbitals','error writing orbital in real space',1)
          deallocate(Psir)
        endif
        !
      enddo
      !  
      deallocate( S, eig )  
      if(allocated(Orbitals)) deallocate(Orbitals)
      if(allocated(Orbs)) deallocate(Orbs)

    else

      if( nspin.ne.1 ) call errore('generate_orbitals','expect nspin==1 (??)',1)
      mixed_basis = (maxnextra > 0)
      lda = norb+maxnextra
      if(okvan .or. okpaw) then
        allocate( spsi(npol*npwx, lda) )
        CALL allocate_bec_type ( nkb, lda, becp )
      endif

      allocate( Orbitals(npwx,lda), S(lda,lda), eigval(lda), weights(lda) )

      do ik=1,nksym
        !
        Orbitals(:,:) = (0.d0,0.d0)
        !
        npw = ngksym(ik)
        !
!        if(okvan .or. okpaw .or. write_psir) &
          CALL gk_sort (xksym (1:3, ik), ngm, g, ecutwfc / tpiba2, &
                npw, igksym(1), g2kin)
        !
        if(okvan .or. okpaw) CALL init_us_2 (npw, igksym(1), xksym (1, ik), vkb)
        !
        call get_orbitals(h5id_input_orbs,ftype,'psig',dfft,Orbitals(:,:), &
                  1,norb,ik,1)
        ! 
        h5id_output_orbs%norbK(ik) = norb
        !
        ! Add states from file if requested
        if(maxnextra > 0) &
            call add_extra_orbitals_from_file(h5id_add_orbs, ik,h5id_output_orbs%norbK(ik),&
                min(norb_extra(ik),nread_from_h5), &
                Orbitals, orb_extra, S(:,1), spsi, dfft, g2g, igksym, nextracut )
        !
        ! check overlaps
        call Overlap(h5id_output_orbs%norbK(ik),h5id_output_orbs%norbK(ik),npwx, &
                (1.d0,0.d0),Orbitals(1,1),npwx, &
                Orbitals(1,1),npwx,(0.d0,0.d0),S(1,1),lda,.true.,spsi)
        do ia=1,h5id_output_orbs%norbK(ik)
          if( abs(S(ia,ia)-1.d0) > 1.d-7 ) then
            write(*,*) 'WARNING: Non-orthogonal orbitals (diag) ik,ia,S:',&
                            ik,ia,S(ia,ia)
          endif
          do ib=ia+1,h5id_output_orbs%norbK(ik)
            if( abs(S(ia,ib)) > 1.d-7 .or. abs(S(ib,ia)) > 1.d-7) then
              write(*,*) 'WARNING: Non-orthogonal orbitals (offdiag) ik,ia,ib,S:',&
                            ik,ia,ib,S(ia,ib),S(ib,ia)
            endif

          enddo
        enddo
        !
        if(gamma_only .AND. gstart == 2 ) &
          Orbitals(1,:) = CMPLX( DBLE( Orbitals(1,:) ), 0.D0 ,kind=DP)
        call esh5_posthf_write(h5id_output_orbs%id, orbsG, 5, ik-1, npwx, &
            h5id_output_orbs%norbK(ik), Orbitals, npwx, error)
        if(error .ne. 0 ) &
          call errore('generate_orbitals','error writing orbital',1)
        if(write_psir) then
          if(okvan .or. okpaw) call errore('missing aug charge!!!',1)
          allocate( Psir(nxxs,h5id_output_orbs%norbK(ik)) )
          call orbitals_g2r(npw,igksym(1),dfft,nxxs,h5id_output_orbs%norbK(ik),&
                Orbitals,Psir)
          ! MAM: this is wrong if nrix != nri, since the meta data has nri as 
          !      the grid dims
          call esh5_posthf_write(h5id_output_orbs%id, orbsR, 5, ik-1, nxxs, &
                h5id_output_orbs%norbK(ik), Psir, nxxs, error)
          if(error .ne. 0 ) &
            call errore('generate_orbitals','error writing orbital in real space',1)
          deallocate(Psir)
        endif
        !
        weights(:) = 0.d0
        eigval(:) = 0.d0
        ! setting values beyond norb to zero
        eigval(1:norb)=et(1:norb,ik)*e2Ha
        scl = 1.d0
        if(abs(wk(ik))>1.d-10) scl = 1.d0/wk(ik)
        weights(1:norb) = wg(1:norb,ik)*scl
        call esh5_posthf_write_et(h5id_output_orbs%id, orbsG, 5, ik-1, & 
            h5id_output_orbs%norbK(ik),eigval, weights, error)
        if(error .ne. 0 ) &
          call errore('generate_orbitals','error writing eigenvalues',1)
        !
      enddo
      !
      if(allocated(Orbitals)) deallocate(Orbitals)
      if(allocated(S)) deallocate(S)
      !
    endif

    !  
    write(*,*) ' ik, # orbitals: '
    do ik=1,nksym  
      write(*,*) ik,h5id_output_orbs%norbK(ik)
    enddo
    !
    ! write gvectors
    allocate(gk(3,npwx))
    write(*,*) ' ik, npw:'
    do ik=1,nksym
      npw = ngksym(ik)
      write(*,*) ik, npw
      call gk_sort(xksym(1:3, ik), ngm, g, ecutwfc/tpiba2, &
                   npw, igksym(1), g2kin)
      gk(:,1:npw) = g(:, igksym(1:npw))
      call esh5_posthf_write_g(h5id_output_orbs%id, gk, npw, ik-1, error)
    enddo
    deallocate(gk)

    if( LEN(TRIM(h5_add_orbs)) > 0 ) then
      if(allocated(norb_extra)) deallocate(norb_extra)
      if(allocated(orb_extra)) deallocate(orb_extra)
      if(allocated(g2g)) deallocate(g2g)
    endif

    if( present(esh5_orbs) ) & 
      call close_esh5_read(h5id_input_orbs)
    if( LEN(TRIM(h5_add_orbs)) > 0 ) & 
      call close_esh5_read(h5id_add_orbs)

    call esh5_posthf_write_norb(h5id_output_orbs%id, orbsG, 5, nksym, &
        h5id_output_orbs%norbK, error)
    if(error .ne. 0 ) &
      call errore('generate_orbitals','error writing norb OrbsG',1)
    if(write_psir) then
      call esh5_posthf_write_norb(h5id_output_orbs%id,orbsR, 5, nksym, &
        h5id_output_orbs%norbK, error)
      if(error .ne. 0 ) &
        call errore('generate_orbitals','error writing norb OrbsR',1)
    endif
    call close_esh5_write(h5id_output_orbs)

    if(allocated(eigval)) deallocate(eigval)
    if(allocated(weights)) deallocate(weights)
    IF( ALLOCATED(spsi) ) DEALLOCATE (spsi)
    if(okvan .or. okpaw) CALL deallocate_bec_type (becp)
  !
  END SUBROUTINE generate_orbitals

  SUBROUTINE generate_full_pw_basis(out_prefix,dfft)
    !
    USE cell_base, ONLY: alat, tpiba, at, bg
    USE fft_types, ONLY: fft_type_descriptor
    ! 
    IMPLICIT NONE
    !
    CHARACTER(len=256), INTENT(IN) :: out_prefix
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    !
    CHARACTER(len=256) :: tmp
    INTEGER :: i, ik, npw, n, error 
    COMPLEX(DP), allocatable :: Orb(:)
    !
    TYPE(h5file_type) :: h5id_output_orbs
    !

    ! only roots of each pool generate orbitals 
    if(me_image .ne. root_image) return

    allocate(Orb(npwx))

    tmp = TRIM( tmp_dir ) // TRIM( out_prefix ) // '.orbitals.h5'
    call open_esh5_write(h5id_output_orbs,dfft,tmp,.false.)
    h5id_output_orbs%norbK(:) = 0
    mixed_basis = .true.   

    if(nspin>1) call errore('generate_full_pw_basis',' nspin>1 ',1)

    do ik=1,nksym
      !
      npw = ngksym(ik) 
      h5id_output_orbs%norbK(ik) = npw
      !
      do n=1,npw
        !   
        Orb(:) = (0.d0,0.d0)
        Orb(n) = (1.0,0.d0) 
        ! not normalized properly with PAW/USPP !!!
        call esh5_posthf_write_band(h5id_output_orbs%id, orbsG, 5, ik-1, &
                    npwx, n-1, Orb, error)
        if(error.ne. 0 ) &
          call errore('generate_full_pw_basis','Error writing orbital.',1)  
        !
      enddo
      !  
    enddo
    !  
    write(*,*) ' ik, # orbitals: '
    do ik=1,nksym
      write(*,*) ik,h5id_output_orbs%norbK(ik)
    enddo
    !
    call esh5_posthf_write_norb(h5id_output_orbs%id, orbsG, 5, nksym, &
        h5id_output_orbs%norbK, error)
    if(error .ne. 0 ) &
      call errore('generate_full_pw_basis','error writing norb OrbsG',1)
    call close_esh5_write(h5id_output_orbs)

    deallocate(Orb)
  !
  END SUBROUTINE generate_full_pw_basis

  ! assuming the same structure as evc for now
  SUBROUTINE get_spanning_basis(ik, n, ntot, A, lda, B, ldb, C, ldc, S, lds, spsi, eig, cutoff)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik, n, lda, ldb, ldc, lds
    INTEGER, INTENT(INOUT) :: ntot
    REAL(DP), INTENT(IN) :: cutoff
    REAL(DP), INTENT(INOUT) :: eig(:) 
    COMPLEX(DP), INTENT(INOUT) :: S(lds,*)
    COMPLEX(DP), INTENT(INOUT) :: spsi(ldb,*)
    COMPLEX(DP), INTENT(INOUT) :: C(ldc,*), A(lda,*), B(ldb,*) 
    !
    INTEGER :: j, ie
    COMPLEX(DP) :: ctemp
    REAL(DP) :: rtemp

    if( lds < 2*n ) &
      call errore('get_spanning_basis','Error: lds < 2*n in get_spanning_basis.',1)
    eig(:) = 0.d0
    S(1:lds,1:lds) = (0.d0,0.d0)
    ! A/A 
    call Overlap(n,n,npwx,(1.d0,0.d0),A,lda,A,lda,(0.d0,0.d0),S,lds,.true.,spsi)
    if(noncolin) then
      ! up/down: 
      call Overlap(n,n,npwx,(1.d0,0.d0),A,lda,A(npwx+1,1),lda,   &
                        (0.d0,0.d0),S(:,n+1:2*n),lds,.true.,spsi)
      ! down/down:
      call Overlap(n,n,npwx,(1.d0,0.d0),A(npwx+1,1),lda,      &
                        A(npwx+1,1),lda,(0.d0,0.d0),S(n+1,n+1),lds,.true.,spsi)
    else
      ! up/down 
      call Overlap(n,n,npwx,(1.d0,0.d0),A,lda,B,ldb,      &
                        (0.d0,0.d0),S(:,n+1:2*n),lds,.true.,spsi)
      ! down/down
      call Overlap(n,n,npwx,(1.d0,0.d0),B,ldb,B,ldb,     &
                        (0.d0,0.d0),S(n+1,n+1),lds,.true.,spsi)
    endif
    call eigsys('V', 'U', .true., 2*n, lds, S(:,1:2*n), eig)

    do ie=1,2*n
      if(verbose) write(*,*) ik, n, ie, eig(ie)
      ctemp = 1.d0/sqrt(dble(eig(ie)))
      S(1:2*n,ie) = S(1:2*n,ie)*ctemp
      if( dble(eig(ie)) < cutoff ) exit
      ntot = ie
    enddo

    call zgemm('N','N',npwx,ntot,n,(1.d0,0.d0),A,lda,S(1,1),lds,&
                              (0.d0,0.d0),C(1,1),ldc)
    if(noncolin) then
      call zgemm('N','N',npwx,ntot,n,(1.d0,0.d0),A(npwx+1,1),lda, &
                        S(n+1,1),lds,(1.d0,0.d0),C(1,1),ldc)
    else
      call zgemm('N','N',npwx,ntot,n,(1.d0,0.d0),B,ldb,S(n+1,1),&
                              lds,(1.d0,0.d0),C(1,1),ldc)
    endif

  END SUBROUTINE get_spanning_basis 

  SUBROUTINE add_extra_orbitals_from_file(h5id_file, ik, norbK, &
                norb_extra, Orbitals,  &
                orb_extra, S, spsi, dfft, g2g, igk, cutoff )
    !
    USE fft_types, ONLY: fft_type_descriptor
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id_file
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: ik, norb_extra
    INTEGER, INTENT(INOUT) :: norbK
    INTEGER, INTENT(IN) :: g2g(:), igk(:)
    REAL(DP), INTENT(IN) :: cutoff
    COMPLEX(DP), INTENT(INOUT) :: S(:) 
    COMPLEX(DP), INTENT(INOUT) :: Orbitals(:,:)
    COMPLEX(DP), INTENT(INOUT) :: orb_extra(:)
    COMPLEX(DP) :: spsi(:,:) 
    !
    INTEGER :: ibnd, nnr_, i, j, error, i0, sz, nread 
    COMPLEX(DP) :: ctemp(1,1)
    REAL(DP) :: rtemp
    COMPLEX(DP) :: gammafac

    gammafac=(1.d0,0.d0)
    if(gamma_only) gammafac=(2.d0,0.d0)

    nread = norb_extra
    !if(read_from_h5 > 0) nread = min(norb_extra,read_from_h5)

    nnr_ = size(orb_extra, 1)
    sz = size(Orbitals, 1)
    i0 = 1 

    if( norbK == 0 ) then
      i0 = 2
      norbK = 1
      call esh5_posthf_read(h5id_file%id,ik-1,0,1,orb_extra,nnr_,error)
      if(error .ne. 0 ) &
        call errore('add_extra_orbitals_from_file','error reading additional orbital',3)
      !
      if(h5id_file%grid_type == 0) then
        psic(:) = (0.d0,0.d0) 
        psic(g2g(1:nnr_)) = orb_extra(1:nnr_)
        orb_extra(1:npwx) = psic(dfft%nl(igk(1:npwx)))
        orb_extra(npwx+1:nnr_) = (0.d0,0.d0)
      elseif(h5id_file%grid_type == 1 ) then
        if(nnr_ .lt. npwx) then
          orb_extra(nnr_+1:npwx) = (0.d0,0.d0)
        else
          orb_extra(npwx+1:nnr_) = (0.d0,0.d0)
        endif
      else
        write(*,*) ' Unknown grid_type (2):',h5id_file%grid_type
        call errore('add_extra_orbitals_from_file',' Unknown grid_type:',1)
      endif
      !  
!      ctemp = (0.d0,0.d0)
!      do j=gstart,npwx
!        ctemp = ctemp + gammafac * orb_extra(j) * conjg(orb_extra(j))
!      enddo
!      if( gstart == 2 ) ctemp = ctemp + orb_extra(1) * conjg(orb_extra(1))
!      if(gamma_only) ctemp = cmplx( real(ctemp), 0.0_dp, kind=dp)  
      call Overlap(1,1,npwx,(1.d0,0.d0),orb_extra,npwx,orb_extra,npwx,&
                            (0.d0,0.d0),ctemp,1,.true.,spsi(:,1))   
      !  
      rtemp = 1.d0 / sqrt(dble(ctemp(1,1)))
      Orbitals(1:npwx,1) = orb_extra(1:npwx) * rtemp
      !  
    endif
    
    do ibnd=i0,nread
      !
      call esh5_posthf_read(h5id_file%id,ik-1,ibnd-1,1,orb_extra,nnr_,error)
      if(error .ne. 0 ) &
        call errore('add_extra_orbitals_from_file','error reading additional orbital',4)
      !
      if(h5id_file%grid_type == 0) then  
        psic(:) = (0.d0,0.d0)
        psic(g2g(1:nnr_)) = orb_extra(1:nnr_)
        orb_extra(1:npwx) = psic(dfft%nl(igk(1:npwx)))
        orb_extra(npwx+1:nnr_) = (0.d0,0.d0)
      elseif(h5id_file%grid_type == 1 ) then
        if(nnr_ .lt. npwx) then
          orb_extra(nnr_+1:npwx) = (0.d0,0.d0)
        else
          orb_extra(npwx+1:nnr_) = (0.d0,0.d0)
        endif
      else
        write(*,*) ' Unknown grid_type (2):',h5id_file%grid_type
        call errore('add_extra_orbitals_from_file',' Unknown grid_type:',1)
      endif  
      ! 
!      ctemp = (0.d0,0.d0)
!      do j=gstart,npwx
!        ctemp = ctemp + gammafac * orb_extra(j) * conjg(orb_extra(j))
!      enddo
!      if( gstart == 2 ) ctemp = ctemp + orb_extra(1) * conjg(orb_extra(1))
!      if(gamma_only) ctemp = cmplx( real(ctemp), 0.0_dp, kind=dp)  
      call Overlap(1,1,npwx,(1.d0,0.d0),orb_extra,npwx,orb_extra,npwx,&
                            (0.d0,0.d0),ctemp,1,.true.,spsi(:,1))   
      !  
      rtemp = 1.d0 / sqrt(dble(ctemp(1,1)))
      orb_extra(1:npwx) = orb_extra(1:npwx) * rtemp
      !  
      call zgemv('C',npwx,norbK,gammafac*(1.d0,0.d0),Orbitals(1,1),npwx, &
                   orb_extra(1),1,(0.d0,0.d0),S,1)
      if(gamma_only .and. gstart == 2) then 
        do j=1,norbK
          S(j) = cmplx( real(S(j) - conjg(Orbitals(1,j)) * orb_extra(1)), &
                            0.0_dp, kind=dp)
        enddo
      endif  
      !
      call zgemv('N',npwx,norbK,(-1.d0,0.d0),Orbitals(1,1),npwx, &
                   S,1,(1.d0,0.d0),orb_extra(1),1)
      !  
!      ctemp = (0.d0,0.d0)
!      do j=gstart,npwx
!        ctemp = ctemp + gammafac * orb_extra(j) * conjg(orb_extra(j)) 
!      enddo
!      if( gstart == 2 ) ctemp = ctemp + orb_extra(1) * conjg(orb_extra(1))
!      if(gamma_only) ctemp = cmplx( real(ctemp), 0.0_dp, kind=dp)  
      call Overlap(1,1,npwx,(1.d0,0.d0),orb_extra,npwx,orb_extra,npwx,&
                            (0.d0,0.d0),ctemp,1,.true.,spsi(:,1))   
      !  
      if( abs(ctemp(1,1)) > cutoff ) then
        norbK = norbK + 1
        rtemp = 1.d0 / sqrt(dble(ctemp(1,1)))
        Orbitals(1:npwx,norbK) = orb_extra(1:npwx) * rtemp 
      endif
      !  
    enddo
    ! 
  END SUBROUTINE add_extra_orbitals_from_file 

  SUBROUTINE add_extra_orbitals(norbK, nvir, Orbitals, S, spsi, cutoff )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nvir
    INTEGER, INTENT(INOUT) :: norbK
    REAL(DP), INTENT(IN) :: cutoff
    COMPLEX(DP), INTENT(INOUT) :: Orbitals(:,:)
    COMPLEX(DP), INTENT(INOUT) :: S(:) 
    COMPLEX(DP) :: spsi(:,:) 
    !
    INTEGER :: ibnd, i, j, error, n0
    COMPLEX(DP) :: ctemp(1,1)
    REAL(DP) :: rtemp
    COMPLEX(DP) :: gammafac

    gammafac=(1.d0,0.d0)
    if(gamma_only) gammafac=(2.d0,0.d0)

    if( norbK == 0 ) &
      call errore('add_extra_orbitals','norbK == 0 in add_extra_orbitals',1) 
    
    n0 = norbK
    do ibnd=n0+1,n0+nvir
      !
!      ctemp = (0.d0,0.d0)
!      do j=gstart,npwx
!        ctemp = ctemp + gammafac * Orbitals(j,ibnd) * conjg(Orbitals(j,ibnd))
!      enddo
!      if( gstart == 2 ) ctemp = ctemp + Orbitals(1,ibnd) * conjg(Orbitals(1,ibnd))
!      if(gamma_only) ctemp = cmplx( real(ctemp), 0.0_dp, kind=dp)  
      call Overlap(1,1,npwx,(1.d0,0.d0),Orbitals(:,ibnd),npwx,Orbitals(:,ibnd),npwx,&
                            (0.d0,0.d0),ctemp,1,.true.,spsi(:,1))
      !  
      rtemp = 1.d0 / sqrt(dble(ctemp(1,1)))
      Orbitals(1:npwx,ibnd) = Orbitals(1:npwx,ibnd) * rtemp
      !  
      call zgemv('C',npwx,norbK,gammafac*(1.d0,0.d0),Orbitals(1,1),npwx, &
                   Orbitals(1,ibnd),1,(0.d0,0.d0),S,1)
      if(gamma_only .and. gstart == 2) then
        do j=1,norbK
          S(j) = cmplx( real(S(j) - conjg(Orbitals(1,j)) * Orbitals(1,ibnd)), &
                            0.0_dp, kind=dp)
        enddo
      endif
      !
      call zgemv('N',npwx,norbK,(-1.d0,0.d0),Orbitals(1,1),npwx, &
                   S,1,(1.d0,0.d0),Orbitals(1,ibnd),1)
      !  
!      ctemp = (0.d0,0.d0)
!      do j=gstart,npwx
!        ctemp = ctemp + gammafac * Orbitals(j,ibnd) * conjg(Orbitals(j,ibnd)) 
!      enddo
!      if( gstart == 2 ) ctemp = ctemp + Orbitals(1,ibnd) * conjg(Orbitals(1,ibnd))
!      if(gamma_only) ctemp = cmplx( real(ctemp), 0.0_dp, kind=dp)  
      call Overlap(1,1,npwx,(1.d0,0.d0),Orbitals(:,ibnd),npwx,Orbitals(:,ibnd),npwx,&
                            (0.d0,0.d0),ctemp,1,.true.,spsi(:,1))
      !  
      if( abs(ctemp(1,1)) > cutoff ) then
        norbK = norbK + 1
        rtemp = 1.d0 / sqrt(dble(ctemp(1,1)))
        Orbitals(1:npwx,norbK) = Orbitals(1:npwx,ibnd) * rtemp 
      endif
      !  
    enddo
    ! 
  END SUBROUTINE add_extra_orbitals

  SUBROUTINE write_trial_wavefunction(h5id_orbs, h5id_hamil, dfft, &
                nelmax,M,ndet,lowcut,highcut, h5wfn) 
    !
    USE klist, ONLY: wk
    USE becmod,   ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE uspp,     ONLY : vkb, nkb
    USE paw_variables, ONLY : okpaw
    USE uspp,       ONLY : okvan
    USE fft_types, ONLY: fft_type_descriptor
    ! 
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id_hamil, h5id_orbs 
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=256), INTENT(IN), OPTIONAL :: h5wfn
    INTEGER, INTENT(IN) :: nelmax,ndet
    REAL(DP), INTENT(IN) :: lowcut,highcut
    COMPLEX(DP), INTENT(OUT) :: M(:,:,:,:)  ! Overlap between basis states and
                                            ! occupied KS states. 
    !
    TYPE(h5file_type) :: h5id_wfn
    CHARACTER(len=10) :: ftype
    INTEGER :: i,j, ik, ispin, ikk, ia, n, mix_, err_,  nel, no
    INTEGER :: non_coll, error, npw
    INTEGER :: maxl(2), Inel(2), maxnorb
    INTEGER, ALLOCATABLE :: nkocc(:,:)
    REAL(DP) :: neltot(2)
    REAL(DP) :: scl, pnorm, ctemp
    REAL(DP), ALLOCATABLE :: wg_(:,:,:)  ! orbital occupations 
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:) 
    COMPLEX(DP), ALLOCATABLE :: Orbs(:,:)  ! old orbitals
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:) 
    COMPLEX(DP), ALLOCATABLE :: Ov(:,:)  ! temporary matrix 
    COMPLEX(DP), ALLOCATABLE :: T1(:,:)  ! temporary matrix 

    nmax_DM=0
    if(present(h5wfn)) then
      ftype='esh5'
      call open_esh5_read(h5id_wfn,TRIM(h5wfn))
      if(noncolin .and. h5id_wfn%npol.ne.2) &
        call errore('write_trial_wavefunction','Error: noncolin and wfn%npol != 2.',1)
    else
      ftype='davcio'
    endif
    maxnorb = maxval(h5id_orbs%norbK(:))

    ! only root works 
    if(me_image .ne. root_image) return

    allocate( nkocc(nksym,numspin) )

    ! generate modified occupation tensor for calculation of trial wfn
!    nelmax = 0
    Inel(:) = 0
    neltot(:) = 0.d0
!    do ispin=1,numspin
!      do ik=1,nksym
!        ikk = ik + nksym*(ispin-1)
!        if(abs(wk(ikk))>1.d-10) then
!            scl = 1.d0/wk(ikk)
!        else 
!            scl = 1.d0
!        endif
!        do ia=1,nbnd
!          ! MAM: for high-T or for highly degenerate states this needs to be reduced
!          ! maybe make it a parameter that defaults to 0.01
!          if( abs(wg(ia,ikk)*scl) > 0.01d0 .and. ia > nelmax )  &
!            nelmax = ia
!        enddo
!      enddo
!    enddo
    allocate( wg_(nelmax,nksym,numspin) )
    wg_(:,:,:) = 0.d0  
    do ispin=1,numspin
      do ik=1,nksym
        ikk = ik + nksym*(ispin-1)
        if(abs(wk(ikk))>1.d-10) then
            scl = 1.d0/wk(ikk)
        else 
            scl = 1.d0
        endif
        do ia=1,nelmax
          neltot(ispin) = neltot(ispin) + wg(ia,ikk)*scl
          if( abs(wg(ia,ikk)*scl) > 0.01d0) wg_(ia,ik,ispin) = wg(ia,ikk)*scl
          ! just checking
          if( wg_(ia,ik,ispin) > (1.d0+1.d-6) ) then
            write(*,*) 'wg: ',wg(ia,ikk),wg(ia,ikk)*scl,wg_(ia,ik,ispin)
            call errore('write_trial_wavefunction','wg>1',1)    
          endif
          if( abs(wg_(ia,ik,ispin)) > 0.95d0 ) then
            wg_(ia,ik,ispin) = 1.d0
            Inel(ispin) = Inel(ispin) + 1 
          endif  
        enddo
      enddo
    enddo

    allocate( Orbs(npol*npwx, nelmax) )

    ! calculate overlap matrix if orbitals are mixed
    mix_=0
    nmax_DM = nelmax
    if(mixed_basis) then
      
      mix_ = 1
  
      ! this assumes that the eigenvalues are energy ordered, which they
      ! should be, otherwise nelmax has no meaning
      ! for n in [1,nelmax]
      ! M( iorb, n, ik, ispin ) = < Orbital(iorb,ik) | Psi_DFT(n, ik) >
      !                         = sum_g=1:npw conjg(Orbitals(g,iorb,ik)) * Orbs(g,n,ik,ispin) 
!      if(noncolin) then
!        allocate( M(maxnorb,nelmax,nksym,2) )
!      else
!        allocate( M(maxnorb,nelmax,nksym,nspin) )
!      endif
      if( (size(M,1) .ne. maxnorb ) .or. &
          (size(M,2) .ne. nelmax ) .or. &
          (size(M,3) .ne. nksym ) .or. &
          (size(M,4) .ne. min(2,nspin) ) ) &
        call errore('write_trial_wavefunction',' Wrong dimensions in OrbMat. ',1)
      M(:,:,:,:) = (0.d0,0.d0)

      if(okvan .or. okpaw) then
        allocate( spsi(npol*npwx, maxnorb) )
        CALL allocate_bec_type ( nkb, maxnorb, becp )
      else
        allocate(spsi(1,1))
      endif
      allocate( Orbitals(npwx, maxnorb) )  

      do ik=1,nksym
        !
        call esh5_posthf_read(h5id_orbs%id,ik-1,0,h5id_orbs%norbK(ik),Orbitals,npwx,error)
        if(error .ne. 0 ) &
          call errore('write_trial_wavefunction','error reading additional orbital',5)
        if(gamma_only .AND. gstart == 2 ) &
          Orbitals(1,:) = CMPLX( DBLE( Orbitals(1,:) ), 0.D0 ,kind=DP)
        !
        do ispin=1,numspin
          !
          ikk = ik + nksym*(ispin-1)
          !  
          npw = ngksym(ikk)
          !
          if(okvan .or. okpaw) then
            CALL gk_sort (xksym (1:3, ikk), ngm, g, ecutwfc / tpiba2, &
                  npw, igksym(1), g2kin)
            CALL init_us_2 (npw, igksym(1), xksym (1, ikk), vkb)
          endif
          !
          call get_orbitals(h5id_wfn,ftype,'psig',dfft,Orbs,&
                1,nelmax,ik,ispin)
          if(gamma_only .AND. gstart == 2 ) &
            Orbs(1,:) = CMPLX( DBLE( Orbs(1,:) ), 0.D0 ,kind=DP)
          ! 
          if(noncolin) then
            CALL Overlap(h5id_orbs%norbK(ik),nelmax,npw,(1.d0,0.d0),Orbitals(1,1),npwx,&
               Orbs(1,1),2*npwx,(0.d0,0.d0),M(1,1,ik,1),maxnorb,.true.,spsi)
            CALL Overlap(h5id_orbs%norbK(ik),nelmax,npw,(1.d0,0.d0),Orbitals(1,1),npwx,&
               Orbs(npwx+1,1),2*npwx,(0.d0,0.d0),M(1,1,ik,2),maxnorb,.true.,spsi)
          else
            CALL Overlap(h5id_orbs%norbK(ik),nelmax,npw,(1.d0,0.d0),Orbitals(1,1),npwx,&
               Orbs(1,1),npwx,(0.d0,0.d0),M(1,1,ik,ispin),maxnorb,.true.,spsi)
          endif
          !
          ! check overlap of DFT state within the Orbital basis
          !
          do ia=1,nelmax
            if( abs(wg_(ia,ik,ispin)) > 0.01d0 ) then
              pnorm = 0.d0
              do n=1,h5id_orbs%norbK(ik)
                ctemp = dble( conjg( M(n,ia,ik,ispin) ) * M(n,ia,ik,ispin) ) 
                pnorm = pnorm + ctemp 
                if( sqrt( ctemp ) > 1.d-6 .and. n .gt. nmax_DM ) & 
                  nmax_DM = n  
              enddo
              if(noncolin) then
                do n=1,h5id_orbs%norbK(ik)
                  ctemp = dble( conjg( M(n,ia,ik,2) ) * M(n,ia,ik,2) )  
                  pnorm = pnorm + ctemp  
                  if( sqrt( ctemp ) > 1.d-6 .and. n .gt. nmax_DM ) &  
                    nmax_DM = n  
                enddo
              endif  
              pnorm = sqrt(pnorm)
              if( abs(pnorm-1.d0) > 1e-4 ) & 
                write(*,*) ' WARNING: Unnormalized orbital in trial wave-function:', &
                             ispin,ik,ia,pnorm 
              pnorm = 1.d0/pnorm
              do n=1,h5id_orbs%norbK(ik)
                M(n,ia,ik,ispin) = M(n,ia,ik,ispin)*pnorm
              enddo  
              if(noncolin) then
                do n=1,h5id_orbs%norbK(ik)
                  M(n,ia,ik,2) = M(n,ia,ik,2)*pnorm
                enddo  
              endif
            else
              M(:,ia,ik,ispin) = (0.d0,0.d0)
              if(noncolin) then
                M(:,ia,ik,2) = (0.d0,0.d0)
              endif  
            endif  
          enddo
        enddo
      enddo
      !
      deallocate( Orbitals )
      !
    else
!      allocate( M(1,1,1,1) )
      M(:,:,:,:) = (0.d0,0.d0)
      do ispin=1,size(M,4)
        do ik=1,size(M,3)
          do ia=1,min(size(M,1),size(M,2))
            M(ia,ia,ik,ispin)=1.d0
          enddo
        enddo
      enddo
    endif 

    neltot(1)=nint(neltot(1))
    neltot(2)=nint(neltot(2))
    write(*,*) ' Number of electrons per spin channel: ',(neltot(i),i=1,numspin)

    ! temporary hack while I write the code to generate ndet>1 in metallic
    ! systems
    if( ndet > 1 ) then
        call errore('write_trial_wavefunction','ndet > 1 not yet implemented.',1)
    else
      ! "fix" occupations for ndet=1 
      do ispin=1,numspin
        do while(Inel(ispin)+0.1d0 < neltot(ispin))
          maxl = maxloc( wg_(:,:,ispin), mask=wg_(:,:,ispin).lt.0.95d0 )
          if(verbose) write(*,*) 'Setting weight to one:', &
                ispin,maxl(1),maxl(2),wg_(maxl(1),maxl(2),ispin) 
          wg_(maxl(1),maxl(2),ispin) = 1.d0
          Inel(ispin) = Inel(ispin) + 1
        enddo 
        do ik=1,nksym 
          do n=1,min(h5id_orbs%norbK(ik),nelmax)
            if(wg_(n,ik,ispin) < 0.95d0) wg_(n,ik,ispin)=0.d0
          enddo 
        enddo 
        write(*,*)'ispin:',neltot(ispin),sum(wg_(:,:,ispin))
        if( abs(neltot(ispin) - sum(wg_(:,:,ispin)) ) > 1e-4) &
          call errore('write_trial_wavefunction','Error with electron count.',1)
      enddo 
    endif

    err_ = 0
    call esh5_posthf_write_wavefunction(h5id_hamil%id,nspin,nksym,h5id_orbs%norbK,maxnorb,&
                nelmax,ndet,wg_,mix_,lowcut,highcut,M,nkocc,err_)
    if(err_ .ne. 0) call errore('write_trial_wavefunction',' Error generating wavefunction.',1)

!MAM...
! since the orbital matrix is stored now, calculate DM on demand, no need to store!
! move this to pp_utilities and call on demand when needed
    ! calculate density matrix for energy evaluation
    write(*,*) 'Largest band index contributing to trial wave-function: ',nelmax
    write(*,*) 'Largest orbital index contributing to trial wave-function: ',nmax_DM
    allocate( DM(nmax_DM, nmax_DM, nksym, nspin), Ov(nelmax, nelmax), &
              T1(nelmax,nmax_DM), DM_mf(nmax_DM, nmax_DM, nksym, nspin) )
    DM(:,:,:,:)=(0.d0,0.d0)
    DM_mf(:,:,:,:)=(0.d0,0.d0)
    ! DM: the spin index encodes the spin sector, 
    ! 1: (up,up)
    ! 2: (dn,dn)  ! lsda and noncolin
    ! 3: (up,dn)  ! only noncolin
    ! 4: (dn,up)  ! only noncolin
    if(noncolin) then  
      ! noncolin is slightly different!  
      do ik=1,nksym

        nel = nkocc(ik,1)
        no = min(nmax_DM,h5id_orbs%norbK(ik))

        ! DM = T( M * inv( H(M) * M ) * H(M) )
        CALL ZGEMM('C','N',nel,nel,h5id_orbs%norbK(ik),(1.d0,0.d0),M(1,1,ik,1),maxnorb,&
                   M(1,1,ik,1),maxnorb,(0.d0,0.d0),Ov,nelmax)
        CALL ZGEMM('C','N',nel,nel,h5id_orbs%norbK(ik),(1.d0,0.d0),M(1,1,ik,2),maxnorb,&
                   M(1,1,ik,2),maxnorb,(1.d0,0.d0),Ov,nelmax)
        CALL zmatinv(nel, Ov, nelmax)
        ! inv( H(M) * M ) * H(M)  up sector
        CALL ZGEMM('N','C',nel,no,nel,(1.d0,0.d0),Ov,nelmax,&
                   M(1,1,ik,1),maxnorb,(0.d0,0.d0),T1,nelmax)
        ! (up,up)
        CALL ZGEMM('T','T',no,no,nel,(1.d0,0.d0),T1,nelmax,&
                   M(1,1,ik,1),maxnorb,(0.d0,0.d0),DM(1,1,ik,1),nmax_DM)
        ! (up,dn)
        CALL ZGEMM('T','T',no,no,nel,(1.d0,0.d0),T1,nelmax,&
                   M(1,1,ik,2),maxnorb,(0.d0,0.d0),DM(1,1,ik,3),nmax_DM)
        ! inv( H(M) * M ) * H(M)  down sector
        CALL ZGEMM('N','C',nel,no,nel,(1.d0,0.d0),Ov,nelmax,&
                   M(1,1,ik,2),maxnorb,(0.d0,0.d0),T1,nelmax)
        ! (dn,up)
        CALL ZGEMM('T','T',no,no,nel,(1.d0,0.d0),T1,nelmax,&
                   M(1,1,ik,1),maxnorb,(0.d0,0.d0),DM(1,1,ik,4),nmax_DM)
        ! (dn,dn)
        CALL ZGEMM('T','T',no,no,nel,(1.d0,0.d0),T1,nelmax,&
                   M(1,1,ik,2),maxnorb,(0.d0,0.d0),DM(1,1,ik,2),nmax_DM)

        !decide how to implement DM_mf 

      enddo
    elseif(mixed_basis) then
      do ispin=1,nspin
        do ik=1,nksym     

          nel = nkocc(ik,ispin)
          no = min(nmax_DM,h5id_orbs%norbK(ik))

          ! DM = T( M * inv( H(M) * M ) * H(M) )
          CALL ZGEMM('C','N',nel,nel,h5id_orbs%norbK(ik),(1.d0,0.d0),M(1,1,ik,ispin),maxnorb,&
                 M(1,1,ik,ispin),maxnorb,(0.d0,0.d0),Ov,nelmax)
          CALL zmatinv(nel, Ov, nelmax)
          CALL ZGEMM('N','C',nel,no,nel,(1.d0,0.d0),Ov,nelmax,&
                 M(1,1,ik,ispin),maxnorb,(0.d0,0.d0),T1,nelmax)
          CALL ZGEMM('T','T',no,no,nel,(1.d0,0.d0),T1,nelmax,&
                 M(1,1,ik,ispin),maxnorb,(0.d0,0.d0),DM(1,1,ik,ispin),nmax_DM)

          ! scale M by wg and form DM_mf
          ikk = ik + nksym*(ispin-1)
          scl = 1.d0
          if(abs(wk(ikk))>1.d-10) scl = 1.d0/wk(ikk)
          do ia=1,nel
            do i=1,maxnorb
              M(i,ia,ik,ispin) = M(i,ia,ik,ispin) * wg(ia,ikk)*scl 
            enddo
          enddo  
          CALL ZGEMM('T','T',no,no,nel,(1.d0,0.d0),T1,nelmax,&
                 M(1,1,ik,ispin),maxnorb,(0.d0,0.d0),DM_mf(1,1,ik,ispin),nmax_DM)

        enddo
      enddo
    else
      if(nspin.gt.1) &
        call errore('write_trial_wavefunction','nspin.gt.1 and mixed_basis=.false.',1)
      do ik=1,nksym
        scl = 1.d0
        if(abs(wk(ik))>1.d-10) scl = 1.d0/wk(ik)
        nel = nkocc(ik,1)
        do ia=1,nel
          DM(ia,ia,ik,1) = (1.d0,0.d0)
        enddo
        do ia=1,nmax_DM
          DM_mf(ia,ia,ik,1) = wg(ia,ik)*scl*(1.d0,0.d0)   ! normalization differs )
        enddo
      enddo
    endif  

    if(present(h5wfn)) call close_esh5_read(h5id_wfn)

    if(allocated(Orbs)) deallocate(Orbs)
    if(allocated(nkocc)) deallocate(nkocc)
    if( allocated(wg_) ) deallocate(wg_)
!    if( allocated(M) )   deallocate( M )
    if( allocated(Ov) )   deallocate( Ov )
    if( allocated(T1) )   deallocate( T1 )
    IF( ALLOCATED(spsi) ) DEALLOCATE (spsi)
    if(okvan .or. okpaw) CALL deallocate_bec_type (becp)
     
  END SUBROUTINE write_trial_wavefunction
  !
  SUBROUTINE diag_hf(dfft,dtype,hamil_file,occ_file,orb_file, &
            new_orb_file,rotate_basis) 
    !
    USE mp,           ONLY: mp_sum, mp_max, mp_bcast, mp_barrier
    USE mp_pools,     ONLY: npool
    USE mp_images, ONLY: intra_image_comm
    USE fft_types, ONLY: fft_type_descriptor
    USE twobody_hamiltonian, ONLY: add_FockM
#if defined(__CUDA)
    USE twobody_hamiltonian, ONLY: add_FockM_gpu
#endif
    USE wvfct, ONLY: et
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    !
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    CHARACTER(len=*), INTENT(IN) :: hamil_file,orb_file,occ_file
    CHARACTER(len=*), INTENT(IN) :: dtype,new_orb_file
    LOGICAL, INTENT(IN), OPTIONAL :: rotate_basis
    !
    LOGICAL :: rotate_basis_
    CHARACTER(len=256) :: h5name
    INTEGER :: n,h5len,oldh5,error,i0,M,ispin,n0
    INTEGER, ALLOCATABLE :: noccK(:,:) 
    COMPLEX(DP) :: ctemp
    COMPLEX(DP), ALLOCATABLE :: FockM(:,:)   ! Fock matrix
    COMPLEX(DP), ALLOCATABLE :: Orbs(:,:,:) ! old orbitals 
    COMPLEX(DP), ALLOCATABLE :: newOrbs(:,:) ! new orbitals 
    REAL(DP), ALLOCATABLE :: eigval(:)     ! eigenvalues 
    REAL(DP), ALLOCATABLE :: weights(:)     
    !
    INTEGER :: ik, ikk, nel(2)
    !
    TYPE(h5file_type) :: h5id_hamil, h5id_input_orbs, h5id_output_orbs, &
                         h5id_occ 
    !
    rotate_basis_ = .true.
    if(present(rotate_basis)) rotate_basis_=rotate_basis

    if(nspin>1 .and. (TRIM(dtype) .ne. 'full')) & 
      call errore('diag_hf','Only dtype=full allowed with nspin>1',1)

    ! open orbital file 
    h5name = TRIM( orb_file )
    call open_esh5_read(h5id_input_orbs,h5name)
    if( h5id_input_orbs%grid_type .ne. 1 ) &
      call errore('diag_hf','grid_type ne 1',1)

    h5name = TRIM( occ_file )
    call open_esh5_read(h5id_occ,h5name)
    if( h5id_occ%grid_type .ne. 1 ) &
      call errore('diag_hf','grid_type ne 1',1)

    !
    ! only insulators now, for metals you need to keep/pass 
    ! noccK(ik) and wg(ik) (or just include based on cutoff from wg(ik))
    !
    allocate( noccK(nksym,numspin), weights(h5id_input_orbs%maxnorb) )
    ! this needs to be improved!
    call get_noccK(noccK,nel,nbnd,nksym,numspin,wg,nbnd,wk)

    if(root_image == me_image) then
      ! open hamiltonian file
      oldh5=1
      h5name = TRIM( hamil_file )
      h5len = LEN_TRIM(h5name)
      CALL esh5_posthf_open_file(h5id_hamil%id,h5name,h5len,oldh5)
      if(oldh5 .ne. 0 ) &
        call errore('diag_hf','error opening hamil file',1)
      ! open new orbital file
      h5name = TRIM( new_orb_file )
      call open_esh5_write(h5id_output_orbs,dfft,h5name)
    endif

    allocate( FockM(h5id_input_orbs%norbK(1),h5id_input_orbs%norbK(1)) ) 
    if(root_image == me_image) &
      allocate( newOrbs(npwx,h5id_input_orbs%maxnorb), &
                eigval(h5id_input_orbs%maxnorb), &
                Orbs(npwx, h5id_input_orbs%maxnorb, 1) )

    do ik=1,nksym

      do ispin = 1,nspin

        ikk = ik + nksym*(ispin-1)
        if( size(FockM,1) .ne. h5id_input_orbs%norbK(ik) ) then 
          deallocate(FockM)
          allocate( FockM(h5id_input_orbs%norbK(ik),h5id_input_orbs%norbK(ik)) )
        endif
        FockM(:,:) = (0.d0,0.d0)

        ! 1. read F(ik) = H1(ik)
        if(root_image == me_image)  then  
          call esh5_posthf_read_h1(h5id_hamil%id,ik,FockM(1,1),error)
          if(error .ne. 0 ) &
            call errore('diag_hf','error reading H1',error)
          ! matrix is written in C ordering, so transpose
          call transpose_mat(h5id_input_orbs%norbK(ik),FockM(:,:))
        endif

        ! 2. add 2-body contribution
#if defined(__CUDA)
        call add_FockM_gpu(h5id_input_orbs,h5id_occ,ik,ispin,dfft,noccK,FockM)          
#else
        call add_FockM(h5id_input_orbs,h5id_occ,ik,ispin,dfft,noccK,FockM)          
#endif

        if(root_image == me_image) then

          if( TRIM(dtype) == 'full' ) then
            M = h5id_input_orbs%norbK(ik)
            i0 = 1    
          elseif( TRIM(dtype) == 'vir_keep_occ_et' ) then
            M = h5id_input_orbs%norbK(ik)-noccK(ik,ispin)
            i0 = noccK(ik,ispin)+1    
            do n=1,i0-1
              eigval(n) = et(n,ikk)*e2Ha !   
            enddo
          elseif( TRIM(dtype) == 'vir_update_occ_et' ) then
            M = h5id_input_orbs%norbK(ik)-noccK(ik,ispin)
            i0 = noccK(ik,ispin)+1    
            do n=1,i0-1
              eigval(n) = FockM(n,n) 
            enddo
          else
            call errore('diag_hf','Unknown dtype.',1)
          endif

          ! 3. diagonalize F (possibly over blocks)
          call eigsys('V','U',.false.,M,h5id_input_orbs%norbK(ik), &
                FockM(i0,i0),eigval(i0))

          if(verbose) then
            write(*,*) 'Eigenvalues for k-point, spin:',ik,ispin
            do n=1,h5id_input_orbs%norbK(ik)
              write(*,*)n,eigval(n)
            enddo
            write(*,*)''
            FLUSH( stdout )
          endif

          if(rotate_basis_) then

            call get_orbitals_set(h5id_input_orbs,'esh5','psig',dfft,ispin,Orbs,&
                  1,h5id_input_orbs%norbK(ik),ik,1)

            ! 4. rotate orbitals and output with eigenvalues
            newOrbs(:,1:i0-1) =  Orbs(:,1:i0-1,1) 
            call zgemm('N','N',npwx,M,M,(1.d0,0.d0), &
                Orbs(1,i0,1),npwx,FockM(i0,i0),h5id_input_orbs%norbK(ik),&
                (0.d0,0.d0),newOrbs(1,i0),npwx)

            call esh5_posthf_write(h5id_output_orbs%id, orbsG, 5, ikk-1, npwx,&
                       h5id_input_orbs%norbK(ik), newOrbs, npwx, error)
            if(error .ne. 0 ) &
              call errore('diag_hf','error writing rotated orbitals',1)

          else

            if( TRIM(dtype) .ne. 'full' ) then
              ! make FockM block diagonal, with identity on occ/occ block 
              FockM(1:i0,:)=(0.d0,0.d0)  
              FockM(:,1:i0)=(0.d0,0.d0)  
              do n=1,i0
                FockM(n,n)=(1.d0,0.d0)  
              enddo  
            endif

            call esh5_posthf_write(h5id_output_orbs%id, 'CanOrbMat', 9, ikk-1, &
                       h5id_input_orbs%norbK(ik), h5id_input_orbs%norbK(ik), &
                       FockM, h5id_input_orbs%norbK(ik), error)
            if(error .ne. 0 ) &
              call errore('diag_hf','error writing CanOrbMat rotation matrix',1)

          endif

          ! fix for metals/finite-T
          weights(:)=0.d0
          do n=1,noccK(ik,ispin)
            weights(n)=1.d0
          enddo

          call esh5_posthf_write_et(h5id_output_orbs%id, orbsG, 5, ikk-1, &
                  h5id_input_orbs%norbK(ik), eigval, weights, error) 
          if(error .ne. 0 ) &
            call errore('diag_hf','error writing eigenvalues',1)
          !
        endif  ! ionode
        !
      enddo ! ispin
      !
    enddo ! ik
    call mp_barrier(intra_image_comm)

    if(allocated(noccK)) deallocate(noccK)
    if(allocated(FockM)) deallocate(FockM)
    if(allocated(eigval)) deallocate(eigval)
    if(allocated(weights)) deallocate(weights)
    if(allocated(Orbs)) deallocate(Orbs)
    if(allocated(newOrbs)) deallocate(newOrbs)

    if(root_image==me_image) then
      h5id_output_orbs%norbK(:) = h5id_input_orbs%norbK(:)
      call esh5_posthf_close_file(h5id_hamil%id)
! can't remember if this is supposed to be size norb or nspin*norb, check!!!
      call esh5_posthf_write_norb(h5id_output_orbs%id,orbsG, 5, nksym, &
                                  h5id_output_orbs%norbK, error)
      if(error .ne. 0 ) &
        call errore('diag_hf','error writing norb OrbsG',1)
      call close_esh5_write(h5id_output_orbs)
    endif
    call close_esh5_read(h5id_input_orbs)
    call close_esh5_read(h5id_occ)

  END SUBROUTINE diag_hf 
  !
  ! collects wavefunctions into esh5 file
  !
  subroutine davcio_to_esh5(h5name,M,dfft)
    !
    USE parallel_include
    USE fft_types, ONLY: fft_type_descriptor
    USE lsda_mod,             ONLY : nspin, isk, lsda
    USE klist,                ONLY : nkstot, wk, nks
    USE wvfct,                ONLY : npwx, et, wg, nbnd
    USE io_files,             ONLY : nwordwfc, iunwfc
    USE mp_pools,             ONLY : intra_pool_comm, me_pool, root_pool, &
                                     inter_pool_comm, npool, my_pool_id
    USE io_global, ONLY: ionode
    USE mp_images, ONLY: intra_image_comm
    USE mp,           ONLY: mp_barrier
    USE wavefunctions, ONLY : evc
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: M
    CHARACTER(len=256), INTENT(IN) :: h5name 
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    !
    INTEGER :: i, j, iks, ik, ikk, n, error, info 
    !
    INTEGER, ALLOCATABLE :: orbs(:)
    INTEGER, EXTERNAL    :: global_kpoint_index
    TYPE(h5file_type) :: h5id_orbs
    REAL(DP), ALLOCATABLE :: occ(:,:,:)
    INTEGER :: idata(2)
    INTEGER :: status(MPI_STATUS_SIZE)
    !

    if(ionode) call print_freemem()
!    if(lsda) write(*,*) ' WARNING: davcio_to_esh5 with lsda, check check check! '
    ! might not be correct in lsda
    !
    ! root_image has full et/wg, no need to communicate that
    ! since pw parallelization is not allowed, npwx is global
    ! 
    idata(1) = nks
    idata(2) = global_kpoint_index (nkstot, 1)
    if( M > nbnd ) &
      call errore('davcio_to_esh5','M > nbnd',1)
    if(me_image == root_image) then
      allocate( occ(nbnd,nks,2) )
      occ(:,:,1) = wg(:,1:nks)  
      occ(:,:,2) = et(:,1:nks)  
      do ik=1,nks
        if(abs(wk(ik))>1.d-10) occ(:,ik,1) = occ(:,ik,1)/wk(ik)
      enddo
      call open_esh5_write(h5id_orbs,dfft,h5name,.false.,npol)
      ! write your own k-points
      iks = global_kpoint_index (nkstot, 1)
      do ik=1,nks
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        ikk = ik + iks - 1
        call esh5_posthf_write(h5id_orbs%id,orbsG,5,ikk-1,npol*npwx, &
            M,evc,npol*npwx,error)
        if(error .ne. 0 ) &
          call errore('davcio_to_esh5','error writing orbital',1)
        call esh5_posthf_write_et(h5id_orbs%id,orbsG,5,ikk-1, &
            M,occ(1,ik,2),occ(1,ik,1),error)
        if(error .ne. 0 ) &
          call errore('davcio_to_esh5','error writing eigenvalues',1)
      enddo
      deallocate(occ)
      !
      ! now receive from other procs
      !
      do n=2,npool
        call mp_barrier(inter_pool_comm)
        CALL MPI_RECV( idata(1), 2, MPI_INT, n-1, 101+n-1, &
                  inter_pool_comm, status, info )
        CALL errore( 'davcio_to_esh5', 'mpi_send', info )
        allocate( occ(nbnd,idata(1),2) )
        CALL MPI_RECV( occ(1,1,1), nbnd*idata(1)*2, MPI_DOUBLE_PRECISION, n-1, 201+n-1, &
                  inter_pool_comm, status, info )
        CALL errore( 'davcio_to_esh5', 'mpi_send', info )
        do ik=1,idata(1)
          CALL MPI_RECV( evc(1,1), 2*npol*npwx*M, MPI_DOUBLE_PRECISION , n-1, 301+n-1, &
                  inter_pool_comm, status, info )
          CALL errore( 'davcio_to_esh5', 'mpi_send', info )
          ikk = ik + idata(2) - 1
          call esh5_posthf_write(h5id_orbs%id,orbsG,5,ikk-1,npol*npwx, &
              M,evc,npol*npwx,error)
          if(error .ne. 0 ) &
            call errore('davcio_to_esh5','error writing orbital',1)
          call esh5_posthf_write_et(h5id_orbs%id,orbsG,5,ikk-1, &
              M,occ(1,ik,2),occ(1,ik,1),error)
          if(error .ne. 0 ) &
            call errore('davcio_to_esh5','error writing eigenvalues',1)
        enddo
        deallocate(occ)
      enddo      
      !  
      allocate(orbs(nkstot))
      orbs(:) = M
      call esh5_posthf_write_norb(h5id_orbs%id,orbsG,5,nkstot,orbs,error)
      if(error .ne. 0 ) &
        call errore('davcio_to_esh5','error writing norb OrbsG',1)        
      deallocate(orbs)
      call close_esh5_write(h5id_orbs)

    else

      if(me_pool == root_pool) then
        allocate( occ(nbnd,nks,2) )
        occ(:,:,1) = wg(:,1:nks)
        occ(:,:,2) = et(:,1:nks)
        do ik=1,nks
          if(abs(wk(ik))>1.d-10) occ(:,ik,1) = occ(:,ik,1)/wk(ik)
        enddo
        do n=2,npool
          call mp_barrier(inter_pool_comm)
          if(my_pool_id+1 .ne. n) cycle 
          CALL MPI_SEND( idata(1), 2, MPI_INT, 0, 101+n-1, &
                    inter_pool_comm, info )
          CALL errore( 'davcio_to_esh5', 'mpi_send', info )
          CALL MPI_SEND( occ(1,1,1), nbnd*nks*2, MPI_DOUBLE_PRECISION, 0, 201+n-1, &
                    inter_pool_comm, info )
          CALL errore( 'davcio_to_esh5', 'mpi_send', info )
          do ik=1,nks
            CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
            CALL MPI_SEND( evc(1,1), 2*npol*npwx*M, MPI_DOUBLE_PRECISION , 0, 301+n-1, &
                    inter_pool_comm, info )
            CALL errore( 'davcio_to_esh5', 'mpi_send', info )
          enddo
        enddo
        deallocate(occ)
      endif  

    endif
    !
    call mp_barrier(intra_image_comm)
    !
  end subroutine davcio_to_esh5
  !
  SUBROUTINE get_noccK(noccK, nel, norb, nk, nsp, wg, nbnd, wk)
    !
    USE KINDS, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: norb,nk,nsp,nbnd
    REAL(DP), INTENT(IN) :: wg(nbnd,*)
    REAL(DP), INTENT(IN), OPTIONAL :: wk(*)
    INTEGER, INTENT(OUT) :: noccK(nk,nsp), nel(nsp)
    !
    INTEGER :: ispin, ik, ikk, ia
    REAL(DP) :: scl, etemp
    !
    nel(:) = 0
    noccK(:,:)=0
    do ispin=1,nsp
      etemp=0.d0
      do ik=1,nk
        ikk = ik + nk*(ispin-1)
        scl = 1.d0
        if(present(wk)) then
          if(abs(wk(ikk))>1.d-10) &
            scl = 1.d0/wk(ikk)
        endif
        do ia=1,norb
          etemp = etemp + wg(ia,ikk)*scl
          if( abs(wg(ia,ikk)*scl) > 0.01d0 ) then
            if( abs(wg(ia,ikk)*scl-1.d0) > 1.d-4 ) &
              call errore('get_noccK',  &
                'Error: Only integer occupations on get_noccK.',1)
            noccK(ik,ispin) = noccK(ik,ispin) + 1
          endif
        enddo
      enddo
      nel(ispin) = nint(etemp)
      if( nel(ispin) .ne. sum(noccK(1:nk,ispin)) ) &
        call errore('get_noccK','Error: Inconsistent number of electrons.',1)
    enddo
    ! 
  END SUBROUTINE get_noccK
!
END MODULE orbital_generators

