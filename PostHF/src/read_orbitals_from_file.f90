!
! Written by Miguel A. Morales, LLNL, 2020
!
MODULE read_orbitals_from_file
  !
  USE KINDS, ONLY : DP
  USE control_flags, ONLY: gamma_only
  USE wavefunctions, ONLY : evc
  USE cell_base, ONLY: tpiba, tpiba2
  USE uspp,     ONLY : vkb, nkb
  USE wvfct, ONLY: npwx, nbnd, g2kin
  USE gvecw, ONLY : ecutwfc
  USE lsda_mod, ONLY: lsda, nspin
  USE gvect, ONLY: ngm, ngm_g, g, ig_l2g, gstart, gg, gcutm
  USE becmod,   ONLY : bec_type,calbec,allocate_bec_type,deallocate_bec_type
  USE fft_types, ONLY: fft_type_descriptor  
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_files, ONLY: nwordwfc, iunwfc
  USE posthf_mod, ONLY: nksym, xksym, norb
  !
  IMPLICIT NONE
  !
  ! can make these specific for read/write for safety
  TYPE h5file_type
    INTEGER*8 id
    INTEGER, ALLOCATABLE :: norbK(:)
    REAL(DP), ALLOCATABLE :: xk(:,:)
    INTEGER :: grid_type, nkpts, nmax_DM, nr1, nr2, nr3, npwx, nspin, npol, maxnorb
  END TYPE h5file_type
  !
  INTEGER, ALLOCATABLE :: igk(:)
  COMPLEX(DP) :: CONE,CZERO,CMINUSONE
  CHARACTER(6), PARAMETER  :: orbsG = 'OrbsG'
  CHARACTER(6), PARAMETER  :: orbsR = 'OrbsR'
  !  
  CONTAINS
  !
  SUBROUTINE open_esh5_read(h5id, fname, check_kpts_)
    !
    USE cell_base, ONLY: tpiba, bg
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(IN) :: fname
    TYPE(h5file_type), INTENT(INOUT) :: h5id
    LOGICAL, OPTIONAL, INTENT(IN) :: check_kpts_
    INTEGER :: h5len, oldh5, n_
    INTEGER :: i, j, ik, error
    REAL(DP), ALLOCATABLE :: xk_(:,:)
    REAL(DP) :: recv(3,3)
    INTEGER, ALLOCATABLE :: norb_(:)
    LOGICAL :: check_kpts
    !
    check_kpts = .true.
    if(present(check_kpts_)) check_kpts = check_kpts_
    ! orbital file open in read_only mode
    allocate(xk_(3,nksym))
    h5len = LEN_TRIM(fname)
    h5id%grid_type = -1
    h5id%nr1 = -1
    h5id%nr2 = -1
    h5id%nr3 = -1
    h5id%npwx=-1
    h5id%nspin=-1
    h5id%npol=-1
    h5id%id = int(-1,kind=8)
    allocate(norb_(4*nksym))
    norb_(:)=0
    call esh5_posthf_open_read(h5id%id,fname,h5len,h5id%nkpts,norb_,&
            h5id%nmax_DM,h5id%nspin,h5id%npol,h5id%npwx,xk_, &
            h5id%grid_type,h5id%nr1,h5id%nr2,h5id%nr3,recv,error)
    ! Note: the number of kpoints in the file can not be larger than nksym right now
    if(h5id%nkpts > nksym) &
      call errore('open_esh5','h5id%nkpts > nksym',1) 
    if(h5id%nspin < 1 .or. h5id%nspin > 2) &
      call errore('open_esh5','h5id%nspin < 1 OR h5id%nspin > 2',1) 
    if(.not.allocated(h5id%norbK)) then
      allocate(h5id%norbK(h5id%nkpts*h5id%nspin))
    else if(size(h5id%norbK,1) .ne. h5id%nkpts*h5id%nspin) then
      deallocate(h5id%norbK)
      allocate(h5id%norbK(h5id%nkpts*h5id%nspin))
    endif    
    h5id%norbK(1:h5id%nkpts*h5id%nspin) = norb_(1:h5id%nkpts*h5id%nspin)
    deallocate(norb_)
    h5id%maxnorb = maxval(h5id%norbK(:))
    if(error .ne. 0 ) &
      call errore('open_esh5','error opening orbital file for read',1)
    if( h5id%grid_type.eq.1 .and. h5id%npwx.gt.npwx) &
      call errore('open_esh5','Inconsistent npwx in esh5 file.',1)
    recv(:,:) = recv(:,:) / tpiba
    do i=1,3
      if(sum( (recv(1:3,i)-bg(1:3,i))**2 ) .gt. 1.d-6) then
        write(*,*) 'rec vec file: ',i,(recv(j,i),j=1,3)
        write(*,*) 'rec vec QE: ',i,(bg(j,i),j=1,3)
        call errore('open_esh5',' error: rec vectors do not agree',1)
      endif
    enddo
    if(check_kpts) then
      if(h5id%nkpts .ne. nksym) &
        call errore('open_esh5','h5id%nkpts != nksym with check_kpts.',1)
      do i=1,nksym
        if(sum( (xk_(1:3,i)/tpiba-xksym(1:3,i))**2 ) .gt. 1.d-6) then
          write(*,*) 'xk file: ',i,(xk_(j,i)/tpiba,j=1,3)
          write(*,*) 'xk  QE: ',i,(xksym(j,i),j=1,3)
          call errore('open_esh5',' error: k-points do not agree',1)
        endif
      enddo
    endif
    if(.not.allocated(h5id%xk)) then
      allocate(h5id%xk(3,h5id%nkpts))
    else if( (size(h5id%xk,1) .ne. 3) .or. &
          (size(h5id%xk,2) .ne. h5id%nkpts) ) then
      deallocate(h5id%xk)
      allocate(h5id%xk(3,h5id%nkpts))
    endif
    h5id%xk(1:3,1:h5id%nkpts) = xk_(1:3,1:h5id%nkpts)/tpiba
    !
    deallocate(xk_)
    !
  END SUBROUTINE open_esh5_read 
  !
  SUBROUTINE close_esh5_read(h5id)
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(INOUT) :: h5id
    !
    call esh5_posthf_close_read(h5id%id)
    !
    if(allocated(h5id%norbK)) deallocate(h5id%norbK)
    if(allocated(h5id%xk)) deallocate(h5id%xk)
    h5id%id = int(-1,kind=8)
    h5id%grid_type = -1
    h5id%nr1 = -1
    h5id%nr2 = -1
    h5id%nr3 = -1
    h5id%npwx=-1
    h5id%npol=-1
    h5id%nmax_DM=-1
    !
  END SUBROUTINE close_esh5_read
  !
  ! Reads wave functions from file and transforms them to real space 
  !
  !-----------------------------------------------------------------------
  SUBROUTINE get_orbitals(h5id,ftype,orbtype,dfft,Psi,b_beg,nbands,ik,ispin, &
                      Q,qkmap,becpsi)
  !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id
    CHARACTER(len=*), INTENT(IN) :: ftype, orbtype
    INTEGER,      INTENT(IN)   :: b_beg,nbands,ik,ispin
    COMPLEX(DP),  INTENT(OUT)  :: Psi(:,:)
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN), OPTIONAL :: Q
    INTEGER, INTENT(IN), OPTIONAL :: qkmap(:,:)
    TYPE(bec_type), OPTIONAL :: becpsi
    !
    CONE  = (1.d0,0.d0)
    CMINUSONE  = (-1.d0,0.d0)
    CZERO = (0.d0,0.d0)
    !
    if(.not.allocated(igk)) allocate(igk(npwx))
    !
    if( (trim(orbtype).ne.'psir') .and. &
        (trim(orbtype).ne.'psig') ) & 
      call errore('get_orbitals','unknown orbital type.',1)
    if( trim(ftype) == 'esh5' ) then
      call psi_from_esh5(h5id,orbtype,dfft,Psi,b_beg,nbands,ik,ispin,Q,qkmap,becpsi)
    else if( trim(ftype) == 'davcio' ) then
      call psi_from_davcio(orbtype,dfft,Psi,b_beg,nbands,ik,ispin,Q,qkmap,becpsi)
    else
      call errore('get_orbitals','unknown file type.',1)
    endif
    !
  END SUBROUTINE get_orbitals
  !
  !-----------------------------------------------------------------------
  SUBROUTINE get_orbitals_set(h5id,ftype,orbtype,dfft,ispin,Psi,b_beg,& 
                nbands,k_beg,nkpts,Q,qkmap,becpsi,nminK)
  !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id
    CHARACTER(len=*), INTENT(IN) :: ftype, orbtype
    INTEGER,      INTENT(IN)   :: b_beg,nbands,k_beg,nkpts,ispin
    COMPLEX(DP),  INTENT(OUT)  :: Psi(:,:,:)
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN), OPTIONAL :: Q
    INTEGER, INTENT(IN), OPTIONAL :: qkmap(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: nminK(:)
    TYPE(bec_type), OPTIONAL :: becpsi(:)
    !
    INTEGER :: ik, ka, b0
    !
    b0 = b_beg
    do ik=1,nkpts
      if(present(nminK)) b0 = nminK(ik)+1 
      ka = ik + k_beg - 1
      if(present(becpsi)) then
        call get_orbitals(h5id,ftype,orbtype,dfft,Psi(:,:,ik),b0,nbands,ka,ispin, &
                      Q,qkmap,becpsi(ik)) 
      else
        call get_orbitals(h5id,ftype,orbtype,dfft,Psi(:,:,ik),b0,nbands,ka,ispin, &
                      Q,qkmap)
      endif
    enddo
    !
  END SUBROUTINE get_orbitals_set
  !
  !-----------------------------------------------------------------------
  SUBROUTINE psi_from_esh5(h5id,orbtype,dfft,Psi,b_beg,nbands,ik,ispin, &
                           Q,qkmap,becpsi)
  !-----------------------------------------------------------------------
    !
    USE fft_interfaces,  ONLY : invfft
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id
    CHARACTER(len=*), INTENT(IN) :: orbtype
    INTEGER,      INTENT(IN)   :: b_beg,nbands,ik,ispin
    COMPLEX(DP),  INTENT(OUT)  :: Psi(:,:)
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN), OPTIONAL :: Q
    INTEGER, INTENT(IN), OPTIONAL :: qkmap(:,:)  
    TYPE(bec_type), OPTIONAL :: becpsi
    !
    TYPE(bec_type) :: becdummy
    INTEGER :: npw,kb,ibnd,error, nktot
    REAL(DP) :: dk(3)
    !
    if(ispin > h5id%nspin) &
      call errore('psi_from_esh5',' ispin > h5id%nspin. ',1)   
    !
    if(PRESENT(becpsi)) call allocate_bec_type(nkb,1,becdummy)
    if( PRESENT(Q) .and. (.not.PRESENT(qkmap))) & 
        call errore('read_orbitals::read_from_esh5 both Q,qkmap. ',1)
    if( PRESENT(qkmap) .and. (.not.PRESENT(Q))) & 
        call errore('read_orbitals::read_from_esh5 both Q,qkmap. ',1)
    !
    Psi(:,:)=(0.d0,0.d0)
    if(PRESENT(becpsi)) then
      if(gamma_only) then
        becpsi%r(:,:) = CZERO  
      elseif(noncolin) then
        becpsi%nc(:,:,:) = CZERO  
      else
        becpsi%k(:,:) = CZERO  
      endif
    endif
    !
    if( PRESENT(Q) ) then
      kb = qkmap(Q,ik)
    else
      kb = ik
    endif  
    !
    CALL gk_sort (h5id%xk (1:3, kb), ngm, g, ecutwfc / tpiba2, &
                  npw, igk(1), g2kin)
    !
    if(PRESENT(becpsi)) then
      CALL init_us_2 (npw, igk(1), h5id%xk(1,kb), vkb)
    endif
    !
    do ibnd=1,nbands
      !
      if( b_beg+ibnd-1 > h5id%norbK(kb+h5id%nkpts*(ispin-1)) ) cycle 
      !
      call esh5_posthf_read(h5id%id,kb-1+h5id%nkpts*(ispin-1),b_beg+ibnd-2,1,evc,npwx,error)
      if(error .ne. 0 ) &
        call errore('psi_from_esh5','error reading orbital',1)
      !
      if( trim(orbtype) .eq. 'psir' ) then
        !
        Psi(dfft%nl(igk(1:npw)),ibnd)=evc(1:npw,1)
        if(gamma_only) &
            Psi(dfft%nlm(igk(1:npw)),ibnd) = CONJG(evc(1:npw,1))
        !
        CALL invfft ('Wave', Psi(1:npwx,ibnd), dfft)
        !
        if(noncolin .and. h5id%npol == 2) then
          !
          Psi(npwx+dfft%nl(igk(1:npw)),ibnd)=evc(npwx+1:npwx+npw,1)
          if(gamma_only) &
              Psi(npwx+dfft%nlm(igk(1:npw)),ibnd) = CONJG(evc(npwx+1:npwx+npw,1))
          !
          CALL invfft ('Wave', Psi(1:npwx,ibnd), dfft)
          !
        endif
        !
      elseif( trim(orbtype) .eq. 'psig' ) then
        !
        Psi(1:npw,ibnd)=evc(1:npw,1)
        if(noncolin .and. h5id%npol == 2) Psi(npwx+1:npwx+npw,ibnd)=evc(npwx+1:npwx+npw,1)
        !
      endif
      !
      if(PRESENT(becpsi)) then
        ! don't know how to do this otherwise  
        CALL  calbec( npw, vkb, evc(:,1:1), becdummy, 1)  
        if(gamma_only) then
          becpsi%r(:,ibnd) = becdummy%r(:,1)
        elseif(noncolin) then
          becpsi%nc(:,:,ibnd) = becdummy%nc(:,:,1)
        else
          becpsi%k(:,ibnd) = becdummy%k(:,1)
        endif
      endif
      !
    enddo
    !
    if(PRESENT(becpsi)) call deallocate_bec_type(becdummy)
    !
  END SUBROUTINE psi_from_esh5
  !  
  !-----------------------------------------------------------------------
  SUBROUTINE psi_from_davcio(orbtype,dfft,Psi,b_beg,nbands,ik,ispin, &
                           Q,qkmap,becpsi)
  !-----------------------------------------------------------------------
    !
    USE fft_interfaces,  ONLY : invfft
    !
    IMPLICIT NONE
    !
    CHARACTER(len=*), INTENT(IN) :: orbtype
    INTEGER,      INTENT(IN)   :: b_beg,nbands,ik,ispin
    COMPLEX(DP),  INTENT(OUT)  :: Psi(:,:)
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    INTEGER, INTENT(IN), OPTIONAL :: Q
    INTEGER, INTENT(IN), OPTIONAL :: qkmap(:,:)  
    TYPE(bec_type), OPTIONAL :: becpsi
    !
    INTEGER :: npw,kb,ibnd,error, nktot
    REAL(DP) :: dk(3)
    !
    if( PRESENT(Q) .and. (.not.PRESENT(qkmap))) & 
        call errore('read_orbitals::read_from_esh5 both Q,qkmap. ',1)
    if( PRESENT(qkmap) .and. (.not.PRESENT(Q))) & 
        call errore('read_orbitals::read_from_esh5 both Q,qkmap. ',1)
    !
    Psi(:,:)=CZERO
    if(present(becpsi)) then
      if(gamma_only) then
        becpsi%r(:,:) = CZERO
      elseif(noncolin) then
       becpsi%nc(:,:,:) = CZERO
      else
        becpsi%k(:,:) = CZERO
      endif
    endif
    !
    if( PRESENT(Q) ) then
      kb = qkmap(Q,ik)
    else
      kb = ik
    endif  
    !
    CALL gk_sort (xksym (1:3, kb), ngm, g, ecutwfc / tpiba2, &
                  npw, igk(1), g2kin)
    !
    if(PRESENT(becpsi)) then
      CALL init_us_2 (npw, igk(1), xksym(1,kb), vkb)
    endif
    !
    CALL davcio (evc, 2*nwordwfc, iunwfc, kb+nksym*(ispin-1), - 1)
    !
    if(PRESENT(becpsi)) then
      CALL  calbec( npw, vkb, evc(:,b_beg:min(b_beg+nbands-1,norb)), &
                    becpsi, min(b_beg+nbands-1,norb)-b_beg)
    endif
    !
    do ibnd=1,nbands
      !
      if( b_beg+ibnd-1 > norb ) cycle 
      !
      if( trim(orbtype) .eq. 'psir' ) then
        !
        Psi(dfft%nl(igk(1:npw)),ibnd)=evc(1:npw,b_beg+ibnd-1)
        if(gamma_only) &
            Psi(dfft%nlm(igk(1:npw)),ibnd) = CONJG(evc(1:npw,b_beg+ibnd-1))
        !
        CALL invfft ('Wave', Psi(1:npwx,ibnd), dfft)
        !
        if(noncolin) then
        !
          Psi(npwx+dfft%nl(igk(1:npw)),ibnd)=evc(npwx+1:npwx+npw,b_beg+ibnd-1)
          if(gamma_only) &
              Psi(npwx+dfft%nlm(igk(1:npw)),ibnd) = CONJG(evc(npwx+1:npwx+npw,b_beg+ibnd-1))
          !
          CALL invfft ('Wave', Psi(npwx+1:2*npwx,ibnd), dfft)
          !
        endif
        !
      elseif( trim(orbtype) .eq. 'psig' ) then
        !
        Psi(1:npw,ibnd)=evc(1:npw,b_beg+ibnd-1)
        if(noncolin) Psi(npwx+1:npwx+npw,ibnd)=evc(npwx+1:npwx+npw,b_beg+ibnd-1)
        !
      endif
      !
    enddo
    !
  END SUBROUTINE psi_from_davcio
  !
  !-----------------------------------------------------------------------
  SUBROUTINE get_psi_esh5(h5id,ibnd,ik,ispin,buff)
  !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id
    INTEGER,      INTENT(IN)   :: ibnd,ik,ispin
    COMPLEX(DP),  INTENT(OUT)  :: buff(:)
    INTEGER :: ikk, error
    !
    if(ispin > h5id%nspin) &
      call errore('get_psi_esh5',' ispin > h5id%nspin. ',1)   
    !
    buff(:)=(0.d0,0.d0)
    ikk = ik + h5id%nkpts*(ispin-1)
    call esh5_posthf_read(h5id%id,ikk-1,ibnd-1,1,buff,size(buff,1),error)
    if(error .ne. 0 ) &
      call errore('get_psi_esh5','error reading orbital',11)
    !
  END SUBROUTINE get_psi_esh5
  !
  !
  ! adding writing routines here too
  !  
  SUBROUTINE open_esh5_write(h5id,dfft,h5name,n_spins,write_psir,npolarization)
  !  
    USE cell_base, ONLY: alat, tpiba, at, bg
    USE fft_types, ONLY: fft_type_descriptor
    ! 
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(INOUT) :: h5id
    CHARACTER(len=*), INTENT(IN) :: h5name 
    INTEGER, INTENT(IN) :: n_spins
    LOGICAL, INTENT(IN), OPTIONAL :: write_psir 
    INTEGER, INTENT(IN), OPTIONAL :: npolarization
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    !
    INTEGER :: ik, h5len, error, npol_
    REAL(DP), ALLOCATABLE :: xkcart_(:,:)
    REAL(DP) :: recv(3,3), at0(3,3)
    ! 

    npol_ = 1
    if(present(npolarization)) npol_ = npolarization

    if(.not.allocated(h5id%norbK)) then
      allocate(h5id%norbK(nksym*n_spins))
    else if(size(h5id%norbK,1) .ne. n_spins*nksym) then
      deallocate(h5id%norbK)
      allocate(h5id%norbK(nksym*n_spins))
    endif
    allocate( h5id%xk(3,nksym) )
    do ik=1,nksym
      h5id%xk(1:3,ik) = xksym(1:3,ik)*tpiba
    enddo
    recv(1:3,1:3) = bg(1:3,1:3) * tpiba
    at0(1:3,1:3) = at(1:3,1:3) * alat
    h5len = LEN_TRIM(h5name)
    CALL esh5_posthf_open_write(h5id%id,h5name,h5len, error)
    if(error .ne. 0 ) &
        call errore('open_esh5_write','error opening orbital file for write',1)
    CALL esh5_posthf_write_meta(h5id%id,orbsG,5,nksym,n_spins,npol_,npwx,h5id%xk, &
                        1,dfft%nr1,dfft%nr2,dfft%nr3,at0,recv,alat,error)
    if(error .ne. 0 ) &
        call errore('open_esh5_write','error writing meta data OrbG',1)
    if(PRESENT(write_psir)) then
      if(write_psir) then
        CALL esh5_posthf_write_meta(h5id%id,orbsR,5,nksym,n_spins,npol_,0,h5id%xk, &
                          0,dfft%nr1,dfft%nr2,dfft%nr3,at0,recv,alat,error)
        if(error .ne. 0 ) &
          call errore('pw2posthf','error writing meta data OrbR',1)
      endif
    endif
    h5id%nspin = n_spins
    h5id%npol = npol_
  !  
  END SUBROUTINE open_esh5_write 
  !
  SUBROUTINE close_esh5_write(h5id)
    ! 
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(INOUT) :: h5id
    !
    call esh5_posthf_close_write(h5id)
    !
    if(allocated(h5id%norbK)) deallocate(h5id%norbK)
    if(allocated(h5id%xk)) deallocate(h5id%xk)
    h5id%id = int(-1,kind=8)
    h5id%grid_type = -1
    h5id%nr1 = -1
    h5id%nr2 = -1
    h5id%nr3 = -1
    h5id%npwx=-1
    h5id%npol=-1
    h5id%nmax_DM=-1
    !
  END SUBROUTINE close_esh5_write
  !
  subroutine update_bands(h5name,dfft)
    !
    USE parallel_include
    USE fft_types, ONLY: fft_type_descriptor
    USE lsda_mod,             ONLY : nspin, isk, lsda
    USE klist,                ONLY : nkstot, wk, nks, xk, ngk, igk_k
    USE wvfct,                ONLY : npwx, et, wg, nbnd
    USE io_files,             ONLY : nwordwfc, iunwfc
    USE mp_pools,             ONLY : intra_pool_comm, me_pool, root_pool, &
                                     inter_pool_comm, npool
    USE mp_images, ONLY: intra_image_comm
    USE mp,           ONLY: mp_barrier
    !
    IMPLICIT NONE
    !
    CHARACTER(len=256), INTENT(IN) :: h5name
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    !
    TYPE(h5file_type) :: h5id
    INTEGER, EXTERNAL    :: global_kpoint_index
    INTEGER :: i, j, ispin, iks, ik, ikk, n, error, info
    !
    if(noncolin) call errore('update_bands','No nonconlin yet.',1)
    !
    if(me_pool == root_pool) then 
      !
      call open_esh5_read(h5id,TRIM(h5name))
      if( h5id%grid_type .ne. 1) &
          call errore('update_bands','grid_type ne 1',1)
      ! 
      ! this is probably wrong with lsda
      !if(lsda) write(*,*) 'WARNING: update_bands with lsda, check check check!'
      iks = global_kpoint_index (nkstot, 1)
      ! it is probably better to loop over kpoints and spins, instead of nks
      ! it is confusing the way QE wraps the spin index in lsda
      do ik=1,nks
        !
        ikk = ik + iks - 1  ! assuming spin independent basis, which so far is
                          ! always the case 
        IF ( lsda ) then
          ispin = isk(ik)
          if(ispin == 2) ikk = ikk - nksym
        endif  
        evc(:,:) = CZERO  
        if(nbnd > h5id%norbK(ikk)) & 
          CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        n = min(nbnd,h5id%norbK(ikk))
        do i=1,n
          call get_psi_esh5(h5id,i,ikk,1,evc(:,i))
        enddo
        !
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, +1)
        !
      enddo
      !
      call close_esh5_read(h5id)
      !
    endif
    !
    call punch('all')
  !
  end subroutine update_bands
  !  
!  subroutine project_qeorbs_on_basis_and_update(h5name,dfft)
!project qe orbitals on given single particle spin-independent basis
!then upload then (call davcio(...,+1); call punch('all') 
!this is meant to be the correct way to do update_bands on uhf and ghf
!calculations
!to allow for PGTO basis, generalize to allow some states to be added "as is"
!or maybe project just up to some value and then write as is
!  end subroutine project_qeorbs_on_basis_and_update
  !
END MODULE read_orbitals_from_file
