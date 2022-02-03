!--------------------------------------------------------------------
! Written by Miguel A. Morales, LLNL, 2020 
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE onebody_hamiltonian
  !----------------------------------------------------------------------
  ! 
  ! Implements MP2 and some other utilities, e.g. pseudo-canonicalization
! MAM: this routine can be parallelized now 
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE    
  !
  CONTAINS
  !
  ! Calculate and print one body hamiltonian
  !
  !-----------------------------------------------------------------------
  SUBROUTINE getH1(h5id_orbs,h5id_hamil,dfft,e1,e1_so,e1_mf,e1_so_mf,tkin)
    USE scf,      ONLY: vltot
    USE wvfct, ONLY: npwx, g2kin
    USE wavefunctions, ONLY : psic
    USE cell_base, ONLY: tpiba2
    USE becmod,   ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
    USE noncollin_module,     ONLY : noncolin, npol
    USE control_flags, ONLY : gamma_only
    USE gvect, ONLY: ngm, g, gstart
    USE gvecw, ONLY : ecutwfc
    USE lsda_mod, ONLY: nspin
    USE uspp,     ONLY : vkb, nkb
    USE realus,   ONLY : real_space
    USE wavefunctions, ONLY : evc
    USE lsda_mod,      ONLY: current_spin
    USE paw_variables, ONLY : okpaw
    USE uspp,       ONLY : okvan    
    use fft_interfaces,       ONLY : invfft, fwfft
    USE fft_types, ONLY: fft_type_descriptor
    USE posthf_mod, ONLY: nksym, xksym, igksym
    USE posthf_mod, ONLY: nmax_DM, DM, DM_mf
    USE read_orbitals_from_file, ONLY: h5file_type
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id_hamil,h5id_orbs
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: e1,e1_so,tkin
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: e1_mf,e1_so_mf
    !
    COMPLEX(DP) :: ctemp
    INTEGER :: ia, ib, i0, no, error, npw
    INTEGER :: ik,ibnd,ispin
    INTEGER :: norb_ik
    COMPLEX(DP) :: CZERO
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:) 
    COMPLEX(DP), ALLOCATABLE :: H1(:,:)
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)
    COMPLEX(DP), ALLOCATABLE :: evc_(:,:)
    !

    CZERO = (0.d0,0.d0)
    if(present(e1)) e1 = CZERO
    if(present(e1_so)) e1_so = CZERO
    if(present(e1_mf)) e1_mf = CZERO
    if(present(e1_so_mf)) e1_so_mf = CZERO
    if(present(tkin)) tkin = CZERO

    ! if paw/uspp, kill Hxc terms 
    ! these (Hartree/EXX) are added in the 2-body part 
    if(okvan .or. okpaw) call reset_deeq()

    allocate( Orbitals(npwx, h5id_orbs%maxnorb) )

    if(noncolin) &
      allocate( evc_(npol*npwx,npol*h5id_orbs%maxnorb) )
    allocate( hpsi(npol*npwx, npol*h5id_orbs%maxnorb) )
    CALL allocate_bec_type ( nkb, npol*h5id_orbs%maxnorb, becp )
    write(*,*) 'USPP nkb:',nkb

    ! set spin, only spin = 1 is used
    current_spin = 1

    do ik=1,nksym
      norb_ik = h5id_orbs%norbK(ik)
      no = min(nmax_DM,norb_ik)
      !
      allocate( H1(npol*norb_ik,npol*norb_ik) )
      !
      Orbitals(:,:) = CZERO  
      !  
      call esh5_posthf_read(h5id_orbs%id,ik-1,0,norb_ik,Orbitals,npwx,error)
      if(error .ne. 0 ) &
        call errore('getH1','error reading additional orbital',2)
      if(gamma_only .AND. gstart == 2 ) &
        Orbitals(1,:) = CMPLX( DBLE( Orbitals(1,:) ), 0.D0 ,kind=DP)
      !
      H1(:,:) = (0.d0,0.d0)
      !
      CALL gk_sort (xksym (1:3, ik), ngm, g, ecutwfc / tpiba2, &
                  npw, igksym(1), g2kin)
      !
      ! kinetic term
      !
      hpsi(:,:) = (0.d0,0.d0)  
      DO ibnd = 1, norb_ik
        hpsi (1:npw, ibnd) = tpiba2 * g2kin (1:npw) * Orbitals(1:npw,ibnd) 
        if(noncolin) hpsi (npwx+1:npwx+npw, ibnd+norb_ik) = tpiba2 * &
                                                g2kin (1:npw) * Orbitals(1:npw,ibnd)
      END DO
      ! kinetic energy
      if(present(tkin)) then
        CALL fillH1(H1, hpsi, Orbitals, norb_ik, h5id_orbs%maxnorb, npw)
        do ispin=1,min(2,nspin)
          i0 = (ispin-1)*(npol-1)*norb_ik
          tkin = tkin + onebody_energy(H1,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
        enddo
      endif
      !
      ! VNL: non local contribution
      !
      IF (nkb > 0 .AND. .NOT. real_space) THEN
        !
        CALL init_us_2 (npw, igksym(1), xksym (1, ik), vkb) 
        if(noncolin) then
          evc_(:,:) = CZERO
          evc_(1:npwx,1:norb_ik) =  Orbitals(1:npwx,1:norb_ik)
          evc_(npwx+1:2*npwx,1+norb_ik:2*norb_ik) = Orbitals(1:npwx,1:norb_ik)
          evc_(npwx+1:2*npwx,1:norb_ik) = CZERO
          evc_(1:npwx,1+norb_ik:2*norb_ik) = CZERO
          call calbec(npw, vkb, evc_(:,1:2*norb_ik), becp, 2*norb_ik)
        else
          CALL calbec(npw, vkb, Orbitals(:,1:norb_ik), becp, norb_ik )
        endif
        CALL add_vuspsi( npwx, npw, npol*norb_ik, hpsi )
        !
      END IF
      !
      ! local potential
      !
      do ibnd=1,norb_ik
        !
        psic (:) = (0.d0,0.d0)
        psic (dfft%nl(igksym(1:npw))) = Orbitals(1:npw,ibnd)
        if(gamma_only) psic (dfft%nlm(igksym(1:npw))) = CONJG(Orbitals(1:npw,ibnd))
        !
        CALL invfft ('Wave', psic, dfft)
        !
        ! vltot lives in dfftp, so doublegrid must be false for now
        psic (1:dfft%nnr) = psic (1:dfft%nnr) * vltot(1:dfft%nnr)
        !
        CALL fwfft ('Wave', psic, dfft)
        !
        hpsi (1:npw, ibnd) = hpsi (1:npw, ibnd) + psic (dfft%nl(igksym(1:npw)))
        if(noncolin) & 
          hpsi (npwx+1:npwx+npw, ibnd+norb_ik) = &  
                                hpsi (npwx+1:npwx+npw, ibnd+norb_ik) + &
                                psic (dfft%nl(igksym(1:npw))) 
        !
      enddo
      if(gamma_only .AND. gstart == 2 ) & 
        hpsi(1,1:norb_ik) = CMPLX( DBLE( hpsi(1,1:norb_ik) ), &
                                              0.D0 ,kind=DP)
      ! H1(a,b) = 1/Ni sum_i conjg(PsiL(i,a)) * PsiR(i, b)
      CALL fillH1(H1, hpsi, Orbitals, norb_ik, h5id_orbs%maxnorb, npw)

      if(present(e1)) then
        do ispin=1,min(2,nspin)
          i0 = (ispin-1)*(npol-1)*norb_ik
          e1 = e1 + onebody_energy(H1,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
        enddo
      endif
      if(present(e1_mf)) then
        do ispin=1,min(2,nspin)
          i0 = (ispin-1)*(npol-1)*norb_ik
          e1_mf = e1_mf + onebody_energy(H1,DM_mf(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
        enddo
      endif
      if(present(e1_so)) then
        if(noncolin) then
          ! (up,down)
          e1_so = e1_so + onebody_energy(H1,DM(:,:,ik,3),0,norb_ik,norb_ik,nmax_DM)
          ! (down,up)
          e1_so = e1_so + onebody_energy(H1,DM(:,:,ik,3),norb_ik,0,norb_ik,nmax_DM)
        endif
      endif
      if(present(e1_so_mf)) then
        if(noncolin) then
          ! (up,down)
          e1_so_mf = e1_so + onebody_energy(H1,DM_mf(:,:,ik,3),0,norb_ik,norb_ik,nmax_DM)
          ! (down,up)
          e1_so_mf = e1_so + onebody_energy(H1,DM_mf(:,:,ik,3),norb_ik,0,norb_ik,nmax_DM)
        endif
      endif

      ! collect on head node!
      ! transpose to account for expected row major format in esh5  
      do ia=1,npol*norb_ik
        do ib=ia+1,npol*norb_ik
          ctemp = H1(ia,ib)
          H1(ia,ib) = H1(ib,ia)
          H1(ib,ia) = ctemp 
        enddo
      enddo
      !
601   CONTINUE
      !
      CALL esh5_posthf_write_h1(h5id_hamil%id,npol*norb_ik,ik,H1)
      !
      deallocate(H1)  
      !
    end do
    
    IF( ALLOCATED(hpsi) ) DEALLOCATE (hpsi)
    IF( ALLOCATED(evc_) ) DEALLOCATE (evc_)
    if(allocated(Orbitals)) deallocate(Orbitals)
    CALL deallocate_bec_type (becp)
 
  END SUBROUTINE getH1

  subroutine fillH1(H1, hpsi, Orbitals, norb, mxorb, npw)
    use noncollin_module, only : noncolin, npol
    use wvfct, only: npwx
    USE posthf_mod, ONLY: e2Ha
    !
    integer, intent(in) :: norb, mxorb, npw
    complex(DP), intent(inout) :: H1(npol*norb,npol*norb)
    complex(DP), intent(in) :: hpsi(npol*npwx, npol*mxorb)
    complex(DP), intent(in) :: Orbitals(npwx, mxorb)
    ! local variables
    COMPLEX(DP) :: CZERO, CONE, CNORM
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
    integer :: ia, ib
    !
    CONE = (1.d0,0.d0)
    CZERO = (0.d0,0.d0)
    CNORM = CONE*e2Ha
    allocate(spsi(1,1))  ! to avoid issues in debugging mode
    !
    CALL Overlap(norb,norb,npw,CNORM,Orbitals(1,1),npwx,&
             hpsi(1,1),npol*npwx,CZERO,H1,npol*norb,.false.,spsi)
    if(noncolin) then
      ! (down,up)
      CALL Overlap(norb,norb,npw,CNORM,Orbitals(1,1),npwx,&
             hpsi(npwx+1,1),npol*npwx,CZERO,H1(norb+1,1), &
             npol*norb,.false.,spsi)
      ! (up,down)
      CALL Overlap(norb,norb,npw,CNORM,Orbitals(1,1),npwx,&
             hpsi(1,norb+1),npol*npwx,CZERO,H1(1,norb+1), &
             npol*norb,.false.,spsi)
      ! (down,down)
      CALL Overlap(norb,norb,npw,CNORM,Orbitals(1,1),npwx,&
             hpsi(npwx+1,norb+1),npol*npwx,CZERO, &
             H1(norb+1,norb+1), &
             npol*norb,.false.,spsi)
      do ia=1,2*norb
        do ib=ia+1,2*norb
          if( abs(H1(ia,ib)-conjg(H1(ib,ia))) > 1e-6 ) &
            write(*,*)ia,ib,H1(ia,ib),H1(ib,ia)
        enddo
      enddo
    endif ! noncollin
    !
    IF( allocated(spsi) ) deallocate(spsi)
  end subroutine fillH1

  pure complex(DP) function onebody_energy(h1, dm, ia0, ib0, norb, nmax)
    use noncollin_module, only : noncolin, npol
    USE lsda_mod, ONLY: nspin
    !
    complex(DP), intent(in) :: h1(npol*norb,npol*norb)
    complex(DP), intent(in) :: dm(nmax,nmax)
    integer, intent(in) :: ia0, ib0, norb, nmax
    ! local variables
    integer :: no, ia, ib
    real(DP) :: fac
    fac=1.d0
    if(nspin==1) fac=2.d0
    no = min(nmax,norb)
    !
    onebody_energy = (0.d0, 0.d0)
    do ib=1,no
      do ia=1,no
        onebody_energy = onebody_energy + fac*h1(ia+ia0,ib+ib0)*dm(ia,ib)
      enddo
    enddo
  end function

END MODULE onebody_hamiltonian 

