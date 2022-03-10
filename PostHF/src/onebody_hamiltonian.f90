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
    USE scf, ONLY: vltot
    USE wvfct, ONLY: npwx, g2kin
    USE wavefunctions, ONLY : psic
    USE cell_base, ONLY: tpiba2
    USE becmod,   ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
    USE noncollin_module, ONLY : noncolin, npol,\
      pointlist, factlist, i_cons
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
    COMPLEX(DP) :: ctemp, e1_pin
    INTEGER :: ia, ib, i0, no, error, npw
    INTEGER :: ik,ibnd,ispin
    INTEGER :: norb_ik, mxorb
    COMPLEX(DP) :: CZERO
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:) 
    COMPLEX(DP), ALLOCATABLE :: H1(:,:), Hpin(:,:), Hpp(:,:), H1copy(:,:)
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)
    COMPLEX(DP), ALLOCATABLE :: evc_(:,:)
    REAL(DP), ALLOCATABLE :: v_of_r(:,:)
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

    mxorb = h5id_orbs%maxnorb
    allocate( Orbitals(npwx, mxorb) )

    if(noncolin) &
      allocate( evc_(npol*npwx,npol*mxorb) )
    allocate( hpsi(npol*npwx, npol*mxorb) )
    CALL allocate_bec_type ( nkb, npol*mxorb, becp )
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
        CALL fillH1(H1, hpsi, Orbitals, norb_ik, mxorb, npw)
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
      CALL fillH1(H1, hpsi, Orbitals, norb_ik, mxorb, npw)

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
      !
      if (i_cons == 1) then
        if ((nspin .ne. 2).and.(nspin .ne. 4)) then
          call errore('onebody', 'nspin b field not implemented', 1)
        endif
        if (gamma_only) then
          call errore('onebody', 'gamma b field not implemented', 1)
        endif
        ! partition real-space FFT grid among lattice sites
        ALLOCATE( pointlist(dfft%nnr) )
        ALLOCATE( factlist(dfft%nnr)  )
        CALL make_pointlists()
        ! add magnetic field
        print*, 'creating external potential'
        !v%of_r(:, :) = 0.d0
        !CALL add_bfield( v%of_r, rho%of_r )
        !CALL report_mag()
        ALLOCATE( v_of_r(dfft%nnr,nspin) )
        CALL add_initial_bfield(v_of_r)
        !
        print*, 'adding external potential'
        allocate( Hpin(norb_ik,norb_ik) )
        Hpin(:,:) = CZERO
        if (noncolin) then
          allocate( Hpp(npol*norb_ik,npol*norb_ik) )
          ia = norb_ik+1
          ib = npol*norb_ik
          ! x
          CALL vlocalH1(Hpin, v_of_r(:,2), Orbitals, dfft, norb_ik, npw)
          Hpp(:,:) = CZERO
          Hpp(ia:ib,1:norb_ik) = Hpin(:,:)
          Hpp(1:norb_ik,ia:ib) = Hpin(:,:)
          do ispin=1,min(2,nspin)
            i0 = (ispin-1)*(npol-1)*norb_ik
            e1_pin = onebody_energy(Hpp,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
            print*, 'ispin, E1 with xpin:', ispin, e1_pin
          enddo
          H1(:,:) = H1(:,:)+Hpp(:,:)
          ! y
          CALL vlocalH1(Hpin, v_of_r(:,3), Orbitals, dfft, norb_ik, npw)
          Hpp(:,:) = CZERO
          Hpp(ia:ib,1:norb_ik) = (0.d0, -1.d0)*Hpin(:,:)
          Hpp(1:norb_ik,ia:ib) = (0.d0,  1.d0)*Hpin(:,:)
          do ispin=1,min(2,nspin)
            i0 = (ispin-1)*(npol-1)*norb_ik
            e1_pin = onebody_energy(Hpp,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
            print*, 'ispin, E1 with ypin:', ispin, e1_pin
          enddo
          H1(:,:) = H1(:,:)+Hpp(:,:)
          ! z
          CALL vlocalH1(Hpin, v_of_r(:, 4), Orbitals, dfft, norb_ik, npw)
          Hpp(:,:) = CZERO
          Hpp(1:norb_ik,1:norb_ik) = Hpin(:,:)
          Hpp(ia:ib,ia:ib) = (-1.d0, 0.d0)*Hpin(:,:)
          do ispin=1,min(2,nspin)
            i0 = (ispin-1)*(npol-1)*norb_ik
            e1_pin = onebody_energy(Hpp,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
            print*, 'ispin, E1 with zpin:', ispin, e1_pin
          enddo
          H1(:,:) = H1(:,:)+Hpp(:,:)
          ! total
          e1_pin = CZERO
          do ispin=1,min(2,nspin)
            i0 = (ispin-1)*(npol-1)*norb_ik
            e1_pin = e1_pin+onebody_energy(H1,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
          enddo
          print*, 'E1 with pin:', e1_pin
          ! transpose to account for expected row major format in esh5
          do ia=1,npol*norb_ik
            do ib=ia+1,npol*norb_ik
              ctemp = H1(ia,ib)
              H1(ia,ib) = H1(ib,ia)
              H1(ib,ia) = ctemp
            enddo
          enddo
          CALL esh5_posthf_write_h1(h5id_hamil%id,npol*norb_ik,ik,H1)
        else ! collinear !!!! HACK: output 2*nkpt H1
          allocate( H1copy(npol*norb_ik, npol*norb_ik) )
          H1copy(:,:) = H1(:,:)
          do ispin=1,2
            CALL vlocalH1(Hpin, v_of_r(:,ispin), Orbitals, dfft, norb_ik, npw)
            H1(:,:) = H1copy(:,:) + Hpin(:,:)
            i0 = 0
            e1_pin = onebody_energy(H1,DM(:,:,ik,ispin),i0,i0,norb_ik,nmax_DM)
            print*, 'E1 up/dn=', ispin, 'with pin:', e1_pin
            ! transpose to account for expected row major format in esh5
            do ia=1,npol*norb_ik
              do ib=ia+1,npol*norb_ik
                ctemp = H1(ia,ib)
                H1(ia,ib) = H1(ib,ia)
                H1(ib,ia) = ctemp
              enddo
            enddo
            CALL esh5_posthf_write_h1(h5id_hamil%id,npol*norb_ik,ik+(ispin-1)*(nksym),H1)
          enddo ! ispin
        endif ! noncolin
      else ! no magnetic constraint
        ! transpose to account for expected row major format in esh5
        do ia=1,npol*norb_ik
          do ib=ia+1,npol*norb_ik
            ctemp = H1(ia,ib)
            H1(ia,ib) = H1(ib,ia)
            H1(ib,ia) = ctemp
          enddo
        enddo
        CALL esh5_posthf_write_h1(h5id_hamil%id,npol*norb_ik,ik,H1)
        !
      endif ! i_cons
      deallocate(H1)  
      !
    end do
    
    IF( ALLOCATED(hpsi) ) DEALLOCATE (hpsi)
    IF( ALLOCATED(evc_) ) DEALLOCATE (evc_)
    if(allocated(Orbitals)) deallocate(Orbitals)
    CALL deallocate_bec_type (becp)
    if( allocated(pointlist) ) deallocate( pointlist )
    if( allocated(factlist) ) deallocate( factlist )
    if( allocated(v_of_r) ) deallocate( v_of_r )
    IF( allocated(Hpin) ) deallocate( Hpin )
    IF( allocated(Hpp) ) deallocate( Hpp )
    IF( allocated(H1copy) ) deallocate( H1copy )
 
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

  subroutine vlocalH1(H1loc, v_of_r, Orbitals, dfft, norb, npw)
    ! construct one-body hamiltonian from local operator on basis
    !
    ! Inputs:
    !   v_of_r (array): (nnr,) local potential
    !   Orbitals (array): (npw, norb) basis functions in PW
    !   dfft (fft_type_descriptor): discrete FFT mesh
    ! Return:
    !   H1loc (array): (norb, norb) one-body hamiltonian
    !
    USE fft_types, ONLY: fft_type_descriptor
    use fft_interfaces, ONLY : invfft, fwfft
    USE wavefunctions, ONLY : psic
    USE wvfct, ONLY: npwx
    USE posthf_mod, ONLY: igksym, e2Ha
    ! inputs and outputs
    complex(DP), intent(out) :: H1loc(norb, norb)
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    real(DP), intent(in) :: v_of_r(dfft%nnr)
    complex(DP), intent(in) :: Orbitals(npw, norb)
    integer, intent(in) :: norb, npw
    ! local variables
    COMPLEX(DP) :: CZERO, CONE, CNORM
    COMPLEX(DP) :: hpsi_loc(npw, norb)
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
    integer :: ibnd
    CZERO = (0.d0,0.d0)
    CONE = (1.d0,0.d0)
    CNORM = CONE*e2Ha
    hpsi_loc(:,:) = CZERO
    allocate(spsi(1,1))  ! to avoid issues in debugging mode
    !
    ! apply local potential in real space
    do ibnd = 1, norb
      psic (:) = (0.d0,0.d0)
      psic (dfft%nl(igksym(1:npw))) = Orbitals(1:npw,ibnd)
      !
      CALL invfft ('Wave', psic, dfft)
      !
      psic (1:dfft%nnr) = psic (1:dfft%nnr) * v_of_r(1:dfft%nnr)
      !
      CALL fwfft ('Wave', psic, dfft)
      !
      hpsi_loc(1:npw, ibnd) = hpsi_loc(1:npw, ibnd) + psic(dfft%nl(igksym(1:npw)))
    end do ! ibnd
    CALL Overlap(norb,norb,npw,CNORM,Orbitals(1,1),npwx,&
            hpsi_loc(1,1),npwx,CZERO,H1loc,norb,.false.,spsi)
    IF( allocated(spsi) ) deallocate(spsi)
  end subroutine vlocalH1

END MODULE onebody_hamiltonian 

