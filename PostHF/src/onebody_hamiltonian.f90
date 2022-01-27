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
  SUBROUTINE getH1(h5id_orbs,h5id_hamil,dfft,e1,e1_so,e1_mf,e1_so_mf)
    USE scf,      ONLY: vltot
    USE wvfct, ONLY: nbnd, npwx, current_k, g2kin
    USE wavefunctions, ONLY : psic
    USE io_global, ONLY: stdout, ionode,  ionode_id
    USE cell_base, ONLY: tpiba2
    USE becmod,   ONLY : bec_type,becp, calbec, allocate_bec_type, deallocate_bec_type
    USE noncollin_module,     ONLY : noncolin, npol
    USE control_flags, ONLY : gamma_only
    USE gvect, ONLY: ngm, g, gstart
    USE gvecw, ONLY : ecutwfc
    USE io_files, ONLY: nwordwfc, iunwfc, tmp_dir, prefix
    USE lsda_mod, ONLY: lsda, nspin
    USE uspp,     ONLY : vkb, nkb
    USE realus,   ONLY : real_space
    USE wavefunctions, ONLY : evc
    USE lsda_mod,      ONLY: current_spin
    USE paw_variables, ONLY : okpaw
    USE uspp,       ONLY : okvan    
    use fft_interfaces,       ONLY : invfft, fwfft
    USE fft_types, ONLY: fft_type_descriptor
    USE posthf_mod, ONLY: nksym, nkfull, xkfull, xksym, &
                        igksym, ngksym,e2Ha, nmax_DM, DM, DM_mf
    USE orbital_generators, ONLY: mixed_basis
    USE read_orbitals_from_file, ONLY: h5file_type
    USE constants, ONLY: BOHR_RADIUS_ANGS
    USE moire, ONLY: lmoire, amoire
    !
    IMPLICIT NONE
    !
    TYPE(h5file_type), INTENT(IN) :: h5id_hamil,h5id_orbs
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: e1,e1_so
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: e1_mf,e1_so_mf
    !
    COMPLEX(DP) :: ctemp
    INTEGER :: ia, ib, i0, no, error, npw,npw2
    INTEGER :: ik,ibnd, ikk, ispin
    REAL(DP) :: fac, wm
    COMPLEX(DP) :: CONE, CZERO, CNORM
    COMPLEX(DP), ALLOCATABLE :: Orbitals(:,:) 
    COMPLEX(DP), ALLOCATABLE :: H1(:,:)
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)
    COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
    COMPLEX(DP), ALLOCATABLE :: evc_(:,:)
    !
    if (lmoire) then
      wm = 1.d0/amoire/amoire
    else ! default band width
      wm = 1.d0
    endif

    CONE = (1.d0,0.d0)
    CZERO = (0.d0,0.d0)
    CNORM = CONE*e2Ha 
    fac=1.d0
    if(nspin==1) fac=2.d0
    if(present(e1)) e1 = CZERO
    if(present(e1_so)) e1_so = CZERO
    if(present(e1_mf)) e1_mf = CZERO
    if(present(e1_so_mf)) e1_so_mf = CZERO

    ! if paw/uspp, kill Hxc terms 
    ! these (Hartree/EXX) are added in the 2-body part 
    if(okvan .or. okpaw) call reset_deeq()

    allocate( Orbitals(npwx, h5id_orbs%maxnorb) )
    allocate(spsi(1,1))  ! to avoid issues in debugging mode

    if(noncolin) &
      allocate( evc_(npol*npwx,npol*h5id_orbs%maxnorb) )
    allocate( hpsi(npol*npwx, npol*h5id_orbs%maxnorb) )
    CALL allocate_bec_type ( nkb, npol*h5id_orbs%maxnorb, becp )
    write(*,*) 'USPP nkb:',nkb

    ! set spin, only spin = 1 is used
    current_spin = 1

    do ik=1,nksym
      !
      allocate( H1(npol*h5id_orbs%norbK(ik),npol*h5id_orbs%norbK(ik) ) )
      !
      Orbitals(:,:) = CZERO  
      !  
      call esh5_posthf_read(h5id_orbs%id,ik-1,0,h5id_orbs%norbK(ik),Orbitals,npwx,error)
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
      DO ibnd = 1, h5id_orbs%norbK(ik)
        hpsi (1:npw, ibnd) = wm*tpiba2 * g2kin (1:npw) * Orbitals(1:npw,ibnd) 
        if(noncolin) hpsi (npwx+1:npwx+npw, ibnd+h5id_orbs%norbK(ik)) = wm*tpiba2 * &
                                                g2kin (1:npw) * Orbitals(1:npw,ibnd)
      END DO
      !
      ! VNL: non local contribution
      !
      IF (nkb > 0 .AND. .NOT. real_space) THEN
        !
        CALL init_us_2 (npw, igksym(1), xksym (1, ik), vkb) 
        if(noncolin) then
          evc_(:,:) = CZERO 
          evc_(1:npwx,1:h5id_orbs%norbK(ik)) =  Orbitals(1:npwx,1:h5id_orbs%norbK(ik))  
          evc_(npwx+1:2*npwx,1+h5id_orbs%norbK(ik):2*h5id_orbs%norbK(ik)) = &
                                                Orbitals(1:npwx,1:h5id_orbs%norbK(ik)) 
          evc_(npwx+1:2*npwx,1:h5id_orbs%norbK(ik)) = CZERO
          evc_(1:npwx,1+h5id_orbs%norbK(ik):2*h5id_orbs%norbK(ik)) = CZERO 
          call calbec(npw, vkb, evc_(:,1:2*h5id_orbs%norbK(ik)), becp, 2*h5id_orbs%norbK(ik))
        else
          CALL calbec(npw, vkb, Orbitals(:,1:h5id_orbs%norbK(ik)), becp, h5id_orbs%norbK(ik) )
        endif
        CALL add_vuspsi( npwx, npw, npol*h5id_orbs%norbK(ik), hpsi )
        !
      END IF
      !
      ! local potential
      !
      do ibnd=1,h5id_orbs%norbK(ik)
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
          hpsi (npwx+1:npwx+npw, ibnd+h5id_orbs%norbK(ik)) = &  
                                hpsi (npwx+1:npwx+npw, ibnd+h5id_orbs%norbK(ik)) + &
                                psic (dfft%nl(igksym(1:npw))) 
        !
      enddo
      if(gamma_only .AND. gstart == 2 ) & 
        hpsi(1,1:h5id_orbs%norbK(ik)) = CMPLX( DBLE( hpsi(1,1:h5id_orbs%norbK(ik)) ), &
                                              0.D0 ,kind=DP)
      ! H1(a,b) = 1/Ni sum_i conjg(PsiL(i,a)) * PsiR(i, b)
      CALL Overlap(h5id_orbs%norbK(ik),h5id_orbs%norbK(ik),npw,CNORM,Orbitals(1,1),npwx,&
               hpsi(1,1),npol*npwx,CZERO,H1,npol*h5id_orbs%norbK(ik),.false.,spsi)
      if(noncolin) then
        ! (down,up)
        CALL Overlap(h5id_orbs%norbK(ik),h5id_orbs%norbK(ik),npw,CNORM,Orbitals(1,1),npwx,&
               hpsi(npwx+1,1),npol*npwx,CZERO,H1(h5id_orbs%norbK(ik)+1,1), &
               npol*h5id_orbs%norbK(ik),.false.,spsi)
        ! (up,down)
        CALL Overlap(h5id_orbs%norbK(ik),h5id_orbs%norbK(ik),npw,CNORM,Orbitals(1,1),npwx,&
               hpsi(1,h5id_orbs%norbK(ik)+1),npol*npwx,CZERO,H1(1,h5id_orbs%norbK(ik)+1), &
               npol*h5id_orbs%norbK(ik),.false.,spsi)
        ! (down,down)
        CALL Overlap(h5id_orbs%norbK(ik),h5id_orbs%norbK(ik),npw,CNORM,Orbitals(1,1),npwx,&
               hpsi(npwx+1,h5id_orbs%norbK(ik)+1),npol*npwx,CZERO, &
               H1(h5id_orbs%norbK(ik)+1,h5id_orbs%norbK(ik)+1), &
               npol*h5id_orbs%norbK(ik),.false.,spsi)
        do ia=1,2*h5id_orbs%norbK(ik)
          do ib=ia+1,2*h5id_orbs%norbK(ik)
            if( abs(H1(ia,ib)-conjg(H1(ib,ia))) > 1e-6 ) &
              write(*,*)ia,ib,H1(ia,ib),H1(ib,ia)
          enddo
        enddo
      endif  

      no = min(nmax_DM,h5id_orbs%norbK(ik))
      if(present(e1)) then
        do ispin=1,min(2,nspin)
          i0 = (ispin-1)*(npol-1)*h5id_orbs%norbK(ik)
          do ib=1,no
            do ia=1,no
              e1 = e1 + fac*H1(ia+i0,ib+i0)*DM(ia,ib,ik,ispin)
            enddo
          enddo
        enddo
      endif
      if(present(e1_mf)) then
        do ispin=1,min(2,nspin)
          i0 = (ispin-1)*(npol-1)*h5id_orbs%norbK(ik)
          do ib=1,no
            do ia=1,no
              e1_mf = e1_mf + fac*H1(ia+i0,ib+i0)*DM_mf(ia,ib,ik,ispin)
            enddo
          enddo
        enddo
      endif
      if(present(e1_so)) then
        if(noncolin) then
          do ib=1,no
            do ia=1,no
              e1_so = e1_so + fac*H1(ia,ib+h5id_orbs%norbK(ik))*DM(ia,ib,ik,3)
            enddo
          enddo
          do ib=1,no
            do ia=1,no
              e1_so = e1_so + fac*H1(ia+h5id_orbs%norbK(ik),ib)*DM(ia,ib,ik,4)
            enddo
          enddo
        endif
      endif
      if(present(e1_so_mf)) then
        if(noncolin) then
          do ib=1,no
            do ia=1,no
              e1_so_mf = e1_so_mf + fac*H1(ia,ib+h5id_orbs%norbK(ik))*DM_mf(ia,ib,ik,3)
            enddo
          enddo
          do ib=1,no
            do ia=1,no
              e1_so_mf = e1_so_mf + fac*H1(ia+h5id_orbs%norbK(ik),ib)*DM_mf(ia,ib,ik,4)
            enddo
          enddo
        endif
      endif

      ! collect on head node!
      ! transpose to account for expected row major format in esh5  
      do ia=1,npol*h5id_orbs%norbK(ik) 
        do ib=ia+1,npol*h5id_orbs%norbK(ik) 
          ctemp = H1(ia,ib)
          H1(ia,ib) = H1(ib,ia)
          H1(ib,ia) = ctemp 
        enddo
      enddo
      !
601   CONTINUE
      !
      CALL esh5_posthf_write_h1(h5id_hamil%id,npol*h5id_orbs%norbK(ik),ik,H1)
      !
      deallocate(H1)  
      !
    end do
    
    IF( ALLOCATED(hpsi) ) DEALLOCATE (hpsi)
    IF( ALLOCATED(spsi) ) DEALLOCATE (spsi)
    IF( ALLOCATED(evc_) ) DEALLOCATE (evc_)
    if(allocated(Orbitals)) deallocate(Orbitals)
    CALL deallocate_bec_type (becp)
 
  END SUBROUTINE getH1

END MODULE onebody_hamiltonian 

