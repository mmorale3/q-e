!
! Written by Miguel A. Morales, LLNL, 2020
!
! Utility functions used by orbital-based pp routines, 
!  - mp2, cholesky, bscorr, etc
!
!
! From exx_base.f90
!-----------------------------------------------------------------------
SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
!-----------------------------------------------------------------------
  ! This routine calculates the 1/|r-r'| part of the exact exchange
  ! expression in reciprocal space (the G^-2 factor).
  ! It then regularizes it according to the specified recipe
  !
  USE KINDS, ONLY : DP
  USE constants, ONLY: tpi, fpi, e2, BOHR_RADIUS_ANGS, eps8
  USE cell_base, ONLY: tpiba2, at, alat, tpiba
  USE coulomb_vcut_module,  ONLY: vcut_get,  vcut_spheric_get
  USE posthf_mod, ONLY:  use_coulomb_vcut_ws, use_coulomb_vcut_spheric, &
                vcut, exxdiv
  USE input_parameters, ONLY: lmoire, amoire_in_ang, epsmoire
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in)    :: ngm   ! Number of G vectors
  REAL(DP), INTENT(in)    :: g(3,ngm) ! Cartesian components of G vectors
  REAL(DP), INTENT(in)    :: xk(3) ! current k vector
  REAL(DP), INTENT(in)    :: xkq(3) ! current q vector
  !
  COMPLEX(DP), INTENT(inout) :: fac(ngm) ! Calculated convolution
  !
  !Local variables
  INTEGER :: ig !Counters
  REAL(DP) :: q(3), qq, x
  LOGICAL :: odg(3)
  REAL(DP)         :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP) :: amoire
  !
  amoire = amoire_in_ang/BOHR_RADIUS_ANGS
  IF( use_coulomb_vcut_ws ) THEN
     DO ig = 1, ngm
        q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
        fac(ig) = vcut_get(vcut,q)
     ENDDO
     RETURN
  ENDIF
  !
  IF ( use_coulomb_vcut_spheric ) THEN
     DO ig = 1, ngm
        q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
        fac(ig) = vcut_spheric_get(vcut,q)
     ENDDO
     RETURN
  ENDIF
  !  
  DO ig=1,ngm
    !
    q(:)= xk(:) - xkq(:) + g(:,ig)
    qq = sum(q(:)**2) * tpiba2
    !
    IF (qq > eps_qdiv) THEN
      !
!       IF ( erfc_scrlen > 0  ) THEN
!          fac(ig)=e2*fpi/qq*(1._DP-exp(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
!       ELSEIF( erf_scrlen > 0 ) THEN
!          fac(ig)=e2*fpi/qq*(exp(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
!       ELSE
        if (lmoire) then
          if (abs(g(3,ig)) > eps8) then
            fac(ig) = 0.d0
          else
            fac(ig) = e2*tpi/sqrt(qq)/(epsmoire*amoire)*at(3,3)*alat
          endif
        else
          fac(ig)=e2*fpi/qq
        endif ! lmoire
!       ENDIF
      !
    ELSE
      !
      fac(ig) = - exxdiv
      !
    ENDIF
    !
  ENDDO
END SUBROUTINE g2_convolution

!-----------------------------------------------------------------------
SUBROUTINE calculate_phase_factor(dfft_, phasefac, dG)
!-----------------------------------------------------------------------
  !
  USE KINDS, ONLY : DP
  USE cell_base, ONLY: at, bg
  USE constants, ONLY: tpi
  USE fft_types, ONLY: fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft_
  COMPLEX(DP),INTENT(OUT)::phasefac(dfft_%nnr,27)
  REAL(DP), INTENT(OUT)   :: dG(3,27)
  INTEGER :: i,j,k,ni,nj,nk,Q,ir
  REAL(DP) :: r(3)

  Q = 1
  do i=-1,1
  do j=-1,1
  do k=-1,1
    dG(1:3,Q) = i * bg(1:3,1) + j * bg(1:3,2) + k * bg(1:3,3)
    phasefac(:,Q) = (0.d0,0.d0)
    do nk=0,dfft_%nr3-1
    do nj=0,dfft_%nr2-1
    do ni=0,dfft_%nr1-1
      ir = 1 + ni + dfft_%nr1x * ( nj + dfft_%nr2x*nk )
      r(1:3) = ni*at(1:3,1)/REAL(dfft_%nr1,DP) + &
                nj*at(1:3,2)/REAL(dfft_%nr2,DP) + &
                nk*at(1:3,3)/REAL(dfft_%nr3,DP)
      phasefac(ir,Q) = EXP((0.d0,1.d0)*sum(tpi*dG(1:3,Q)*r(1:3)))
    enddo
    enddo
    enddo
    Q = Q+1
  enddo
  enddo
  enddo

END SUBROUTINE calculate_phase_factor

  ! setup Q-point maps, QKtoK2 and kminus
  ! Q(ka,kb) = (ka-kb)+G, where nk pais of {ka,kb} map to a given Q 
  ! Inverse mapping (QK lists the specific pair that maps to Q):
  !  {Q, K} -> { ka = k,  kb = K - Q + G'  }
  !  for a given {Q, K}, the corresponding pairs is then given by: {K,
  !  QKtoK2(Q,K)} generate Qpts grid
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_Qpt_map(nk,xk,Qpts,QKtoK2,kminus,nq1,nq2,nq3,use_regularization)
  !-----------------------------------------------------------------------
    !
    USE KINDS, ONLY : DP
    USE cell_base, ONLY: bg, at
    !
    IMPLICIT NONE
    !
    INTEGER,   INTENT(IN)   :: nk
    REAL(DP),  INTENT(IN)   :: xk(3,nk)
    LOGICAL,   INTENT(IN)   :: use_regularization
    REAL(DP),  INTENT(OUT)  :: Qpts(3,nk)
    INTEGER,   INTENT(OUT)  :: QKtoK2(nk,nk)
    INTEGER,   INTENT(OUT)  :: kminus(nk)
    INTEGER, INTENT(OUT) :: nq1,nq2,nq3
    REAL(DP) :: xk_cryst(3,nk)
    INTEGER :: ik,ka,kb,i,j,k
    REAL(DP) :: dk(3)
    !

    do ik=1,nk
      Qpts(1:3,ik) = xk(1:3,ik) - xk(1:3,1)
    enddo
    ! Can be parallelized if it is slow!!!
    QKtoK2(:,:)=0
    do ik=1,nk
      do ka=1,nk
        do kb=1,nk
          do i=-1,1
          do j=-1,1
          do k=-1,1
            dk(1:3) =  xk(1:3,ka) - xk(1:3,kb) - Qpts(1:3,ik) - &
                       i * bg(1:3,1) - j * bg(1:3,2) - k * bg(1:3,3)
            if( sum(dk(1:3)**2) .lt. 1.0E-8_DP ) then
              if( QKtoK2(ik,ka) .gt. 0 ) &
                call errore('pw2qmcpack',' Problems setting up QKtoK2 map.',10)
              QKtoK2(ik,ka) = kb
              GOTO 301
            endif
          enddo
          enddo
          enddo
        enddo
301     if( QKtoK2(ik,ka) .lt. 1 ) &
          call errore('pw2qmcpack',' Could not find QKtoK2 map.',11)
      enddo
    enddo
    kminus(:)=0
    do ik=1,nk
      do ka=1,nk
        do i=-1,1
        do j=-1,1
        do k=-1,1
          dk(1:3) =  Qpts(1:3,ik) + Qpts(1:3,ka) -  &
                     i * bg(1:3,1) - j * bg(1:3,2) - k * bg(1:3,3)
          if( sum(dk(1:3)**2) .lt. 1.0E-8_DP ) then
            if( kminus(ik) .gt. 0 ) &
              call errore('pw2qmcpack',' Problems setting up kminus map.',12)
            kminus(ik) = ka
            GOTO 302
          endif
        enddo
        enddo
        enddo
      enddo
302   if( kminus(ik) .lt. 1 ) &
        call errore('pw2qmcpack',' Could not find kminus map.',13)
    enddo
    do ik=1,nk
      if( kminus( kminus(ik) ) .ne. ik ) &
        call errore('pw2qmcpack',' Error in Q(-Q) map.',14)
    enddo
    if (use_regularization) then! find dimensions of Q-grid
    ! MAM: is this generic?
    if( (nq1==0) .or. (nq2==0) .or. (nq3==0) ) then
      nq1=0
      nq2=0
      nq3=0
      xk_cryst(:,:) = Qpts(:,:)
      CALL cryst_to_cart( nk, xk_cryst, at, -1 )
      do ik=1,nk
        if( (abs(xk_cryst(2,ik)) < 1e-6) .and. (abs(xk_cryst(3,ik)) < 1e-6) ) nq1=nq1+1
        if( (abs(xk_cryst(1,ik)) < 1e-6) .and. (abs(xk_cryst(2,ik)) < 1e-6) ) nq3=nq3+1
        if( (abs(xk_cryst(3,ik)) < 1e-6) .and. (abs(xk_cryst(1,ik)) < 1e-6) ) nq2=nq2+1
      enddo
    endif
    write(*,*) ' Dimensions of the Q-point grid: ',nq1,nq2,nq3
    if( nq1*nq2*nq3 .ne. nk ) & 
      call errore('calculate_Qpt_map','Problems with Qpoint grid. nq1*nq2*nq3 ne nk',1)
    endif ! dimension of Q-grid
  END SUBROUTINE calculate_Qpt_map
 !
  !-----------------------------------------------------------------------
  SUBROUTINE orbitals_g2r(npw, igk, dfft_, nnr, m, PsiG, PsiR)
  !-----------------------------------------------------------------------
  !
    USE KINDS, ONLY : DP
    USE wavefunctions, ONLY : evc
    USE control_flags, ONLY : gamma_only
    USE fft_types, ONLY: fft_type_descriptor
    use fft_interfaces,       ONLY : invfft
    USE wvfct, ONLY: npwx
    !
    IMPLICIT NONE
    !
    INTEGER,      INTENT(IN)   :: m, npw, nnr
    INTEGER,      INTENT(IN)   :: igk(npwx) 
    COMPLEX(DP),  INTENT(IN)   :: PsiG(npwx,m)
    COMPLEX(DP),  INTENT(OUT)  :: PsiR(nnr,m)
    TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft_
    !
    INTEGER :: ibnd
    !
    PsiR(:,:)=(0.d0,0.d0)
    do ibnd=1,m
      !
      PsiR(dfft_%nl(igk(1:npw)),ibnd)=PsiG(1:npw,ibnd)
      if(gamma_only) &
          PsiR(dfft_%nlm(igk(1:npw)),ibnd) = CONJG(PsiG(1:npw,ibnd))
      !
      CALL invfft ('Wave', PsiR(:,ibnd), dfft_)
      !
    enddo
    !
  END SUBROUTINE orbitals_g2r
  !
  ! pos in range [1,npart]
  subroutine fair_divide(e1,e2,pos,npart,nelms)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: pos,npart,nelms
    INTEGER, INTENT(OUT) :: e1,e2
    INTEGER :: d,r

    d = nelms / npart
    r = mod(nelms, npart)
    if( pos <= r ) then
      e1 = (pos-1) * (d + 1) + 1
      e2 = pos * (d + 1)
    else
      e1 = (pos-1) * d + r + 1
      e2 = pos * d + r
    endif

  end subroutine fair_divide
  !
  SUBROUTINE eigsys(jobz, uplo, reverse, n, lda, A, W)
    !
    USE KINDS, ONLY : DP
    !
    IMPLICIT NONE
    !
    CHARACTER(len=1),INTENT(IN) :: uplo, jobz
    INTEGER, INTENT(IN) :: n, lda
    COMPLEX(DP), INTENT(INOUT) :: A(lda,n)
    REAL(DP), INTENT(OUT) :: W(n)
    LOGICAL, INTENT(IN) :: reverse
    !
    INTEGER :: i,j,IL,IU,M,lwork,lrwork,liwork,info
    REAL(DP) :: VL,VU,ABSTOL
    INTEGER, ALLOCATABLE :: isuppz(:), iwork(:)
    REAL(DP), ALLOCATABLE :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: Z(:,:), work(:)

    allocate( Z(n,n), isuppz(2*n), work(1), rwork(1), iwork(1) )

    lwork = -1
    lrwork = -1
    liwork = -1
    call zheevr(jobz,'A',uplo,n,A,lda,VL,VU,IL,IU,ABSTOL,M,W,Z,n,isuppz,&
                work,lwork,rwork,lrwork,iwork,liwork,info)

    lwork = INT(dble(work(1)))
    lrwork = INT(rwork(1))
    liwork = INT(iwork(1))
    deallocate(work, rwork, iwork)
    allocate( work(lwork), rwork(lrwork), iwork(liwork) )

    call zheevr(jobz,'A',uplo,n,A,lda,VL,VU,IL,IU,ABSTOL,M,W,Z,n,isuppz,&
                work,lwork,rwork,lrwork,iwork,liwork,info)

    if(info.ne.0) call errore('pw2qmcpack','info != 0',info)
    if(M.ne.n) call errore('pw2qmcpack','Not enough eigenvalues in eigsys.',1)

    ! copy eigenvalues to A
    if(reverse) then
      do j=1,n/2
        VL = W(j)
        W(j) = W(n-j+1)
        W(n-j+1) = VL
      enddo
      if(jobz=='V') then
        do j=1,n
          do i=1,n
            A(i,n-j+1) = Z(i,j)
          enddo
        enddo
      endif
    else
      if(jobz=='V') then
        do j=1,n
          do i=1,n
            A(i,j) = Z(i,j)
          enddo
        enddo
      endif
    endif

    deallocate(Z,isuppz,work, rwork, iwork)

  END SUBROUTINE eigsys

  ! lazy for now
  SUBROUTINE eigsysD(jobz, uplo, reverse, n, lda, A, W)
    !
    USE KINDS, ONLY : DP
    !
    IMPLICIT NONE
    !
    CHARACTER(len=1),INTENT(IN) :: uplo, jobz
    INTEGER, INTENT(IN) :: n, lda
    REAL(DP), INTENT(INOUT) :: A(lda,n)
    REAL(DP), INTENT(OUT) :: W(n)
    LOGICAL, INTENT(IN) :: reverse
    !
    INTEGER :: i,j,IL,IU,M,lwork,liwork,info
    REAL(DP) :: VL,VU,ABSTOL
    INTEGER, ALLOCATABLE :: isuppz(:), iwork(:)
    REAL(DP), ALLOCATABLE :: Z(:,:), work(:)

    allocate( Z(n,n), isuppz(2*n), work(1), iwork(1) )

    lwork = -1
    liwork = -1
    call dsyevr(jobz,'A',uplo,n,A,lda,VL,VU,IL,IU,ABSTOL,M,W,Z,n,isuppz,&
                work,lwork,iwork,liwork,info)

    lwork = INT(dble(work(1)))
    liwork = INT(iwork(1))
    deallocate(work, iwork)
    allocate( work(lwork), iwork(liwork) )

    call dsyevr(jobz,'A',uplo,n,A,lda,VL,VU,IL,IU,ABSTOL,M,W,Z,n,isuppz,&
                work,lwork,iwork,liwork,info)

    if(info.ne.0) call errore('pw2qmcpack','info != 0',info)
    if(M.ne.n) call errore('pw2qmcpack','Not enough eigenvalues in eigsys.',1)

    ! copy eigenvalues to A
    if(reverse) then
      do j=1,n/2
        VL = W(j)
        W(j) = W(n-j+1)
        W(n-j+1) = VL
      enddo
      if(jobz=='V') then
        do j=1,n
          do i=1,n
            A(i,n-j+1) = Z(i,j)
          enddo
        enddo
      endif
    else
      if(jobz=='V') then
        do j=1,n
          do i=1,n
            A(i,j) = Z(i,j)
          enddo
        enddo
      endif
    endif

    deallocate(Z,isuppz,work, iwork)

  END SUBROUTINE eigsysD
  !
  !-----------------------------------------------------------------------
  SUBROUTINE zmatinv (n, a, lda)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, lda
  COMPLEX (DP), INTENT(INOUT)  :: a(lda,n)
  !
  INTEGER :: info, lwork
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  INTEGER, ALLOCATABLE :: ipiv (:)
  ! ipiv  : work space for pivoting
  COMPLEX(DP), ALLOCATABLE :: work (:)
  !
  lwork=64*n
  ALLOCATE(ipiv(n), work(lwork) )
  !
  CALL zgetrf (n, n, a, lda, ipiv, info)
  CALL errore ('matinv', 'error in ZGETRF', abs (info) )
  CALL zgetri (n, a, lda, ipiv, work, lwork, info)
  CALL errore ('matinv', 'error in ZGETRI', abs (info) )
  !
  DEALLOCATE ( work, ipiv )
  !
  END SUBROUTINE zmatinv
  !-----------------------------------------------------------------------
  !
  SUBROUTINE find_2d_partition( N, rk, np, i0, i1, j0, j1)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, rk, np
    INTEGER, INTENT(OUT) :: i0,i1,j0,j1
    INTEGER :: facts(2), i, imax, k, r, c
    !
    i0 = 1
    j0 = 1
    i1 = N
    j1 = N
    if(np == 1) return
    facts(1) = 1
    facts(2) = np 
    imax = nint(sqrt(np*1.d0))+1
    do i=2,imax
      if( MOD ( np, i ) == 0 ) then
        k = np/i
        if( abs(k-i) < abs(facts(1)-facts(2)) ) then
          facts(1)=i
          facts(2)=k
        endif
      endif
    enddo
    if(facts(1)*facts(2) .ne. np ) &
        call errore('pw2qmcpack','Error in find_2d_partition.',1)
    i = (rk-1)/facts(2) + 1
    k = MOD( (rk-1), facts(2) ) + 1
    call fair_divide(i0,i1,i,facts(1),N)
    call fair_divide(j0,j1,k,facts(2),N)
    !
  END SUBROUTINE find_2d_partition
  !
  SUBROUTINE transpose_mat(n, A)
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    COMPLEX (DP), INTENT(INOUT)  :: a(n,n)
    !
    INTEGER :: ia, ib
    COMPLEX (DP) :: ctemp
    !
    do ia=1,n
      do ib=ia+1,n
        ctemp = A(ia,ib)
        A(ia,ib) = A(ib,ia)
        A(ib,ia) = ctemp
      enddo
    enddo
    ! 
  END SUBROUTINE transpose_mat

  ! make sure init_us_2 and allocate_bec_type have been 
  ! called for the current kpoint
  subroutine Overlap(m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,addAug,spsi)
    !
    USE kinds, ONLY : DP
    USE control_flags, ONLY : gamma_only
    USE gvect, ONLY: gstart
    USE paw_variables,        ONLY : okpaw
    USE uspp,       ONLY : okvan
    USE becmod,   ONLY : becp, calbec
    USE uspp,     ONLY : vkb, nkb
    USE realus,   ONLY : real_space
    !
    IMPLICIT NONE
    !
    INTEGER, intent(IN) :: m,n,k,lda,ldb,ldc
    COMPLEX(DP), INTENT(IN) :: alpha,beta
    LOGICAL, INTENT(IN) :: addAug
    COMPLEX(DP) :: A(lda,*),B(ldb,*)
    COMPLEX(DP) ::  C(ldc,*)
    COMPLEX(DP) ::  spsi(ldb,*)
    COMPLEX(DP) :: gammafac
    INTEGER :: i,j
    !
    gammafac=(1.d0,0.d0)
    if(gamma_only) gammafac=(2.d0,0.d0)

    if(addAug .and. (okpaw .or. okvan)) then
      if(gamma_only) call errore('pw2qmcpack: gamma_only and addAug!',1)
      call calbec(k, vkb, B(:,1:n), becp, n)
      call s_psi(ldb, k, n, B, spsi)
      call zgemm('C','N',m,n,k,gammafac*alpha,A,lda,spsi,ldb,beta,C,ldc)
      if(gamma_only) then
        if(gstart==2) call zgerc(m,n,(-1.d0,0.d0)*alpha,A,lda,spsi,ldb,C,ldc)
        do j=1,n
          do i=1,m
            C(i,j) = cmplx( real(C(i,j)), 0.0_dp, kind=dp )
          enddo
        enddo
      endif
    else
      call zgemm('C','N',m,n,k,gammafac*alpha,A,lda,B,ldb,beta,C,ldc)
      if(gamma_only) then
        if(gstart==2) call zgerc(m,n,(-1.d0,0.d0)*alpha,A,lda,B,ldb,C,ldc)
        do j=1,n
          do i=1,m
            C(i,j) = cmplx( real(C(i,j)), 0.0_dp, kind=dp )
          enddo
        enddo
      endif
    endif
  end subroutine Overlap

  ! raken from newd()
  ! this routine resets deeq to their bare values, removing onecenter PAW
  ! Hartree, Exc terms and Veff*Qnm terms coming from the 
  ! soft augmentation part these terms are added in the 2-electron hamiltonian
  subroutine reset_deeq()
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : deeq, dvan, deeq_nc, dvan_so, okvan
  USE paw_variables,        ONLY : okpaw
  USE uspp_param,           ONLY : upf, nh
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE scf,                  ONLY : v
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : newq_r
  USE dfunct,        ONLY : newq
  !
  IMPLICIT NONE
  !
  INTEGER :: nt, ih, jh, na, is
  REAL(kind=dp), ALLOCATABLE  :: v_(:,:)
  !
  if(.not.(okvan .or. okpaw)) return
  !
  deeq(:,:,:,:) = 0.D0
  !
  ! add only vltot part!
  allocate(v_(dfftp%nnr,nspin))
  v_(:,:)=0.d0
  IF (tqr) THEN
     CALL newq_r(v_,deeq,.false.)
  ELSE
     CALL newq(v_,deeq,.false.)
  END IF
  deallocate(v_)
  !
  atoms : &
  DO na = 1, nat
    !
    nt  = ityp(na)
    if_noncolin:&
    IF ( noncolin ) THEN
      call errore('noncollinear not yet implemented with paw/uspp',1)
      !
      IF (upf(nt)%has_so) THEN
         !
!         CALL newd_so(na)
         !
      ELSE
         !
!         CALL newd_nc(na)
         !
      END IF
      !
    ELSE if_noncolin
      !
      DO is = 1, nspin
         !
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
               deeq(ih,jh,na,is) = deeq(ih,jh,na,is) + dvan(ih,jh,nt)
               deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
            END DO
         END DO
         !
      END DO
      !
    END IF if_noncolin
    !
  END DO atoms
  end subroutine reset_deeq
  !
  SUBROUTINE calculate_factorized_paw_one_center( ke_fact, nke )
  !=----
  ! factorize paw's ke%k for each species  
  !=----
    USE kinds, ONLY : DP
    USE constants, ONLY: e2
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp_param,         ONLY : nh, upf
    USE paw_variables,        ONLY : okpaw
    USE mp_images, ONLY: intra_image_comm, me_image, root_image
    USE mp, ONLY: mp_bcast
    USE paw_exx,            ONLY : ke, PAW_init_fock_kernel, &
                                       PAW_clean_fock_kernel
    USE posthf_mod, ONLY: ke_factorization
    !
    IMPLICIT NONE
    !
    type(ke_factorization), INTENT(INOUT) :: ke_fact(ntyp)
    INTEGER, INTENT(OUT) :: nke
    !
    INTEGER :: ns, ij, ou, na,n, lda
    INTEGER :: i,ih, jh, oh, uh
    REAL(DP), ALLOCATABLE :: k(:,:)
    REAL(DP), ALLOCATABLE :: ev(:)
    !
    ! on the first call, symmetrize ke tensor
    nke=0
    if(.not.okpaw) return
    if(me_image==root_image) then
      CALL PAW_init_fock_kernel()
      do ns = 1,ntyp
        ij=0
        ! not using symmetry to keep it simple
        allocate(k(nh(ns)*nh(ns),nh(ns)*nh(ns)), ev(nh(ns)*nh(ns)))
        do jh = 1, nh(ns)
        do ih = 1, nh(ns)
        ij=ij+1
        ou=0
        do uh = 1, nh(ns)
        do oh = 1, nh(ns)
          !
          ou=ou+1
          ! factor of 1/e2 is to convert to Hartree, factor of 0.5
          ! comes from the QE implementation. Not sure I understand why!
          ! Keeping it for consistency with QE
          k(ij,ou) = (0.5d0/e2) * 0.125_dp * ( &
              ke(ns)%k(ih,jh,oh,uh) + ke(ns)%k(oh,uh,ih,jh) + &
              ke(ns)%k(jh,ih,oh,uh) + ke(ns)%k(uh,oh,ih,jh) + &
              ke(ns)%k(ih,jh,uh,oh) + ke(ns)%k(oh,uh,jh,ih) + &
              ke(ns)%k(jh,ih,uh,oh) + ke(ns)%k(uh,oh,jh,ih))

!          k(ij,ou) = (0.5d0/e2) * 0.25_dp * ( &
!              ke(ns)%k(ih,jh,oh,uh) + ke(ns)%k(oh,uh,ih,jh) + &
!              ke(ns)%k(jh,ih,uh,oh) + ke(ns)%k(uh,oh,jh,ih))

!          k(ij,ou) = (0.5d0/e2) * ke(ns)%k(ih,jh,oh,uh)

          !
        enddo
        enddo
        enddo
        enddo

        lda = nh(ns)*nh(ns)
        call eigsysD('V', 'U', .true.,lda,lda,k,ev)
        n=0
        do i=1, nh(ns)* nh(ns)
          if( abs(ev(i)) > 1.e-14 ) n = n + 1
        enddo
        allocate(ke_fact(ns)%L(nh(ns)*nh(ns),n))
        n=0
        do i=1,nh(ns)* nh(ns)
          if( abs(ev(i)) > 1.e-14 ) then
            n = n + 1
            ke_fact(ns)%L(:,n) = CMPLX(k(:,i),0._dp,kind=DP) * &
                                 sqrt(CMPLX(ev(i),0._dp,kind=DP))
          endif
        enddo
        call mp_bcast(n, root_image, intra_image_comm)
        call mp_bcast(ke_fact(ns)%L(:,:), root_image, intra_image_comm)
        ! if truncating, modify nke
        if ( upf(ns)%tpawp ) then
          do na = 1, nat
            if (ityp(na)==ns) then
              nke = nke + size( ke_fact(ns)%L, 2)
            endif
          enddo
        endif
        deallocate(k,ev)
      enddo
      call PAW_clean_fock_kernel()
      ! bcast nke
      call mp_bcast(nke, root_image, intra_image_comm)
    else
      do ns = 1,ntyp
        call mp_bcast(n, root_image, intra_image_comm)
        allocate(ke_fact(ns)%L(nh(ns)*nh(ns),n))
        call mp_bcast(ke_fact(ns)%L(:,:), root_image, intra_image_comm)
      enddo
      call mp_bcast(nke, root_image, intra_image_comm)
    endif

  END SUBROUTINE calculate_factorized_paw_one_center

  ! From PAW_xx_energy
  !=----------------------------------------------------------------------------=!
  SUBROUTINE contract_paw_one_center(ke,A,beca,becb)
  !=----------------------------------------------------------------------------=!
  !
    USE kinds, ONLY : DP
    USE constants, ONLY: e2
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp_param,         ONLY : nh, upf
    USE uspp,               ONLY : nkb, ofsbeta 
    USE paw_variables,        ONLY : okpaw
    USE posthf_mod, ONLY: ke_factorization
    !
    IMPLICIT NONE
    !
    type(ke_factorization), INTENT(IN) :: ke(ntyp)
    COMPLEX(DP), INTENT(OUT) :: A(nkb)
    COMPLEX(DP),INTENT(in) :: beca(nkb), becb(nkb)
    !
    INTEGER :: np, na, ns, nke, n
    INTEGER :: ih, jh, oh, uh, ij
    INTEGER :: ikb, jkb, okb, ukb, ijkb0

    A(:) = (0.d0,0.d0) 
    if(.not.okpaw) return
    nke = 1
    DO np = 1, ntyp
      IF ( upf(np)%tpawp ) THEN
        DO na = 1, nat
          IF (ityp(na)==np) THEN
            ijkb0 = ofsbeta(na)
            !
            DO n = 1,size(ke(np)%L, 2)
              ij=1
              DO jh = 1, nh(np)
                jkb = ijkb0 + jh
                DO ih = 1, nh(np)
                  ikb = ijkb0 + ih
                  A(nke) = A(nke) + ke(np)%L(ij,n)* &
                              CONJG(beca(ikb)) * becb(jkb)
                  ij = ij + 1
                ENDDO
                !
              ENDDO
              nke = nke + 1
              !
            ENDDO
            !
          END IF
        ENDDO ! nat
      END IF
    ENDDO

  END SUBROUTINE contract_paw_one_center

  ! From PAW_xx_energy
  !=----------------------------------------------------------------------------=!
  FUNCTION PAW_J_energy(beca, becb, becc, becd)
    !=----------------------------------------------------------------------------=!
    ! Compute onsite contribution to 2-electron integral: (ab|cd) 
    !
    USE kinds, ONLY : DP
    USE constants, ONLY: e2
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp
    USE uspp_param,         ONLY : nh, upf
    USE uspp,               ONLY : nkb, ofsbeta 
    USE io_global,          ONLY : ionode
    USE paw_exx,        ONLY: ke
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_fock_kernel
    USE posthf_mod, ONLY: ke_factorization
    !
    IMPLICIT NONE
    !
    COMPLEX(DP),INTENT(in) :: beca(nkb), becb(nkb), becc(nkb), becd(nkb)
    !
    COMPLEX(DP) :: PAW_J_energy
    !
    INTEGER :: np, na, ns
    INTEGER :: ih, jh, oh, uh
    INTEGER :: ikb, jkb, okb, ukb, ijkb0
    REAL(DP), ALLOCATABLE :: k(:,:,:,:)
    LOGICAL,SAVE :: first = .true.
    !
    PAW_J_energy = (0._dp,0._dp)
    !
    if(.not.okpaw) return
    !
    ! on the first call, symmetrize ke tensor
    if(first) then
      first=.false.
      CALL PAW_init_fock_kernel()
      do ns = 1,ntyp
        allocate(k(nh(ns),nh(ns),nh(ns),nh(ns)))
        k(:,:,:,:) = ke(ns)%k(:,:,:,:)
        do uh = 1, nh(ns)
        do oh = 1, nh(ns)
        do jh = 1, nh(ns)
        do ih = 1, nh(ns)
          !
          ke(ns)%k(ih,jh,oh,uh) = 0.125_dp * ( &
              k(ih,jh,oh,uh) + k(oh,uh,ih,jh) + &
              k(jh,ih,oh,uh) + k(uh,oh,ih,jh) + &
              k(ih,jh,uh,oh) + k(oh,uh,jh,ih) + &
              k(jh,ih,uh,oh) + k(uh,oh,jh,ih))
          !
        enddo
        enddo
        enddo
        enddo
        deallocate(k)
      enddo
    endif
    ! 
    DO np = 1, ntyp
      ONLY_FOR_PAW : &
      IF ( upf(np)%tpawp ) THEN
        DO na = 1, nat
        IF (ityp(na)==np) THEN
          ijkb0 = ofsbeta(na)
          !
          DO uh = 1, nh(np)
            ukb = ijkb0 + uh
            DO oh = 1, nh(np)
              okb = ijkb0 + oh
              DO jh = 1, nh(np)
                jkb = ijkb0 + jh
                DO ih = 1, nh(np)
                  ikb = ijkb0 + ih
                  ! Eq. 32 and 42 Ref. 1 :
                  ! k(i,j,o,u) in Chemistry notation (ij|ou)  
                  ! you can calculate the cholesky decomposition of k(i,j,o,u)
                  ! to speed up this part, reducing scaling to nh(np)^2 instead
                  ! of nh(np)^4   
                  PAW_J_energy = PAW_J_energy + ke(np)%k(ih,jh,oh,uh) &
                                * CONJG(beca(ikb)) * becb(jkb) &
                                * CONJG(becc(okb)) * becd(ukb)
                  !
                ENDDO !ih, ukb
              ENDDO !jh, okb
            ENDDO !oh, jkb
          ENDDO !uh, ikb
          !
        END IF
        ENDDO ! nat
      END IF &
      ONLY_FOR_PAW
    ENDDO
    ! MAM: I don't understand the need for the factor of 0.5 here,
    ! but I can't reproduce QE energies without it!!!
    ! There is a comment about a factor of 0.5 coming from Kresse's paper, eq 32
    ! but this is the factor of 1/2 in the hamiltonian which I think is already
    ! applied in PW/src/electron.f90.  Keeping it here until I understand
    ! better.
    PAW_J_energy = 0.5d0 * PAW_J_energy / e2
    !
    RETURN
    !=----------------------------------------------------------------------------=!
  END FUNCTION PAW_J_energy
  !=----------------------------------------------------------------------------=!

  !--------------------------------------------------------------------------
  FUNCTION count_free_unit()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: count_free_unit
    INTEGER :: iunit
    LOGICAL :: opnd
    !
    !
    count_free_unit = 0 
    unit_loop: DO iunit = 999, 1, -1
       !
       INQUIRE( UNIT = iunit, OPENED = opnd )
       !
       IF ( opnd ) count_free_unit = count_free_unit + 1 
       !
    END DO unit_loop
    !
    RETURN
    !
  END FUNCTION count_free_unit

!-----------------------------------------------------------------------
SUBROUTINE find_Q_symm(nk, xk, Qpts, QKtoK2, nQ, xQ, wQ)
!-----------------------------------------------------------------------
! finds the set of irreducble Q points
! only symmetries that leave the kpoint grid invariant are considered
! the routine assumes a full (unfolded) set of kpoints 
! Note: For now, will not return relations between Q points and symmetry ops, 
!       can add later if needed
  !
  USE KINDS, ONLY : DP
  USE io_global,          ONLY : ionode
  USE cell_base, ONLY: bg, at, tpiba2
  USE symm_base,         ONLY : nsym, s
  USE constants, ONLY: tpi, fpi, e2
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nk
  REAL(DP), INTENT(IN)  :: xk(3,nk)      ! list of k-points 
  REAL(DP), INTENT(IN)  :: Qpts(3,nk)     ! list of Q-points
  INTEGER,  INTENT(IN)  :: QKtoK2(nk,nk)
  INTEGER,  INTENT(OUT) :: nQ 
  REAL(DP), INTENT(OUT) :: wQ(nk) 
  INTEGER,  INTENT(OUT) :: xQ(nk)
  !
  !
  INTEGER :: ik,ikk,isym, ngoodsym
  REAL(DP) :: xk_cryst(3,nk), sk(3), dk(3)
  REAL(DP) :: eps = 1.d-4 
  REAL(DP) :: one  
  LOGICAL :: good_syms(nsym), unique(nk), found 

  xk_cryst(:,:) = xk(:,:)
  CALL cryst_to_cart( nk, xk_cryst, at, -1 )
  ngoodsym=0
  good_syms(:) = .true. ! assume that all are "good"
  ! find subset of symmetries that leave the kpoint grid invariant  
  isymloop: do isym=1,nsym
    ikloop: do ik=1,nk
      ! check that s(i)*xk(k) = xk(k'), for all k. 
      sk(:) = s(:,1,isym)*xk_cryst(1,ik) + &
              s(:,2,isym)*xk_cryst(2,ik) + &
              s(:,3,isym)*xk_cryst(3,ik) 
      found = .false.
      ikkloop: do ikk = 1, nk
        dk(:) = sk(:)-xk_cryst(:,ikk) - NINT(sk(:)-xk_cryst(:,ikk))
        if ( abs(dk(1))<=eps .and. &
             abs(dk(2))<=eps .and. &
             abs(dk(3))<=eps ) then
          found = .true.
          exit ikkloop  
        endif 
      end do ikkloop
      if( .not. found ) then
        good_syms(isym) = .false.
        exit ikloop
      endif
    end do ikloop
    if( good_syms(isym) ) ngoodsym = ngoodsym + 1
  end do isymloop

  if(ionode) write(*,*) ' nsym, nQsym: ',nsym,ngoodsym

  ! now look for irreducible set of Q within subgroup of "good" symmetries
  ! reuse xk_cryst for Q vectors in crystal axis
  xk_cryst(:,:) = Qpts(:,:)
  CALL cryst_to_cart( nk, xk_cryst, at, -1 )
  nQ = 0  
  wQ(:) = 0.d0
  xQ(:) = 0 
  unique(:) = .true.
  do ik=1,nk
    if(unique(ik)) then
      nQ=nQ+1
      xQ(nQ) = ik
      wQ(nQ) = wQ(nQ) + 1.d0
      isymloop2: do isym=1,nsym 
        if(.not.good_syms(isym)) cycle isymloop2
        sk(:) = s(:,1,isym)*xk_cryst(1,ik) + &
                s(:,2,isym)*xk_cryst(2,ik) + &
                s(:,3,isym)*xk_cryst(3,ik)
        ikkloop2: do ikk=ik+1,nk
          if(.not.unique(ikk)) cycle ikkloop2
          dk(:) = sk(:)-xk_cryst(:,ikk) - NINT(sk(:)-xk_cryst(:,ikk))      
          if ( abs(dk(1))<=eps .and. &
             abs(dk(2))<=eps .and. &
             abs(dk(3))<=eps ) then
            ! found match, tag as not unique and add weight of nQ
            wQ(nQ) = wQ(nQ) + 1.d0
            unique(ikk) = .false.
            exit ikkloop2
          endif        
        end do ikkloop2    
      end do isymloop2   
    endif
  enddo    

  if(abs(sum(wq(:))-dble(nk)) > 1e-8 ) then
    write(*,*) 'Problems in find_Q_symm, sum(wq):',sum(wq),dble(nk)
    call errore('find_Q_symm',' Problems in find_Q_symm.',10)
  endif  
  wq(:) = wq(:) / dble(nk) 

!-----------------------------------------------------------------------
END SUBROUTINE find_Q_symm
!----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
! Finds number of occupied orbitals from vector of Fermi weights  
!-----------------------------------------------------------------------
SUBROUTINE get_nelmax(nbnd,nk,ns,wk,n1,wg,nelmax)
  !
  USE KINDS, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nbnd,nk,ns,n1
  REAL(DP), INTENT(IN) :: wg(n1,*), wk(*) 
  INTEGER, INTENT(OUT) :: nelmax    
  !
  INTEGER :: ispin,ik,ikk,ia
  REAL(DP) :: scl  
  !  
  nelmax = 0
  do ispin=1,ns
    do ik=1,nk
      ikk = ik + nk*(ispin-1)
      if(abs(wk(ikk))>1.d-10) then
          scl = 1.d0/wk(ikk)
      else
          scl = 1.d0
      endif
      do ia=1,nbnd
        ! MAM: for high-T or for highly degenerate states this needs to be reduced
        if( abs(wg(ia,ikk)*scl) > 0.01d0 .and. ia > nelmax )  &
          nelmax = ia
      enddo
    enddo
  enddo
END SUBROUTINE get_nelmax
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Generates an N-point Gauss-Legendre grid (with scaling function f(x)=1)
! Not efficient, but simple!
!-----------------------------------------------------------------------
SUBROUTINE gaussleg_quad(N,ix,iw)
  !
  USE KINDS, ONLY : DP
  USE constants, ONLY: pi
  !
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(OUT) :: ix(N), iw(N)
  !
  INTEGER :: i, cnt, N2
  REAL(DP) :: x,x0,pn(2),dpnx
  
  N2 = N/2
  if( MOD(N,2) == 1 ) N2=N2+1

  do i = 1, N2
    x0 = cos(pi*(i-0.25d0)/(N+0.5d0))
    pn = LegPoly(N,x0)
    dpnx = (x0*pn(1) - pn(2)) * (1.d0*N)/(x0*x0-1.d0)
    x = x0 - pn(1)/dpnx
    cnt=1
    do while( abs(x-x0) > 1e-12 )
      x0 = x 
      pn = LegPoly(N,x0)
      dpnx = (x0*pn(1) - pn(2)) * (1.d0*N)/(x0*x0-1.d0)
      x = x0 - pn(1)/dpnx
      cnt=cnt+1
    end do
    ix(i) = x 
    iw(i) = 2.d0/((1.d0-x*x)*dpnx*dpnx)
    ! using symmetry
    if( i <= N/2 ) then
      ix(N-i+1) = -1.d0*ix(i)
      iw(N-i+1) = iw(i)
    endif
  end do
!  do i = 1, N
!    write(*,*) 'xi,wi: ',ix(i),iw(i)
!  end do
 
  CONTAINS 

  ! could evaluate both Pn and Pn-1 simultaneously, but keeping it simple
  FUNCTION LegPoly(n,x) result(y)
  implicit none
  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: y(2)
  INTEGER :: i
  REAL(DP) :: y1,y2  ! y(n), y(n-1)
  !
  if(n==0) then
    y(1)=1.d0
    y(2)=0.d0
  elseif(n==1) then
    y(1)=x
    y(2)=1.d0
  else
    y1=x     
    y2=1.d0
    do i=2,n
      y(1)= ( (2.d0*i-1.d0)*x*y1 - (1.d0*i-1.d0)*y2) / (1.d0*i)
      y(2)= y1
      y1=y(1)
      y2=y(2)  
    enddo
  endif 
  !
  END FUNCTION LegPoly
  !
END SUBROUTINE gaussleg_quad

SUBROUTINE gaussleg_quad2(N,ix,iw)
  !
  USE KINDS, ONLY : DP
  USE constants, ONLY: pi
  !
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(OUT) :: ix(N), iw(N)
  !
  INTEGER :: i, cnt, N2
  REAL(DP) :: x,x0,pnx,pnx1,dpnx
  
  N2 = N/2
  if( MOD(N,2) == 1 ) N2=N2+1

  do i = 1, N2
    x0 = cos(pi*(i-0.25d0)/(n+0.5d0))
    pnx = Pn(N,x0)
    pnx1 = Pn(N-1,x0)
    dpnx = (x0*pnx - pnx1) * (1.d0*n)/(x0*x0-1.d0)
    x = x0 - pnx/dpnx
    cnt=1
    do while( abs(x-x0) > 1e-12 )
      x0 = x 
      pnx = Pn(N,x0)
      pnx1 = Pn(N-1,x0)
      dpnx = (x0*pnx - pnx1) * (1.d0*n)/(x0*x0-1.d0)
      x = x0 - pnx/dpnx
      cnt=cnt+1
    end do
    ix(i) = x 
    iw(i) = 2.d0/((1.d0-x*x)*dpnx*dpnx)
    ! using symmetry
    if( i <= N/2 ) then
      ix(N-i+1) = -1.d0*ix(i)
      iw(N-i+1) = iw(i)
    endif
  end do
  do i = 1, N
    write(*,*) 'xi,wi: ',ix(i),iw(i)
  end do
 
  CONTAINS 

  ! could evaluate both Pn and Pn-1 simultaneously, but keeping it simple
  RECURSIVE FUNCTION Pn(n,x) result(y)
  implicit none
  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: y
  !
  if(n==0) then
    y=1.d0
  elseif(n==1) then
    y=x
  else
    y = ( (2.d0*n-1.d0)*x*Pn(n-1,x) - (1.d0*n-1.d0)*Pn(n-2,x)) / (1.d0*n)
  endif 
  !
  END FUNCTION Pn
  !
END SUBROUTINE gaussleg_quad2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transforms an N-point Gauss-Legendre integration grid from [-1:1] to [0:inf]
!  Follows: X. Ren, et al.,  New Journal of Physics, Volume 14, 053020 (2012).
!-----------------------------------------------------------------------
  SUBROUTINE transform_gaussleg(N,x0,ix,iw)
  !
  USE KINDS, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(INOUT) :: ix(N), iw(N)
  !
  INTEGER :: i
  REAL(DP) :: y
  do i=1,N
    y = 1.d0/(1.d0-ix(i))
    iw(i) = 2.d0*x0*iw(i)*y*y
    ix(i) = x0*(1.d0+ix(i))*y
  enddo
  !
  END SUBROUTINE transform_gaussleg
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Test Gauss-Legendre integration on f(x) = a / (a*a + x*x) = pi/2  
!-----------------------------------------------------------------------
  SUBROUTINE test_gaussleg_grid(N,a,ix,iw)
  !
  USE KINDS, ONLY : DP
  USE constants, ONLY: pi
  !
  implicit none
  INTEGER, INTENT(IN) :: N
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: ix(N), iw(N)
  !
  INTEGER :: i
  REAL(DP) :: err
  !
  err = 0.0d0
  do i=1,N
    err = err + iw(i) * a / (a*a + ix(i)*ix(i)) 
  enddo
  err = err - pi*0.5d0
  write(*,*) ' Testing GL integration grid. N, a, error: ',N,a,err
  !
  END SUBROUTINE test_gaussleg_grid
!-----------------------------------------------------------------------

