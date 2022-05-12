module madelung
  use kinds, only: dp
  use constants, only: pi
  implicit none
  real(dp) :: volume, alpha
  integer :: ndim
contains

  subroutine madelung_init(alat, at, ndim_in)
    real(dp), intent(in) :: alat, at(3,3)
    integer, intent(in) :: ndim_in
    alpha = 2.0d0/alat
    ndim = ndim_in
    if ((ndim .ne. 2).and.(ndim .ne. 3)) then
      call errore('madelung_sum', 'unknown ndim', ndim)
    endif
    volume = alat**ndim*det(at(:ndim,:ndim),ndim)
  end subroutine

  pure function madelung_sum(alat, at, bg) result(vmad)
    real(dp), intent(in) :: alat, at(3,3), bg(3,3)
    real(dp) :: vmad
    ! local
    real(dp) :: r(3), rmag, blat
    real(dp) :: esr, elr, vlr_r0, vsr_k0
    integer :: i, j, k, mx, my, mz
    integer :: nr, nk

    !!!! hard-code number of lattices to sum over
    nr = 5
    nk = 5

    ! constants
    vlr_r0 = 2.d0*alpha/sqrt(pi)
    if (ndim.eq.2) then
      vsr_k0 = 2.d0*sqrt(pi)/alpha/volume
    else if (ndim.eq.3) then
      vsr_k0 = pi/(alpha*alpha)/volume
    end if

    ! sums
    ! short-range
    esr = 0.d0
    mx = nr; my = nr; mz = nr
    if (ndim.eq.2) mz=0
    do concurrent (i=-mx:mx, j=-my:my, k=-mz:mz)
      if ((i.eq.0).and.(j.eq.0).and.(k.eq.0)) cycle
      r = i*at(:,1)+j*at(:,2)+k*at(:,3)
      rmag = alat*sqrt(dot_product(r, r))
      esr = esr + vsr_r(rmag)
    enddo
    esr = 0.5*esr
    ! long-range
    blat = 2*pi/alat
    elr = 0.d0
    mx = nk; my = nk; mz = nk
    if (ndim.eq.2) mz=0
    do concurrent (i=-mx:mx, j=-my:my, k=-mz:mz)
      if ((i.eq.0).and.(j.eq.0).and.(k.eq.0)) cycle
      r = i*bg(:,1)+j*bg(:,2)+k*bg(:,3)
      rmag = blat*sqrt(dot_product(r, r))
      if (ndim.eq.2) then
        elr = elr + vlr_k_2d(rmag)
      else if (ndim.eq.3) then
        elr = elr + vlr_k_3d(rmag)
      endif
    enddo
    elr = 0.5*elr

    ! assemble
    vmad = esr + elr -0.5*(vlr_r0+vsr_k0)
  end function madelung_sum

  pure function vlr_k_2d(k) result(val)
    real(dp), intent(in) :: k
    real(dp) :: val
    val = 2*pi/k*erfc(k/(2*alpha))/volume
  end function vlr_k_2d

  pure function vlr_k_3d(k) result(val)
    real(dp), intent(in) :: k
    real(dp) :: val
    ! local
    real(dp) :: k2, b2
    k2 = k*k
    b2 = (2*alpha)*(2*alpha)
    val = 4*pi/k2*exp(-k2/b2)/volume
  end function vlr_k_3d

  pure function vsr_r(r) result(val)
    real (dp), intent(in) :: r
    real(dp) :: val
    val = erfc(alpha*r)/r
  end function vsr_r

  recursive pure function det(a, n) result(acc)
    ! https://rosettacode.org/wiki/Determinant_and_permanent#Fortran
    integer, intent(in) :: n
    real(dp), dimension(n,n), intent(in) :: a(:,:)
    real(dp), dimension(n-1,n-1) :: b
    real(dp) :: acc
    ! local variables
    integer :: i, sgn
    if (n .eq. 1) then
      acc = a(1,1)
    else
      acc = 0
      sgn = 1
      do i=1, n
        b(:, :(i-1)) = a(2:, :i-1)
        b(:, i:) = a(2:, i+1:)
        acc = acc + sgn * a(1, i) * det(b, n-1)
        sgn = -1 * sgn
      enddo
    endif
  end function det
end module madelung
