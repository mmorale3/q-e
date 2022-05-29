module iocc
  use kinds, only : dp
  USE wvfct, ONLY: et
  implicit none
  private

  public :: fill_bands, n_of_mu
  public :: unique

contains

  pure integer function fill_bands(et, mu) result(nmu)
    real(dp), intent(in) :: et(:,:), mu
    ! local
    integer :: ib, ik, nb, nk
    nb = size(et, 1)
    nk = size(et, 2)
    nmu = 0
    do ik=1,nk
      do ib=1,nb
        if (et(ib,ik) <= mu) nmu = nmu+1
      enddo
    enddo
  end function fill_bands

  subroutine unique(n, a, m, au)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n)
    integer, intent(out) :: m
    real(dp), intent(out) :: au(n)
    ! local
    integer :: i, j
    real(dp) :: tol, ai
    logical :: found
    tol = 1e-6 ! hard-code tolerance
    m = 0
    au(:) = 0
    do i=1,n
      ai = a(i)
      found = .false.
      do j=1,m
        if (abs(ai-au(j)) < tol) then
          found = .true.
          exit
        endif
      enddo
      if (.not.found) then
        m = m+1
        au(m) = ai
      endif
    enddo
    call quicksort(au, 1, m)
  end subroutine unique

  recursive subroutine quicksort(a, first, last)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: first, last
    ! local
    integer i, j
    real(dp) :: x, t
  
    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

  function n_of_mu(mxmu, mus, etv) result(nmu)
    real(dp) :: nmu(mxmu)
    integer, intent(in) :: mxmu
    real(dp), intent(in) :: mus(mxmu)
    real(dp), intent(in) :: etv(:)
    ! local
    integer :: ib, nb, imu
    real(dp) :: e1
    nb = size(etv)
    nmu(:) = 0
    do ib=1,nb
      e1 = etv(ib)
      do imu=1,mxmu
        if (e1 <= mus(imu)) nmu(imu) = nmu(imu) + 1
      enddo
    enddo
  end function n_of_mu

end module iocc
