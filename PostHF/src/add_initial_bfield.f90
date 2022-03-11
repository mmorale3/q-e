subroutine add_initial_bfield(v_of_r)
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE cell_base,        ONLY : omega
  USE fft_base,         ONLY : dfftp
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : noncolin, i_cons, mcons, lambda, &
                               pointlist, factlist
  !
  implicit none
  ! inputs and outputs
  real(DP), intent(out) :: v_of_r(dfftp%nnr,nspin)
  ! local variables
  real(DP) :: fact
  integer :: nt, na, ipol, ir
  !
  if (i_cons .ne. 1) call errore('add_initial_bfield', 'i_cons != 1', 1)
  v_of_r(:,:) = 0.d0
  ! compute local potential due to constraints
  do ir = 1, dfftp%nnr
    if (pointlist(ir) == 0 ) cycle
    na = pointlist(ir)
    nt = ityp(na)
    fact = -1.d0*lambda*factlist(ir)*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
    if (noncolin) then
      do ipol = 1,3
        v_of_r(ir,ipol+1) = v_of_r(ir,ipol+1) + fact*mcons(ipol,nt)
      enddo ! ipol
    else
      v_of_r(ir,1) = v_of_r(ir,1) + fact*mcons(1,nt)
      v_of_r(ir,2) = v_of_r(ir,2) - fact*mcons(1,nt)
    endif ! noncolin
  enddo ! ir
end subroutine add_initial_bfield
