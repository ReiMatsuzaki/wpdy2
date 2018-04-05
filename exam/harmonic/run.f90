#include "macros.fpp"
module Mod_Harmo0
  double precision :: x0 = 1.0d0
  double precision :: p0 = 0.0d0
  double precision :: a0 = 1.0d0
  double precision :: m = 1.0d0
  double precision :: k = 1.0d0
contains
  subroutine calc_psi0(x, res, ierr)
    use Mod_const, only : II
    double precision, intent(in) :: x
    complex(kind(0d0)), intent(out) :: res(:)
    integer, intent(out) :: ierr
    ierr = 0
    res(:) = 0
    res(1) = exp(-a0*(x-x0)**2 + II*p0*(x-x0))
  end subroutine calc_psi0
  subroutine calc_psit(x, t, res, ierr)
    use Mod_const, only : PI, II
    double precision, intent(in) :: x, t
    complex(kind(0d0)), intent(out) :: res
    integer, intent(out) :: ierr
    double precision xt, pt, w, gt
    complex(kind(0d0)) at
    ierr = 0
    w = sqrt(k/m)
    xt = x0*cos(w*t) + p0/(w*m)*sin(w*t)
    pt = (-x0*sin(w*t) + p0/(w*m)*cos(w*t)) *w*m
    at = a0
    gt = 0.0d0
    res = (2*real(at)/PI)**0.25d0 * exp(-at*(x-xt)**2 + II*pt*(x-xt) + II*gt)
  end subroutine calc_psit
  subroutine calc_vs(x, res, ierr)
    double precision, intent(in) :: x
    double precision, intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    ierr = 0
    res = x**2
  end subroutine calc_vs
end module Mod_Harmo0

module Mod_Main
  implicit none
contains
  subroutine run
    use Mod_WPDySOp
    use Mod_Harmo0
    integer, parameter :: nstate =  1    
    integer :: ierr, it
    
    call WPDySOp_new(nstate, 256, ierr); CHK_ERR(ierr)
    dt_ = 0.1d0
    call WPDy_set_xs(-5.0d0, 0.1d0, ierr); CHK_ERR(ierr)
    call WPDy_set_psi0(calc_psi0, ierr); CHK_ERR(ierr)
    call WPDy_set_vs(calc_vs, ierr); CHK_ERR(ierr)
    call WPDySOp_setup(ierr); CHK_ERR(ierr)
     
    do it = 0, 2
       call WPDy_con(it, ierr)
       call WPDySOp_inte(ierr)
    end do
 
    call WPDySOp_delete(ierr)  
  end subroutine run
end module Mod_Main

program main
  use Mod_Main
  call run
end program main
