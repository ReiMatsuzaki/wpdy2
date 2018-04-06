#include "macros.fpp"
module Mod_Harmo0
  double precision :: x0 = 1.0d0
  double precision :: p0 = 0.0d0
  double precision :: a0 = 5.0d0
  double precision :: m = 10.0d0
  double precision :: w = 1.0d0  
contains
  subroutine calc_psi0(x, res, ierr)
    use Mod_const, only : II, PI
    double precision, intent(in) :: x
    complex(kind(0d0)), intent(out) :: res(:)
    integer, intent(out) :: ierr
    ierr = 0
    res(:) = 0
    res(1) = (2*a0/PI)**0.25d0*exp(-a0*(x-x0)**2 + II*p0*(x-x0))
  end subroutine calc_psi0
  subroutine calc_psit(x, t, res, ierr)
    use Mod_const, only : PI, II
    double precision, intent(in) :: x, t
    complex(kind(0d0)), intent(out) :: res
    integer, intent(out) :: ierr
    double precision xt, pt, a
    complex(kind(0d0)) at, gt
    ierr = 0
    a = m*w/2
    xt = x0*cos(w*t) + p0/(w*m)*sin(w*t)
    pt = p0*cos(w*t) - m*x0*w  *sin(w*t)
    at = a * (a0*cos(w*t) + II*a*sin(w*t)) / (II*a0*sin(w*t) + a*cos(w*t))
    gt = (pt*xt-p0*x0)/2 + II/2*log(II*a0/a*sin(w*t) + cos(w*t))
    res = (2*real(at)/PI)**0.25d0 * exp(-at*(x-xt)**2 + II*pt*(x-xt) + II*gt)
  end subroutine calc_psit
  subroutine calc_vs(x, res, ierr)
    double precision, intent(in) :: x
    double precision, intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    double precision :: k 
    ierr = 0
    k = m*w*w
    res = k/2*x**2
  end subroutine calc_vs
end module Mod_Harmo0

module Mod_Main
  implicit none
contains
  subroutine run
    use Mod_WPDySOp
    use Mod_Harmo0
    use con
    integer, parameter :: nstate =  1    
    integer :: ierr, it, ix
    integer, parameter :: nx = 256
    double precision :: L = 10.0d0
    complex(kind(0d0)) :: f(nx)
    
    call WPDySOp_new(nstate, nx, ierr); CHK_ERR(ierr)
    dt_ = 0.1d0
    m_ = m
    call WPDy_set_xs(-L/2, L/nx, ierr); CHK_ERR(ierr)
    call WPDy_set_psi0(calc_psi0, ierr); CHK_ERR(ierr)
    call WPDy_set_vs(calc_vs, ierr); CHK_ERR(ierr)
    call WPDySOp_setup(ierr); CHK_ERR(ierr)
     
    do it = 0, 100
       call WPDy_con(it, ierr)
       call WPDySOp_inte(ierr)
       do ix = 1, nx
          call calc_psit(xs_(ix), abs(dt_)*it, f(ix), ierr)
       end do
       call con_wf1("fa_re", it, real(f))
       call con_wf1("fa_im", it, aimag(f))
    end do
 
    call WPDySOp_delete(ierr)  
  end subroutine run
end module Mod_Main

program main
  use Mod_Main
  call run
end program main
