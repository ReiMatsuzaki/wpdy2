#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_UTestWpdy
  use Mod_UTest
  use Mod_WPDySOp
  implicit none
contains
  subroutine calc_psit(x0, p0, a0, m, w, x, t, res, ierr)
    use Mod_const, only : II, PI
    double precision, intent(in) :: x0, p0, a0, m, w, x, t
    complex(kind(0d0)), intent(out) :: res
    integer, intent(out) :: ierr
    double precision :: a, xt, pt
    complex(kind(0d0)) :: at, gt
    ierr = 0
    a = m*w/2
    xt = x0*cos(w*t) + p0/(w*m)*sin(w*t)
    pt = p0*cos(w*t) - m*x0*w  *sin(w*t)
    at = a * (a0*cos(w*t) + II*a*sin(w*t)) / (II*a0*sin(w*t) + a*cos(w*t))
    gt = (pt*xt-p0*x0)/2 + II/2*log(II*a0/a*sin(w*t) + cos(w*t))
    res = (2*real(at)/PI)**0.25d0 * exp(-at*(x-xt)**2 + II*pt*(x-xt) + II*gt)
  end subroutine calc_psit
  subroutine UTestWPDy_run
    use Mod_Timer
    type(Obj_Timer) :: timer
    integer ierr
    call Timer_new(timer, "UTestWPDy", .true., ierr)
    write(*,*) 
    write(*,*) "UTestWPDy begin"
    write(*,*)

    call Timer_begin(timer, "harmonic", ierr)
    call test_harmonic
    call Timer_end(  timer, "harmonic", ierr)
    
    write(*,*) 
    write(*,*) "UTestWPDy end"
    write(*,*) 

    call Timer_delete(timer, ierr)
    
  end subroutine UTestWPDy_run
  subroutine test_harmonic
    use Mod_LambdaHarmonic
    use Mod_LambdaGWP
    integer it, ierr
    integer, parameter :: nstate=1
    double precision :: norm2, r, p
    integer, parameter :: nx = 256
    double precision :: L = 10.0d0, w
    complex(kind(0d0)) :: y0, y1

    call WPDySOp_new(nstate, nx, ierr); CHK_ERR(ierr)
    dt_ = 0.1d0
    m_ = 10.0d0

    ! Harmonics
    w = 1.2d0
    k = m_*w*w

    ! GWP
    x0 = 1.0d0
    p0 = 0.2d0
    a0 = 3.0d0
    
    call WPDy_set_xs(-L/2, L/nx, ierr); CHK_ERR(ierr)    
    call WPDy_set_psi0(calc_psi0, ierr); CHK_ERR(ierr)
    call WPDy_set_vs(calc_vs, ierr); CHK_ERR(ierr)
    call WPDySOp_setup(ierr); CHK_ERR(ierr)
    call WPDy_psi(1, 120, y0, ierr)    
    call calc_psit(x0, p0, a0, m_, w, xs_(120), 0.0d0, y1, ierr)
    EXPECT_NEAR_C(y0, y1, 1.0d-14, ierr)
    
    call WPDy_rn(0, norm2, ierr)
    EXPECT_EQ_D(1.0d0, norm2, ierr)
    call WPDySOp_pn(0, norm2, ierr)
    EXPECT_EQ_D(1.0d0, norm2, ierr)

    call WPDy_rn(1, r, ierr)
    call WPDySOp_pn(1, p, ierr)
    EXPECT_NEAR_D(x0, r, 1.0d-14, ierr)
    EXPECT_NEAR_D(p0, p, 1.0d-14, ierr)
    
    do it = 1, 4
       call WPDySOp_inte(ierr)
    end do
    call WPDy_rn(0, norm2, ierr)
    EXPECT_EQ_D(1.0d0, norm2, ierr)

    call WPDy_psi(1, 120, y0, ierr)
    call calc_psit(x0, p0, a0, m_, w, xs_(120), 4*abs(dt_), y1, ierr)
    EXPECT_NEAR_C(y0, y1, 1.0d-3, ierr)
    
    call WPDySOp_delete(ierr)
    
  end subroutine Test_Harmonic
end module Mod_UTestWpdy

program main
  use Mod_UTestWpdy  
  call UTestWPDy_run
end program main
