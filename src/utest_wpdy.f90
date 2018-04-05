#include "macros.fpp"
#include "macros_utest.fpp"
module Mod_Harmo0
contains
  subroutine calc_psi0(x, res, ierr)
    use Mod_const, only : II
    double precision, intent(in) :: x
    complex(kind(0d0)), intent(out) :: res(:)
    integer, intent(out) :: ierr
    double precision x0, p0
    x0 = 0.0d0
    p0 = 0.0d0
    ierr = 0
    res(:) = 0
    res(1) = exp(-5*(x-x0)**2 + II*p0*(x-x0))
  end subroutine calc_psi0
  subroutine calc_vs(x, res, ierr)
    double precision, intent(in) :: x
    double precision, intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    ierr = 0
    res = x**2
  end subroutine calc_vs
end module Mod_Harmo0
module Mod_UTestWpdy
  use Mod_UTest
  use Mod_WPDySOp
  implicit none
contains
  subroutine UTestWPDy_run
    use Mod_Timer
    type(Obj_Timer) :: timer
    integer ierr
    call Timer_new(timer, "UTestWPDy", .true., ierr)
    write(*,*) 
    write(*,*) "UTestWPDy begin"
    write(*,*)

    call Timer_begin(timer, "test1", ierr)
    call Test1
    call Timer_end(  timer, "test1", ierr)

    call Timer_begin(timer, "test2", ierr)
    call Test2
    call Timer_end(  timer, "test2", ierr)    
    
    write(*,*) 
    write(*,*) "UTestWPDy end"
    write(*,*) 

    call Timer_delete(timer, ierr)
    
  end subroutine UTestWPDy_run
  subroutine Test1
    use Mod_const, only : ii
    use Mod_Harmo0
    integer it, ierr
    integer, parameter :: nstate=1
    double precision :: x0, p0, norm2
    
    x0 = 0.0d0
    p0 = 0.0d0
    call WPDySOp_new(nstate, 256, ierr); CHK_ERR(ierr)    
    call WPDy_set_xs(-5.0d0, 0.1d0, ierr); CHK_ERR(ierr)
    call WPDy_set_psi0(calc_psi0, ierr); CHK_ERR(ierr)
    call WPDy_set_vs(calc_vs, ierr); CHK_ERR(ierr)
    call WPDySOp_setup(ierr); CHK_ERR(ierr)
    
    call WPDy_rn(0, norm2, ierr)
    EXPECT_EQ_D(1.0d0, norm2, ierr)
    
    do it = 1, 2
       call WPDySOp_inte((0.1d0, 0.0d0), ierr)       
    end do
    call WPDy_rn(0, norm2, ierr)
    EXPECT_EQ_D(1.0d0, norm2, ierr)
    call WPDySOp_delete(ierr)
  end subroutine Test1
  subroutine Test2
    use Mod_const, only : ii
    integer ix, it, ierr
    integer, parameter :: nstate=2
    double precision :: x0, p0, norm2
    complex(kind(0d0)) :: f
    
    x0 = 1.0d0
    p0 = 1.0d0
    call WPDySOp_new(nstate, 256, ierr); CHK_ERR(ierr)
    call WPDy_set_xs(-5.0d0, 0.1d0, ierr); CHK_ERR(ierr)
    do ix = 1, nx_
       f = exp(-(xs_(ix)-x0)**2 + ii*p0*(xs_(ix)-x0))
       frs_(1,2*(ix-1))   = real(f)
       frs_(1,2*(ix-1)+1) = aimag(f)
       vs_(1,1,ix) = xs_(ix)**2
       vs_(2,2,ix) = xs_(ix)**2 + 1.0d0
       vs_(1,2,ix) = exp(-xs_(ix)**2)
       vs_(2,1,ix) = exp(-xs_(ix)**2)
    end do
    call WPDySOp_setup(ierr); CHK_ERR(ierr)
    
    do it = 1, 10
       call WPDySOp_inte((0.1d0, 0.0d0), ierr)
    end do
    call WPDy_rn(0, norm2, ierr)
    EXPECT_EQ_D(1.0d0, norm2, ierr)    
  end subroutine Test2
end module Mod_UTestWpdy

program main
  use Mod_UTestWpdy  
  call UTestWPDy_run
end program main
