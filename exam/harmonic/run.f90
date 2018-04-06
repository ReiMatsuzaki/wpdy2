#include "macros.fpp"

module Mod_Main
  implicit none
contains
  subroutine run
    use Mod_WPDySOp
    use Mod_LambdaHarmonic
    use Mod_LambdaGWP
    use con
    integer, parameter :: nstate =  1    
    integer :: ierr, it
    integer, parameter :: nx = 256
    double precision :: L = 10.0d0, w
    
    call WPDySOp_new(nstate, nx, ierr); CHK_ERR(ierr)
    dt_ = 0.1d0
    m_ = 10.0d0

    w = 1.2d0
    k = m_*w*w

    x0 = 1.0d0
    p0 = 0.2d0
    a0 = 3.0d0
    
    call WPDy_set_xs(-L/2, L/nx, ierr); CHK_ERR(ierr)
    call WPDy_set_psi0(calc_psi0, ierr); CHK_ERR(ierr)
    call WPDy_set_vs(calc_vs, ierr); CHK_ERR(ierr)
    call WPDySOp_setup(ierr); CHK_ERR(ierr)
     
    do it = 0, 100
       call WPDy_con(it, ierr)
       call WPDySOp_inte(ierr)
       !do ix = 1, nx
       !   call calc_psit(xs_(ix), abs(dt_)*it, f(ix), ierr)
       !end do
       !call con_wf1("fa_re", it, real(f))
       !call con_wf1("fa_im", it, aimag(f))
    end do
 
    call WPDySOp_delete(ierr)  
  end subroutine run
end module Mod_Main

program main
  use Mod_Main
  call run
end program main
