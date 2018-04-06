#include "macros.fpp"

module Mod_WPDy
  implicit none
  integer :: nx_, nstate_
  double precision :: dx_
  double precision dt_
  integer nt_, ntskip_    
  double precision, allocatable :: xs_(:), vs_(:,:,:), cvs_(:), frs_(:,:), frs0_(:,:)
  double precision :: m_
  complex(kind(0d0)) :: dtdy_
contains
  ! -- main --
  subroutine WPDy_new(nstate, nx, ierr)
    integer, intent(in) :: nstate, nx
    integer, intent(out) :: ierr
    nx_ = nx
    nstate_ = nstate
    dx_ = 0.1d0
    dt_ = 1.0d0
    nt_ = 100
    ntskip_ = 1
    m_ = 1.0d0
    ierr = 0
    allocate(xs_(nx))                     ! grid points
    allocate(vs_(nstate_, nstate_, nx))     ! potential matrix
    allocate(cvs_(nx))    ! imaginary part of potential matrix
    allocate(frs_(nstate_, 0:2*nx-1))      ! real and imaginary part of WP
    allocate(frs0_(nstate_, 0:2*nx-1))     ! WP at t = 0
    
    xs_(:)     = 0.0d0
    vs_(:,:,:) = 0.0d0
    cvs_(:)    = 0.0d0
    frs_(:,:)  = 0.0d0
    frs0_(:,:) = 0.0d0
    
  end subroutine WPDy_new
  subroutine WPDy_set_xs(x0, dx, ierr)
    double precision, intent(in) :: x0, dx
    integer, intent(out) :: ierr
    integer i
    ierr = 0
    do i = 1, nx_
       xs_(i) = x0 + (i-1)*dx       
    end do
  end subroutine WPDy_set_xs
  subroutine WPDy_set_vs(calc_v, ierr)
    interface
       subroutine calc_v(x, v, ierr)
         double precision, intent(in) :: x
         double precision, intent(out) :: v(:,:)
         integer, intent(out) :: ierr
       end subroutine calc_v
    end interface
    integer, intent(out) :: ierr
    integer ix
    ierr = 0
    do ix = 1, nx_
       call calc_v(xs_(ix), vs_(:,:,ix), ierr); CHK_ERR(ierr)
    end do
  end subroutine WPDy_set_vs
  subroutine WPDy_set_psi0(calc_psi0, ierr)
    interface
       subroutine calc_psi0(x, y, ierr)
         double precision, intent(in) :: x
         complex(kind(0d0)), intent(out) :: y(:)
         integer, intent(out) :: ierr
       end subroutine calc_psi0
    end interface
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: y(nstate_)
    integer ix, n
    ierr = 0
    do ix = 1, nx_       
       call calc_psi0(xs_(ix), y(:), ierr)       
       do n = 1, nstate_
          frs_(n,2*(ix-1))   = real(y(n))
          frs_(n,2*(ix-1)+1) = aimag(y(n))
       end do
    end do    
  end subroutine WPDy_set_psi0
  subroutine WPDy_setup(ierr)
    integer, intent(out) :: ierr
    double precision norm
    double precision, parameter :: tol = 1.0d-10
    ierr = 0
    dx_ = xs_(2)-xs_(1)
    call WPDy_rn(0, norm, ierr); CHK_ERR(ierr)
    norm = sqrt(norm)
    if(norm < tol) then
       MSG_ERR("norm is too small")
       write(0,*) "norm:", norm
       ierr = 1; return
    end if    
    frs_(:,:) = frs_(:,:) / norm
    frs0_(:,:) = frs_(:,:)
    dtdy_ = dt_/ntskip_
  end subroutine WPDy_setup
  subroutine WPDy_con(it, ierr)
    use Mod_const, only : II
    use con ! mangan4
    integer, intent(in) :: it
    integer, intent(out) :: ierr
    double precision :: d, prob(nstate_)
    integer n
    ierr = 0

    if(it.eq.0) then
       call con_wi("_nx", 0, nx_)
       call con_wi("_nstate", 0, nstate_)
       call con_wf("_dx", 0, dx_)
       call con_wf("_dt", 0, dt_)
       call con_wf("_dtdy", 0, abs(dtdy_))
       call con_wf1("_x", 0, xs_)
       call con_wf3("_v", 0, vs_)
       call con_wf("_m",  0, m_)
    end if
    
    call con_wf("t", it,  it*dt_)
    call con_wf2("fr", it, frs_)

    call WPDy_rn(0, d, ierr)
    call con_wf("norm", it, d)
    
    call WPDy_rn(1, d, ierr)
    call con_wf("r1", it, d)

    do n = 1, nstate_
       call WPDy_prob(n, prob(n), ierr)
    end do
    call con_wf1("prob", it, prob(:))

  end subroutine WPDy_con
  subroutine WPDy_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(xs_, vs_, cvs_, frs_, frs0_)
  end subroutine WPDy_delete
  ! -- calc --
  subroutine WPDy_psi(n, i, res, ierr)
    complex(kind(0d0)), intent(out) :: res
    integer, intent(in)  :: i, n
    integer, intent(out) :: ierr
    
    ierr = 0
    res = dcmplx(frs_(n,2*(i-1)), frs_(n,2*(i-1)+1))
    
  end subroutine WPDy_psi
  subroutine WPDy_rn(n, res, ierr)
    integer, intent(in) :: n
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    integer ii, i
    ierr = 0
    res = 0.0d0
    do ii = 1, nstate_
       do i = 0, nx_-1
          res = res + (frs_(ii,2*i)**2 + frs_(ii,2*i+1)**2) * xs_(i+1)**n
       end do
    end do
    res = res*dx_
  end subroutine WPDy_rn
  subroutine WPDy_ac(res, ierr)
    complex(kind(0d0)), intent(out) :: res
    integer, intent(out) :: ierr
    double precision :: reac, imac
    integer i, ii
    ierr = 0
    reac = 0.0d0
    imac = 0.0d0
    do i = 0, nx_-1
       do ii = 1, nstate_
          reac = reac + frs0_(ii,2*i)*frs_(ii,2*i)   + frs0_(ii,2*i+1)*frs_(ii,2*i+1)
          imac = imac + frs0_(ii,2*i)*frs_(ii,2*i+1) - frs0_(ii,2*i+1)*frs_(ii,2*i)
       end do
    end do
    reac = reac * dx_
    imac = imac * dx_
    res = dcmplx(reac, imac)

  end subroutine WPDy_ac
  subroutine WPDy_prob(ii, res, ierr)
    integer, intent(in) :: ii
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    integer ix
    res = 0
    ierr = 0
    do ix = 0, nx_-1
       res = res + frs_(ii,2*ix)**2 + frs_(ii,2*ix+1)**2
    end do
    res = res * dx_
  end subroutine WPDy_prob
end module Mod_WPDy
