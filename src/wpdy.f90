#include "macros.fpp"

module Mod_WPDy
  implicit none
  integer :: nx_, nstate_
  double precision :: dx_
  double precision, allocatable :: xs_(:), vs_(:,:,:), cvs_(:), frs_(:,:), frs0_(:,:)
  double precision :: m_
  logical :: setupq_=.false.
contains
  ! -- main --
  subroutine WPDy_new(nstate, nx, ierr)
    integer, intent(in) :: nstate, nx
    integer, intent(out) :: ierr
    nx_ = nx
    nstate_ = nstate
    dx_ = 0.1d0
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
  subroutine WPDy_set_vs(ifile, ierr)
    integer, intent(in) :: ifile
    integer, intent(out) :: ierr
    integer ii, jj, ix
    double precision val
    ierr = 0
    read(ifile, *)
    do
       read(ifile, *, end=100) ii, jj, ix, val
       vs_(ii,jj,ix) = val
    end do
100 continue
  end subroutine WPDy_set_vs
  subroutine WPDy_set_frs(ifile, ierr)
    integer, intent(in) :: ifile
    integer, intent(out) :: ierr
    integer ix, n
    double precision re, im
    ierr = 0
    read(ifile, *)
    do
       read(ifile, *, end=100) ix, n, re, im
       frs_(n,2*(ix-1))   = re
       frs_(n,2*(ix-1)+1) = im
    end do
100 continue
  end subroutine WPDy_set_frs
  subroutine WPDy_setup(ierr)
    integer, intent(out) :: ierr
    double precision norm
    double precision, parameter :: tol = 1.0d-10
    ierr = 0
    call WPDy_rn(0, norm, ierr); CHK_ERR(ierr)
    norm = sqrt(norm)
    if(norm < tol) then
       MSG_ERR("norm is too small")
       write(0,*) "norm:", norm
       ierr = 1; return
    end if
    frs_(:,:) = frs_(:,:) / norm
    frs0_(:,:) = frs_(:,:)
    setupq_ = .true.
  end subroutine WPDy_setup
  subroutine WPDy_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(xs_, vs_, cvs_, frs_, frs0_)
  end subroutine WPDy_delete
  ! -- calc --
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
