#include "macros.fpp"

module Mod_WPDySOp
  use Mod_WPDy
  use Mod_fft
  implicit none
  double precision :: dk_, K_
  double precision, allocatable :: ks_(:)
contains
  ! -- main --
  subroutine WPDySOp_new(nstate, nx, ierr)
    use Mod_FFT
    use Mod_Const, only : pi
    integer, intent(in) :: nstate
    integer, intent(in) :: nx
    integer, intent(out) :: ierr
    ierr = 0

    if(nstate.ne.2 .and. nstate.ne.1) then
       MSG_ERR("illegal value: nstate_")
       write(0,*) "nstate:", nstate
       ierr = 1; return
    end if

    call FFT_new(2*nx, ierr); CHK_ERR(ierr)
    call WPDy_new(nstate, nx, ierr); CHK_ERR(ierr)
    allocate(ks_(nx))
    
  end subroutine WPDySOp_new
  subroutine WPDySOp_setup(ierr)
    use Mod_const, only : PI
    integer, intent(out) :: ierr
    integer i
    
    if(nstate_.ne.2 .and. nstate_.ne.1) then
       write(0,*) "only nstate=2 is supported"
       stop
    end if

    call WPDy_setup(ierr); CHK_ERR(ierr)

    K_ = pi/dx_
    dk_ = 2*pi/(nx_*dx_)
    do i = 1, nx_
       ks_(i) = (i-nx_/2)*dk_
    end do
    
  end subroutine WPDySOp_setup
  subroutine WPDySOp_delete(ierr)
    integer, intent(out) :: ierr
    call WPDy_delete(ierr); CHK_ERR(ierr)
    deallocate(ks_)
    call FFT_delete(ierr); CHK_ERR(ierr)
  end subroutine WPDySOp_delete
  ! ==== calc ====
  subroutine WPDySOp_inte(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    if(nstate_.eq.1) then
       call inte_1(dtdy_, ierr); CHK_ERR(ierr)
    else if(nstate_.eq.2) then
       call inte_2(dtdy_, ierr); CHK_ERR(ierr)
    else
       MSG_ERR("unsupported nstate")
       ierr = 1; return
    end if
  end subroutine WPDySOp_inte
  subroutine inte_1(dt, ierr)
    use Mod_const, only : ii
    complex(kind(0d0)), intent(in) :: dt
    integer, intent(out) :: ierr
    double precision :: sqrt_nx, ene
    integer :: i
    complex(kind(0d0)) :: fi    
    ! dx = X/N
    ! dk = 2pi/(N dx)
    ! dx dk = 2pi /N
    ! exp( i x_n k_m ) = exp( i (n dx) (m dk-K/2) )
    !                  = exp( 2pi i nm/N) exp(-i n dx K/2)
    !                  = exp( 2pi i nm/N) exp(-i xK/2)
    
    ierr = 0
    if(nstate_ .ne. 1) then
       MSG_ERR("nstate must be 1")
       ierr = 1; return
    end if
    sqrt_nx = sqrt(1.0d0 * nx_)
    do i = 0, nx_-1
       fi = dcmplx(frs_(1,2*i), frs_(1,2*i+1))
       fi = exp(-II*vs_(1,1,i+1)*0.5d0 * dt)*fi
       frs_(1,2*i)   = real(fi)
       frs_(1,2*i+1) = aimag(fi)
    end do
    call fft_backward(frs_(1,:), int(xs_(1)/dx_), -nx_/2+1)    
    frs_(1,:)   = frs_(1,:)/sqrt_nx
    do i = 0, nx_-1
       fi = dcmplx(frs_(1,2*i), frs_(1,2*i+1))
       ene = 1.0d0/(2.0d0*m_) * ks_(i+1)**2
       fi = exp(-II*ene*dt) * fi
       frs_(1,2*i)   = real(fi)
       frs_(1,2*i+1) = aimag(fi)
    end do
    call fft_forward(frs_, -nx_/2+1, int(xs_(1)/dx_))
    frs_(1,:)   = frs_(1,:)/sqrt_nx
    do i = 0, nx_-1
       fi = dcmplx(frs_(1,2*i), frs_(1,2*i+1))
       fi = exp(-II*vs_(1,1,i+1)*0.5d0 * dt)*fi
       frs_(1,2*i)   = real(fi)
       frs_(1,2*i+1) = aimag(fi)
    end do

    !do i = 0, 2*nx_-1
    !   if(abs(frs_(1,i)) < 1.0d-14) then
    !      frs_(1,i) = 0.0d0
    !   end if
    !end do
    
  end subroutine inte_1
  subroutine inte_2(dt, ierr)
    use Mod_const, only : ii
    complex(kind(0d0)), intent(in) :: dt
    integer, intent(out) :: ierr
    double precision :: sqrt_nx, ene
    complex(kind(0d0)) :: sqD, v1, v2, v12, v21
    integer          :: i, is
    complex(kind(0d0)) :: fi, f1i, f2i, tmp, ex
    ierr = 0
    if(nstate_ .ne. 2) then
       MSG_ERR("nstate must be 2")
       ierr = 1; return
    end if
    sqrt_nx = sqrt(1.0d0 * nx_)

    ! -- operate exp[-iVdt/2] --
    do i = 0, nx_-1
       v1 = dcmplx(vs_(1,1,i+1), cvs_(i+1)); v2 = dcmplx(vs_(2,2,i+1), cvs_(i+1))
       v12 = vs_(1,2,i+1); v21 = vs_(2,1,i+1)
       sqD = sqrt((v2-v1)**2 + 4.0d0*v12*v21)
       if(abs(sqD)<1.0d-10) then
          MSG_ERR("sqD is too small")
          ierr = 1; return
       end if
       ex = exp(-ii*(v1+v2)*dt/4.0d0)
       f1i = dcmplx(frs_(1,2*i), frs_(1,2*i+1))
       f2i = dcmplx(frs_(2,2*i), frs_(2,2*i+1))
       tmp = ex*(cos(sqD*dt/4.0d0)     * f1i + &
            ii*  sin(sqD*dt/4.0d0)/sqD * ((v2-v1)*f1i -2.0d0*v12*f2i))       
       frs_(1,2*i)   = real(tmp)
       frs_(1,2*i+1) = aimag(tmp)
       tmp = ex*(cos(sqD*dt/4.0d0)     * f2i + &
            ii*  sin(sqD*dt/4.0d0)/sqD * (-2.0d0*v21*f1i +(v1-v2)*f2i))
       frs_(2,2*i)   = real(tmp)
       frs_(2,2*i+1) = aimag(tmp)
    end do
    call fft_backward(frs_(1,:), int(xs_(1)/dx_), -nx_/2+1)
    call fft_backward(frs_(2,:), int(xs_(1)/dx_), -nx_/2+1)
    frs_(:,:) = frs_(:,:)/sqrt_nx

    ! -- operate exp[-iTdt] --
    do is = 1, nstate_
       do i = 0, nx_-1
          fi = dcmplx(frs_(is,2*i), frs_(is,2*i+1))
          ene = 1.0d0/(2.0d0*m_) * ks_(i+1)**2
          fi = exp(-II*ene*dt)*fi
          frs_(is,2*i) = real(fi)
          frs_(is,2*i+1) = aimag(fi)
       end do
    end do

    call fft_forward(frs_(1,:), -nx_/2+1, int(xs_(1)/dx_))
    call fft_forward(frs_(2,:), -nx_/2+1, int(xs_(1)/dx_))
    frs_(:, :) = frs_(:, :)/sqrt_nx

    ! -- operate exp[-iVdt/2] --
    do i = 0, nx_-1
       v1 = dcmplx(vs_(1,1,i+1), cvs_(i+1)); v2 = dcmplx(vs_(2,2,i+1), cvs_(i+1))
       v12 = vs_(1,2,i+1); v21 = vs_(2,1,i+1)
       sqD = sqrt((v2-v1)**2 + 4.0d0*v12*v21)
       ex = exp(-ii*(v1+v2)*dt/4.0d0)
       f1i = dcmplx(frs_(1,2*i), frs_(1,2*i+1))
       f2i = dcmplx(frs_(2,2*i), frs_(2,2*i+1))
       tmp = ex*(cos(sqD*dt/4.0d0)     * f1i + &
            ii*  sin(sqD*dt/4.0d0)/sqD * ((v2-v1)*f1i -2.0d0*v12*f2i))
       frs_(1,2*i)   = real(tmp)
       frs_(1,2*i+1) = aimag(tmp)
       tmp = ex*(cos(sqD*dt/4.0d0)     * f2i + &
            ii*  sin(sqD*dt/4.0d0)/sqD * (-2.0d0*v21*f1i +(v1-v2)*f2i))
       frs_(2,2*i)   = real(tmp)
       frs_(2,2*i+1) = aimag(tmp)
    end do
    
  end subroutine Inte_2
  subroutine WPDySOp_pn(n, res, ierr)
    use Mod_const, only : II
    integer, intent(in) :: n
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    double precision :: fr(0:2*nx_-1), sqrt_nx
    complex(kind(0d0)) :: ff
    integer i
    
    ierr = 0
    if(nstate_.ne.1) then
       MSG_ERR("nstate!=1")
       ierr = 1; return
    end if

    sqrt_nx = sqrt(1.0d0*nx_)
    fr(0:2*nx_-1) = frs_(1, 0:2*nx_-1)
    ! call fft_backward(fr(:), int(xs_(1)/dx_), -nx_/2+1)
    call fft_backward(fr(:))
    fr(:) = fr(:)/sqrt_nx
    do i = 0, nx_-1
       fr(2*i)   = fr(2*i)   * ks_(i+1)**n
       fr(2*i+1) = fr(2*i+1) * ks_(i+1)**n
    end do
    !call fft_forward(fr(:),  -nx_/2+1, int(xs_(1)/dx_))
    call fft_forward(fr(:))
    fr(:) = fr(:)/sqrt_nx
    res = 0
    do i = 0, nx_-1
       ff = dcmplx(frs_(1,2*i), -frs_(1,2*i+1)) * dcmplx(fr(2*i), fr(2*i+1))
       res = res + real(ff)
    end do
    res = res*dx_
    
  end subroutine WPDySOp_pn
end module Mod_WPDySOp
