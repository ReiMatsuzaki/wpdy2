#include "macros.fpp"
  
module Mod_FFT
  implicit none
  private
  integer                       :: n2
  double precision, allocatable :: w(:)
  integer, allocatable          :: ip(:)
  public :: FFT_new, FFT_delete, FFT_forward, FFT_backward  
contains
  ! -- constructors --
  subroutine FFT_new(in_n2, ierr)
    integer, intent(in) :: in_n2
    integer, intent(out) :: ierr
    integer n

    ierr = 0
    n = in_n2/2
    do
       if(n.eq.1) goto 100
       if(mod(n,2).ne.0) then
          MSG_ERR("illegal in_n2")
          write(0,*) "in_n2:", n2
          ierr = 1; return
       end if
       n = n/2
    end do
100 continue
    
    
    n2 = in_n2
    allocate(w(0: n2/2-1))
    allocate(ip(0: 2+n2))
    ip(0) = 0

  end subroutine FFT_new
  subroutine fft_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(w)
    deallocate(ip)
    
  end subroutine fft_delete
  ! -- Calc --  
  subroutine FFT_forward(a, j0, k0)
    double precision    :: a(0:*)
    integer, intent(in), optional :: j0, k0
    integer val_j0, val_k0

    if(present(j0)) then
       val_j0 = j0
    else
       val_j0 = 0
    end if
    if(present(k0)) then
       val_k0 = k0
    else
       val_k0 = 0
    end if

    call shift_cdft_impl(a, val_j0, val_k0)
    
  end subroutine fft_forward
  subroutine FFT_backward(a, j0, k0)
    double precision    :: a(0:*)
    integer, intent(in), optional :: j0, k0
    integer val_j0, val_k0
    if(present(j0)) then
       val_j0 = j0
    else
       val_j0 = 0
    end if
    if(present(k0)) then
       val_k0 = k0
    else
       val_k0 = 0
    end if
    call shift_cdft_inv_impl(a, val_j0, val_k0)
    
  end subroutine fft_backward
  ! -- private --
  subroutine shift_cdft_impl(a, j0, k0)
    implicit none
    double precision ::  a(0: *)
    double precision :: re, im
    integer j, k, j0, k0
    double precision pi, t
    integer n

    n = n2
    pi = acos(-1.0d0)
    
    do j = 0, n/2-1
       re = a(2*j)
       im = a(2*j+1)
       t = 2.0d0*pi*j*k0/(n/2)
       a(2*j)   = re*cos(t) - im*sin(t)
       a(2*j+1) = re*sin(t) + im*cos(t)
    end do
  
    call cdft(n, 1, a, ip, w)
    
    do k = 0, n/2-1
       re = a(2*k)
       im = a(2*k+1)
       t = 2.0d0*pi*(j0*k0+j0*k)/(n/2)
       a(2*k)   = re*cos(t) - im*sin(t)
       a(2*k+1) = re*sin(t) + im*cos(t)
    end do
  end subroutine shift_cdft_impl
  subroutine shift_cdft_inv_impl(a, j0, k0)
    implicit none
    
    double precision ::  a(0: *)
    double precision :: re, im
    integer j, k, j0, k0
    double precision pi, t
    integer n
    pi = acos(-1.0d0)    

    n = n2
    
    do j = 0, n/2-1
       re = a(2*j)
       im = a(2*j+1)
       t = -2.0d0*pi*j*k0/(n/2)
       a(2*j)   = re*cos(t) - im*sin(t)
       a(2*j+1) = re*sin(t) + im*cos(t)
    end do
    
    call cdft(n, -1, a, ip, w)
    
    do k = 0, n/2-1
       re = a(2*k)
       im = a(2*k+1)
       t = -2.0d0*pi*(j0*k0+j0*k)/(n/2)
       a(2*k)   = re*cos(t) - im*sin(t)
       a(2*k+1) = re*sin(t) + im*cos(t)
    end do
  end subroutine shift_cdft_inv_impl
end module mod_fft

module mod_fft2d
  implicit none
  private
  integer :: n1, n2
  double precision, allocatable :: t(:)
  integer, allocatable          :: ip(:)
  double precision, allocatable :: w(:)  
  public :: fft2d_begin, fft2d_end, &
       fft2d_forward, fft2d_backward, &
       fft2d_forward_slow
contains
  ! -- constructors --
  subroutine fft2d_begin(in_n1, in_n2)
    integer, intent(in) :: in_n1, in_n2
    n1 = in_n1
    n2 = in_n2

    allocate(t(0:8*n2-1))
    allocate(ip(0:2+max(n1,n2)))
    allocate(w(0:max(n1/2,n2/2)))

    ip(0) = 0
    
  end subroutine fft2d_begin
  subroutine fft2d_end
    deallocate(t)
    deallocate(ip)
    deallocate(w)    
  end subroutine fft2d_end
  ! -- main --
  subroutine fft2d_forward(a, j1_0, k1_0, j2_0, k2_0)
    double precision :: a(0:2*n1-1, 0:n2-1)
    integer :: j1_0, k1_0, j2_0, k2_0
    integer :: j1, k1, j2, k2
    double precision :: pi, re, im, theta

    pi = acos(-1.0d0)
    do j1 = 0, n1-1
       do j2 = 0, n2-1
          re = a(2*j1,   j2)
          im = a(2*j1+1, j2)
          theta = 2.0d0*pi*(1.0d0*j1*k1_0/n1 + 1.0d0*j2*k2_0/n2)
          a(2*j1, j2)   = (re*cos(theta) - im*sin(theta))
          a(2*j1+1, j2) = re*sin(theta) + im*cos(theta)
       end do
    end do
    
    call cdft2d(n1*2, 2*n1, n2, +1, a, t, ip, w)

    do k1 = 0, n1-1
       do k2 = 0, n2-1
          re = a(2*k1,   k2)
          im = a(2*k1+1, k2)
          theta = 2.0d0*pi*((j1_0*k1_0 + j1_0*k1)/(1.0d0*n1) + (j2_0*k2_0 + j2_0*k2)/(1.0d0*n2))
          a(2*k1,   k2) = re*cos(theta) - im*sin(theta)
          a(2*k1+1, k2) = re*sin(theta) + im*cos(theta)
       end do
    end do
    
  end subroutine fft2d_forward
  subroutine fft2d_backward(a, j1_0, k1_0, j2_0, k2_0)

    double precision :: a(0:2*n1-1, 0:n2-1)
    integer :: j1_0, k1_0, j2_0, k2_0
    integer :: j1, k1, j2, k2
    double precision :: pi, re, im, theta

    pi = acos(-1.0d0)
    do j1 = 0, n1-1
       do j2 = 0, n2-1
          re = a(2*j1,   j2)
          im = a(2*j1+1, j2)
          theta = -2.0d0*pi*(1.0d0*j1*k1_0/n1 + 1.0d0*j2*k2_0/n2)
          a(2*j1, j2)   = re*cos(theta) - im*sin(theta)
          a(2*j1+1, j2) = re*sin(theta) + im*cos(theta)
       end do
    end do

    call cdft2d(n1*2, 2*n1, n2, -1, a, t, ip, w)

    do k1 = 0, n1-1
       do k2 = 0, n2-1
          re = a(2*k1,   k2)
          im = a(2*k1+1, k2)
          theta = -2.0d0*pi*((j1_0*k1_0 + j1_0*k1)/(1.0d0*n1) + (j2_0*k2_0 + j2_0*k2)/(1.0d0*n2))
          a(2*k1,   k2) = re*cos(theta) - im*sin(theta)
          a(2*k1+1, k2) = re*sin(theta) + im*cos(theta)
       end do
    end do
    
  end subroutine fft2d_backward
  subroutine fft2d_forward_slow(a, j1_0, k1_0, j2_0, k2_0)
    !double precision    :: a(0:, 0:)
    double precision    :: a(0:2*n1-1, 0:n2-1)
    integer, intent(in) :: j1_0, k1_0, j2_0, k2_0
    double precision    :: pi    
    integer j1, j2, k1, k2
    double precision    :: re, im, rea, ima, t
    double precision    :: b(0:100, 0:100)
    pi = acos(-1.0d0)

    do k1 = k1_0, n1-1+k1_0
       do k2 = k2_0, n2-1+k2_0
          re = 0.0d0
          im = 0.0d0
          do j1 = j1_0, n1-1+j1_0
             do j2 = j2_0, n2-1+j2_0
                rea = a(2*(j1-j1_0),   j2-j2_0)
                ima = a(2*(j1-j1_0)+1, j2-j2_0)
                t = 2.0d0*pi*((1.0d0*j1*k1)/n1 + (1.0d0*j2*k2)/n2)
                re = re + rea*cos(t) - ima*sin(t)
                im = im + rea*sin(t) + ima*cos(t)
             end do
          end do
          b(2*(k1-k1_0),   k2-k2_0) = re
          b(2*(k1-k1_0)+1, k2-k2_0) = im
       end do
    end do

    do k1=0, n1-1
       do k2=0, n2-1
          a(2*k1,   k2) = b(2*k1,   k2)
          a(2*k1+1, k2) = b(2*k1+1, k2)
       end do
    end do
  end subroutine fft2d_forward_slow
end module mod_fft2d
