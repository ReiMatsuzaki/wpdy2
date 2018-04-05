#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_UTestFFT
  implicit none
contains
  subroutine UTestFFT_run
    call test_fft2
    call test_fft3
  end subroutine UTestFFT_run
  subroutine slow_cdft(a, n, b)
    implicit none
    integer :: n
    double precision :: a(0:2*n-1)
    double precision :: b(0:2*n-1)
    complex(kind(0d0)) :: v
    integer :: i, j
    double precision :: pi
    pi = acos(-1.0d0)
    do i = 0, n-1
       v = 0.0d0
       do j = 0, n-1
          v = v + dcmplx(a(2*j), a(2*j+1)) *&
               exp(dcmplx(0.0d0, 2.0d0*pi*i*j/n))
       end do
       b(2*i)   = real(v)
       b(2*i+1) = aimag(v)
    end do
    
  end subroutine slow_cdft
  subroutine test_fft2
    use mod_utest
    use mod_fft
    implicit none

    integer, parameter :: n = 2
    double precision   :: a(0:2*n*n-1)
    integer :: i, j

    write(*,*) "start test fft2"

    do i = 0, n-1
       do j = 0, n-1
          a(2*n*i+2*j)   = 10.0d0*i + j
          a(2*n*i+2*j+1) = 1000.0d0*i + 100.0d0*j
       end do
    end do

    write(*, *) a

    call fft_begin(2*n)
    call transpose_as_cmat(a)
    write(*,*) 

    write(*,*) a

    call fft_end

  end subroutine test_fft2
  subroutine test_fft3
    use mod_utest
    use mod_fft
    implicit none

    integer, parameter :: n = 4
    double precision   :: a(0:2*n*n-1), b(0:2*n*n-1)
    integer :: i, j

    write(*,*) "start test fft3"

    do i = 0, n-1
       do j = 0, n-1
          a(2*n*i + 2*j)     = sin(1.12d0*i)+0.1d0*j
          a(2*n*i + 2*j + 1) = cos(1.2d0*i) + 0.01d0*j
       end do
    end do

    b(:) = a(:)

    call fft_begin(2*n)
    call fft_2d_forward(a, 0, 0, 0,0)
    call fft_2d_forward_slow(b,0,0,0,0)
    write(*,*) a(4), b(4)

    call fft_end

  end subroutine test_fft3
end module UTest_FFT
program main
  use Mod_UTestFFT
  call UTestFFT_run
end program main
