! functions used for lambda in WPDY.
module Mod_LambdaEchart
  implicit none
  double precision :: a, xb, V0
contains
  subroutine calc_vs(x, res, ierr)
    double precision, intent(in) :: x
    double precision, intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    double precision :: sech, z
    ierr = 0
    z = a*(x-xb)
    sech = 2/(exp(z)+exp(-z))
    res = V0*sech**2
  end subroutine calc_vs
end module Mod_LambdaEchart

module Mod_LambdaHarmonic
  implicit none
  double precision :: k
contains
  subroutine calc_vs(x, res, ierr)
    double precision, intent(in) :: x
    double precision, intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    ierr = 0
    res = k/2*x**2
  end subroutine calc_vs
end module Mod_LambdaHarmonic

module Mod_LambdaGWP
  implicit none
  double precision :: x0, p0, a0
contains
  subroutine calc_psi0(x, res, ierr)
    use Mod_const, only : II, PI
    double precision, intent(in) :: x
    complex(kind(0d0)), intent(out) :: res(:)
    integer, intent(out) :: ierr
    ierr = 0
    res(:) = 0
    res(1) = (2*a0/PI)**0.25d0*exp(-a0*(x-x0)**2 + II*p0*(x-x0))
    if(abs(res(1)) < 1.0d-14) then
       res(1) = 0
    end if
  end subroutine calc_psi0
end module Mod_LambdaGWP
