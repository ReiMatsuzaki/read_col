module Mod_Utils
  implicit none
  real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
contains
  subroutine LA_GEES_GEN(A, B, W, VL)
    ! Solve generalized eigen values problem
    !    A UL = diag(W) B UL
    complex*16, intent(in) :: A(:, :)
    complex*16, intent(in) :: B(:, :)
    complex*16, intent(out) :: W(:, :)
    complex*16, intent(out) :: VL(:, :)

    complex*16, allocatable :: X(:, :)
    complex*16, allocatable :: YB(:, :)

    integer :: n

    n = size(A(:, 1))
    allocate(X(n, n))
    allocate(YB(n, n))

    ! diagonalize B matrix
    ! B = VLB^H diag(W) VLB
    call LA_GEES(B, X, YB)

    ! A VL = diag(W) UL^H diag(WB) UL VL
    ! (UL A UL^H) (UL VL) = diag(W) 
    
  end subroutine LA_GEES_GEN
end module Mod_Utils
 
