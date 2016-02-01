module Mod_Utils
  implicit none
  real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
contains
  subroutine LA_GEEV_GEN(A, B, Wa, Va)
    use f95_lapack, only : LA_GEEV
    ! Solve generalized eigen values problem
    !    A UL = diag(W) B UL
    complex*16, intent(inout) :: A(:, :)
    complex*16, intent(inout) :: B(:, :)
    complex*16, intent(out) :: Wa(:)
    complex*16, intent(out) :: Va(:, :)

    complex*16, allocatable :: Bph(:, :)   ! B^(+1/2)
    complex*16, allocatable :: Bmh(:, :)   ! B^(-1/2)
    complex*16, allocatable :: Ap(:, :)
    complex*16, allocatable :: Vap(:, :)
    complex*16, allocatable :: Wb(:)
    complex*16, allocatable :: Vb(:, :)

    integer :: n, i, j, k, l

    n = size(A(:, 1))
    allocate(Bph(n, n))
    allocate(Bmh(n, n))
    allocate(Ap(n, n))
    allocate(Vap(n, n))
    allocate(Wb(n))
    allocate(Vb(n, n))
    

    ! diagonalize B matrix
    ! B.Vb(:, i) = Vb(:, i) Wb(i)
    ! B.Vb = Vb.diag(Wb)
    ! Vb^T B.Vb = Vb^T.Vb.diag(Wb) = diag(Wb)
    ! B^(1/2) = Vb diag(Sqrt(Wb)) Vb^T

    ! A.Va(:,i) = Wa(i) B.Va(:,i)
    ! B^(-1/2) A B^(-1/2) B^(+1/2) Va(:,i) = Wa(i) B^(-1/2) B B^(-1/2) B^(+1/2) Va(:,i)
    ! A' Va'(:,i) = Wa(i) Va'(:,i)
    ! A' = B^(-1/2) A B^(-1/2)
    ! Va'(:,i) = B^(+1/2) Va(:,i)

    call LA_GEEV(B, Wb, Vb)

    do i = 1, n
       do j = 1, n
          Bph(i, j) = (0.0d0, 0.0d0)
          Bmh(i, j) = (0.0d0, 0.0d0)
          do k = 1, n
             Bph(i,j) = Bph(i,j) + Vb(i, k) * sqrt(Wb(k)) * Vb(j, k)
             Bmh(i,j) = Bmh(i,j) + Vb(i, k) * (1.0d0/sqrt(Wb(k))) * Vb(j, k)
          end do
       end do
    end do

    do i = 1, n
       do j = 1, n
          Ap(i, j) = (0.0d0, 0.0d0)
          do k = 1, n
             do l = 1, n
                Ap(i, j) = Ap(i, j) + Bmh(i, k) * A(k, l) * Bmh(l, j)
             end do
          end do
       end do
    end do

    call LA_GEEV(Ap, Wa, Vap)

    do i = 1, n
       do j = 1, n
          Va(i, j) = (0.0d0, 0.0d0)
          do k = 1, n
             Va(i,j) = Va(i,j) + Bmh(i, k) * Vap(k, j)
          end do
       end do
    end do

    deallocate(Bph)
    deallocate(Bmh)
    deallocate(Ap)
    deallocate(Vap)
    deallocate(Wb)
    deallocate(Vb)
    
    
  end subroutine LA_GEEV_GEN
end module Mod_Utils
 
