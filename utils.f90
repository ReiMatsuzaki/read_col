module Mod_Utils
  
  implicit none
  real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
contains
  function conj_transpose(A)
    complex*16, intent(in) :: A(:, :)
    complex*16 :: conj_transpose(size(A, 2), size(A, 1))
    integer i, j
    do i = 1, size(A, 1)
       do j = 1, size(A, 2)
          conj_transpose(j, i) = conjg(A(i, j))
       end do
    end do
  end function conj_transpose
  subroutine solve_eigen_sym(A, W, V)
    use f95_lapack, only : LA_GEEV
    complex*16, intent(inout) :: A(:, :)
    complex*16, intent(out)   :: W(size(A, 1))
    complex*16, intent(out)   :: V(size(A, 1), size(A, 1))
    complex*16  :: Vl(size(A, 1), size(A, 1))
    complex*16  :: Vr(size(A, 1), size(A, 1))
    complex*16 :: VV(size(A, 1), size(A, 1))
    integer i, j, n
    
    n = size(A, 1)
    call LA_GEEV(A, W, Vl, Vr)
    
    VV = matmul(transpose(Vr), Vr)

    do i = 1, n
       do j = 1, n
          V(i, j) = Vr(i, j) / sqrt(VV(j, j))
       end do
    end do
        
  end subroutine solve_eigen_sym
  subroutine solve_eigen_sym_gen(A, B, W, V)
    use f95_lapack, only : LA_GEEV
    complex*16, intent(inout) :: A(:, :)
    complex*16, intent(inout) :: B(:, :)
    complex*16, intent(out)   :: W(size(A, 1))
    complex*16, intent(out)   :: V(size(A, 1), size(A, 1))

    complex*16 :: Ap(size(A, 1), size(A, 1))
    complex*16 :: Bph(size(A, 1), size(A, 1))
    complex*16 :: Vp(size(A, 1), size(A, 1))

    Bph = mat_half_inv(B)
    Ap = matmul(matmul(Bph, A), Bph)
    call solve_eigen_sym(Ap, W, Vp)
    V = matmul(Bph, Vp)
    
  end subroutine solve_eigen_sym_gen
  function id_mat(n)
    integer, intent(in) :: n
    complex*16 :: id_mat(n, n)
    integer i, j
    do i = 1, n
       do j = 1, n
          id_mat(i, j) = (0.0d0, 0.0d0)
       end do
       id_mat(i, i) = (1.0d0, 0.0d0)
    end do
  end function id_mat
  function diag_mat(v)
    complex*16, intent(in) :: v(:)
    complex*16 :: diag_mat(size(v), size(v))
    integer    :: i, j
    do i = 1, size(v)
       do j = 1, size(v)
          diag_mat(i, j) = (0.0d0, 0.0d0)
       end do
       diag_mat(i, i) = v(i)
    end do
  end function diag_mat
  function diag_mat_mult(V, A)
    complex*16, intent(in) :: V(:)
    complex*16 :: diag_mat_mult(size(V), size(V))    
    complex*16, intent(in) :: A(:, :)
    integer n, i, j
    
    n = size(V)
    if(n .ne. size(A, 1) .or. n .ne. size(A, 2)) then
       write(*, *) "diag_mat_mul. size mismath"
       stop
    end if

    ! {diag(V) B}_ij = sum_k delta(i,k) V(i) * B(k, j) = V(i)B(i,j)
    
    do i = 1, n
       diag_mat_mult(i, :) = (/ (V(i) * A(i, j), j = 1, n) /)
    end do
    
  end function diag_mat_mult

  ! B Vr(:,i) = Vr(:,i) W(i)
  ! VrT(i,:) B = W(i) VrT(i,:)
  ! VlH(i, :) B = W(i) VlH(i, :)
  !
  ! 0 = {W(i) - W(j)} VlH(i, :). Vr(:, j)
  !
  ! B = Vr W VlH^T = (Vr) W (VlH)
  function mat_half_inv(A)
    use f95_lapack, only : LA_GEEV
    complex*16, intent(inout) :: A(:, :)
    complex*16 :: mat_half_inv(size(A, 1), size(A, 1))
    complex*16 :: W(size(A, 1))
    complex*16 :: Vr(size(A, 1), size(A, 1))
    complex*16 :: y(size(A, 1))
    integer :: k, n

    n = size(A, 1)
    
    call solve_eigen_sym(A, W, Vr)

    y = (/ (1.0d0/sqrt(W(k)), k = 1, n) /)
    
    mat_half_inv = matmul(Vr, diag_mat_mult(y, transpose(Vr)))
    
  end function mat_half_inv
  subroutine LA_GEEV_GEN(A, B, Wa, Va)
    use f95_lapack, only : LA_GEEV
    ! Solve generalized eigen values problem
    !    A UL = diag(W) B UL
    !
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
    complex*16, intent(inout) :: A(:, :)
    complex*16, intent(inout) :: B(:, :)
    complex*16, intent(out) :: Wa(:)
    complex*16, intent(out) :: Va(:, :)

    complex*16, allocatable :: Bmh(:, :)   ! B^(-1/2)
    complex*16, allocatable :: Ap(:, :)
    complex*16, allocatable :: Vap(:, :)

    integer :: n

    n = size(A(:, 1))
    allocate(Bmh(n, n))
    allocate(Ap(n, n))
    allocate(Vap(n, n))

    Bmh(:, :) = mat_half_inv(B)
    
    Ap = matmul(matmul(Bmh, A), Bmh)

    call LA_GEEV(Ap, Wa, Vap)

    Va = matmul(Bmh, Vap)

    deallocate(Bmh)
    deallocate(Ap)
    deallocate(Vap)
    
  end subroutine LA_GEEV_GEN
end module Mod_Utils
 
