module Mod_LinearAlgebra
  implicit none
  interface eigen_sym
     module procedure eigen_sym_id, eigen_sym_gen
  end interface eigen_sym
contains
  ! ==== basic linear algebra ====
  function matmul3(A, B, C)
    complex*16, intent(in) :: A(:, :), B(:, :), C(:, :)
    complex*16 :: matmul3(size(A, 1), size(C, 2))
    complex*16 :: BC(size(B, 1), size(C, 2))
    integer i, j, k
    
    if(size(A, 2) .ne. size(B, 1) .or. size(B, 2) .ne. size(C, 1)) then
       write(*, *) "size mismatch in matmul3"
       stop
    end if

    BC(:, :) = (0.0d0, 0.0d0)
    do i = 1, size(B, 1)
       do j = 1, size(C, 2)
          BC(i, j) = (0.0d0, 0.0d0)
          do k = 1, size(C, 1)
             BC(i, j) = BC(i, j)  + B(i, k) * C(k, j)
          end do
       end do
    end do

    matmul3(:, :) = (0.0d0, 0.0d0)
    do i = 1, size(A, 1)
       do j = 1, size(C, 2)
          matmul3(i, j) = (0.0d0, 0.0d0)
          do k = 1, size(A, 2)
             matmul3(i, j) = matmul3(i, j) + A(i, k) * BC(k, j)
          end do
       end do
    end do
    
  end function matmul3
  function conjg_transpose(A)
    complex*16, intent(in) :: A(:, :)
    complex*16 :: conjg_transpose(size(A, 2), size(A, 1))
    integer i, j
    do i = 1, size(A, 1)
       do j = 1, size(A, 2)
          conjg_transpose(j, i) = conjg(A(i, j))
       end do
    end do
  end function conjg_transpose
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
  ! ==== Eigen values ====
  subroutine sort_re(W, V)
    complex*16, intent(inout) :: W(:)
    complex*16, intent(inout) :: V(:, :)
    integer n, i, j
    complex*16 tmp
    complex*16 tmpV(size(V,1))
    n = size(W)

    do i = 1, n-1
       do j = n, i+1, -1
          if(real(W(j)) < real(W(j-1))) then
             tmp = W(j)
             W(j) = W(j-1)
             W(j-1) = tmp

             tmpV = V(:,j)
             V(:,j) = V(:,i)
             V(:,i) = tmpV
          end if
       end do
    end do
    
  end subroutine sort_re
  subroutine eigen_sym_id(A, W, V)
    use f95_lapack, only : LA_GEEV
    complex*16, intent(inout) :: A(:, :)
    complex*16, intent(out)   :: W(size(A, 1))
    complex*16, intent(out)   :: V(size(A, 1), size(A, 1))
    complex*16  :: Vl(size(A, 1), size(A, 1))
    complex*16  :: Vr(size(A, 1), size(A, 1))
    complex*16 :: VV(size(A, 1), size(A, 1))
    integer j, n
    
    n = size(A, 1)
    call LA_GEEV(A, W, Vl, Vr)
    VV = matmul(transpose(Vr),Vr)
    
    do j = 1, n
       V(:, j) = Vr(:, j) / sqrt(VV(j, j))
    end do
        
  end subroutine eigen_sym_id
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
  subroutine eigen_sym_gen(A, B, W, V)
    
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
    call eigen_sym(Ap, W, Vp)
    V = matmul(Bph, Vp)

    call sort_re(W, V)
    
  end subroutine eigen_sym_gen
  !  function diag_mat_mult(V, A)
  !    complex*16, intent(in) :: V(:)
  !    complex*16 :: diag_mat_mult(size(V), size(V))    
  !    complex*16, intent(in) :: A(:, :)
  !    integer n, i, j
  !    
  !    n = size(V)
  !    if(n .ne. size(A, 1) .or. n .ne. size(A, 2)) then
  !       write(*, *) "diag_mat_mul. size mismath"
  !       stop
  !    end if
  !
  !    ! {diag(V) B}_ij = sum_k delta(i,k) V(i) * B(k, j) = V(i)B(i,j)
  !    
  !    do i = 1, n
  !       diag_mat_mult(i, :) = (/ (V(i) * A(i, j), j = 1, n) /)
  !    end do
  !    
  !  end function diag_mat_mult
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
    
    call eigen_sym(A, W, Vr)

    y = (/ (1.0d0/sqrt(W(k)), k = 1, n) /)
    
    mat_half_inv = matmul(Vr, matmul(diag_mat(y), transpose(Vr)))
    
  end function mat_half_inv
end module Mod_LinearAlgebra
