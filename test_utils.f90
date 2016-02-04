program main
  use Mod_Utils
  use Mod_LinearAlgebra
  use Mod_UnitTest
  call expect_eq("Pi", 3.1415d0, pi, 0.0001d0)
  call run_test("sort", test_sort_re)
  call run_test("eigen_sym", test_eigen_sym)
  call run_test("mat_half_inv", test_mat_half_inv)
  call run_test("eigen_sym_gen", test_eigen_sym_gen)
  call run_test("cond", test_cond)
contains
  subroutine test_sort_re()
    complex*16 :: w(3) = (/(0.1d0, 0.3d0), (-0.45d0, 1.1d0), (2.1d0, 0.3d0)/)
    complex*16 :: w0(3)
    complex*16 :: V(3,3), V0(3,3)

    V(1,:)=(/(1.0d0, 1.0d0), (1.0d0, 2.0d0), (1.0d0, 3.0d0)/)
    V(2,:)=(/(2.0d0, 1.0d0), (2.0d0, 2.0d0), (2.0d0, 3.0d0)/)
    V(3,:)=(/(3.0d0, 1.0d0), (3.0d0, 2.0d0), (3.0d0, 3.0d0)/)

    w0 = w
    V0 = V
    call sort_re(w, V)
    
    call expect_eq("w(1)", w0(2), w(1))
    call expect_eq("w(2)", w0(1), w(2))
    call expect_eq("w(3)", w0(3), w(3))

    call expect_eq("V(:,1)", V0(:,2), V(:,1))
    call expect_eq("V(:,2)", V0(:,1), V(:,2))
    call expect_eq("V(:,3)", V0(:,3), V(:,3))
    
  end subroutine test_sort_re
  subroutine test_eigen_sym()
    integer,parameter :: n=3
    complex*16 :: A(n,n), AA(n, n), one(n, n)
    complex*16:: w(n), v(n, n)

    A(:, :) = (0.0d0, 0.0d0)
    A(1, 1) = (1.0d0, 0.4d0)
    A(1, 2) = (2.0d0, 0.2d0)
    A(2, 1) = (2.0d0, 0.2d0)
    A(2, 3) = (2.0d0, -0.2d0)
    A(3, 2) = (2.0d0, -0.2d0)    
    A(2 ,2) = (1.0d0, 0.0d0)
    A(3, 3) = (1.1d0, 0.4d0)
    AA(:, :) = A(:, :)

    call eigen_sym(AA, w, v)

    call expect_eq("A=AT", A, transpose(A))
    
    one = id_mat(n)
    call expect_eq("Av=v.diag(w)", matmul(A, v), matmul(v, diag_mat(w)))
    call expect_eq("vT.v=1", one, matmul(transpose(v), v))
    call expect_eq("v.vT=1", one, matmul(v, transpose(v)))
    
  end subroutine test_eigen_sym
  subroutine test_mat_half_inv()
    integer,parameter :: n=3
    complex*16 :: B(n,n), BB(n, n), Bmh(n, n)

    B(1, 1) = (2.0d0, 0.2d0) 
    B(1, 2) = (-2.0d0, 0.0d0) 
    B(2, 1) = (-2.0d0, 0.0d0) 
    B(2 ,2) = (3.0d0, 0.0d0) 
    B(3, 3) = (2.1d0, -0.4d0) 
    BB(:, :) = B(:, :)

    Bmh =  mat_half_inv(BB)

    call expect_eq("B^(-1/2)", id_mat(n), matmul(B, matmul(Bmh, Bmh)))
    
  end subroutine test_mat_half_inv
  subroutine test_eigen_sym_gen()
    integer,parameter :: n=3
    complex*16 :: A(n,n), AA(n, n), B(n,n), BB(n,n)
    complex*16:: w(n), vs(n, n)
    
    A(:, :) = (0.0d0, 0.0d0); B(:, :) = (0.0d0, 0.0d0)
    A(1, 1) = (1.0d0, 0.0d0); B(1, 1) = (2.0d0, 0.2d0) 
    A(1, 2) = (2.0d0, 0.4d0); B(1, 2) = (-2.0d0, 0.0d0) 
    A(2, 1) = (2.0d0, 0.4d0); B(2, 1) = (-2.0d0, 0.0d0) 
    A(2 ,2) = (1.0d0, 0.0d0); B(2 ,2) = (3.0d0, 0.0d0) 
    A(3, 3) = (1.1d0, 1.0d0); B(3, 3) = (2.1d0, -0.4d0) 
    AA(:, :) = A(:, :);       BB(:, :) = B(:, :)
    
    call eigen_sym(AA, BB, w, vs)

    call expect_eq("eigen_sym_gen", &
         matmul(A, vs), &
         matmul(matmul(B, vs), diag_mat(w)))
    call expect_eq("V^BT V=1", id_mat(n), matmul3(transpose(vs), B, vs))

    call expect_true("Re part order", all((/ (real(w(i)) < real(w(i+1)), i=1,n-1 )/)))
    call expect_eq("w(1)", w(3), (6.43614712d0, -0.8792356d0), 10.0d0**(-9.0))
    
  end subroutine test_eigen_sym_gen
  subroutine test_cond()
    call expect_eq("cond_int", 3, cond(.true., 3, 1))
    call expect_eq("cond_int", 1, cond(.false., 3, 1))
    call expect_eq("cond_double", 1.1d0, cond(.true., 1.1d0, 1.2d0))
    call expect_eq("cond_complex", (1.1d0, 0.1d0), &
         cond(.false., (1.1d0, 0.0d0), (1.1d0, 0.1d0)))
  end subroutine test_cond
end program main
