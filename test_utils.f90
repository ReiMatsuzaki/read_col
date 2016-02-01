program main
  use Mod_Utils
  use Mod_UnitTest
  call DNear("Pi", 3.1415d0, pi, 0.0001d0)
  call test_la_geev()
  call test_mat_half_inv()
  call test_eigen_sym_gen()
contains
  subroutine test_la_geev()
    use f95_lapack
    integer,parameter :: n=3
    complex*16 :: A(n,n), AA(n, n), one(n, n)
    complex*16:: w(n), v(n, n)

    A(1, 1) = (1.0d0, 0.4d0)
    A(1, 2) = (2.0d0, 0.2d0)
    A(2, 1) = (2.0d0, 0.2d0)
    A(2 ,2) = (1.0d0, 0.0d0)
    A(3, 3) = (1.1d0, 0.4d0)
    AA(:, :) = A(:, :)

    write(*, *)
    write(*, *) "solve_eigen_sym"
    call solve_eigen_sym(AA, w, v)

    one = id_mat(n)
    call CMatEq("vT.v=1", one, matmul(transpose(v), v))
    call CMatEq("v.vT=1", one, matmul(v, transpose(v)))
    call CMatEq("Av=v.diag(w)", matmul(A, v), matmul(v, diag_mat(w)))
    
  end subroutine test_la_geev
  subroutine test_mat_half_inv()
    integer,parameter :: n=3
    complex*16 :: B(n,n), BB(n, n), Bmh(n, n)

    write(*, *)
    write(*, *) "mat_half_inv"    
    
    B(1, 1) = (2.0d0, 0.2d0) 
    B(1, 2) = (-2.0d0, 0.0d0) 
    B(2, 1) = (-2.0d0, 0.0d0) 
    B(2 ,2) = (3.0d0, 0.0d0) 
    B(3, 3) = (2.1d0, -0.4d0) 
    BB(:, :) = B(:, :)

    Bmh =  mat_half_inv(BB)

    call CMatEq("B^(-1/2)", id_mat(n), matmul(B, matmul(Bmh, Bmh)))

!    write(*, *) matmul(B, matmul(Bmh, Bmh))

!    write(*, *)
!    write(*, *) matmul(Bmh, matmul(B, Bmh))
    
  end subroutine test_mat_half_inv
  subroutine test_eigen_sym_gen
    integer,parameter :: n=3
    complex*16 :: A(n,n), AA(n, n), B(n,n), BB(n,n)
    complex*16:: w(n), vs(n, n)
    A(1, 1) = (1.0d0, 0.0d0); B(1, 1) = (2.0d0, 0.2d0) 
    A(1, 2) = (2.0d0, 0.4d0); B(1, 2) = (-2.0d0, 0.0d0) 
    A(2, 1) = (2.0d0, 0.4d0); B(2, 1) = (-2.0d0, 0.0d0) 
    A(2 ,2) = (1.0d0, 0.0d0); B(2 ,2) = (3.0d0, 0.0d0) 
    A(3, 3) = (1.1d0, 1.0d0); B(3, 3) = (2.1d0, -0.4d0) 
    AA(:, :) = A(:, :);       BB(:, :) = B(:, :)

    write(*, *) 
    write(*, *) "test_eigen_sym_gen"
    
    call solve_eigen_sym_gen(AA, BB, w, vs)

    call CNEar("w(1)", w(1), (6.43614712d0, -0.8792356d0), 10.0d0**(-9.0))

    call CMatEq("eigen_sym_gen", &
         matmul(A, vs), &
         matmul(matmul(B, vs), diag_mat(w)))
    
  end subroutine test_eigen_sym_gen
end program main
