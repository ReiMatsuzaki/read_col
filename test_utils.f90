program main
  use Mod_Utils
  use Mod_UnitTest
  call DNear("Pi", 3.1415d0, pi, 0.0001d0)
  !call test_la_geev()
  call test_la_geev_gen()
contains
  subroutine test_la_geev()
    use f95_lapack, only : LA_GEEV
    integer,parameter :: n=3
    integer i, j, k, l
    complex*16 :: A(n,n), AA(n, n)
    complex*16:: w(n), vs(n, n)
    complex*16 :: x = (0.0d0, 0.0d0)
    A(1, 1) = (1.0d0, 0.0d0)
    A(1, 2) = (2.0d0, 0.0d0)
    A(2, 1) = (2.0d0, 0.0d0)
    A(2 ,2) = (1.0d0, 0.0d0)
    A(3, 3) = (1.1d0, 0.0d0)
    AA(:, :) = A(:, :)

    call LA_GEEV(AA, w, vs)

    write(*, *) "vs^T A vs"
    do i = 1, n
       do j = 1, n
          x = (0.0d0, 0.0d0)
          do k = 1, n
             do l = 1, n
                x = x + vs(k, i) * A(k, l) * vs(l, j)
             end do
          end do
          write(*, *) i, j, x
       end do
    end do
    
    write(*, *) "w"
    do i = 1, n
       write(*, *) i, w(i)
    end do
    
  end subroutine test_la_geev
  subroutine test_la_geev_gen
    integer,parameter :: n=3
    integer i, j, k
    complex*16 :: A(n,n), AA(n, n), B(n,n), BB(n,n)
    complex*16:: w(n), vs(n, n)
    complex*16 :: x, y
    A(1, 1) = (1.0d0, 0.0d0); B(1, 1) = (2.0d0, 0.0d0) 
    A(1, 2) = (2.0d0, 0.0d0); B(1, 2) = (-2.0d0, 0.0d0) 
    A(2, 1) = (2.0d0, 0.0d0); B(2, 1) = (-2.0d0, 0.0d0) 
    A(2 ,2) = (1.0d0, 0.0d0); B(2 ,2) = (3.0d0, 0.0d0) 
    A(3, 3) = (1.1d0, 0.0d0); B(3, 3) = (2.1d0, 0.0d0) 
    AA(:, :) = A(:, :);       BB(:, :) = B(:, :)

    call LA_GEEV_GEN(AA, BB, w, vs)

    write(*, *) "A vs"
    do j = 1, n
       write(*, *) "j:", j
       write(*, *) "A vs(:, j)"
       do i = 1, n
          x = (0.0d0, 0.0d0)
          y = (0.0d0, 0.0d0)
          do k = 1, n
             x = x + A(i, k) * vs(k, j)
             y = y + B(i, k) * vs(k, j) * w(j)
          end do
          call CEq("la_geev_gen", x, y)
       end do
    end do
    
  end subroutine test_la_geev_gen
end program main
