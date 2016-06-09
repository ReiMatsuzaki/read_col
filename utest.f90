module Mod_UnitTest
  implicit none
  real*8, parameter :: eps_double_prec = 10.0d0**(-14.0d0)
  interface expect_eq
     module procedure &
          expect_eq_int, &
          expect_eq_double, &
          expect_eq_complex, &
          expect_eq_complex_array, &
          expect_eq_complex_matrix
  end interface expect_eq
contains
  ! ==== Utilities ====
  subroutine WriteSuccess(label)
    character(*) label
    print *, "[  PASSED  ] ", label
  end subroutine WriteSuccess
  subroutine run_test(label, test_sub)
    character(*), intent(in) :: label
    interface
       subroutine test_sub()
       end subroutine test_sub
    end interface
    print *, 
    print *, "[----------] ", label
    call test_sub()
    print *, "[----------] "
    
  end subroutine run_test
  subroutine WriteFailed(label)
    character(*) label
    write(*, *) "[  FAILED  ] ", label    
  end subroutine WriteFailed
  ! ==== Boolean ====
  subroutine expect_true(label, a)
    character(*) label
    logical a
    if(a) then
       call WriteSuccess(label)
    else
       call WriteFailed(label)
    end if
  end subroutine Expect_True
  subroutine expect_false(label, a)
    character(*) label
    logical a
    call expect_true(label, .not. a)
  end subroutine Expect_False
  ! ==== Equality ====
  subroutine expect_eq_int(label, a, b)
    character(*) label
    integer a, b
    if(a .eq. b) then
       call WriteSuccess(label)
    else
       call WriteFailed(label)
       write(*, *) "   a = ", a
       write(*, *) "   b = ", b
       call exit(1)
    end if
  end subroutine expect_eq_int
  subroutine expect_eq_double(label, a, b, epsin)
    real*8, intent(in)   :: a, b
    character(*), intent(in) ::  label
    real*8, optional, intent(in) :: epsin
    real*8 :: eps
    if(present(epsin)) then
       eps = epsin
    else
       eps =  eps_double_prec
    end if

    if(abs(a-b) > eps) then
       call WriteFailed(label)
       write(*, *) "a = ", a
       write(*, *) "b = ", b
       write(*, *) "a-b=", a-b
       write(*, *) "|a-b|=", abs(a-b)
    else
       call WriteSuccess(label)
    end if
  end subroutine expect_eq_double
  subroutine expect_eq_complex(label, a, b, epsin)
    complex*16, intent(in)   :: a, b
    character(*), intent(in) :: label
    real*8, optional, intent(in) :: epsin
    real*8 :: eps
    if(present(epsin)) then
       eps = epsin
    else
       eps = eps_double_prec
    end if

    if(abs(a-b) > eps) then
       call WriteFailed(label)
       write(*, *) "a = ", a
       write(*, *) "b = ", b
       write(*, *) "a-b=", a-b
       write(*, *) "|a-b|=", abs(a-b)
    else
       call WriteSuccess(label)
    end if
  end subroutine Expect_Eq_Complex
  subroutine expect_eq_complex_array(label, a, b, epsin)
    character(*), intent(in) :: label
    complex*16, intent(in) :: a(:), b(:)
    real*8, optional, intent(in) :: epsin
    integer :: n
    real*8, allocatable :: c(:)
    integer :: i
    real*8 :: eps
    if(present(epsin)) then
       eps = epsin
    else
       eps = eps_double_prec
    end if


    if(size(a) .ne. size(b)) then
       write(*, *) "CArrayEq. size(a) != size(b)"
       stop
    end if
    
    n = size(a)
    
    allocate(c(n))

    do i = 1, n
       c(i) = abs(a(i) - b(i))
    end do

    if(maxval(c) > eps) then
       call WriteFailed(label)
       write(*, *) "a = ", a
       write(*, *) "b = ", b
       write(*, *) "|a-b|_oo = ", maxval(c)
    else
       call WriteSuccess(label)
    end if

    deallocate(c)
    
  end subroutine Expect_Eq_Complex_Array
  subroutine expect_eq_complex_matrix(label, a, b, epsin)
    character(*), intent(in) :: label
    complex*16, intent(in) :: a(:, :), b(:, :)
    real*8, optional, intent(in) :: epsin
    complex*16 :: c(size(a, 1), size(a, 2))!, absc(size(a, 1), size(a, 2))
    integer :: i, j
    real*8 :: eps
    if(present(epsin)) then
       eps = epsin
    else
       eps = eps_double_prec
    end if
    
    if(all(shape(a) .ne. shape(b))) then
       write(*, *) "CArrayEq. shape(a) != shape(b)"
       stop
    end if
    
    do i = 1, size(a, 1)
       do j = 1, size(b, 1)
          c(i, j) = a(i, j) - b(i, j)
       end do
    end do

    if(maxval(abs(c)) > eps) then
       call WriteFailed(label)
       do i = 1, size(a, 1)
          do j = 1, size(a, 2)
             if(abs(c(i,j))>eps) then
                write(*, *) "i , j = ", i, j
                write(*, *) "a", a(i,j)
                write(*, *) "b", b(i,j)
             end if
          end do
       end do
       write(*, *) "|a-b|_oo = ", maxval(abs(c))
    else
       call WriteSuccess(label)
    end if
    
  end subroutine Expect_Eq_Complex_Matrix

end module Mod_UnitTest
