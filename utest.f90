module Mod_UnitTest
  implicit none
contains
  subroutine WriteSuccess(label)
    character(*) label
    write(*, *) "Success: ", label
  end subroutine WriteSuccess
  subroutine WriteFailed(label)
    character(*) label
    write(*, *) "FAILED: ", label    
  end subroutine WriteFailed
  subroutine BTrue(label, a)
    character(*) label
    logical a
    if(a) then
       write(*, *) "Success: ", label
    else
       write(*, *) "FAIED: ", label
    end if
  end subroutine BTrue
  subroutine BFalse(label, a)
    character(*) label
    logical a
    call BTrue(label, .not. a)
  end subroutine BFalse
  subroutine IEq(label, a, b)
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
  end subroutine IEq
  subroutine DNear(label, a, b, eps)
    real*8 a, b, eps
    character(*) label
    if(abs(a-b) > eps) then
       write(*, *) "Failed: ", label
       write(*, *) "   a = ", a
       write(*, *) "   b = ", b
       call exit(1)
       call exit(1)
    else
       write(*, *) "Success: ", label
    end if
  end subroutine DNear
  subroutine CNear(label, a, b, eps)
    complex*16 a, b
    character(*) label
    real*8 eps
    if(abs(a-b) > eps) then
       call WriteFailed(label)
       write(*, *) "   a = ", a
       write(*, *) "   b = ", b
       write(*, *) "   a-b = ", a-b
       write(*, *) "  |a-b|= ", abs(a-b)
       write(*, *) "   eps = ", eps
       call exit(1)
    else
       call WriteSuccess(label)
    end if
  end subroutine CNear
  subroutine CEq(label, a, b)
    complex*16 a, b
    character(*) label
    real*8 eps
    eps = 10.0 ** (-14.0)
    if(abs(a-b) > eps) then
       call WriteFailed(label)
       write(*, *) "a = ", a
       write(*, *) "b = ", b
       write(*, *) "a-b=", a-b
       write(*, *) "|a-b|=", abs(a-b)
    else
       call WriteSuccess(label)
    end if
  end subroutine CEq
  subroutine CArrayEq(label, a, b)
    character(*), intent(in) :: label
    complex*16, intent(in) :: a(:), b(:)
    real*8, parameter      :: eps = 10.0 ** (-14.0)
    integer :: n
    real*8, allocatable :: c(:)
    integer :: i

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
    
  end subroutine CArrayEq
  subroutine CMatEq(label, a, b)
    character(*), intent(in) :: label
    complex*16, intent(in) :: a(:, :), b(:, :)
    real*8 :: c(size(a, 1), size(a, 2))
    real*8, parameter      :: eps = 10.0 ** (-14.0)
    integer :: i, j

    if(all(shape(a) .ne. shape(b))) then
       write(*, *) "CArrayEq. shape(a) != shape(b)"
       stop
    end if
    
    do i = 1, size(a, 1)
       do j = 1, size(b, 1)
          c(i, j) = abs(a(i, j) - b(i, j))
       end do
    end do

    if(maxval(c) > eps) then
       call WriteFailed(label)
       write(*, *) "a = ", a
       write(*, *) "b = ", b
       write(*, *) "|a-b|_oo = ", maxval(c)
    else
       call WriteSuccess(label)
    end if
    
  end subroutine CMatEq
end module Mod_UnitTest
