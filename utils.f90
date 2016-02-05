module Mod_Utils
  implicit none
  
  real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
  
  interface cond
     module procedure &
          cond_int, &
          cond_double, &
          cond_complex
  end interface cond

  interface c_random
     module procedure &
          c_random_scalar, &
          c_random_vec, &
          c_random_mat
  end interface c_random
contains
  ! ==== condition operator ====
  function cond_int(q, a, b)
    logical, intent(in) :: q
    integer, intent(in) :: a, b
    integer :: cond_int
    if(q) then
       cond_int = a
    else
       cond_int = b
    end if
  end function cond_int
  function cond_double(q, a, b)
    logical, intent(in) :: q
    real*8, intent(in) :: a, b
    real*8 :: cond_double
    if(q) then
       cond_double = a
    else
       cond_double = b
    end if
  end function cond_double
  function cond_complex(q, a, b)
    logical, intent(in) :: q
    complex*16, intent(in) :: a, b
    complex*16 :: cond_complex
    if(q) then
       cond_complex = a
    else
       cond_complex = b
    end if
  end function cond_complex

  ! ==== complex random ====
  subroutine c_random_scalar(x)
    complex*16, intent(out) :: x
    real*8 :: a, b
    call random_number(a)
    call random_number(b)
    x = a + (0.0d0, 1.0d0) * b
  end subroutine c_random_scalar
  subroutine c_random_vec(x)
    complex*16, intent(out) :: x(:)
    integer i
    do i = 1, size(x)
       call c_random_scalar(x(i))
    end do
  end subroutine c_random_vec
  subroutine c_random_mat(x)
    complex*16, intent(out) :: x(:, :)
    integer i, j
    do i = 1, size(x, 1)
       do j = 1, size(x, 2)
          call c_random_scalar(x(i, j))
       end do
    end do
    
  end subroutine c_random_mat
end module Mod_Utils
 
