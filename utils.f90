module Mod_Utils
  implicit none
  real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
  interface cond
     module procedure &
          cond_int, &
          cond_double, &
          cond_complex
  end interface cond
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

  ! ====  ====
end module Mod_Utils
 
