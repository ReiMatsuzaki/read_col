program main
  use Mod_IntIn
  implicit none
  integer, parameter :: ifile = 3
  type(IntIn) int_in
!  open(unit = ifile, file = 'int.in', status='old')
  call IntIn_new_read(int_in)
end program main
