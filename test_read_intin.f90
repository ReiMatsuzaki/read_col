program main
  use Mod_IntIn
  implicit none
  integer, parameter :: ifile = 3
  type(IntIn) int_in

  call IntIn_new_read(int_in)
  call IntIn_show_basis_symmetry(int_in)

end program main
!subroutine practice_proj_wave(int_in)
!  type(int_in), intent(in) :: int_in
!end subroutine practice_proj_wave
