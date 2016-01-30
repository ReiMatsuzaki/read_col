program main
  use Mod_IntIn
  use Mod_SymVec
  use Mod_UnitTest
  implicit none
!  integer, parameter :: ifile = 3
  type(IntIn) int_in

  call IntIn_new_read(int_in, 5, .false.)
  call IntIn_show_basis_symmetry(int_in)
contains
  subroutine test_SymVec()
    type(SymVec) vec, vec2
    integer, parameter :: ifile = 13
    complex*16 :: v
    
    call SymVec_new(vec, (/3, 2, 4/))
    call SymVec_set(vec, 1, 1, (1.0d0, 0.0d0))
    call SymVec_set(vec, 1, 2, (1.0d0, 1.0d0))
    call SymVec_set(vec, 1, 3, (1.0d0, 2.0d0))

    call SymVec_set(vec, 2, 1, (2.0d0, 0.0d0))
    call SymVec_set(vec, 2, 2, (2.0d0, 1.0d0))

    call SymVec_set(vec, 3, 1, (3.0d0, 0.0d0))
    call SymVec_set(vec, 3, 2, (3.0d0, 1.0d0))
    call SymVec_set(vec, 3, 3, (3.0d0, 2.0d0))
    call SymVec_set(vec, 3, 4, (3.0d0, 3.0d0))

    open(unit=ifile, file='test_symvec.dat')
    call SymVec_write(vec, ifile)
    close(ifile)
    open(unit=ifile, file='test_symvec.dat')
    call SymVec_new_read(vec2, ifile)
    close(ifile)


    call SymVec_get(vec2, 1, 1, v)
    call CEq("vec2(1,1)", v, (1.0d0, 0.0d0))

    call SymVec_get(vec2, 3, 3, v)
    call CEq("vec2(1,1)", v, (3.0d0, 2.0d0))    
    
  end subroutine test_SymVec

end program main
!subroutine practice_proj_wave(int_in)
!  type(int_in), intent(in) :: int_in
!end subroutine practice_proj_wave
