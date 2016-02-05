program main
  use Mod_AoInts
  use Mod_UnitTest  
  implicit none

  call run_test("AOINTS", test_AoInts)
  !call run_test("DIPINTS", test_DIPINTS)
  
contains
  subroutine test_AoInts()
    type(AoInts) :: ao_ints
    integer, parameter :: ifile = 13
    open(unit=ifile, file='AOINTS', status='old', form='unformatted')
    call AoInts_new_read(ao_ints, ifile)
    write(*, *) "Show>>>"
    call AoInts_show(ao_ints)
    write(*, *) "Show<<<"
    call AoInts_delete(ao_ints)
    close(ifile)
  end subroutine test_AoInts
  subroutine test_DIPINTS()
    type(AoInts) :: ao_ints
    integer, parameter :: ifile = 13
    open(unit=ifile, file='DIPINTS', status='old', form='unformatted')
    call AoInts_new_read(ao_ints, ifile)
    call AoInts_show(ao_ints)
    call AoInts_delete(ao_ints)
    close(ifile)
  end subroutine test_DIPINTS
 
end program main

