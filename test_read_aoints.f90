program main
  use Mod_AoInts
  implicit none
  type(AoInts) ao_ints
  call AoInts_new_read(ao_ints)
  call AoInts_show(ao_ints)
  call AoInts_delete(ao_ints)
end program main

