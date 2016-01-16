program main
  use Mod_AoInts
  implicit none
      
  integer, parameter ::  ifile = 3
  integer :: ist, is, iso
  type(AoInts) ao_ints
  open(unit=ifile, file='AOINTS', status='old', form='unformatted')
  call AoInts_new(ao_ints, ifile)
  close(unit=ifile)
end program main

