program main
  use Mod_AoInts
  implicit none
      
  integer, parameter ::  ifile = 3
  integer :: ist, is, iso
  type(AoInts) ao_ints
  open(unit=ifile, file='AOINTS', status='old', form='unformatted')
  call AoInts_new_read(ao_ints, ifile)
  !  call AoInts_show(ao_ints)
  !  call SymMat_show(ao_ints % s_mat)
  ! call SymERI_show(ao_ints % eri)
  close(unit=ifile)
  call AoInts_delete(ao_ints)
end program main

