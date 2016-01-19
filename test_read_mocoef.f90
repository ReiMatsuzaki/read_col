program main
  use Mod_MOCoef
  use Mod_AoInts
  implicit none

  type(AoInts) ao
  type(MOCoef) mo
  integer, allocatable :: n_ist(:)
  integer i
  
  call AoInts_new_read(ao)

  allocate(n_ist(ao % nst))
  do i = 1, ao % nst
     n_ist(i) = ao % nso(i)
  end do

  print *, "mo coef new"
  call MOCoef_new(mo, n_ist)
  print *, "mo coef set read"
  call MOCoef_set_read(mo)

  print *, "finalize"
  call AoInts_delete(ao)
  call MOCoef_delete(mo)
  
end program main
