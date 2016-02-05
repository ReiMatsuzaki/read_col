program main
  use Mod_MOCoef
  use Mod_AoInts
  use Mod_UnitTest
  implicit none

  type(AoInts) ao
  type(MOCoef) mo
  integer, allocatable :: n_ist(:)
  integer i
  integer, parameter :: ifile = 13
  real*8 :: eps = 0.000002d0

  open(unit = ifile, file = "AOINTS", form='unformatted')
  call AoInts_new_read(ao, ifile)
  close(unit = ifile)

  allocate(n_ist(ao % nst))
  do i = 1, ao % nst
     n_ist(i) = ao % nso(i)
  end do
  
  open(unit = ifile, file='MOCOEF')
  call MOCoef_new(mo, n_ist)
  call MOCoef_set_read(mo, ifile)
  close(unit = ifile)

  call expect_eq("MO(1,1)(1,1), ", (0.702767d0, 0.0d0), &
       BlockMat_val(mo % mo_coef, 1, 1, 1, 1), eps)
  call expect_eq("MO(1,1)(3,1), ", (0.837580d0, 0.0d0), &
       BlockMat_val(mo % mo_coef, 1, 1, 3, 1), eps)
  call expect_eq("MO(1,1)(3,2), ", (-3.679934d0, 0.0d0), &
       BlockMat_val(mo % mo_coef, 1, 1, 3, 2), eps)
  call expect_eq("MO(3,3)(2,3), ", (-2.623901d0, -1.318718d0), &
       BlockMat_val(mo % mo_coef, 3, 3, 2, 3), eps)
  call expect_eq("Eig(1)(2)", (0.385171d0, 0.000000d0), &
       BlockVec_val(mo % mo_eig, 1, 2), eps)
  call expect_eq("Eig(3)(3)", (5.750098d0, -0.823218d0), &
       BlockVec_val(mo % mo_eig, 3, 3), eps)
  write(*, *) 

  print *, "finalize"
  call AoInts_delete(ao)
  call MOCoef_delete(mo)
  
end program main
