program main
  use Mod_ci
  implicit none
  type(CI) wave_func
!  allocate(wave_func % eig(10))
  !  write(*, *) wave_func % eig(1)
  call CI_new_read_civec(wave_func)
end program main
