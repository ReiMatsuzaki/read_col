program main
  use Mod_ci
  use Mod_UnitTest
  implicit none
  type(CI) wave_func
  type(CSF) csf3
  integer, parameter :: ifile_civec = 13
  integer, parameter :: ifile_csf = 21
  integer, parameter :: ifile_ciin = 22
  integer, parameter :: ifile_phoxsec = 23
  integer :: num_mo_sym(8)
  real*8 eps
  data num_mo_sym/8,6,6,2,6,2,2,1/
  eps = 10.0 ** (-10.0)
  
  !
  ! ==== new ====
  !
  open(unit=ifile_csf, file='CSF', status='old', err=999)
  open(unit=ifile_phoxsec, file='phoxsec', status='old', err=999)
  call CI_new_read_ciin(wave_func, 5)
  call CI_new_read_civec(wave_func)
  call CI_new_read_csf(wave_func, ifile_csf)
  rewind ifile_csf
  call CI_set_read_csf(wave_func, ifile_csf)
  call CI_set_read_phoxsec(wave_func, ifile_phoxsec)
  
  !
  ! ==== Unit test ====
  !
  call expect_eq("size(num_mo_sym)", 8, size(wave_func % num_mo_sym))
  call expect_eq("num_mo_sym(1)", 8, wave_func % num_mo_sym(1))
  call expect_eq("num_mo_sym(2)", 6, wave_func % num_mo_sym(2))
  call expect_eq("num_mo_sym(3)", 6, wave_func % num_mo_sym(3))
  call expect_eq("num_mo_sym(4)", 2, wave_func % num_mo_sym(4))
  call expect_eq("mo_offset(3)", 8+6, wave_func % mo_offset(3))
  call expect_eq("mo_offset(4)", 8+6+6, wave_func % mo_offset(4))
  call expect_eq("eig_1", wave_func % eig(1), (-1.61242973860195d0, -0.01813005550149d0), eps)
  
  ! copied from result of bin_civec
  ! 3   5  0.13058947001191537       0.35021669340222594
  call expect_eq("ci_3_5", wave_func % coef(3, 5), &
       (0.13058947001191537d0, 0.35021669340222594d0), eps)

  call expect_eq("d_ene(3)", (0.9164897d0, 0.01810293d0),&
       wave_func % d_ene(3))
  
  call expect_eq("tdm_l(2)", &
       (0.12d0, 0.23d0), &
       wave_func % tdm_l(2))
  call expect_eq("tdm_v(2)", &
       (0.23d0, 0.34d0), &
       wave_func % tdm_l(2))
  
  csf3 = wave_func % csf(3)
  call expect_eq("csf(3).n_sd", csf3 % n_sd, 2)
  call expect_eq("csf(3).n_ele", csf3 % n_ele, 2)
  call Expect_Eq("csf(3).coef(1)", +0.7071067811865475d0 , csf3 % coef(1), eps)
  call Expect_Eq("csf(3).coef(2)", -0.7071067811865475d0 , csf3 % coef(2), eps)
  call expect_eq("csf(3).mo(1,1)",     csf3 % mo(1, 1), 3)     ! 11 (3Pz)
  call expect_eq("csf(3).mo_sym(1,1)", csf3 % mo_sym(1, 1), 2) !
  call expect_eq("csf(3).mo_sym(1,1)", csf3 % spin(1, 1), 1)   ! 
  call expect_eq("csf(3).mo(1,2)", csf3 % mo(1, 2), 1)         ! -1 (1S)
  call expect_eq("csf(3).mo(1,2)", csf3 % mo_sym(1, 2), 1)     ! 
  call expect_eq("csf(3).mo_sym(1,1)", csf3 % spin(1, 2), -1)  ! 
  call expect_eq("csf(3).mo(1,2)", csf3 % mo(2, 1),  3)        ! -11 (3Pz)
  call expect_eq("csf(3).mo(1,2)", csf3 % mo_sym(2, 1), 2)     !
  call expect_eq("csf(3).mo_sym(1,1)", csf3 % spin(2, 1), -1)  !  
  call expect_eq("csf(3).mo(1,2)", csf3 % mo(2, 2), 1)         ! 1  (1S)
  call expect_eq("csf(3).mo(1,2)", csf3 % mo_sym(2, 2), 1)     !
  call expect_eq("csf(3).mo_sym(1,1)", csf3 % spin(2, 2), 1) 

  
  STOP
999 write(*, *) "failed to find CSF file"
end program main
