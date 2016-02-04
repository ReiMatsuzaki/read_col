! compute projection of N-1 electron wave function to
! each N electron determinant to obtain one particle orbital

module Mod_ProjOne
  use Mod_CI
  use Mod_AoInts
  use Mod_MOCoef
  implicit none
  type ProjOne
     type(AoInts) :: ao
     type(MOCoef) :: mo
     type(CI)     :: ci_wf
  end type ProjOne
contains
  subroutine new_read(this)
    type(ProjOne) :: this
    integer, parameter :: ifile_in = 13
    integer, parameter :: ifile_aoints = 14
    integer, parameter :: ifile_mocoef = 15
    integer, parameter :: ifile_csf = 16
    integer, parameter :: ifile_ciin = 17
    character(32)      :: aoints_path  = "no_file"
    character(32)      :: mocoef_path = "no_file"
    character(32)      :: csf_path  = "no_file"
    character(32)      :: ciin_path = "no_file"
    character(32)      :: in_path
    integer n_ist(8)

    namelist/proj_one/aoints_path, mocoef_path, csf_path, ciin_path

    if(iargc() .ne. 1) then
       call print_usage
    end if
    
    call getarg(1, in_path)
    open(unit = ifile_in, file = in_path, status='old')
    read(unit = ifile_in, nml = proj_one)
    close(unit = ifile_in)

    write(*, *) "aoints_path: ", aoints_path
    write(*, *) "mocoef_path: ", mocoef_path
    write(*, *) "ciin_path: ", ciin_path
    write(*, *) "csf_path: ", csf_path

    open(unit = ifile_aoints, file = aoints_path, status = 'old', form='unformatted')
    call AoInts_new_read(this % ao, ifile_aoints)
    close(unit = ifile_aoints)

    n_ist(:) = this % ao % nso(1 : this % ao % nst)
    call MOCoef_new(this % mo, n_ist)
!    open(unit = ifile_mocoef, file = mocoef_path, status = 'old')
    call MOCoef_set_read(this % mo)
 !   close(unit = ifile_aoints)    

    open(unit = ifile_ciin, file = ciin_path, status = 'old')
    call CI_new_read_ciin(this % ci_wf, ifile_ciin)
    call CI_new_read_civec(this % ci_wf)
    close(unit = ifile_ciin)
  
    open(unit = ifile_csf, file = csf_path, status = 'old')
    call CI_new_read_csf(this % ci_wf, ifile_csf)
    close(unit = ifile_csf)
    
  end subroutine new_read
  subroutine print_usage()
    write(*, *) "proj_one INPUT"
    write(*, *) "------------------"
    write(*, *) "&proj_one"
    write(*, *) '"csf_path="CSF"'
    write(*, *) '"ciin_path="ci2.in"""'
    write(*, *) '"num_e=2' 
    write(*, *) "/"
    write(*, *) "------------------"
  end subroutine print_usage  
end module Mod_ProjOne

program main
  use Mod_ProjOne
  implicit none

  type(ProjOne) this
  call new_read(this)

end program main
