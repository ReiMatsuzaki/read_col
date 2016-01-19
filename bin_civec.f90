! This program is copied from cCOLUMBUS/tools/bin_civec/
! 2015/1/20

program main
  implicit none
  integer, parameter :: nfl12 = 12
  integer            :: ncidim, ndum, i, nstate, n, n_arg
  character          :: blabel(1), name1(10)
  complex*16         :: eig
  complex*16, allocatable :: coef(:)
  character          :: file_name*50, c_nci*10

  n_arg = iargc()

  if(n_arg .eq. 0) then
     call show_usage()
     stop
  end if

  call getarg(1, file_name)

  open(unit=nfl12,file=file_name,status='old',form='unformatted')
  READ(nfl12) ncidim,blabel,name1, ndum
  write(*,*) "# n cidim= ", ncidim
  write(*,*) "# blabel= ", blabel
  write(*,*) "# name1= ", name1
  write(*,*) "# ndum= ", ndum

  allocate(coef(ncidim))
  
  if (iargc() .gt. 1) then
     call getarg(2, c_nci)
     read (c_nci,*) nstate
  else
     nstate = ncidim
  end if

  write (*,*) "# n state= ", nstate
  write (*,*) "# n, Re[E_n], Im[E_n] ; eigen energy"
  write (*,*) "# n, i, Re[C_ni], Re[C_ni] ; coefficient"

  do n = 1, nstate
     READ(nfl12) eig, coef(:)
     write (*,*) n, real(eig), aimag(eig)
     do i = 1, ncidim
        write (*,*) n, i, real(coef(i)), aimag(coef(i))
     end do
  end do
  
  close(nfl12)

end program main

subroutine show_usage()
  write(*,*) "Usage : bin_civec file_name"
  write(*,*) "  or  : bin_civec file_name #_of_state"
  write(*,*) ""
  write(*,*) "This program converts a binary file of CI vector to a text file."
  write(*,*) "The binary file is produced by cCOLUMBUS and usually called CIVEC."
end subroutine show_usage

