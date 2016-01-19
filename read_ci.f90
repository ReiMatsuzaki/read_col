  ! Read CSF and CIVEC files produced by cCOLMBUS program and
  ! store CI coefficient, its eigene nergies and configurations.

module Mod_CI
  implicit none
  type CI
     integer   :: ncidim
     integer   :: nstate
     complex*16, allocatable :: eig(:)
     complex*16, allocatable :: coef(:,:)
  end type CI
contains
  subroutine CI_new_read_civec(this)
    implicit none
    type(CI) this
    integer, parameter :: ifile = 12
    real*8 blabel(10) ! see ciprop.f
    integer name1(20) ! see ciprop.f
    integer ncidim, nstate, n

    open(unit=ifile, file='CIVEC', status='old', form='unformatted')

    read(ifile) ncidim, blabel, name1, nstate
    this % ncidim = ncidim
    this % nstate = nstate
    write(*, *) "ncidim: ", this % ncidim
    write(*, *) "nstate:   ", this % nstate

    allocate(this % eig(nstate))
    allocate(this % coef(nstate,  ncidim))
        
    do n = 1, this % nstate
       read(ifile) this % eig(n), this % coef(n,:)
       write(*, *) this % eig(n)
    end do
    
    close(unit=ifile)
    
  end subroutine CI_new_read_civec
  subroutine CI_delete(this)
    type(CI) this
    deallocate(this % coef)
    deallocate(this % eig)
  end subroutine CI_delete
end module Mod_CI

module Mod_CSF
  implicit none
  type CSF
  end type CSF
contains
  subroutine CSF_new_read()
end module Mod_CSF

