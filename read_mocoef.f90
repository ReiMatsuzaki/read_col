module Mod_MOCoef
  use Mod_BlockMat
  use Mod_SymVec
  implicit none
  type MOCoef
     integer nst         ! number of symmetry
     integer, allocatable    :: n_ist(:)    ! # of orbital of symmetry ist
     complex*16, allocatable :: coef(:)  ! MO coefficient
     complex*16, allocatable :: eig(:)   ! orbital energy
     type(BlockMat) :: mo_coef
     type(SymVec)   :: mo_eig
  end type MOCoef
contains
  subroutine MOCoef_new(this, n_ist)
    type(MOCoef), intent(inout) :: this
    integer,      intent(in)    :: n_ist(:)
    integer :: n_coef, n_eig
    integer :: ist
    integer :: isym_iblock(size(n_ist))
    
    this % nst = size(n_ist)
    allocate(this % n_ist(this % nst))

    print *, this % nst
    this % n_ist = n_ist
    print *, this % n_ist

    n_coef = 0
    n_eig = 0
    do ist = 1, this % nst
       n_coef = n_coef + this % n_ist(ist) ** 2
       n_eig = n_eig + this % n_ist(ist)
    end do

    print *, "Allocation:", n_coef, n_eig
    allocate(this % coef(n_coef))
    allocate(this % eig(n_eig))

    isym_iblock(:) = (/(ist, ist = 1, size(n_ist))/)
    call BlockMat_new(this % mo_coef, n_ist, isym_iblock, isym_iblock)
    write(*, *) ""
    write(*, *) "symvec_new>>>"
    call SymVec_new(this % mo_eig, n_ist)
    write(*, *) "symvec_new<<<"
    
  end subroutine MOCoef_new
  subroutine MOCoef_delete(this)
    type(MOCoef) this
    deallocate(this % n_ist)
    deallocate(this % coef)
    deallocate(this % eig)

    call BlockMat_delete(this % mo_coef)
  end subroutine MOCoef_delete
  subroutine MOCoef_set_read(this)
    type(MOCoef) this
    integer ifile
    real*8 fmt(4)
    integer ist, k, j
    integer nn
    integer ic_max, ic_min, i_start, i_fin

    ifile = 13
    open(unit = ifile, file='MOCOEF')
    read(ifile, '(10A8)') fmt

    ic_max = 0
    do ist = 1, this % nst
       
       nn = this % n_ist(ist)
       ic_min = ic_max + 1
       ic_max = ic_max + nn * nn
       
       i_start = ic_min
       i_fin   = ic_min + nn - 1
       do k = 1, nn
!          write(*,'(4I)') ist, k, i_start, i_fin
          read(ifile, fmt) (this % coef(j), j = i_start, i_fin)
          !          write(*,*)  (this % coef(j), j = i_start, i_fin)
          do j = 1, nn
             call BlockMat_set_element(this % mo_coef, ist, ist, j, k, &
                  this % coef(j+i_start-1))
          end do          
          i_start = i_fin + 1
          i_fin = i_fin + nn

       end do
    end do

    ic_max = 0
    do ist = 1, this % nst
       nn = this % n_ist(ist)
       ic_min = ic_max + 1
       ic_max = ic_max + nn
!       print *,  ist, ic_min, ic_max
       read(ifile, fmt) (this % eig(j), j = ic_min, ic_max)
       do j = 1, nn
          call SymVec_set(this % mo_eig, ist, j, this % eig(ic_min+j-1))
       end do
!       write(*,*)  (this % eig(j), j = ic_min, ic_max)
    end do

    close(unit = ifile)
    
  end subroutine MOCoef_set_read
end module Mod_MOCoef
