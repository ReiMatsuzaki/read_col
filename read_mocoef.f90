module Mod_MOCoef
  use Mod_BlockMat
  use Mod_BlockVec
  implicit none
  type MOCoef
!     integer nst         ! number of symmetry
!     integer, allocatable    :: n_ist(:)    ! # of orbital of symmetry ist
!     complex*16, allocatable :: coef(:)  ! MO coefficient
!     complex*16, allocatable :: eig(:)   ! orbital energy
     type(BlockMat) :: mo_coef
     type(BlockVec)   :: mo_eig
  end type MOCoef
contains
  ! ==== Constructors ====
  subroutine MOCoef_new(this, n_ist)
    type(MOCoef), intent(inout) :: this
    integer,      intent(in)    :: n_ist(:)
    integer :: ist
    integer :: isym_iblock(size(n_ist))
    
    isym_iblock(:) = (/(ist, ist = 1, size(n_ist))/)
    call BlockMat_new(this % mo_coef, n_ist, isym_iblock, isym_iblock)
    call BlockVec_new(this % mo_eig, n_ist)
  end subroutine MOCoef_new
  subroutine MOCoef_delete(this)
    type(MOCoef) this

    call BlockMat_delete(this % mo_coef)
    call BlockVec_delete(this % mo_eig)
    
  end subroutine MOCoef_delete

  ! ==== Accessor ====
  function MOCoef_num_sym(this)
    type(MOCoef), intent(in) :: this
    integer :: MOCoef_num_sym
    MOCoef_num_sym = BlockVec_num_sym(this % mo_eig)
  end function MOCoef_num_sym
  function MOCoef_size_isym(this, isym)
    type(MOCoef), intent(in) :: this
    integer, intent(in)      :: isym
    integer :: MOCoef_size_isym
    MOCoef_size_isym = BlockVec_size(this % mo_eig, isym)
  end function MOCoef_size_isym
  
  ! ==== I/O ====
  subroutine MOCoef_set_read(this, ifile)
    type(MOCoef) this
    integer, intent(in) ::  ifile
    real*8 fmt(4)
    integer ist, k, j
    integer nn
    complex*16, allocatable :: vals(:)

    read(ifile, '(10A8)') fmt
    do ist = 1, MOCoef_num_sym(this)
       
       nn = MOCoef_size_isym(this, ist)
       allocate(vals(nn))
       do k = 1, nn
          read(ifile, fmt) (vals(j), j = 1, nn)
          do j = 1, nn
             call BlockMat_set_val(this % mo_coef, ist, ist, j, k, vals(j))
          end do
       end do
       deallocate(vals)
    end do


    do ist = 1, MOCoef_num_sym(this)
       nn = MOCoef_size_isym(this, ist)
       allocate(vals(nn))
       read(ifile, fmt) (vals(j), j = 1, nn)
       do j = 1, nn
          call BlockVec_set_val(this % mo_eig, ist, j, vals(j))
       end do
       deallocate(vals)
    end do
    
  end subroutine MOCoef_set_read
end module Mod_MOCoef
