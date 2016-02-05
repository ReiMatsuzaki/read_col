module Mod_BlockVec
  implicit none
  type BlockVec
     integer :: nst
     integer, allocatable    :: offset_ist(:)
     integer, allocatable    :: num_ist(:)
     complex*16, allocatable :: val(:)
  end type BlockVec
contains
  ! ==== Constructors ====
  subroutine BlockVec_new(this, num_ist)
    type(BlockVec), intent(out) :: this
    integer, intent(in)       :: num_ist(:)

    integer ist ! index for symmetry
    integer nst ! number of symmetry

    nst = size(num_ist)
    this % nst = nst
    allocate(this % offset_ist(nst))
    allocate(this % num_ist(nst))
    allocate(this % val(sum(num_ist)))
    
    this % offset_ist(1) = 0
    do ist = 2, nst
       this % offset_ist(ist) = &
            this % offset_ist(ist-1) + num_ist(ist-1)
    end do
    this % num_ist(:) = num_ist(:)
    
  end subroutine BlockVec_new
  subroutine BlockVec_new_copy(this, from)
    type(BlockVec), intent(out)  :: this
    type(BlockVec), intent(in) :: from
    call BlockVec_new(this, from % num_ist)

    this % nst           = from % nst
    this % offset_ist(:) = from % offset_ist(:)
    this % num_ist(:)    = from % num_ist(:)
    this % val(:)        = from % val(:)
    
  end subroutine BlockVec_new_copy
  subroutine BlockVec_delete(this)
    type(BlockVec), intent(inout) :: this
    deallocate(this % offset_ist)
    deallocate(this % num_ist)
    deallocate(this % val)
  end subroutine BlockVec_delete

  ! ==== Accessor ====
  function BlockVec_num_sym(this)
    type(BlockVec), intent(in) :: this
    integer :: BlockVec_num_sym
    BlockVec_num_sym = size(this % offset_ist)
  end function BlockVec_num_sym
  function BlockVec_index(this, ist, i)
    type(BlockVec), intent(in) :: this
    integer, intent(in)      :: ist, i
    integer :: BlockVec_index
    BlockVec_index = this % offset_ist(ist) + i
  end function BlockVec_index
  function BlockVec_size(this, ist)
    type(BlockVec), intent(in) :: this
    integer, intent(in) :: ist
    integer :: BlockVec_size
    ! BlockVec_size = this % offset_ist(ist+1) - this % offset_ist(ist)
    BlockVec_size = this % num_ist(ist)
  end function BlockVec_size
  subroutine BlockVec_set_val(this, ist, i, ele)
    type(BlockVec), intent(inout) :: this
    integer, intent(in)         :: ist, i
    complex*16, intent(in)      :: ele
    integer :: idx
    idx = BlockVec_index(this, ist, i)
    this % val(idx) = ele
  end subroutine BlockVec_set_val
  subroutine BlockVec_set_block(this, ist, vals)
    type(BlockVec), intent(inout) :: this
    integer, intent(in)         :: ist
    complex*16, intent(in)      :: vals(:)
    integer :: n0, n1

    if(BlockVec_size(this, ist) /= size(vals)) then
       write(*, *) "BlockVec_set_block"
       write(*, *) "size mismatch"
       stop
    end if

    n0 = this % offset_ist(ist)+1
    n1 = this % offset_ist(ist+1)
    this % val(n0:n1) = vals(:)
    
  end subroutine BlockVec_set_block
  function BlockVec_val(this, ist, i)
    type(BlockVec), intent(in) :: this
    integer, intent(in)         :: ist, i
    integer :: idx
    complex*16 :: BlockVec_val
    idx =  BlockVec_index(this, ist, i)
    BlockVec_val = this % val(idx)
  end function BlockVec_val
  subroutine BlockVec_get_block(this, ist, vals)
    type(BlockVec), intent(inout) :: this
    integer, intent(in)         :: ist
    complex*16, intent(out)      :: vals(:)
    integer :: n0, n1

    if(BlockVec_size(this, ist) /= size(vals)) then
       write(*, *) "BlockVec_get_block"
       write(*, *) "size mismatch"
       stop
    end if
    
    n0 = this % offset_ist(ist)+1
    n1 = this % offset_ist(ist+1)
    vals(:) = this % val(n0:n1) 
    
  end subroutine BlockVec_get_block
  
  ! ==== I/O ====
  subroutine BlockVec_show(this)
    type(BlockVec), intent(in) :: this
    integer    :: ist, i
    do ist = 1, this % nst
       do i = 1, this % num_ist(ist)
          write(*, *) ist, BlockVec_val(this, ist, i)
       end do
    end do
  end subroutine BlockVec_show
  subroutine BlockVec_write(this, ifile)
    type(BlockVec), intent(in) :: this
    integer, intent(in)      :: ifile
    integer :: nst, ist, i
    nst = this % nst
    write(ifile, *) this % nst
    write(ifile, *) (this % num_ist(ist), ist = 1, nst)
    do i = 1, size(this % val)
       write(ifile, *) this % val(i)
    end do
    
  end subroutine BlockVec_write
  subroutine BlockVec_new_read(this, ifile)
    type(BlockVec), intent(out) :: this
    integer, intent(in)      :: ifile
    integer :: nst, ist, i
    integer, allocatable :: num_ist(:)
    
    read(ifile, *) nst
    allocate(num_ist(nst))
    read(ifile, *) (num_ist(ist), ist = 1, nst)
    call BlockVec_new(this, num_ist)

    do i = 1, size(this % val)
       read(ifile, *) this % val(i)
    end do
    
  end subroutine BlockVec_new_read

end module Mod_BlockVec

module Mod_BlockMat  
  ! Mod_BlockMat
  ! read from AOINTS file and store store symmetry distinct block matrix
  !
  ! Variables
  ! ----------
  ! val           : [complex] : store matrix values
  ! offset_iblock : [integer] : offset_iblock(i) gives offset in val array for 
  !                             ith block
  ! isym_iblock   : [integer] : isym_iblock(i) gives symmetry index for bra function
  !                             for ith block
  ! jsym_iblock   : [integer] : jsym_iblock(i) gives symmetry index for ket function
  !                             for ith block
  ! iblock_ijsym : [integer, integer] : iblock_ijsym(i, j) gives block index of
  !                                     ith symmetry and jth symmetry if exist and
  !                                     gives 0 if not exist
  !
  ! Subroutines
  ! ------------
  ! new : constructor
  ! delete : destructor
  ! index : compute index in val array
  ! set : set matrix element value (Warning: given index is not checked)
  ! get : get matrix element value (Warning: given index is not checked)
  ! get_block_size : get given symmetry block size if exist
  ! new_read : read AOINTS and allocate memories
  ! set_read : read AOINTS and set values
  ! show : shows internal datas
  implicit none
  type BlockMat
     complex*16, allocatable :: val(:)
     integer, allocatable    :: offset_iblock(:)
     integer, allocatable    :: isym_iblock(:)
     integer, allocatable    :: jsym_iblock(:)
     integer, allocatable    :: num_isym(:)
     integer, allocatable    :: iblock_ijsym(:, :)
  end type BlockMat
contains
  ! ==== Constructors ====
  subroutine BlockMat_new(this, num_isym, isym_iblock, jsym_iblock)
    type(BlockMat) this
    integer num_isym(:), isym_iblock(:), jsym_iblock(:)
    integer num_sym, num_val, num_block
    integer isym, jsym, iblock

    ! input validation
    if(size(isym_iblock) .ne. size(jsym_iblock)) then
       write(*, *) "Error on BlockMat_new."
       write(*, *) " size(isym_iblock) .ne. size(jsym_iblock)"
       stop
    end if

    ! compute numbers
    num_sym = size(num_isym)
    num_block = size(isym_iblock)    
    num_val = 0
    do iblock = 1, num_block
       num_val = num_val + &
            num_isym(isym_iblock(iblock)) * &
            num_isym(jsym_iblock(iblock))
    end do

    ! memory allocation
    allocate(this % val(num_val))
    allocate(this % offset_iblock(num_block+1))
    allocate(this % isym_iblock(num_block))
    allocate(this % jsym_iblock(num_block))
    allocate(this % num_isym(num_sym))
    allocate(this % iblock_ijsym(num_sym, num_sym))
    this % val = (7.77d0, 7.77d0)

    ! offset_iblock
    this % offset_iblock(1) = 0
    do iblock = 2, num_block + 1
       this % offset_iblock(iblock) = this % offset_iblock(iblock-1) + &
            num_isym(isym_iblock(iblock-1)) * &
            num_isym(jsym_iblock(iblock-1))       
    end do

    ! isym_iblock, jsym_iblock
    this % isym_iblock(:) = isym_iblock(:)
    this % jsym_iblock(:) = isym_iblock(:)

    ! num_isym
    this % num_isym(:) = num_isym(:)

    ! iblock_ijsym
    this % iblock_ijsym(:, :) = 0
    do iblock = 1, num_block
       isym = isym_iblock(iblock)
       jsym = jsym_iblock(iblock)
       this % iblock_ijsym(isym, jsym) = iblock
    end do
    
  end subroutine BlockMat_new
  subroutine BlockMat_delete(this)
    type(BlockMat) this
    deallocate(this % val)
    deallocate(this % offset_iblock)
    deallocate(this % isym_iblock)
    deallocate(this % jsym_iblock)
    deallocate(this % num_isym)
    deallocate(this % iblock_ijsym)
  end subroutine BlockMat_delete

  ! ==== Basic ====
  function BlockMat_index(this, isym, jsym, i, j)
    type(BlockMat), intent(in) ::this
    integer, intent(in)           :: isym, jsym, i, j
    integer :: BlockMat_index
    BlockMat_index = this % offset_iblock(this % iblock_ijsym(isym, jsym)) + &
         (j-1) * this % num_isym(jsym) + (i-1) + 1
  end function BlockMat_index
  function BlockMat_exist_block(this, isym, jsym)
    type(BlockMat), intent(in) :: this
    integer, intent(in) :: isym, jsym
    logical BlockMat_exist_block, res
    res = (this % iblock_ijsym(isym, jsym) .ne. 0)
    BlockMat_exist_block = res
  end function BlockMat_exist_block
  subroutine BlockMat_check_block(this, isym, jsym, label)
    type(BlockMat), intent(in) :: this
    integer, intent(in)           :: isym, jsym
    character(*), intent(in)      :: label
    
    if(.not. BlockMat_exist_block(this, isym, jsym)) then
       write(*, *) "Error"
       write(*, *) "Failed to find block matrix", isym, jsym
       write(*, *) label
       stop
    end if

  end subroutine BlockMat_check_block
  
  ! ==== Accessor ====
  function BlockMat_num_sym(this)
    type(BlockMat), intent(in) :: this
    integer :: BlockMat_num_sym
    BlockMat_num_sym = size(this % num_isym(:))
  end function BlockMat_num_sym
  function BlockMat_block_size(this, isym, jsym)
    type(BlockMat), intent(in) :: this
    integer, intent(in) :: isym, jsym
    integer :: BlockMat_block_size(2)

    call BlockMat_check_block(this, isym, jsym, "BlockMat_block_size")

    BlockMat_block_size(:) = (/this % num_isym(isym), this % num_isym(jsym) /)
    
  end function BlockMat_block_size
  subroutine BlockMat_set_val(this, isym, jsym, i, j, val)
    type(BlockMat)      :: this
    integer, intent(in)    :: isym, jsym, i, j
    complex*16, intent(in) :: val
    this % val(BlockMat_index(this, isym, jsym, i, j)) = val
  end subroutine BlockMat_set_val
  subroutine BlockMat_set_block(this, isym, jsym, mat)
    type(BlockMat)      :: this
    integer, intent(in)    :: isym, jsym
    complex*16, intent(in) :: mat(:, :)
    integer :: iblock, idx1, idx2, num, numij(2)

    call BlockMat_check_block(this, isym, jsym, "BlockMat_set_block")

    iblock = this % iblock_ijsym(isym, jsym)
    idx1 = this % offset_iblock(iblock)
    idx2 = this % offset_iblock(iblock + 1)
    numij = BlockMat_block_size(this, isym, jsym)
    num = numij(1) * numij(2)
    this % val(idx1+1:idx2) = reshape(mat,(/num/))
    
  end subroutine BlockMat_set_block
  function BlockMat_val(this, isym, jsym, i, j)
    type(BlockMat), intent(in) :: this
    integer, intent(in) :: isym, jsym, i, j
    complex*16 :: BlockMat_val
    BlockMat_val = this % val(BlockMat_index(this, isym, jsym, i, j))
  end function BlockMat_Val
  subroutine BlockMat_block(this, isym, jsym, res)
    type(BlockMat) this
    integer, intent(in) :: isym, jsym
    integer :: iblock, idx1, idx2, num(2)
    complex*16 :: res(:, :)
    
    iblock = this % iblock_ijsym(isym, jsym)
    idx1 = this % offset_iblock(iblock)
    idx2 = this % offset_iblock(iblock + 1)
    num = BlockMat_block_size(this, isym, jsym)
    res(:, :) = reshape(this % val(idx1+1 : idx2), num)
         
  end subroutine BlockMat_block

  ! ==== I/O ====
  subroutine BlockMat_show(this)
    type(BlockMat) this
    integer i, j, isym, jsym, iblock, idx
    write(*, *) "BlockMat"
    write(*, *) "size(val)    : ", size(this % val)
    write(*, *) "# of block   : ", size(this % offset_iblock)
    write(*, *) "# of symmetry: ", size(this % num_isym)
    do iblock = 1, size(this % offset_iblock)
       isym = this % isym_iblock(iblock)
       jsym = this % jsym_iblock(iblock)
       write(*, *) "(isym, jsym): ", isym, jsym
       do i = 1, this % num_isym(isym)
          do j = 1, this % num_isym(jsym)
             idx = BlockMat_index(this, isym, jsym, i, j)
             write(*, *) i, j, this % val(idx)
          end do
       end do
    end do
  end subroutine BlockMat_show
  subroutine BlockMat_new_read(this, ifile)
    type(BlockMat)   :: this
    integer, intent(in) :: ifile
    
    integer iblk, ibuf
    integer lbli(1080)
    complex*16 spi(1080)

    integer i, j, i_sym, int, iblock
    integer, parameter :: num_sym_u = 10
    integer num_sym, num_block
    integer num_isym(num_sym_u)
    integer, allocatable ::  isym_iblock(:)
    integer, parameter ::  mask1 =  "000007FF"X

    num_isym(:) = 0

    iblk = 0
    num_sym = 0
    do while(iblk .eq. 0)
       read(ifile) iblk, ibuf, lbli, spi
       do int = 1, ibuf
          j = iand(lbli(int), mask1)
          i = iand(ishft(lbli(int), -15), mask1)
          i_sym = ishft(lbli(int), -26)
          if(num_sym < i_sym) then
             num_sym = i_sym
          end if
          if(num_isym(i_sym) < i) then
             num_isym(i_sym) = i
          end if
          if(num_isym(i_sym) < j) then
             num_isym(i_sym) = j
          end if
       end do
    end do

    num_block = num_sym
    allocate(isym_iblock(num_block))
    do iblock = 1, num_block
       isym_iblock(iblock) = iblock
    end do
    
    call BlockMat_new(this, num_isym(1:num_sym), isym_iblock, isym_iblock)
    deallocate(isym_iblock)
    
  end subroutine BlockMat_new_read
  subroutine BlockMat_set_read(this, ifile)
    type(BlockMat) this
    integer ifile
    integer iblk, ibuf
    integer int
    integer lbli(1080)
    integer, parameter ::  mask1 = "000007FF"X
    complex*16 spi(1080), val
    integer i, j, isym

    iblk = 0
    do while(iblk .eq. 0)
       read(ifile) iblk, ibuf, lbli, spi
       do int = 1, ibuf
          j = iand(lbli(int), mask1)
          i = iand(ishft(lbli(int), -15), mask1)
          isym = ishft(lbli(int), -26)
          val = spi(int)
          call BlockMat_set_val(this, isym, isym, i, j, val)
          call BlockMat_set_val(this, isym, isym, j, i, val)
       end do
    end do
    
  end subroutine BlockMat_set_read

end module Mod_BlockMat

module Mod_BlockLinearAlgebra
  use Mod_BlockVec
  use Mod_BlockMat
  implicit none
contains
  function block_matmul(M, v, t)
    type(BlockMat), intent(in) :: M
    type(BlockVec), intent(in) :: v
    character, optional :: t
    type(BlockVec) :: block_matmul

    if(BlockVec_num_sym(v) /= BlockMat_num_sym(M)) then
       write(*, *) "block_matmul"
       write(*, *) "number of symmetry is not match"
       stop
    end if

    if((t .eq. 'T') .or. (t .eq. 't')) then
       call block_matmul_normal(M, v, block_matmul)
    end if
    
  end function block_matmul
  subroutine block_matmul_normal(M, v, res)
    type(BlockMat), intent(in) :: M
    type(BlockVec), intent(in) :: v
    type(BlockVec), intent(out) :: res
    integer iblock, isym, jsym, num_sym, idx0, idx1, ni, nj, idx2, idx3

    num_sym = BlockVec_num_sym(v)
    call BlockVec_new_copy(res, v)

    res % val(:) = (0.0d0, 0.0d0)

    do iblock = 1, size(M % isym_iblock)
       isym = M % isym_iblock(iblock)
       jsym = M % jsym_iblock(iblock)
       ni = M % num_isym(isym)
       nj = M % num_isym(jsym)
       idx0 = M % offset_iblock(iblock) + 1
       idx1 = M % offset_iblock(iblock+1)
       idx2 = v % offset_ist(iblock) + 1
       idx3 = v % offset_ist(iblock+1)
!       do i = 1, ni
!          res % val(res % offset_ist(isym)+i) = res % val(res % offset_ist(isym)+i) &
!               + M % val(M % offset_iblock(iblock) + )
!          res % offset_ist(isym) + 1, res % offset_ist(isym+1)
!          do j = 1, nj
!             res % val(i) = res % val(i)
!          end do
!       end do
       call BlockVec_set_block(res, isym, &
            matmul(reshape(M % val(idx0:idx1), (/ni, nj/)), v % val(idx2:idx3)))
    end do
  end subroutine block_matmul_normal
!  subroutine block_matmul_transpose(M, v, res)
!    type(BlockMat), intent(in) :: M
!    type(BlockVec), intent(in) :: v
!    type(BlockVec), intent(out) :: res    
!  end subroutine block_matmul_transpose
end module Mod_BlockLinearAlgebra
