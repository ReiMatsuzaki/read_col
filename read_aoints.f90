module Mod_SymMat
  type SymMat
     integer  num
     complex*16, allocatable :: val(:, :)
  end type SymMat
contains
  subroutine SymMat_new(this, num)
    type(SymMat) this
    integer num
    integer i, j
    this % num = num
    allocate(this % val(num, num))
    do i = 1, num
       do j = 1, num
          this % val(i, j) = 0.0d0
       end do
    end do
  end subroutine SymMat_new
  subroutine SymMat_delete(this)
    type(SymMat) this
    deallocate(this % val)
  end subroutine SymMat_delete
  subroutine SymMat_set(this, i, j, val)
    type(SymMat) this
    integer i, j
    complex*16 val
    this % val(i, j) = val
    this % val(j, i) = val
  end subroutine SymMat_set
  subroutine SymMat_show(this)
    type(SymMat) this
    integer i,j,num
    num = this % num
    write(*, *) num
    do i = 1, num
       write(*, *) (this % val(i, j), j =1, num)
    end do
  end subroutine SymMat_show
end module Mod_SymMat


module Mod_SymBlockMat
  use Mod_SymMat
  implicit none
  type SymBlockMat
     complex*16, allocatable :: val(:)
     integer, allocatable    :: offset_iblock(:)
     integer, allocatable    :: isym_iblock(:)
     integer, allocatable    :: jsym_iblock(:)
     integer, allocatable    :: num_isym(:)
     integer, allocatable    :: iblock_ijsym(:, :)
  end type SymBlockMat
contains
  subroutine SymBlockMat_new(this, num_isym, isym_iblock, jsym_iblock)
    type(SymBlockMat) this
    integer num_isym(:), isym_iblock(:), jsym_iblock(:)
    integer num_sym, num_val, num_block
    integer isym, jsym, iblock

    ! input validation
    if(size(isym_iblock) .ne. size(jsym_iblock)) then
       write(*, *) "Error on SymBlockMat_new."
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
    allocate(this % offset_iblock(num_block))
    allocate(this % isym_iblock(num_block))
    allocate(this % jsym_iblock(num_block))
    allocate(this % num_isym(num_sym))
    allocate(this % iblock_ijsym(num_sym, num_sym))

    ! offset_iblock
    this % offset_iblock(1) = 0
    do iblock = 2, num_block
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
    do isym = 1, num_sym
       do jsym = 1, num_sym
          this % iblock_ijsym(isym, jsym) = 0
       end do
    end do
    do iblock = 1, num_block
       isym = isym_iblock(iblock)
       jsym = jsym_iblock(iblock)
       this % iblock_ijsym(isym, jsym) = iblock
    end do
    
  end subroutine SymBlockMat_new
  subroutine SymBlockMat_delete(this)
    type(SymBlockMat) this
    deallocate(this % val)
    deallocate(this % offset_iblock)
    deallocate(this % isym_iblock)
    deallocate(this % jsym_iblock)
    deallocate(this % num_isym)
    deallocate(this % iblock_ijsym)
  end subroutine SymBlockMat_delete
  integer function SymBlockMat_index(this, isym, jsym, i, j) result(res)
    type(SymBlockMat) this
    integer, intent(in)    :: isym, jsym, i, j
    res = this % offset_iblock(this % iblock_ijsym(isym, jsym)) + &
         (i-1) * this % num_isym(jsym) + (j-1) + 1
  end function SymBlockMat_index
  subroutine SymBlockMat_set(this, isym, jsym, i, j, val)
    type(SymBlockMat)      :: this
    integer, intent(in)    :: isym, jsym, i, j
    complex*16, intent(in) :: val
    this % val(SymBlockMat_index(this, isym, jsym, i, j)) = val
  end subroutine SymBlockMat_set
  subroutine SymBlockMat_get(this, isym, jsym, i, j, res)
    type(SymBlockMat) this
    integer, intent(in) :: isym, jsym, i, j
    complex*16, intent(out) :: res
    res = this % val(SymBlockMat_index(this, isym, jsym, i, j))
  end subroutine SymBlockMat_get
  subroutine SymBlockMat_get_block_size(this, isym, jsym, num_i, num_j, have_val)
    type(SymBlockMat) :: this
    integer, intent(in) :: isym
    integer, intent(in) :: jsym
    integer, intent(out) :: num_i
    integer, intent(out) :: num_j
    logical, intent(out) :: have_val
    integer iblock

    iblock = this % iblock_ijsym(isym, jsym)
    if(iblock .eq. 0) then
       num_i = 0
       num_j = 0
       have_val = .false.
    else
       num_i = this % num_isym(isym)
       num_j = this % num_isym(jsym)
       have_val = .true.
    end if
    
  end subroutine SymBlockMat_get_block_size
  subroutine SymBlockMat_new_read(this, ifile)
    type(SymBlockMat)   :: this
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

    ! call SymBlockMat_show(this)

    do i_sym = 1, num_sym_u
       num_isym(i_sym) = 0
    end do

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
    
    call SymBlockMat_new(this, num_isym(1:num_sym), isym_iblock, isym_iblock)
    deallocate(isym_iblock)
    
  end subroutine SymBlockMat_new_read
  subroutine SymBlockMat_set_read(this, ifile)
    type(SymBlockMat) this
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
          call SymBlockMat_set(this, isym, isym, i, j, val)
          call SymBlockMat_set(this, isym, isym, j, i, val)
       end do
    end do
    
  end subroutine SymBlockMat_set_read
  subroutine SymBlockMat_show(this)
    type(SymBlockMat) this
    integer i, j, isym, jsym, iblock, idx
    write(*, *) "SymBlockMat"
    write(*, *) "size(val)    : ", size(this % val)
    write(*, *) "# of block   : ", size(this % offset_iblock)
    write(*, *) "# of symmetry: ", size(this % num_isym)
    do iblock = 1, size(this % offset_iblock)
       isym = this % isym_iblock(iblock)
       jsym = this % jsym_iblock(iblock)
       write(*, *) "(isym, jsym): ", isym, jsym
       do i = 1, this % num_isym(isym)
          do j = 1, this % num_isym(jsym)
             idx = SymBlockMat_index(this, isym, jsym, i, j)
             write(*, *) i, j, this % val(idx)
          end do
       end do
    end do
  end subroutine SymBlockMat_show
end module Mod_SymBlockMat


module Mod_SymERI
  implicit none
  type SymERI
     integer num
     integer, allocatable :: i(:)
     integer, allocatable :: j(:)
     integer, allocatable :: k(:)
     integer, allocatable :: l(:)
     integer, allocatable :: isym(:)
     integer, allocatable :: jsym(:)
     integer, allocatable :: ksym(:)
     integer, allocatable :: lsym(:)
     complex*16, allocatable :: val(:)     
  end type SymERI
contains
  subroutine SymERI_new(this, num)
    type(SymERI) this
    integer num
    this % num = num
    allocate(this % i(num))
    allocate(this % j(num))
    allocate(this % k(num))
    allocate(this % l(num))
    allocate(this % isym(num))
    allocate(this % jsym(num))
    allocate(this % ksym(num))
    allocate(this % lsym(num))
    allocate(this % val(num))
  end subroutine SymERI_new
  subroutine SymERI_delete(this)
    type(SymERI) this
    deallocate(this % i)
    deallocate(this % j)
    deallocate(this % k)
    deallocate(this % l)
    deallocate(this % isym)
    deallocate(this % jsym)
    deallocate(this % ksym)
    deallocate(this % lsym)
  end subroutine SymERI_delete
  subroutine SymERI_new_read(this, ifile)
    type(SymERI) this
    integer ifile
    integer num, iblk, ibuf
    integer lbli(1620)
    complex*16 spi(810)

    num = 0
    iblk = 0
    do while(iblk .eq. 0)
       read(ifile) iblk, ibuf, lbli, spi
       num = num + ibuf
    end do

    call SymERI_new(this, num)
    
  end subroutine SymERI_new_read
  subroutine SymERI_set_read(this, ifile)
    type(SymERI) this
    integer ifile
    integer iblk, ibuf, idx
    integer int
    integer lbli(1620)
    integer ijlbl, kllbl, m1, m2
    complex*16 spi(810)

    m1 = "000007FF"X
    m2 = 15
    iblk = 0
    idx = 0
    do while(iblk .eq. 0)
       read(ifile) iblk, ibuf, lbli, spi
       do int = 1, ibuf
          idx = idx + 1
          ijlbl = lbli(int + int -1)
          kllbl = lbli(int + int)
          this % ISYM(idx) =ISHFT(IJLBL,-26)
          this % I(idx)=IAND(ISHFT(IJLBL,-15),M1)
          this % JSYM(idx)=IAND(ISHFT(IJLBL,-11),M2)
          this % J(idx)=IAND(IJLBL,M1)
          this % KSYM(idx)=ISHFT(KLLBL,-28)
          this % K(idx)=IAND(ISHFT(KLLBL,-17),M1)
          this % LSYM(idx)=IAND(ISHFT(KLLBL,-13),M2)
          this % L(idx)=IAND(ISHFT(KLLBL,-2),M1)
          this % VAL(idx)=spi(INT)
       end do
    end do
    
  end subroutine SymERI_set_read
  subroutine SymERI_show(this)
    type(SymERI) this
    integer i
    write(*, *) "Symmetric ERI"
    write(*, *) "num = ", this % num
    do i = 1, this % num
       write(*, "(8I3, 2f15.5)") this%i(i), this%j(i), this%k(i),&
            this%l(i), this%isym(i), this%jsym(i), this%ksym(i), this%lsym(i),this%val(i)
    end do
  end subroutine SymERI_show
end module Mod_SymERI


module Mod_AoInts
  use Mod_SymBlockMat
  use Mod_SymERI
  implicit none
  type AoInts
     character blabel(80)
     character*4 ityp(8), mtyp(10)
     real*8 repfunc
     complex*16 zscale
     integer nst, ns, isfr
     integer*2 nd(8), nso(8), ms(142), mnl(142), kstar(142)
     type(SymBlockMat) s_mat, t_mat, v_mat
     type(SymERI) eri
  end type AoInts
contains
  subroutine AoInts_read_header(this, ifile)
    type(AoInts)        :: this
    integer, intent(in) :: ifile
    integer ist, is, iso
    read(ifile) &
         this % blabel, &
         this % repfunc, &
         this % nst, &
         (this % nd(ist), ist=1, this % nst), &
         (this % ityp(ist), ist=1, this % nst), &
         (this % nso(ist), ist=1, this % nst), &
         this % ns, &
         (this % mtyp(is), is=1, this % ns), &
         this % isfr, &
         (this % ms(iso), iso=1, this % isfr), &
         (this % mnl(iso), iso=1, this % isfr), &
         (this % kstar(iso), iso=1, this % isfr), &
         this % zscale

  end subroutine AoInts_read_header
  subroutine AoInts_new_read(this)
    type(AoInts)        :: this
    integer, parameter :: ifile = 3
    open(unit=ifile, file='AOINTS', status='old', form='unformatted')

    call AoInts_read_header(this, ifile)
    call SymBlockMat_new_read(this % s_mat, ifile)
    call SymBlockMat_new_read(this % t_mat, ifile)
    call SymBlockMat_new_read(this % v_mat, ifile)
    call SymERI_new_read(this % eri, ifile)

    rewind ifile

    call AoInts_read_header(this, ifile)
    call SymBlockMat_set_read(this % s_mat, ifile)
    call SymBlockMat_set_read(this % t_mat, ifile)
    call SymBlockMat_set_read(this % v_mat, ifile)
    call SymERI_set_read(this % eri, ifile)

    close(unit = ifile)
    
  end subroutine AoInts_new_read
  subroutine AoInts_delete(this)
    type(AoInts) this
    call SymBlockMat_delete(this % s_mat)
    call SymBlockMat_delete(this % t_mat)
    call SymBlockMat_delete(this % v_mat)
    call SymERI_delete(this % eri)
  end subroutine AoInts_delete
  subroutine AoInts_show(this)
    type(AoInts) this
    integer ist, is, iso
    write(*, *) "label: ", this % blabel
    write(*, *)
    write(*, *) "Symmetry: ", this % nst
    write(*, *) (this % nd(ist), ist=1, this % nst)
    write(*, *)(this % ityp(ist), ist=1, this % nst)
    write(*, *) (this % nso(ist), ist=1, this % nst)
    write(*, *)
    write(*, *) "Symmetry distinct atoms"
    write(*, *) this % ns
    write(*, *) (this % mtyp(is), is=1, this % ns)
    write(*, *)
    write(*, *) "basis"
    write(*, *) this % isfr
    write(*, *) (this % ms(iso), iso=1, this % isfr)
    write(*, *) (this % mnl(iso), iso=1, this % isfr)
    write(*, *) (this % kstar(iso), iso=1, this % isfr)
    write(*, *) this % zscale
  end subroutine AoInts_show
end module Mod_AoInts

