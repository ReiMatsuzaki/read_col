
!

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
  use Mod_BlockMat
  use Mod_SymERI
  implicit none
  type AoInts
     character blabel(80)
     character*4 ityp(8), mtyp(10)
     real*8 repfunc
     complex*16 zscale
     integer nst, ns, isfr
     integer*2 nd(8), nso(8), ms(142), mnl(142), kstar(142)
     type(BlockMat) s_mat, t_mat, v_mat
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

    write(*, *) &
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

    write(*, *) "End"

  end subroutine AoInts_read_header
  subroutine AoInts_new_read(this, ifile)
    type(AoInts)        :: this
    integer, intent(in) :: ifile

    call AoInts_read_header(this, ifile)
    call BlockMat_new_read(this % s_mat, ifile)
    call BlockMat_new_read(this % t_mat, ifile)
    call BlockMat_new_read(this % v_mat, ifile)
    call SymERI_new_read(this % eri, ifile)

    rewind ifile

    call AoInts_read_header(this, ifile)
    call BlockMat_set_read(this % s_mat, ifile)
    call BlockMat_set_read(this % t_mat, ifile)
    call BlockMat_set_read(this % v_mat, ifile)
    call SymERI_set_read(this % eri, ifile)

    close(unit = ifile)
    
  end subroutine AoInts_new_read
  subroutine AoInts_delete(this)
    type(AoInts) this
    call BlockMat_delete(this % s_mat)
    call BlockMat_delete(this % t_mat)
    call BlockMat_delete(this % v_mat)
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
    ! write(*, *) this % zscale
    write(*, *) this % isfr
    do iso = 1, this % isfr
       write(*, *) this % ms(iso), this % mnl(iso)
    end do
    !write(*, *) (this % kstar(iso), iso=1, this % isfr)
    
  end subroutine AoInts_show
end module Mod_AoInts

