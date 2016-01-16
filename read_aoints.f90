module Mod_SymERI
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
end module Mod_SymERI

module Mod_SymMat
    type SymMatrix
     integer num
     integer, allocatable :: i(:)
     integer, allocatable :: j(:)
     integer, allocatable :: isym(:)
     complex*16, allocatable :: val(:)
  end type SymMatrix
end module Mod_SymMat

module Mod_AoInts
  implicit none
  type AoInts
     character blabel(80)
     character*4 ityp(8), mtyp(10)
     real*8 repfunc
     complex*16 zscale
     integer nst, ns, isfr
     integer*2 nd(8), nso(8), ms(142), mnl(142), kstar(142)
  end type AoInts
contains
  subroutine AoInts_new_read(this, ifile)
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

    write(*, *) this % nst, this % ns
    write(*, *) this % blabel
    
  end subroutine AoInts_new_read
  subroutine AoInts_show(this)
    type(AoInts) this
    integer ist, is, iso
    write(*, *) "label: ", this % blabel
    write(*, *) "Symmetry: ", this % nst
    write(*, *) (this % nd(ist), ist=1, this % nst)
    write(*, *)(this % ityp(ist), ist=1, this % nst)
    write(*, *) (this % nso(ist), ist=1, this % nst)
    write(*, *) this % ns
    write(*, *) (this % mtyp(is), is=1, this % ns)
    write(*, *) this % isfr
    write(*, *) (this % ms(iso), iso=1, this % isfr)
    write(*, *) (this % mnl(iso), iso=1, this % isfr)
    write(*, *) (this % kstar(iso), iso=1, this % isfr)
    write(*, *) this % zscale
  end subroutine AoInts_show
end module Mod_AoInts

