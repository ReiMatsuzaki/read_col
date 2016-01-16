module Mod_SymERI
  type SymERI
     integer num
     integer, allocatable i(:)
     integer, allocatable j(:)
     integer, allocatable k(:)
     integer, allocatable l(:)
     integer, allocatable isym(:)
     integer, allocatable jsym(:)
     integer, allocatable ksym(:)
     integer, allocatable lsym(:)
     complex*16, allocatable val(:)     
  end type SymERI
end module Mod_SymERI

module Mod_SymMat
    type SymMatrix
     integer num
     integer, allocatable i(:)
     integer, allocatable j(:)
     integer, allocatable isym(:)
     complex*16, allocatable val(:)
  end type SymMatrix
end module Mod_SymMat

module Mod_AoInts
  type AoInts
     character blabel(80)
     real*8 repfunc
     integer nst, ns, isfr
     integer*2 nd(8), nso(8), ms(142), mnl(142), kstar(142)
  end type AoInts
end module Mod_AoInts

program read_aoints
  implicit none

  ! This program read AOINTS binary file produced by molint
  ! and int.in which is input file of molint and store the values in memory.
  ! The program molint is part of ccolumbus.
  ! 2015/1/15 R.Matsuzaki
  
  integer, parameter ::  ifile = 3
  character :: blabel(10), ityp(8)
  integer :: i, nst, nso(8), ns, mtype(10), isfr, ms(142), mnl(142), kstar(142)
  real :: repnuc, nd(8), zscale

  open(unit=ifile, file='AOINTS', status='old', form='unformatted')
  
  read(ifile, *) blabel, repnuc, &
       nst, (nd(i),i=1,nst), (ityp(i),i=1,nst), (nso(i),i=1,nst), &
       ns, (mtype(i),i=1,ns), &
       isfr, (ms(i),i=1,isfr), (mnl(i),i=1,isfr), (kstar(i),i=1,isfr), zscale

  close(unit=ifile)
  
end program read_aoints
