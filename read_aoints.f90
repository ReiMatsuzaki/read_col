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
