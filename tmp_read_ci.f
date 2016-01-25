C      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 BLABEL(10)
      integer NAME1(20), N, NSTATE
      NFL12=12
      write(*, *) "MAMAMM"
      open(unit=nfl12,file='CIVEC',status='old',form='unformatted')
      REWIND NFL12
      READ(NFL12) N,BLABEL,NAME1,NSTATE
      write(*, *) "N=", N
      write(*, *) "NSTATE=", NSTATE
      STOP
      END
