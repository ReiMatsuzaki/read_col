C      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 REPNUC
      character*8 BLABEL
C      INTEGER*2 ND(8),NSO(8),MS(142),MNL(142),KSTAR(142)
C      COMPLEX*16 PA,PB,AI,ZSCALE
      DIMENSION BLABEL(10)
C ALBL(810),AI(810),ITYP(8),INTYPE(3),NSOFF(8)
C      DIMENSION LBLI(1080),ISYMOF(8),IOFFST(256),MTYPE(10)
C      EQUIVALENCE(ALBL(1),LBLI(1))
C      COMMON/PASS/PA(1296),PB(1296),LBLINT(1296),NINT
C      DATA INTYPE /4H  S ,4H  T ,4H  V /
      IFILE=3
C 905  FORMAT(//' INTEGRALS ARE READ IN FROM   FILE ',I2,' (IFILE).'//
C     1' INTEGRALS ARE WRITTEN OUT ON FILE ',I2,' (NFILE).'//)
      open (unit=IFILE, file='AOINTS',status='OLD',form='unformatted')
C      READ(IFILE) BLABEL,REPNUC,NST,(ND(IST),IST=1,NST),
C     1   (ITYP(IST),IST=1,NST),(NSO(IST),IST=1,NST),NS,
C     2   (MTYPE(IS),IS=1,NS),ISFR,(MS(ISO),ISO=1,ISFR),(MNL(ISO),
C     3   ISO=1,ISFR),(KSTAR(ISO),ISO=1,ISFR),ZSCALE
      READ(IFILE) BLABEL,REPNUC
      WRITE(6,910) BLABEL,REPNUC
 910  FORMAT(' ',10A8//' NUCLEAR REPULSION = ',D25.15//)


C      character blabel
C      DIMENSION blabel(10)
C      open(unit=3, file='AOINTS', status='old', form='unformatted')
C      read(3) blabel, rep_func
C      write(*,910) blabel, rep_func
C 910  format(' ',10A8//' nuclear repulsion = ', D25.15//)
      stop
      end
