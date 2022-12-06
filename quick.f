      PROGRAM inandout
      IMPLICIT NONE
      INTEGER i,j,n,kounter
      PARAMETER(n=1000)
      DOUBLE PRECISION freq(n),field(n),auxx

      OPEN(unit=20,file='SUSZI.dat',STATUS='unknown')
      kounter=0
      DO i=1,n
      READ(20,*,END=95) freq(i),field(i)
      kounter=kounter+1
      END DO
 95   CONTINUE
c      DO i=1,kounter
c         field(i)=1.0*freq(i)**(-0.9d0)
c         write(77,*) freq(I),field(i)
c      END DO
c      STOP
      auxx=dlog(10.0d0)
      OPEN(UNIT=88,FILE='exponent.dat',STATUS='unknown')
      DO i=1,kounter
         IF(freq(i).GT.0.0d0) THEN
            WRITE(88,*) dlog(freq(i))/auxx,dlog(dabs(field(i)))/auxx
         END IF
      END DO
      CLOSE(88)
      STOP
      END
