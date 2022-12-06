CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Programm liest Spektralfunktionen und berechnet die Selfenergien
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program selfplot
      IMPLICIT NONE
      INTEGER n,n2,i,j,m
      PARAMETER(n=282)
      PARAMETER(n2=64)
      DOUBLE PRECISION omega(n),af(n),ab(n),gf(n),templist(n2)
      DOUBLE PRECISION pi,gb(n),sigmaf(n),sigmab(n)
      DOUBLE PRECISION temp,dummy,beta,temp2,xxx(n2),yyy(n2),RG
      DOUBLE PRECISION zzz(n2),vvv(n2)
      PARAMETER(pi=3.1415926535898d0)
      CHARACTER*7   handle
      CHARACTER*10  name1,name2
      CHARACTER*14  name3,name4

      OPEN(UNIT=20,FILE='TEMPLIST.dat',STATUS='UNKNOWN')
      DO i=1,n2
         READ(20,*,END=99) temp
         templist(i)=temp
      END DO
 99   CONTINUE

      DO j=1,n2
         temp=templist(j)
         beta=1.0/temp
         WRITE(handle,'(E7.2)') temp
         name2='ab_'//handle
         name1='af_'//handle
         name3='sigmab_'//handle
         name4='sigmaf_'//handle
         OPEN(UNIT=91,FILE=name1,STATUS='UNKNOWN')
         OPEN(UNIT=92,FILE=name2,STATUS='UNKNOWN')
         OPEN(UNIT=93,FILE=name3,STATUS='UNKNOWN')
         OPEN(UNIT=94,FILE=name4,STATUS='UNKNOWN')
         DO i=1,n
            READ(91,*) omega(i),af(i),temp2,RG
            READ(92,*) dummy,ab(i)
         END DO
         CALL KRAKRO(n,omega,af,gf)
         CALL KRAKRO(n,omega,ab,gb)
         DO i=1,n
            sigmaf(i)=af(i)/(af(i)**2+gf(i)**2)
            sigmab(i)=ab(i)/(ab(i)**2+gf(i)**2)
         END DO
         DO i=1,n
            WRITE(94,*) omega(i),sigmaf(i),temp,RG
            WRITE(93,*) omega(i),sigmab(i),temp,RG
         END DO
         CLOSE(91)
         CLOSE(92)
         CLOSE(93)
         CLOSE(94)

      END DO


      STOP
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL*8 x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
