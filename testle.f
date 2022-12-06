      PROGRAM sublead
      IMPLICIT NONE
      INTEGER n,n2,i,j,m,choice
      PARAMETER(n=55)
      DOUBLE PRECISION pi,alpha,a(n),temp(n),rg,amplitude
      PARAMETER(pi=3.1415926535898d0)
      CHARACTER*17  name1,name2
      CHARACTER*6   griff
*************************************
      WRITE(*,*) 'which g?'
      READ(*,*) rg
      WRITE(griff,'(F6.4)') RG

      WRITE(*,*) 'which file?'
      PRINT*
      PRINT*,'A_F ---> press 1'
      PRINT*,'A_B ---> press 2'
      READ(*,*) choice
      WRITE(*,*) 'leading exponent?'
      READ(*,*) alpha
      WRITE(*,*) 'amplitude?'
      READ(*,*) amplitude
      name2='ab_T.dat_'//griff
      name1='af_T.dat_'//griff
      print*,name2
      print*,name1
      IF(choice.EQ.1) THEN
         OPEN(UNIT=91,FILE=name1,STATUS='UNKNOWN')
         DO i=1,n
            READ(91,*) temp(i),a(i)
         END DO
         CLOSE(91)
      ELSEIF(choice.EQ.2) THEN
         OPEN(UNIT=92,FILE=name2,STATUS='UNKNOWN')
         DO i=1,n
            READ(92,*) temp(i),a(i)
         END DO
         CLOSE(92)
      ELSE
         STOP
      ENDIF

      DO i=1,n
c         write(77,*) temp(i),temp(i)**alpha*
c     &        (temp(i)**(-alpha)*a(i)-amplitude)
         write(66,*) temp(i),a(i)-amplitude*temp(i)**alpha
      END DO

      STOP
      END
