CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     LOESEN DES KONDO_MODELL mit Kopplung an bosonischen Bad Aphi im 
C     Grenzfall N,M -> unendlich und ratio=M/N endlich
C     Programm liest Spektralfunktionen ein und berechnet Faltung
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program faltung
      IMPLICIT NONE
      INTEGER n,n2,i,j,m
      PARAMETER(n=242)
      PARAMETER(n2=55)
      DOUBLE PRECISION omega(n),af(n),ab(n),tmat(n),templist(n2)
      DOUBLE PRECISION suszi(n),resus(n),pi,gb(n)
      DOUBLE PRECISION temp,dummy,beta,temp2,xxx(n2),yyy(n2),RG
      DOUBLE PRECISION zzz(n2),vvv(n2)
      PARAMETER(pi=3.1415926535898d0)
      CHARACTER*7   handle
      CHARACTER*10  name1,name2
      CHARACTER*17  name3
      CHARACTER*6   griff
      CHARACTER*16  name4
      CHARACTER*18  name5
*************************************
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
         OPEN(UNIT=91,FILE=name1,STATUS='UNKNOWN')
         OPEN(UNIT=92,FILE=name2,STATUS='UNKNOWN')
         DO i=1,n
            READ(91,*) omega(i),af(i),temp2,RG
            READ(92,*) dummy,ab(i)
         END DO
*************************************
         CALL FERMIFUNC0
         CALL BOSEFUNC0(omega,temp)
         
***************************T-MARTIX*************************************
         CALL tmatrix(omega,n,beta,af,ab,tmat)
         CALL DSCAL(n,1.0d0/pi,tmat,1)
         CALL LOCATE(omega,n,0.0d0,m)
         xxx(j)=tmat(m)+(tmat(m+1)-tmat(m))/(omega(m+1)-omega(m))
     &        *(0.0d0-omega(m))
*********************LOCAL SUSCEPTIBILITY*******************************
      CALL SUSCEPTIBILITY(beta,af,omega,n,suszi)
      CALL DSCAL(n,1.0d0/pi,suszi,1)
      CALL KRAKRO(n,omega,suszi,resus)
      CALL LOCATE(omega,n,0.0d0,m)
         yyy(j)=-resus(m)-(resus(m+1)-resus(m))/(omega(m+1)-omega(m))
     &           *(0.0d0-omega(m))
************************************************************************
         zzz(j)=af(m)+(af(m+1)-af(m))/(omega(m+1)-omega(m))
     &           *(0.0d0-omega(m))
         CALL KRAKRO(n,omega,ab,gb)
         vvv(j)=gb(m)+(gb(m+1)-gb(m))/(omega(m+1)-omega(m))
     &           *(0.0d0-omega(m))
************************************************************************
      END DO
      WRITE(griff,'(F6.4)') RG
      name4='TMATT.dat_'//griff
      name5='STATSUS.dat_'//griff
      OPEN(UNIT=99,FILE=name4,STATUS='UNKNOWN')
      OPEN(UNIT=98,FILE=name5,STATUS='UNKNOWN')
      DO j=1,n2
         write(99,*) templist(j),xxx(j)
         write(98,*) templist(j),yyy(j)
      END DO
      CLOSE(99)
      CLOSE(98)
      name3='tmat.dat_'//griff
      OPEN(UNIT=77,FILE=name3,STATUS='UNKNOWN')
      name3='susi.dat_'//griff
      OPEN(UNIT=78,FILE=name3,STATUS='UNKNOWN')
      DO i=1,n
         write(77,*) omega(i),tmat(i)
         write(78,*) omega(i),suszi(i)
      END DO
      CLOSE(77)
      CLOSE(78)
      name3='af_T.dat_'//griff
      OPEN(UNIT=79,FILE=name3,STATUS='UNKNOWN')
      name3='ab_T.dat_'//griff
      OPEN(UNIT=76,FILE=name3,STATUS='UNKNOWN')
      DO j=1,n2
         write(79,*) templist(j),zzz(j)
         write(76,*) templist(j),vvv(j)
      END DO
      CLOSE(79)
      CLOSE(76)
      STOP
      END

      


************************************************************************
************************************************************************
*************************CALCULATE T-MATRIX*****************************
************************************************************************
      SUBROUTINE tmatrix(omega,n,beta,afp,abp,tmat)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION omega(n),afp(n),abp(n),beta,tmat(n)
      INTEGER i,j,k,kk,l,m,mm,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION ferm0,bose0
      DOUBLE PRECISION wfreq1,wfreq2,xxx,yyy,www,D,pi
      EXTERNAL ferm0,bose0

      pi=4.0d0*datan(1.0d0)
      D=omega(n)+1.0d-10
      DO i=1,n
         exom=omega(i)
*****************Integrand bei Null (ohne Bose+Fermi-Funktion)
         CALL LOCATE(omega,n,0.0d0,m)
         xxx=abp(m)+(abp(m+1)-abp(m))/(omega(m+1)-omega(m))
     &           *(0.0d0-omega(m))
         xxx=afp(i)*xxx
***Integral (bose+fermi)
         yyy=-exom+1.0/beta*(dlog(1.0d0+dexp(-beta*(D-exom)))
     &        -1.0/beta*dlog(1.0d0+dexp(-beta*(D+exom))))

         www=xxx*yyy            !to be added to integral
*******************************************************
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (abp(j+1)*
     &          (bose0(omega(j+1))+ferm0(beta*(omega(j+1)+exom)))*
     &          wert2+abp(j)*
     &          (bose0(omega(j))+ferm0(beta*(omega(j)+exom)))*
     &              wert1)
***********************************************************integration
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               wfreq1=abp(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  CALL LOCATE(omega,n,freq2,mm)
                  wfreq2=abp(mm)+(abp(mm+1)-abp(mm))/
     &              (omega(mm+1)-omega(mm))
     &           *(freq2-omega(mm)) !interpolate A_f
                  wertn=afp(m+k)
***********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (wfreq2*
     &             (bose0(freq2)+ferm0(beta*(freq2+exom)))*
     &             wertn+wfreq1*
     &                 (bose0(freq1)+ferm0(beta*(freq1+exom)))*
     &                 wert1)
***********************************************************integration
                  freq1=freq2
                  wfreq1=wfreq2
                  wert1=wertn
               END DO
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (abp(j+1)*
     &          (bose0(omega(j+1))+ferm0(beta*(omega(j+1)+exom)))*
     &              wert2+wfreq1*
     &              (bose0(freq1)+ferm0(beta*(freq1+exom)))*
     &              wert1)
***********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         tmat(i)=-(aux+www)
      END DO
    
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
      SUBROUTINE SUSCEPTIBILITY(beta,afp,omega,n,suszi)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,mm,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION ferm0
      DOUBLE PRECISION omega(n),afp(n),beta
      DOUBLE PRECISION suszi(n),wfreq1,wfreq2
      EXTERNAL ferm0

      DO i=1,n
         exom=omega(i)
         DO j=1,n
            omex=exom+omega(j)  !frequency argument of abm
            CALL locate(omega,n,omex,m)
            iscr(j)=m           !index of epsilon on eps+exom
         END DO
         aux=0.0d0
         m=iscr(1)
         IF(m.EQ.0) THEN
            wert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
            ENDIF
            IF(kk.eq.0) THEN
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (afp(j+1)*
     &          (ferm0(beta*omega(j+1))-ferm0(beta*(omega(j+1)+exom)))*
     &          wert2+afp(j)*
     &          (ferm0(beta*omega(j))-ferm0(beta*(omega(j)+exom)))*
     &              wert1)
***********************************************************integration
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               wfreq1=afp(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  CALL LOCATE(omega,n,freq2,mm)
                  wfreq2=afp(mm)+(afp(mm+1)-afp(mm))/
     &              (omega(mm+1)-omega(mm))
     &           *(freq2-omega(mm)) !interpolate A_f
                  wertn=afp(m+k)
***********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (wfreq2*
     &             (ferm0(beta*freq2)-ferm0(beta*(freq2+exom)))*
     &             wertn+wfreq1*
     &                 (ferm0(beta*freq1)-ferm0(beta*(freq1+exom)))*
     &                 wert1)
***********************************************************integration
                  freq1=freq2
                  wfreq1=wfreq2
                  wert1=wertn
               END DO
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (afp(j+1)*
     &          (ferm0(beta*omega(j+1))-ferm0(beta*(omega(j+1)+exom)))*
     &              wert2+wfreq1*
     &              (ferm0(beta*freq1)-ferm0(beta*(freq1+exom)))*
     &              wert1)
***********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         suszi(i)=aux
      END DO
    
      RETURN
      END
