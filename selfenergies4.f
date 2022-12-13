************************************************************************
*****************************SIGMAF*************************************
************************************************************************
      SUBROUTINE SIGMF1m(beta,afp,afm,omega,n,sefpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n),next
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION vert1,vert2,vertn
      DOUBLE PRECISION xxx,yyy,www,D
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION aphi,ferm1,bose0
      DOUBLE PRECISION omega(n),afp(n),afm(n),beta
      DOUBLE PRECISION sefpnew(n)
      EXTERNAL aphi,bose0,ferm1


      D=omega(n)+1.0d-10
c      DO i=1,n
         DO i=1,n/2
         exom=omega(i)
*****************Integrand bei Null (ohne Bose+Fermi-Funktion)
c         CALL LOCATE(omega,n,0.0d0,m)
c         xxx=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
c     &           *(0.0d0-omega(m))
         xxx=afp(i)*aphi(0.0d0)
***Integral (bose+fermi)
         yyy=exom-1.0/beta*dlog(1.0d0+dexp(beta*(exom-D)))
     &        +1.0/beta*dlog(1.0d0+dexp(-beta*(exom+D)))

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
            vert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
            vert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
            vert1=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_bm
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               vert2=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_bm
            ENDIF
            IF(kk.eq.0) THEN
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (aphi(omega(j+1))*vert2-bose0(omega(j+1))
     &              *aphi(omega(j+1))*wert2+
     &              aphi(omega(j))*vert1-bose0(omega(j))
     &              *aphi(omega(j))*wert1)
**********************************************************integration
cc               write(86,*) omega(j+1),wert2
            ELSE
cc            IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
cc                  write(86,*) omega(m+k)-exom,afp(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afp(m+k)
                  vertn=afm(m+k)
**********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (aphi(freq2)*vertn-bose0(freq2)
     &                 *aphi(freq2)*wertn+
     &                 aphi(freq1)*vert1-bose0(freq1)
     &                 *aphi(freq1)*wert1)
**********************************************************integration
                  freq1=freq2
                  wert1=wertn
                  vert1=vertn
               END DO
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (aphi(omega(j+1))*vert2-bose0(omega(j+1))
     &              *aphi(omega(j+1))*wert2+
     &              aphi(freq1)*vert1-bose0(freq1)
     &              *aphi(freq1)*wert1)
**********************************************************integration
            ENDIF
 30         CONTINUE
            wert1=wert2
            vert1=vert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
c         sefpnew(i)=-(aux+www)/2.0
         sefpnew(i)=aux/2.0
      END DO
      DO i=1,n/2
         sefpnew(n+1-i)=sefpnew(i)
      END DO


      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC SIGMFb CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SIGMF2m(beta,abp,abm,omega,n,sefpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION vert1,vert2,vertn
      DOUBLE PRECISION xxx,yyy,www,D
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION cspec,ferm1,bose0
      DOUBLE PRECISION omega(n),abp(n),abm(n),beta
      DOUBLE PRECISION sefpnew(n)
      EXTERNAL cspec,bose0,ferm1
      
        
      D=omega(n)+1.0d-10
c      DO i=1,n
       DO i=1,n/2
         exom=omega(i)
*****************Integrand bei Null (ohne Bose+Fermi-Funktion)
         CALL LOCATE(omega,n,0.0d0,m)
         xxx=abp(m)+(abp(m+1)-abp(m))/(omega(m+1)-omega(m))
     &           *(0.0d0-omega(m))
         xxx=xxx*cspec(exom)
***Integral (bose+fermi)
         yyy=exom-1.0/beta*dlog(1.0d0-dexp(beta*(exom-D)))
     &        +1.0/beta*dlog(1.0d0-dexp(-beta*(exom+D)))

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
            vert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
            vert1=0.0d0
         ELSE
            wert1=abp(m)+(abp(m+1)-abp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
            vert1=abm(m)+(abm(m+1)-abm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_bm
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSE
               wert2=abp(l)+(abp(l+1)-abp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               vert2=abm(l)+(abm(l+1)-abm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_bm
            ENDIF
            IF(kk.eq.0) THEN
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (cspec(-omega(j+1))*vert2+
     &              ferm1(-omega(j+1))*cspec(-omega(j+1))*wert2+
     &              cspec(-omega(j))*vert1+
     &              ferm1(-omega(j))*cspec(-omega(j))*wert1)
**********************************************************integration
cc               write(86,*) omega(j+1),wert2
            ELSE
cc            IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
cc                  write(86,*) omega(m+k)-exom,afp(m+k)
                  freq2=omega(m+k)-exom
                  wertn=abp(m+k)
                  vertn=abm(m+k)
**********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (cspec(-freq2)*vertn+
     &                 ferm1(-freq2)*cspec(-freq2)*wertn+
     &                 cspec(-freq1)*vert1+ferm1(-freq1)
     &                 *cspec(-freq1)*wert1)
**********************************************************integration
                  freq1=freq2
                  wert1=wertn
                  vert1=vertn
               END DO
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (cspec(-omega(j+1))*vert2+
     &              ferm1(-omega(j+1))*cspec(-omega(j+1))*wert2+
     &              cspec(-freq1)*vert1+ferm1(-freq1)
     &              *cspec(-freq1)*wert1)
**********************************************************integration
            ENDIF
 30         CONTINUE
            wert1=wert2
            vert1=vert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
c         sefpnew(i)=(aux+www)/2.0
         sefpnew(i)=-aux/2.0d0
      END DO
        DO i=1,n/2
         sefpnew(n+1-i)=sefpnew(i)
      END DO


      RETURN
      END



************************************************************************
*****************************SIGMAB*************************************
************************************************************************
      SUBROUTINE SIGMBm(beta,afp,afm,omega,n,sebpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION vert1,vert2,vertn
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION cspec,ferm1
      DOUBLE PRECISION omega(n),afp(n),afm(n),beta
      DOUBLE PRECISION sebpnew(n)
      EXTERNAL cspec,ferm1
      
      
      DO i=1,n/2
cc         exom=-omega(i)
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
            vert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
            vert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
            vert1=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_bm
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               vert2=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_bm
            ENDIF
            IF(kk.eq.0) THEN
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (cspec(omega(j+1))*vert2-
     &          ferm1(omega(j+1))*cspec(omega(j+1))*wert2+
     &          cspec(omega(j))*vert1-
     &          ferm1(omega(j))*cspec(omega(j))*wert1)
***********************************************************integration
!               write(87,*) omega(j+1),vert2,wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               DO k=1,kk
!                  write(86,*) omega(m+k)-exom,afm(m+k),afp(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afp(m+k)
                  vertn=afm(m+k)
***********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (cspec(freq2)*vertn-
     &             ferm1(freq2)*cspec(freq2)*wertn+
     &             cspec(freq1)*vert1-
     &             ferm1(freq1)*cspec(freq1)*wert1)
***********************************************************integration
                  freq1=freq2
                  wert1=wertn
                  vert1=vertn
               END DO
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (cspec(omega(j+1))*vert2-
     &              ferm1(omega(j+1))*cspec(omega(j+1))*wert2+
     &              cspec(freq1)*vert1-
     &              ferm1(freq1)*cspec(freq1)*wert1)
***********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
            vert1=vert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sebpnew(i)=-aux/2.0
      END DO
      DO i=1,n/2
         sebpnew(n+1-i)=-sebpnew(i)
      END DO
      

      RETURN
      END
*********************************************************************************
*********************************************************************************
*********************************************************************************
*****************************SIGMAFP1********************************************
*********************************************************************************
      SUBROUTINE SIGMF1p(beta,afp,afm,omega,n,sefpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n),next
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION vert1,vert2,vertn
      DOUBLE PRECISION xxx,yyy,www,D
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION aphi,ferm1,bose0
      DOUBLE PRECISION omega(n),afp(n),afm(n),beta
      DOUBLE PRECISION sefpnew(n)
      EXTERNAL aphi,bose0,ferm1

c      CALL DZERO(afm,n)
      D=omega(n)+1.0d-10
cc      DO i=1,n
         DO i=1,n/2
         exom=omega(i)
*****************Integrand bei Null (ohne Bose+Fermi-Funktion)
c         CALL LOCATE(omega,n,0.0d0,m)
c         xxx=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
c     &           *(0.0d0-omega(m))
         xxx=afp(i)*aphi(0.0d0)
***Integral (bose+fermi)
         yyy=exom-1.0/beta*dlog(1.0d0+dexp(beta*(exom-D)))
     &        +1.0/beta*dlog(1.0d0+dexp(-beta*(exom+D)))

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
            vert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
            vert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
            vert1=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_bm
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
            m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               vert2=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_bm
            ENDIF
            IF(kk.eq.0) THEN
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (bose0(omega(j+1))*aphi(omega(j+1))*vert2-
     &              aphi(omega(j))*wert2+
     &              bose0(omega(j))*aphi(omega(j))*vert1-
     &              aphi(omega(j))*wert1)
**********************************************************integration
cc               write(86,*) omega(j+1),wert2
            ELSE
cc            IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
cc                  write(86,*) omega(m+k)-exom,afp(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afp(m+k)
                  vertn=afm(m+k)
**********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (bose0(freq2)*aphi(freq2)*vertn-
     &                 aphi(freq2)*wertn+
     &                 bose0(freq1)*aphi(freq1)*vert1-
     &                 aphi(freq1)*wert1)
**********************************************************integration
                  freq1=freq2
                  wert1=wertn
                  vert1=vertn
               END DO
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (bose0(omega(j+1))*aphi(omega(j+1))*vert2-
     &              aphi(omega(j+1))*wert2+
     &              bose0(freq1)*aphi(freq1)*vert1-
     &              aphi(freq1)*wert1)
**********************************************************integration
            ENDIF
 30         CONTINUE
            wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sefpnew(i)=-(aux+www)/2.0
      END DO
      DO i=1,n/2
         sefpnew(n+1-i)=-sefpnew(i)
      END DO



      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC SIGMAFp2 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SIGMF2p(beta,abp,abm,omega,n,sefpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION vert1,vert2,vertn
      DOUBLE PRECISION xxx,yyy,www,D
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION cspec,ferm1,bose0
      DOUBLE PRECISION omega(n),abp(n),abm(n),beta
      DOUBLE PRECISION sefpnew(n)
      EXTERNAL cspec,bose0,ferm1
      
        
      D=omega(n)+1.0d-10
cc      DO i=1,n
       DO i=1,n/2
         exom=omega(i)
*****************Integrand bei Null (ohne Bose+Fermi-Funktion)
         CALL LOCATE(omega,n,0.0d0,m)
         xxx=abp(m)+(abp(m+1)-abp(m))/(omega(m+1)-omega(m))
     &           *(0.0d0-omega(m))
         xxx=xxx*cspec(exom)
***Integral (bose+fermi)
         yyy=exom-1.0/beta*dlog(1.0d0-dexp(beta*(exom-D)))
     &        +1.0/beta*dlog(1.0d0-dexp(-beta*(exom+D)))

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
            vert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
            vert1=0.0d0
         ELSE
            wert1=abp(m)+(abp(m+1)-abp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
            vert1=abm(m)+(abm(m+1)-abm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_bm
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSE
               wert2=abp(l)+(abp(l+1)-abp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               vert2=abm(l)+(abm(l+1)-abm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_bm
            ENDIF
            IF(kk.eq.0) THEN
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (ferm1(-omega(j+1))*cspec(-omega(j+1))*vert2+
     &              cspec(-omega(j+1))*wert2+
     &              ferm1(-omega(j))*cspec(-omega(j))*vert1+
     &              cspec(-omega(j))*wert1)
**********************************************************integration
cc               write(86,*) omega(j+1),wert2
            ELSE
cc            IF((m.eq.0).OR.(m.eq.n)) GOTO 30
               freq1=omega(j)
               DO k=1,kk
cc                  write(86,*) omega(m+k)-exom,afp(m+k)
                  freq2=omega(m+k)-exom
                  wertn=abp(m+k)
                  vertn=abm(m+k)
**********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (ferm1(-freq2)*cspec(freq2)*vertn+
     &                 cspec(-freq2)*wertn+
     &                 ferm1(-freq1)*cspec(-freq1)*vert1+
     &                 cspec(-freq1)*wert1)
**********************************************************integration
                  freq1=freq2
                  wert1=wertn
                  vert1=vertn
               END DO
**********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &              (ferm1(-omega(j+1))*cspec(-omega(j+1))*vert2+
     &              cspec(-omega(j+1))*wert2+
     &              ferm1(-freq1)*cspec(-freq1)*vert1+
     &              cspec(-freq1)*wert1)
**********************************************************integration
            ENDIF
 30         CONTINUE
            wert1=wert2
            vert1=vert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sefpnew(i)=-(aux+www)/2.0
      END DO
        DO i=1,n/2
         sefpnew(n+1-i)=-sefpnew(i)
      END DO

      RETURN
      END



*************************************************************************
*****************************SIGMABp*************************************
*************************************************************************
      SUBROUTINE SIGMBp(beta,afp,afm,omega,n,sebpnew)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,n,m,iscr(n)
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION vert1,vert2,vertn
      DOUBLE PRECISION freq1,freq2,wertn
      DOUBLE PRECISION exom,omex
      DOUBLE PRECISION cspec,ferm1
      DOUBLE PRECISION omega(n),afp(n),afm(n),beta
      DOUBLE PRECISION sebpnew(n)
      EXTERNAL cspec,ferm1
      
      
      DO i=1,n/2
cc         exom=-omega(i)
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
            vert1=0.0d0
         ELSEIF(m.EQ.n) THEN
            wert1=0.0d0
            vert1=0.0d0
         ELSE
            wert1=afp(m)+(afp(m+1)-afp(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_b
             vert1=afm(m)+(afm(m+1)-afm(m))/(omega(m+1)-omega(m))
     &           *(omega(1)+exom-omega(m)) !interpolate A_bm
         ENDIF
         DO j=1,n-1             !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
            l=iscr(j+1)
             m=iscr(j)
            kk=l-m
            IF(l.EQ.0) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSEIF(l.EQ.n) THEN
               wert2=0.0d0
               vert2=0.0d0
            ELSE
               wert2=afp(l)+(afp(l+1)-afp(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_b
               vert2=afm(l)+(afm(l+1)-afm(l))/(omega(l+1)-omega(l))
     &              *(omega(j+1)+exom-omega(l)) !interpolate A_bm
            ENDIF
            IF(kk.eq.0) THEN
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-omega(j))*
     &              (ferm1(omega(j+1))*cspec(omega(j+1))*vert2-
     &              cspec(omega(j+1))*wert2+
     &              ferm1(omega(j))*cspec(omega(j))*vert1-
     &              cspec(omega(j))*wert1)
***********************************************************integration
c               write(86,*) omega(j),wert2
            ELSE
c            IF((m.eq.0).OR.(m.eq.n)) GOTO 20
               freq1=omega(j)
               DO k=1,kk
c                  write(86,*) omega(m+k)-exom,abm(m+k)
                  freq2=omega(m+k)-exom
                  wertn=afp(m+k)
                  vertn=afm(m+k)
***********************************************************integration
                  aux=aux+0.5*(freq2-freq1)*
     &                 (ferm1(freq2)*cspec(freq2)*vertn-
     &                 cspec(freq2)*wertn+
     &                 ferm1(freq1)*cspec(freq1)*vert1-
     &                 cspec(freq1)*wert1)
***********************************************************integration
                  freq1=freq2
                  wert1=wertn
                  vert1=vertn
               END DO
***********************************************************integration
               aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &             (ferm1(omega(j+1))*cspec(omega(j+1))*vert2-
     &              cspec(omega(j+1))*wert2+
     &              ferm1(freq1)*cspec(freq1)*vert1-
     &              cspec(freq1)*wert1)
***********************************************************integration
            ENDIF
 20         CONTINUE
            wert1=wert2
            vert1=vert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
         END DO
         sebpnew(i)=aux/2.0
      END DO
      DO i=1,n/2
         sebpnew(n+1-i)=sebpnew(i)
      END DO


      RETURN
      END



