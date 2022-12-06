CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     LOESEN DES KONDO_MODELL mit Kopplung an bosonischen Bad Aphi im 
C     Grenzfall N,M -> unendlich mit ratio=M/N endlich
C     in the presence of leads kept at finite bias voltage and temperature
C     m indicates the Keldysh component: afm,sefm,...
C     p indicates the retarded-advanced component
C     KRAKRO is g(y)=1/Pi P \int dx f(x)/(y-x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM NCA
      IMPLICIT REAL*8 (a-h,o-y)
      CHARACTER(len=4) check 
      INTEGER nspin,nNeu
      INTEGER N,nadd,naux,ngrid,nnn,mnew,nwv,nvolt !Grid parameters, ngrid=#gridpoints
      INCLUDE 'PARAMETER' 
cc      PARAMETER (naux=N+(N-1)*nadd,ngrid=2*naux)
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER (nspin=2)
      INTEGER i,j,k,l,numb0,number
      INTEGER iterate,maxiterate,kounter
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535898d0)
      DOUBLE PRECISION temp,beta,D,RJ0,RJ,RG,xupdt
      DOUBLE PRECISION Gamma,aphi    !grid parameter
      DOUBLE PRECISION delta,x,xxx,xx0,xtol,yyy
      DOUBLE PRECISION omega(ngrid),gam_arr(n)
      DOUBLE PRECISION sefp(ngrid),sefm(ngrid),afp(ngrid),afm(ngrid),
     &     gfr(ngrid),sefr(ngrid),sefmnew(ngrid),sefpnew(ngrid),
     &     sefmold(ngrid),sefpold(ngrid),afpold(ngrid),afmold(ngrid)
      DOUBLE PRECISION sebp(ngrid),sebm(ngrid),abp(ngrid),abm(ngrid),
     &     gbr(ngrid),sebr(ngrid),sebmnew(ngrid),sebpnew(ngrid),
     &     sebmold(ngrid),sebpold(ngrid),abpold(ngrid),abmold(ngrid)
      DOUBLE PRECISION suszi(ngrid),resus(ngrid),tmat(ngrid) 
      DOUBLE PRECISION scr(ngrid),out
      DOUBLE PRECISION omegan(nNeu),omega1(nnn),omega2(nwv)
      DOUBLE PRECISION templist(100)
      DOUBLE PRECISION V0,T1,T2
******
      INTEGER nBroyden
      PARAMETER (nBroyden=100)
      DOUBLE PRECISION Fmatrix(nBroyden,4*ngrid)
      DOUBLE PRECISION Vmatrix(nBroyden,4*ngrid)
      DOUBLE PRECISION Vout(4*ngrid)
****** 
      EXTERNAL cspec,ferm0,bose0,ferm1,aphi
      CHARACTER*7 handle
      CHARACTER*10 name1,name2,name3,name4
      COMMON /versuch/badgamB,badcutB,badbetaB,pref,xhlpB
**********************READ IN PARAMETERS********************************
      OPEN(UNIT=14,FILE='PARAMETER.dat',STATUS='UNKNOWN')
      READ(14,*) temp0,temp1,numb0,D
      READ(14,*) RJ0,RG,ratio
      READ(14,*) Gamma
      READ(14,*) maxiterate,xtol,xupdt
      READ(14,*) badgamB,badcutB,badbetaB
      READ(14,*) badgam0,badcut0,cut20,badbeta0
      READ(14,*) V0,T1,T2
      CLOSE(14)

      bdgm=badgam0              !this is the DOS exponent
      bdct=badcut0
      bdbet=badbeta0
      cut2=cut20
      xhlpB=badcutB**badgamB
      pref=1.0d0/gammaf(1.0d0+bdgm)
      IF (DABS(RJ0).LT.1.0d-14) RJ0=1.0d-14
*************************************
      OPEN(UNIT=20,FILE='TEMPLIST.dat',STATUS='UNKNOWN')
       kounter=0
      DO i=1,100
         READ(20,*,END=99) temp
         templist(i)=temp
          kounter=kounter+1
      END DO
 99   CONTINUE
      number=kounter
      CLOSE(20)
*************************************
      OPEN(UNIT=97,FILE='STATUS.dat',STATUS='UNKNOWN')

      DO lkm=1,number
      check='nein'
      temp=templist(lkm)
*************************************
      beta=1.0d0/temp
*************************************
      WRITE(handle,'(E7.2)') temp
      name2='ab_'//handle
      name1='af_'//handle
      name3='Kf_'//handle
      name4='Kb_'//handle
      OPEN(UNIT=91,FILE=name1,STATUS='UNKNOWN')
      OPEN(UNIT=92,FILE=name2,STATUS='UNKNOWN')
      OPEN(UNIT=93,FILE=name3,STATUS='UNKNOWN')
      OPEN(UNIT=94,FILE=name4,STATUS='UNKNOWN')
*************************************
      OPEN(UNIT=15,FILE='RESULTS.dat',STATUS='UNKNOWN')
      READ(15,*) rlambd0
      CLOSE(15)
********************build frequency grid omega**************************
      CALL GRID(Gamma,N,Nadd,naux,nNeu,omegan)
      CALL RENORM(D,nNeu,omegan)

      cutoff=bdct
      width=bdct/10.0d0
      CALL XRAPOINTS(omegan,nNeu,new,cutoff,width,omega1,nnn)

      cutoff=max(V0,5.0d-9)
      width=min(2.0*temp,cutoff)
      CALL XRAPOINTS(omega1,nnn,nvolt,cutoff,width,omega2,nwv)

      cutoff=cut2
      width=cut2/5.0d0
      print*,cut2
      print*,width
      CALL XRAPOINTS(omega2,nwv,mnew,cutoff,width,omega,ngrid)



************************************************************************      
*************************initialize functions***************************
      CALL READIN(sefpold,sebpold,sefmold,sebmold,omega,ngrid)

      DO i=1,ngrid
         sefp(i)=sefpold(i)/pi
         sebp(i)=sebpold(i)/pi
         sefm(i)=sefmold(i)/pi
         sebm(i)=sebmold(i)/pi
      END DO

      CALL KRAKRO(ngrid,omega,sefp,sefr) !real part of selfenergies
      CALL DSCAL(ngrid,-1.0d0,sefr,1)
      CALL KRAKRO(ngrid,omega,sebp,sebr) !real part of selfenergies
      CALL DSCAL(ngrid,-1.0d0,sebr,1)
      
      CALL FERMIFUNC0
      CALL BOSEFUNC0(omega,temp)
      DO i=1,ngrid
         sefm(i)=(1.0-2.0*ferm0(beta*omega(i)))*sefm(i)
         sebm(i)=bose0(omega(i))*sebm(i)
         WRITE(71,*) omega(i),sebm(i),sefm(i)
         WRITE(72,*) omega(i),bose0(omega(i))
      END DO
      RG=0.0d0
************************************************************************
C      RJ=RJ0*xxx
      RJ=T1+T2
      IF (DABS(RJ).LT.1.0d-14) RJ=1.0d-14

      DO i=1,ngrid
         afp(i)=sefp(i)/
     &          ((omega(i)-rlambd0+pi*sefr(i))**2+pi**2*sefp(i)**2)
         abp(i)=sebp(i)/
     &          ((-1.0/RJ+pi*sebr(i))**2+pi**2*sebp(i)**2)
C
         afm(i)=sefm(i)/
     &          ((omega(i)-rlambd0+pi*sefr(i))**2+pi**2*sefp(i)**2)
         abm(i)=sebm(i)/
     &          ((-1.0/RJ+pi*sebr(i))**2+pi**2*sebp(i)**2)
      END DO
      
************************************************************************
******************get DoS and Fermi function****************************
      CALL FERMIFUNC0
      CALL FERMFUNCV(omega,temp,T1,T2,V0)
      CALL BOSEFUNC0(omega,temp)
      CALL DOSINTERPOL(omega,ngrid,bdgm,bdct,cut2,bdbet,xxx)
      

      print*,'weight: ',xxx
      
************************************************************************
     
      OPEN(unit=18,FILE='bosBath.dat',STATUS='unknown')
      OPEN(unit=15,FILE='fermBath.dat',STATUS='unknown')
      DO i=1,ngrid
         WRITE(18,*) omega(i),aphi(omega(i))
         WRITE(15,*) omega(i),cspec(omega(i))
         WRITE(17,*) omega(i),afp(i),afm(i)
c         WRITE(16,*) omega(i),abp(i),abm(i)
         WRITE(27,*) omega(i),sebp(i)
         WRITE(26,*) omega(i),sefp(i)
         WRITE(28,*) omega(i),ferm1(omega(i))
      END DO
      CLOSE(18)
      CLOSE(15)
      
************************************************************************
************************************************************************
**********************START ITERATION***********************************
************************************************************************
************************************************************************
C      OPEN(UNIT=98,FILE='STATSUS.dat',STATUS='UNKNOWN')
C      OPEN(UNIT=99,FILE='TMATRIX.dat',STATUS='UNKNOWN')
 
      DO iterate=1,maxiterate
         
         CALL SIGMF1m(beta,afp,afm,omega,ngrid,scr)
         CALL SIGMF2m(beta,abp,abm,omega,ngrid,sefpnew)
         CALL SIGMBm(beta,afp,afm,omega,ngrid,sebpnew)

      DO i=1,ngrid
         WRITE(79,*) omega(i),sefpnew(i)
         WRITE(78,*) omega(i),sebpnew(i)
         WRITE(77,*) omega(i),scr(i)
      END DO
     
      

      DO i=1,ngrid
         sefpnew(i)=ratio*sefpnew(i)+RG**2*scr(i)
      END DO

      CALL SIGMF1p(beta,afp,afm,omega,ngrid,scr)
      CALL SIGMF2p(beta,abp,abm,omega,ngrid,sefmnew)
      CALL SIGMBp(beta,afp,afm,omega,ngrid,sebmnew)

      DO i=1,ngrid
         WRITE(81,*) omega(i),sefmnew(i)
         WRITE(82,*) omega(i),sebmnew(i)
         WRITE(83,*) omega(i),scr(i)
      END DO
      STOP

      DO i=1,ngrid
         sefmnew(i)=ratio*sefmnew(i)+RG**2*scr(i)
      END DO

***********************update selfenergies******************************
c         DO i=1,ngrid
c            sefpold(i)=sefp(i)
c            sefp(i)=(xupdt*sefpnew(i)+sefpold(i))/(1.0d0+xupdt)
c            sebpold(i)=sebp(i)
c            sebp(i)=(xupdt*sebpnew(i)+sebpold(i))/(1.0d0+xupdt)
c         END D
         IF(iterate.LE.nBroyden) THEN
            kount=iterate
            DO i=1,ngrid
               Vmatrix(kount,i)=sefp(i)
               Vmatrix(kount,i+ngrid)=sebp(i)
               Vmatrix(kount,i+2*ngrid)=sefm(i)
               Vmatrix(kount,i+3*ngrid)=sebm(i)
               Fmatrix(kount,i)=sefpnew(i)-sefp(i)
               Fmatrix(kount,i+ngrid)=sebpnew(i)-sebp(i)
               Fmatrix(kount,i+2*ngrid)=sefmnew(i)-sefm(i)
               Fmatrix(kount,i+3*ngrid)=sebmnew(i)-sebm(i)
            END DO
         ELSE
            kount=nBroyden
            DO l=2,kount
               DO i=1,ngrid
                  Vmatrix(l-1,i)=Vmatrix(l,i)
                  Vmatrix(l-1,i+ngrid)=Vmatrix(l,i+ngrid)
                  Vmatrix(l-1,i+2*ngrid)=Vmatrix(l,i+2*ngrid)
                  Vmatrix(l-1,i+3*ngrid)=Vmatrix(l,i+3*ngrid)
                  Fmatrix(l-1,i)=Fmatrix(l,i)
                  Fmatrix(l-1,i+ngrid)=Fmatrix(l,i+ngrid)
                  Fmatrix(l-1,i+2*ngrid)=Fmatrix(l,i+2*ngrid)
                  Fmatrix(l-1,i+3*ngrid)=Fmatrix(l,i+3*ngrid)
               END DO

            END DO
            DO i=1,ngrid
               Vmatrix(kount,i)=sefp(i)
               Vmatrix(kount,i+ngrid)=sebp(i)
               Vmatrix(kount,i+2*ngrid)=sefm(i)
               Vmatrix(kount,i+3*ngrid)=sebm(i)
               Fmatrix(kount,i)=sefpnew(i)-sefp(i)
               Fmatrix(kount,i+ngrid)=sebpnew(i)-sebp(i)
               Fmatrix(kount,i+2*ngrid)=sefmnew(i)-sefm(i)
               Fmatrix(kount,i+3*ngrid)=sebmnew(i)-sebm(i)
            END DO
         ENDIF


             
         CALL  BROYDEN(kount,4*ngrid,nBroyden,Fmatrix,Vmatrix,Vout)


         DO i=1,ngrid
            sefp(i)=Vout(i)
            sebp(i)=Vout(i+ngrid)
            sefm(i)=Vout(i+2*ngrid)
            sebm(i)=Vout(i+3*ngrid)
         END DO
         
*******************************************
         CALL KRAKRO(ngrid,omega,sefp,sefr) !real part of selfenergies
         CALL DSCAL(ngrid,-1.0d0,sefr,1)
         CALL KRAKRO(ngrid,omega,sebp,sebr) !real part of selfenergies
         CALL DSCAL(ngrid,-1.0d0,sebr,1)
******************calculate rlambda0 update*****************************
         rlambold=rlambd0

ctest         CALL RUPDATE(sefp,sefr,omega,beta,rlambd0,delta)
        rlambd0=0.0d0

c         CALL SHIFT(delta,omega,ngrid,sefp,sefm,sefr,sebp,sebm,sebr)
        
************************************************************************
****************calculate new spectral functions************************
         DO i=1,ngrid
            afpold(i)=afp(i)
            afmold(i)=afm(i)
            afp(i)=sefp(i)/
     &           ((omega(i)-rlambd0+pi*sefr(i))**2+pi**2*sefp(i)**2)
            afm(i)=sefm(i)/
     &           ((omega(i)-rlambd0+pi*sefr(i))**2+pi**2*sefp(i)**2)
            abpold(i)=abp(i)
            abmold(i)=abm(i)
            abp(i)=sebp(i)/
     &          ((-1.0/RJ+pi*sebr(i))**2+pi**2*sebp(i)**2)
            abm(i)=sebm(i)/
     &          ((-1.0/RJ+pi*sebr(i))**2+pi**2*sebp(i)**2)
         END DO

         DO i=1,ngrid
            WRITE(13,*) omega(i),afm(i),afmold(i)
            WRITE(14,*) omega(i),afp(i),afpold(i)
            WRITE(18,*) omega(i),abp(i),abpold(i)
            WRITE(16,*) omega(i),abm(i),abmold(i)
         END DO
         

         CALL ROUT('afp.dat',ngrid,omega,afp)
         CALL ROUT('abp.dat',ngrid,omega,abp)

************************************************************************
         CALL PNNESS(ngrid,omega,afm,afp,rnf)
         print*,'occupation: ',2.0*rnf,' lambda: ',rlambd0
         
********************CHECK for CONVERGENCE*******************************
         xx0=0.0d0
         xx1=0.0d0
         xx2=0.0d0
         xx3=0.0d0
         i0=0
         DO i=1,ngrid
            xxx=dabs(afp(i)-afpold(i))/(dabs(afp(i))+1.0d-70)
            IF(xxx.GT.xx0) THEN
            xx0=xxx
            i0=i
            out=omega(i0)
            ENDIF
            xxx=dabs(abp(i)-abpold(i))/(dabs(abp(i))+1.0d-70)
            IF(xxx.GT.xx1) THEN
            xx1=xxx
            i0=ngrid+i
            out=omega(i)
            ENDIF
            xxx=dabs(afm(i)-afmold(i))/(dabs(afm(i))+1.0d-70)
            IF(xxx.GT.xx2) THEN
            xx2=xxx
            i0=i
            out=omega(i0)
            ENDIF
            xxx=dabs(abm(i)-abmold(i))/(dabs(abm(i))+1.0d-70)
            IF(xxx.GT.xx3) THEN
            xx3=xxx
            i0=ngrid+i
            out=omega(i)
            ENDIF
         END DO
         print*,'deltax: ',xx0,xx1,xx2,xx3,i0,out
         print*,'Iteration: ',iterate
         
         IF(xx0.LT.xtol) THEN
         check='ja'
         GOTO 300
         ENDIF
      END DO
 300  CONTINUE

      write(88,*) 'converged',iterate,xxx,xtol


************************************************************************
         IF(check.EQ.'ja') THEN
            OPEN(unit=51,FILE='sefp.dat'//handle,STATUS='unknown')
            OPEN(unit=52,FILE='sebp.dat'//handle,STATUS='unknown')
            OPEN(unit=53,FILE='sefm.dat'//handle,STATUS='unknown')
            OPEN(unit=54,FILE='sebm.dat'//handle,STATUS='unknown')
           DO i=1,ngrid
                 WRITE(91,*) omega(i),afp(i),temp
                 WRITE(92,*) omega(i),abp(i),temp
                 WRITE(93,*) omega(i),afm(i),temp
                 WRITE(94,*) omega(i),abm(i),temp
                 WRITE(51,*) omega(i),sefp(i),temp
                 WRITE(52,*) omega(i),sebp(i),temp
                 WRITE(53,*) omega(i),sefm(i),temp
                 WRITE(54,*) omega(i),sebm(i),temp
           END DO

           CLOSE(51)
           CLOSE(52)
           CLOSE(53)
           CLOSE(54)
           CLOSE(91)
           CLOSE(92)
           CLOSE(93)
           CLOSE(94)

            CALL ROUT('SEFP0.dat',ngrid,omega,sefp)
            CALL ROUT('SEBP0.dat',ngrid,omega,sebp)
            OPEN(UNIT=15,FILE='RESULTS0.dat',STATUS='UNKNOWN')
            WRITE(15,*) rlambd0
            CLOSE(15)
         END IF
      CALL ROUT('SEFP.dat',ngrid,omega,sefp)
      CALL ROUT('SEBP.dat',ngrid,omega,sebp)
**************************T-MATRIX**************************************
      CALL tmatrix(omega,ngrid,beta,afp,abp,tmat)
      CALL DSCAL(ngrid,1.0d0/pi,suszi,1)
      OPEN(UNIT=15,FILE='tmat.dat'//handle,STATUS='UNKNOWN')
      DO i=1,ngrid
         WRITE(15,*) omega(i),tmat(i)
      END DO
      CLOSE(15)
      yyy=tmat(m)+(tmat(m+1)-tmat(m))/(omega(m+1)-omega(m))
     &     *(0.0d0-omega(m))
      WRITE(99,*) temp,yyy
************************************************************************
      OPEN(UNIT=15,FILE='RESULTS.dat',STATUS='UNKNOWN')
      WRITE(15,*) rlambd0
      CLOSE(15)

******************************
      write(97,*) temp,' ',check
      END DO
      CLOSE(97)
      CLOSE(98)
      CLOSE(99)
******************************
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
