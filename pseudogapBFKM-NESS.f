CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     LOESEN DES KONDO_MODELL mit Kopplung an bosonischen Bad Aphi im 
C     Grenzfall N,M -> unendlich mit ratio=M/N endlich
C     in the presence of leads kept at finite bias voltage and temperature
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
         sefp(i)=sefpold(i)
         sebp(i)=sebpold(i)
         sefm(i)=sefmold(i)
         sebm(i)=sebmold(i)
      END DO

      CALL KRAKRO(ngrid,omega,sefp,sefr)
      CALL KRAKRO(ngrid,omega,sebp,sebr)


************************************************************************
      RJ=T1+T2
      IF (DABS(RJ).LT.1.0d-14) RJ=1.0d-14

      DO i=1,ngrid
         afp(i)=sefp(i)/
     &          ((omega(i)+rlambd0-sefr(i))**2+(sefp(i))**2)
         abp(i)=sebp(i)/
     &          ((1.0/RJ-sebr(i))**2+(sebp(i))**2)
C
         afm(i)=sefm(i)/
     &          ((omega(i)+rlambd0-sefr(i))**2+(sefp(i))**2)
         abm(i)=sebm(i)/
     &          ((1.0/RJ-sebr(i))**2+(sebp(i))**2)
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
      END DO
      CLOSE(18)
      CLOSE(15)
      
************************************************************************
************************************************************************
**********************START ITERATION***********************************
************************************************************************
************************************************************************
 
      DO iterate=1,maxiterate

         CALL SIGMF1m(beta,afp,afm,omega,ngrid,scr)
         CALL SIGMF2m(beta,abp,abm,omega,ngrid,sefpnew)
         CALL SIGMBm(beta,afp,afm,omega,ngrid,sebpnew)
         CALL DSCAL(ngrid,1.0d0/pi,scr,1)
         CALL DSCAL(ngrid,1.0d0/pi,sefpnew,1)
         CALL DSCAL(ngrid,1.0d0/pi,sebpnew,1)

 
      

      DO i=1,ngrid
         sefpnew(i)=ratio*sefpnew(i)+RG**2*scr(i)
      END DO

      CALL SIGMF1p(beta,afp,afm,omega,ngrid,scr)
      CALL SIGMF2p(beta,abp,abm,omega,ngrid,sefmnew)
      CALL SIGMBp(beta,afp,afm,omega,ngrid,sebmnew)
      CALL DSCAL(ngrid,1.0d0/pi,scr,1)
      CALL DSCAL(ngrid,1.0d0/pi,sefmnew,1)
      CALL DSCAL(ngrid,1.0d0/pi,sebmnew,1)

      DO i=1,ngrid
         sefmnew(i)=ratio*sefmnew(i)+RG**2*scr(i)
      END DO
      
     


         
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
         CALL KRAKRO(ngrid,omega,sefp,sefr)
         CALL KRAKRO(ngrid,omega,sebp,sebr)
******************calculate rlambda0 update*****************************
         rlambold=rlambd0

        rlambd0=0.0d0
        
************************************************************************
****************calculate new spectral functions************************
         DO i=1,ngrid
            afpold(i)=afp(i)
            afmold(i)=afm(i)
            afp(i)=sefp(i)/
     &           ((omega(i)+rlambd0-sefr(i))**2+(sefp(i))**2)
            afm(i)=sefm(i)/
     &           ((omega(i)+rlambd0-sefr(i))**2+(sefp(i))**2)
            abpold(i)=abp(i)
            abmold(i)=abm(i)
            abp(i)=sebp(i)/
     &           ((1.0/RJ-sebr(i))**2+(sebp(i))**2)
            abm(i)=sebm(i)/   
     &           ((1.0/RJ-sebr(i))**2+(sebp(i))**2)
         END DO

        

         CALL ROUT('afp.dat',ngrid,omega,afp)
         CALL ROUT('abp.dat',ngrid,omega,abp)

        
************************************************************************
         CALL PNNESS(ngrid,omega,afm,rnf)
         print*,'occupation: ',2.0*rnf,' lambda: ',rlambd0
         
********************CHECK for CONVERGENCE*******************************
         xx0=0.0d0
         xx1=0.0d0
         xx2=0.0d0
         xx3=0.0d0
         i0=0
         i1=0
         i2=0
         i3=0
         DO i=1,ngrid
            xxx=dabs(afp(i)-afpold(i))/(dabs(afp(i))+1.0d-70)
            IF(xxx.GT.xx0) THEN
            xx0=xxx
            i0=i
            out=omega(i)
            ENDIF
            xxx=dabs(abp(i)-abpold(i))/(dabs(abp(i))+1.0d-70)
            IF(xxx.GT.xx1) THEN
            xx1=xxx
            i1=ngrid+i
            out=omega(i)
            ENDIF
            xxx=dabs(afm(i)-afmold(i))/(dabs(afm(i))+1.0d-70)
            IF(xxx.GT.xx2) THEN
            xx2=xxx
            i2=2*ngrid+i
            out=omega(i)
            ENDIF
            xxx=dabs(abm(i)-abmold(i))/(dabs(abm(i))+1.0d-70)
            IF(xxx.GT.xx3) THEN
            xx3=xxx
            i3=3*ngrid+i
            out=omega(i)
            ENDIF
         END DO
         xx0=max(xx0,xx1)
         xx0=max(xx0,xx2)
         xx0=max(xx0,xx3)
         print*,'deltax: ',xx0,xx1,xx2,xx3
         print*,i0,i1,i2,i3
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
      CALL ROUT('SEFM.dat',ngrid,omega,sefm)
      CALL ROUT('SEBM.dat',ngrid,omega,sebm)
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
      
