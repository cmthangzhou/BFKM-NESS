***********************************************************************
***********************************************************************
      SUBROUTINE GRID(Gamma,N,Nadd,naux,ngrid,gridvec)
      IMPLICIT NONE
      INTEGER N,Nadd,ngrid,naux
      INTEGER i,j
      DOUBLE PRECISION D,GAMMA
      DOUBLE PRECISION aux,delta,VECTOR(1:naux)
      DOUBLE PRECISION gam_arr(N) !contains GAMMA**(-i)
      DOUBLE PRECISION gridvec(1:ngrid)

      DO i=1,naux
      vector(i)=0.0d0
      END DO
      
      aux=Gamma**(-N+1)
      DO i=1,N
         IF(i.ne.N) THEN
            DO j=0,Nadd
               delta=aux*(Gamma-1.0)/DFLOAT(Nadd+1)
               IF(j.eq.0) THEN
                  vector((i-1)*(Nadd+1)+1)=aux
               ELSE
                  vector((i-1)*(Nadd+1)+1+j)=
     &              vector((i-1)*(Nadd+1)+j)+delta
               ENDIF
            END DO
            aux=aux*Gamma
         ELSE
            vector(naux)=aux
         END IF
      END DO
      
      DO i=1,naux
         gridvec(i)=-vector(naux+1-i)
         gridvec(i+naux)=vector(i)
      END DO

      aux=1.0d0
      DO i=1,N
         aux=aux*Gamma
         gam_arr(i)=aux
      END DO

      RETURN
      END
**************neues GRID**************
***********************************************************************
***********************************************************************
***********************************************************************
****add 2*new points to grid at lambda in a window of width width
***********************************************************************
      SUBROUTINE XRAPOINTS(gridVec,ngrid,new,lambda,width,omega,nNeu)
      IMPLICIT NONE
      INTEGER nOld,nNeu,i,j,m,m1,m2,new,mRechts,mLinks
      INTEGER N,nAdd,naux,ngrid
c      INCLUDE 'PARAMETER'       !should define parameter new
c      PARAMETER (naux=N+(N-1)*nadd,ngrid=2*naux)
c      PARAMETER(nNeu=ngrid+4*new)
      DOUBLE PRECISION gridVec(ngrid),lambda,width,omega(nNeu)
      DOUBLE PRECISION xxx,xxxR,xxxL,tst,wst1,wst2
      
      lambda=dabs(lambda) !to ensure symmetric grid
      IF (lambda-width.LT.0.0d0) THEN
         PRINT*, 'readjust width'
         STOP
      ENDIF
      xxx=lambda-width
      CALL locate(gridVec,ngrid,xxx,m1)
      xxx=lambda+2.0*width
      CALL locate(gridVec,ngrid,xxx,m2)
      
      m=m2-m1
      mRechts=m/2+new
      mLinks=m-m/2+new

      DO i=1,m1-ngrid/2
         omega(ngrid/2+2*new+i)=gridvec(ngrid/2+i)
      END DO

      tst=(width/100.0d0)
      wst1=datan((width)/1.0d0)
      wst2=datan((2.0*width)/1.0d0)
      xxxR=(wst2-tst)/(mRechts*1.0d0)
      xxxL=(wst1-tst)/(mLinks*1.0d0)

      DO i=1,mLinks
         omega(m1+2*new+i)=-dtan(tst+xxxL*(mLinks-i))+lambda
      END DO
      DO i=1,mRechts
         omega(2*new+mLinks+m1+i)=dtan(tst+xxxR*(i-1))+lambda
      END DO
      DO i=m2+1,ngrid
         omega(i+4*new)=gridvec(i)
      END DO

      DO i=1,ngrid/2+2*new
         omega(i)=-omega(ngrid+4*new-i+1)
      END DO

c      DO i=1,ngrid+4*new
c         write(70,*) i,omega(i)
c      END DO
      RETURN
      END
***********************************************************************
***********************************************************************
      SUBROUTINE RENORM(D,ngrid,freq)
      IMPLICIT NONE
      INTEGER ngrid
      INTEGER i
      DOUBLE PRECISION D,freq(1:ngrid)

      DO i=1,ngrid
         freq(i)=D*freq(i)
      END DO
      END
***********************************************************************
***********************************************************************
      DOUBLE PRECISION FUNCTION cspec(xxx)
      IMPLICIT NONE
      INTEGER N,naux,nadd,nNeu,nnn,new,ngrid,nspin,mnew,m,nvolt,nwv
      INCLUDE 'PARAMETER'
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER (nspin=2)
      DOUBLE PRECISION ferm0,xxx,r
      DOUBLE PRECISION pi,weight,aux,beta,cut
      DOUBLE PRECISION dos(ngrid),xw(ngrid)
      parameter (pi=3.1415926535898d0) 
      EXTERNAL ferm0
      COMMON /doskram/dos,xw
      
      CALL LOCATE(xw,ngrid,xxx,m)
      IF(m.LT.1) THEN
         cspec=0.0d0
      ELSEIF(m.GE.ngrid) THEN
         cspec=0.0d0
      ELSE
         cspec=dos(m)+(dos(m+1)-dos(m))/
     &        (xw(m+1)-xw(m))*(xxx-xw(m))
      ENDIF

      RETURN
      END
***********************************************************************
***********************************************************************
      SUBROUTINE DOSINTERPOL(omega,m,r0,aux0,cut0,beta0,xxx)
      IMPLICIT NONE
      INTEGER i,m1,m2,m
      INTEGER N,naux,nadd,nNeu,nnn,new,ngrid,nspin,mnew,nwv,nvolt
      INCLUDE 'PARAMETER'
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER (nspin=2)
      DOUBLE PRECISION pi,xxx
      DOUBLE PRECISION ferm0,beta,cut,aux,r,r0,weight,beta0,cut0,aux0
      DOUBLE PRECISION baux,fac,omega(m)
      DOUBLE PRECISION dos(ngrid),xw(ngrid),gfr(ngrid)
      PARAMETER (pi=3.1415926535898d0)
      EXTERNAL ferm0
      COMMON /doskram/dos,xw
      
      
c      beta=30.0d0
c      cut=10.0d0
c      aux=5.0d0
      DO i=1,ngrid
         xw(i)=omega(i)
      END DO

      beta=beta0
      cut=cut0
      aux=aux0
      baux=-aux0
      r=r0
      CALL LOCATE(xw,ngrid,baux,m1)
      CALL LOCATE(xw,ngrid,aux,m2)
      DO i=1,m1
         dos(i)=exp(-xw(i)**2/pi)
      END DO
      fac=exp(-xw(m1)**2/pi)/dabs(xw(m1))**r
      DO i=m1+1,m2
         dos(i)=fac*dabs(xw(i))**r
      END DO
      DO i=m2+1,ngrid
         dos(i)=exp(-xw(i)**2/pi)
      END DO
c      DO i=1,ngrid
c         write(19,*) xw(i),dos(i)
c      END DO

c      DO i=1,n
c         IF(dabs(xw(i)).LT.aux) THEN
c            dos(i)=dabs(xw(i))**r
c         ELSE
c            dos(i)=dabs(aux)**r*ferm0(beta*(xw(i)-cut-15.0/beta))*
c     &           ferm0(beta*(-xw(i)-cut-15.0/beta))
c         ENDIF
c      END DO
c      
c      DO i=1,n
c         dos(i)=exp(-xw(i)**2/pi)
c      END DO
**integration
      weight=0.0d0
      DO i=1,ngrid-1
         weight=weight+0.5*(xw(i+1)-xw(i))*(dos(i+1)+dos(i))
      END DO
      weight=weight/pi

      DO i=1,ngrid
         dos(i)=dos(i)/weight
      END DO
      xxx=weight
***********************************************

      RETURN
      END

***********************************************************************
***********************************************************************
      SUBROUTINE KRAKRO(n,x,im,re)
c kramers kronig-relation
       implicit double precision (a-h,o-z), integer (i-n)
       double precision x,im,re,lg
       dimension x(n),im(n),re(n),lg(n),dx(n),st(n)

       pi=4d0*datan(1.0d0)
       do i=1,n
        dx(i)=x(min(n,i+1))-x(max(1,i-1))
        st(i)=(im(min(n,i+1))-im(max(1,i-1)))/dx(i)
        if ((i.gt.1).and.(i.lt.n)) then
         lg(i)=dlog(dabs((x(i)-x(1))/(x(i)-x(n))))*im(i)
        endif
       enddo

       do j=1,n
        re(j)=0d0
        if ((j.gt.1).and.(j.lt.n)) re(j)=-2.0d0*lg(j)
        do i=1,n
         if (i.ne.j) then
          re(j)=re(j)+(im(i)-im(j))/(x(i)-x(j))*dx(i)
         else
          re(j)=re(j)+st(i)*dx(i)
         endif
        enddo
       enddo
       call dscal(n,-0.5d0/pi,re,1)
      end
***********************************************************************
***********************************************************************
      SUBROUTINE FERMIFUNC0
      IMPLICIT NONE
      INTEGER i,nyyy
      DOUBLE PRECISION x,dx,width
      PARAMETER (nyyy=4000)
      PARAMETER (width=20.0d0)
      DOUBLE PRECISION  fermiown(nyyy+1)
      COMMON /fermikram/ fermiown,dx
      
      dx=width/dfloat(nyyy)
      DO i=1,nyyy+1
         x=(i-1)*dx-width/2.0d0
         IF(x.LE.0.0d0) THEN
            fermiown(i)=1.0d0/(dexp(x)+1.0d0)
         ELSE
            fermiown(i)=dexp(-x)/(dexp(-x)+1.0d0)
         ENDIF
      END DO

      RETURN
      END
***********************************************************************
***********************************************************************
      DOUBLE PRECISION FUNCTION ferm0(xxx)
      IMPLICIT NONE
      INTEGER i,nyyy
      DOUBLE PRECISION width,dx,xxx,xyz,a
      PARAMETER (nyyy=4000)
      PARAMETER (width=20.0d0)
      DOUBLE PRECISION  fermiown(nyyy+1)
      COMMON /fermikram/ fermiown,dx

      a=width/2.0d0
      IF (xxx.LT.-a) THEN
         ferm0=1.0d0
      ELSEIF(xxx.GE.a) THEN
         ferm0=0.0d0
      ELSE
         xyz=mod(xxx+a,dx)
         i=(xxx+a)/dx+1
         ferm0=fermiown(i)+xyz/dx*(fermiown(i+1)-fermiown(i))
      ENDIF
      RETURN
      END
***********************************************************************
***********************************************************************
      SUBROUTINE BOSEFUNC0(omega,temp)
      IMPLICIT NONE
      INTEGER i,nn,n0,n,nadd,naux,ngrid,new,nNeu,mnew,nnn,nvolt,nwv
      INCLUDE 'PARAMETER' 
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER(n0=ngrid)
      DOUBLE PRECISION  boseown(n0),omega(n0)
      DOUBLE PRECISION omega0(n0),temp,beta
      COMMON /bosekram/ boseown,omega0

      beta=1.0d0/temp

      DO i=1,n0
         omega0(i)=omega(i)
      END DO
      DO i=1,n0
         IF(omega(i).LT.-1.0d-12) THEN
            boseown(i)=1.0d0/(dexp(beta*omega(i))-1.0d0)
         ELSEIF(omega(i).GT.1.0d-12) THEN
            boseown(i)=-1.0d0*dexp(-beta*omega(i))/
     &           (dexp(-beta*omega(i))-1.0d0)
         ELSE
            boseown(i)=0.d0
         ENDIF
         boseown(i)=2.0*boseown(i)+1.0
      END DO
      
      
      RETURN
      END
***********************************************************************
***********************************************************************
      DOUBLE PRECISION FUNCTION bose0(xxx)
      IMPLICIT NONE
      INTEGER i,n,nadd,naux,ngrid,n0,new,Nneu,mnew,nnn,nvolt,nwv
      INTEGER m
      INCLUDE 'PARAMETER' 
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER(n0=ngrid)
      DOUBLE PRECISION xxx,a
      DOUBLE PRECISION  boseown(n0),omega(n0)
      COMMON /bosekram/ boseown,omega

      CALL LOCATE(omega,n0,xxx,m)
      IF(m.LT.1) THEN
         bose0=-1.0d0
      ELSEIF(m.GE.n0) THEN
         bose0=1.0d0
      ELSE
         bose0=boseown(m)+(boseown(m+1)-boseown(m))/
     &        (omega(m+1)-omega(m))*(xxx-omega(m))
      ENDIF
      RETURN
      END
***********************************************************************
***********************************************************************

***********************************************************************
***********************************************************************
      subroutine dzero(darray,ndim)
      implicit real*8 (a-h,o-z)
      
c---------------------------------------------------------------------
c     initialize array <darray(1,...,ndim)> of floats to zero
c---------------------------------------------------------------------
      
      dimension darray(ndim)
      
      m=mod(ndim,7)
      if  (m.ne.0) then
         do 10 i = 1,m
            darray(i) = 0.0d0
 10      continue
      endif
      mp1 = m + 1
      do 20 i = mp1,ndim,7
         darray(i) = 0.0d0
         darray(i+1) = 0.0d0
         darray(i+2) = 0.0d0
         darray(i+3) = 0.0d0
         darray(i+4) = 0.0d0
         darray(i+5) = 0.0d0
         darray(i+6) = 0.0d0
 20   continue
      return
      end

      subroutine izero(iarray,ndim)
      implicit real*8 (a-h,o-z)
      
c---------------------------------------------------------------------
c     initialize array <darray(1,...,ndim)> of floats to zero
c---------------------------------------------------------------------
      
      dimension iarray(ndim)
      
      m=mod(ndim,7)
      if  (m.ne.0) then
         do 10 i = 1,m
            iarray(i) = 0.0d0
 10      continue
      endif
      mp1 = m + 1
      do 20 i = mp1,ndim,7
         iarray(i) = 0.0d0
         iarray(i+1) = 0.0d0
         iarray(i+2) = 0.0d0
         iarray(i+3) = 0.0d0
         iarray(i+4) = 0.0d0
         iarray(i+5) = 0.0d0
         iarray(i+6) = 0.0d0
 20   continue
      return
      end
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
      DOUBLE PRECISION FUNCTION ferm1(xxx)
      IMPLICIT NONE
      INTEGER i,n,nadd,naux,ngrid,n0,new,Nneu,mnew,nnn,nvolt,nwv
      INTEGER m
      INCLUDE 'PARAMETER' 
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER(n0=ngrid)

      DOUBLE PRECISION xxx
      DOUBLE PRECISION  fermown(n0),omega(n0)
      COMMON /fermkram/ fermown,omega

      CALL LOCATE(omega,n0,xxx,m)
      IF(m.LE.1) THEN
         ferm1=-1.0d0
      ELSEIF(m.GE.n0) THEN
         ferm1=1.0d0
      ELSE
         ferm1=fermown(m)+(fermown(m+1)-fermown(m))/
     &        (omega(m+1)-omega(m))*(xxx-omega(m))
      ENDIF
      RETURN
      END
***********************************************************************
***********************************************************************
      DOUBLE PRECISION FUNCTION aphi(xxx)
      IMPLICIT NONE
      INTEGER nxxx,j0
      DOUBLE PRECISION pi,bandwidth,dx,xxx,beta,cut,ferm0,gamma,dummy
      parameter (pi=3.1415926535898d0) 
      DOUBLE PRECISION pref,xhelp
      EXTERNAL ferm0
      COMMON /versuch/gamma,cut,beta,pref,xhelp
   

      IF(dabs(xxx).LT.(1.0d-12)) THEN
         aphi=0.0d0
      ELSEIF ((xxx.LT.0.0d0).AND.(xxx.GT.-cut)) THEN
         aphi=pi*pref*dabs(xxx)**gamma

      ELSEIF ((xxx.GT.0.0d0).AND.(xxx.LT.cut)) THEN
         aphi=-pi*pref*dabs(xxx)**gamma

      ELSE

         aphi=-sign(1.0d0,xxx)*pi*pref*xhelp*
     &   ferm0(beta*(xxx-cut-4.0/beta))*ferm0(beta*(-xxx-cut-4.0/beta))   
      ENDIF

      RETURN
      END
***********************************************************************
***********************************************************************
      DOUBLE PRECISION FUNCTION gammaf(xx)
      DOUBLE PRECISION xx,gammln
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+dlog(stp*ser/x)
      gammaf=dexp(gammln)
      return
      END
*********************************************************************
*********************************************************************
*********************************************************************
      SUBROUTINE FERMFUNCV(omega,temp,J1,J2,V)
      IMPLICIT NONE
      INTEGER i,nn,n0,n,nadd,naux,ngrid,new,nNeu,mnew,nnn,nwv,nvolt
      INCLUDE 'PARAMETER'
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      PARAMETER (n0=ngrid)
      DOUBLE PRECISION J1,J2,x1,x2,V
      DOUBLE PRECISION fermown(n0),omega(n0)
      DOUBLE PRECISION omega0(n0),temp,beta,xxx
      COMMON /fermkram/ fermown,omega0

      beta=1.0d0/temp
      DO i=1,n0
         omega0(i)=omega(i)
         x1=beta*(omega(i)-V)
         x2=beta*(omega(i)+V)
         IF(x1.LE.0.0d0) THEN
            xxx=J1-2.0*J1/(dexp(x1)+1.0d0)
            IF(x2.LE.0.0d0) THEN
               fermown(i)=xxx+J2-2.0*J2/(dexp(x2)+1.0d0)
            ELSE
               fermown(i)=xxx+J2-2.0*J2*dexp(-x2)/(dexp(-x2)+1.0d0)
            ENDIF
         ELSE
            xxx=J1-2.0*J1*dexp(-x1)/(dexp(-x1)+1.0d0)
            IF(x2.LE.0.0d0) THEN
               fermown(i)=xxx+J2-2.0*J2/(dexp(x2)+1.0d0)
            ELSE
               fermown(i)=xxx+J2-2.0*J2*dexp(-x2)/(dexp(-x2)+1.0d0)
            ENDIF
         ENDIF
         fermown(i)=fermown(i)/(J1+J2)
      END DO

      
      RETURN
      END
************************************************************************
************************************************************************
*********************CALCULATES PARTICLE NUMBERS************************
      SUBROUTINE PNNESS(n,omega,field1,rnf)
      IMPLICIT NONE
      INTEGER i,j,n
      DOUBLE PRECISION exom,rnf,beta,field1(n),omega(n)
      DOUBLE PRECISION wert1,wert2,wertn,freq1,freq2
      DOUBLE PRECISION omex,aux
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535898d0)

      aux=0.0d0
     
      DO j=1,n-1
         aux=aux+0.5*(omega(j+1)-omega(j))*(field1(j+1)+field1(j))
      END DO

      rnf=0.5*(1.0-aux)
      
      RETURN
      END

