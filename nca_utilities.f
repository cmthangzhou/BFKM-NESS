************************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE OWNBUBBLE(beta,omegf,nf,omegu,nu,af,au,ad)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      IMPLICIT NONE
      INTEGER nf,nu,nca
      INTEGER nnn,kkk,numbre,mumber,i,l,k,izaehler
      DOUBLE PRECISION pi,beta
      PARAMETER(pi=3.1415926535898d0) 
      DOUBLE PRECISION omegf(nf),omegu(nu),af(nf),au(nu)
      DOUBLE PRECISION cernel(nu+nf,nf),aun(nu+nf,nf),afn(nu+nf,nf)
      DOUBLE PRECISION omegfa(nu+nf,nf),ad(nf)
      DOUBLE PRECISION omegfv(nu+nf,nf),aux(nu+nf),aux2(nu+nf),
     &     aux3(nu+nf)
      INTEGER indx(nu+nf),inaux(nu+nf),kndex(nu+nf,nf),izeiger(2)
      DOUBLE PRECISION help(2),depsilon(nu+nf,nf)
      DOUBLE PRECISION ferm0
      DOUBLE PRECISION aaa,dummy,peak,xx1
      EXTERNAL ferm0
      nca=nu+nf

      do 1111 l=1,nf            !all the omegas
         peak=omegf(l)
         do k=1,nu
            omegfv(k,l)=omegu(k)+peak
         end do                 
         
         do i=1,nu              !indx is an kndex field
            indx(i)=i
            aux(i)=omegfv(i,l)            
         end do
         
c!==================================!dummy is flag for second grid
         dummy=-666666.0d0
         do i=1,nf
            aux2(i+nu)=dummy
            aux3(i+nu)=af(i)
         end do
c!==================================
         do i=1,nu
            aux2(i)=au(i)
            aux3(i)=dummy
         end do
         do i=1,nf
            indx(i+nu)=i+nu
            aux(i+nu)=omegf(i)
         end do
c----------------------------
         call ownsort(nca,nu,aux,indx) !sort the frequencies
c----------------------------
         do k=1,nca
            inaux(k)=indx(k)
         end do
         
         do i=1,nca-1
            if ((aux(i+1)-aux(i)).LT.1.0d-9) then
               if(i.eq.(nca-1)) then
                  aux(i+1)=2.0d0*aux(i)-aux(i-1)
               else
                  aux(i+1)=0.5d0*(aux(i)+aux(i+2))
                  
                  if(indx(i+1).lt.indx(i)) then !indx(i+1) from omegu
                     numbre=indx(i+1) !indx for interpol of omegu
                     indx(i+1)=indx(i)
                     indx(i)=numbre
                  else
                     numbre=inaux(i+1)
                     inaux(i+1)=inaux(i)
                     inaux(i)=numbre
                  endif
                  
               end if
            endif
         end do
         
         
         do i=1,nca
            kndex(i,l)=indx(i)  !store kndex table
         end do
         do i=1,nca
            omegfa(i,l)=aux(i)  !store frequency vectors
         end do
         
c-----------------------
         do i=1,nca             !sorting while copying
            aun(i,l)=aux2(indx(i))
            afn(i,l)=aux3(inaux(i))
         end do
*----------------------------
c!        loop over ordered vector; picks out dummy fields and interpolates
         help(1)=0.0d0
         izeiger(1)=0
         izaehler=0
         do i=1,nca
            IF (aun(i,l).NE.dummy) THEN
               izaehler=izaehler+1
               help(2)=aun(i,l)
               izeiger(2)=i
               IF (izeiger(1).NE.(i-1)) THEN
                  mumber=i-izeiger(1)-1
                  nnn=izeiger(1)
                  DO k=1,mumber
                     kkk=nnn+k
                     if(nnn.eq.0) then
                        aun(kkk,l)=0.0d0
                     else
                        xx1=aun(i,l)-aun(nnn,l)
                        aaa=(omegfa(kkk,l)-omegfa(nnn,l))/
     &                       (omegfa(i,l)-omegfa(nnn,l))
                        aun(kkk,l)=aun(nnn,l)+xx1*aaa
                     endif
                  END DO
               END IF
               IF (izaehler.EQ.nu) THEN
                  DO kkk=i+1,nca
                     aun(kkk,l)=0.0d0
                  END DO
                  GOTO 111
               END IF
               help(1)=help(2)
               izeiger(1)=izeiger(2)
            END IF
         END DO
         
 111     CONTINUE
         
         
*----------------------------
c!        loop over ordered vector; picks out dummy fields and interpolates
         help(1)=0.0d0
         izeiger(1)=0
         izaehler=0
         do i=1,nca
            IF (afn(i,l).NE.dummy) THEN
               izaehler=izaehler+1
               help(2)=afn(i,l)
               izeiger(2)=i
               IF (izeiger(1).NE.(i-1)) THEN
                  mumber=i-izeiger(1)-1
                  nnn=izeiger(1)
                  DO k=1,mumber
                     kkk=nnn+k
                     if(nnn.eq.0) then
                        afn(kkk,l)=0.0d0
                     else
                        xx1=afn(i,l)-afn(nnn,l)
                        aaa=(omegfa(kkk,l)-omegfa(nnn,l))/
     &                       (omegfa(i,l)-omegfa(nnn,l))
                        afn(kkk,l)=afn(nnn,l)+xx1*aaa
                     endif
                  END DO
               END IF
               IF (izaehler.EQ.nf) THEN
                  DO kkk=i+1,nca
                     afn(kkk,l)=0.0d0
                  END DO
                  GOTO 222
               END IF
               help(1)=help(2)
               izeiger(1)=izeiger(2)
            END IF
         END DO
         
 222     CONTINUE
         
c------------------------
         DO i=1,nca
            if (i.eq.1) then
               depsilon(i,l)=(omegfa(2,l)-omegfa(1,l))/2.0d0
            else if (i.eq.nca) then
               depsilon(i,l)=(omegfa(nca,l)-omegfa(nca-1,l))/2.0d0
            else
               depsilon(i,l)=(omegfa(i+1,l)-omegfa(i-1,l))/2.0d0
            endif
            cernel(i,l)=ferm0(beta*omegfa(i,l))-
     &           ferm0(beta*(omegfa(i,l)+omegf(l)))
         END DO
c-----------------------------
         ad(l)=0.0d0
         DO i=1,nca
C            ad(l)=ad(l)+aun(i,l)*afn(i,l)*depsilon(i,l)
            ad(l)=ad(l)+cernel(i,l)*aun(i,l)*afn(i,l)*depsilon(i,l)
         END DO
         ad(l)=ad(l)*pi
 1111 CONTINUE
      
      RETURN
      END
      
************************************************************************
************************************************************************
      subroutine ownsort(nges,n1,w,indx)
      implicit integer (i,j,k,l,m,n)
      implicit real*8 (a-h,o-z)
      dimension w(nges),aux(nges),indx(nges)
      
      
      n2=nges-n1
      l=1
      k=n1+1
      do i=1,nges
         mem=i
         if(w(l).lt.w(k)) then
            indx(i)=l
            aux(i)=w(l)
            if(l.eq.n1) then
               goto 30
            end if
            l=l+1
         else
            indx(i)=k
            aux(i)=w(k)
            if(k.eq.nges) then
               goto 40
            end if
            k=k+1
         end if
      end do

 30   continue
      do i=k,nges
         aux(i)=w(i)
         indx(i)=i
      end do
      goto 50
 40   continue
      do i=l,n1
         aux(i+n2)=w(i)
         indx(i+n2)=i
      end do

 50   continue
      do i=1,nges
         w(i)=aux(i)
      end do
      return
      end
************************************************************************
************************************************************************
*********************CALCULATES PARTICLE NUMBERS************************
      SUBROUTINE PNNESS(n,omega,field1,field2,nf)
      IMPLICIT NONE
      INTEGER i,j,n,m,iscr(n),kk,k,l
      DOUBLE PRECISION exom,nf,beta,field1(n),field2(n),omega(n)
      DOUBLE PRECISION wert1,wert2,wertn,freq1,freq2
      DOUBLE PRECISION omex,aux
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535898d0)

      aux=0.0d0
      DO i=1,n
         field1(i)=field1(i)-field2(i)
      END DO
      DO j=1,n-1
         aux=aux+0.5*(omega(j+1)-omega(j))*
     &        (field1(j+1)+field1(j))
      END DO

      nf=-aux/2.0
      
      RETURN
      END
************************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE OWNBUB(beta,omegf,nf,af,ad)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      IMPLICIT NONE
      INTEGER nf
      INTEGER nnn,kkk,numbre,mumber,i,l,k,izaehler
      DOUBLE PRECISION pi,beta
      PARAMETER(pi=3.1415926535898d0) 
      DOUBLE PRECISION omegf(nf),af(nf),au(nf)
      DOUBLE PRECISION cernel(2*nf,nf),aun(2*nf,nf),afn(2*nf,nf)
      DOUBLE PRECISION omegfa(2*nf,nf),ad(nf)
      DOUBLE PRECISION omegfv(2*nf,nf),aux(2*nf),aux2(2*nf),
     &     aux3(2*nf)
      INTEGER indx(2*nf),inaux(2*nf),kndex(2*nf,nf),izeiger(2)
      DOUBLE PRECISION help(2),depsilon(2*nf,nf)
      DOUBLE PRECISION ferm0
      DOUBLE PRECISION aaa,dummy,peak,xx1
      EXTERNAL ferm0

      DO i=1,nf
         au(i)=af(i)
      END DO
      do 1111 l=1,nf            !all the omegas
         peak=omegf(l)
         do k=1,nf
            omegfv(k,l)=omegf(k)+peak
         end do                 
         
         do i=1,nf              !indx is an kndex field
            indx(i)=i
            aux(i)=omegfv(i,l)            
         end do
         
c!==================================!dummy is flag for second grid
         dummy=-666666.0d0
         do i=1,nf
            aux2(i+nf)=dummy
            aux3(i+nf)=af(i)
         end do
c!==================================
         do i=1,nf
            aux2(i)=au(i)
            aux3(i)=dummy
         end do
         do i=1,nf
            indx(i+nf)=i+nf
            aux(i+nf)=omegf(i)
         end do
c----------------------------
         call ownsort(2*nf,nf,aux,indx) !sort the frequencies
c----------------------------
         do k=1,2*nf
            inaux(k)=indx(k)
         end do
         
         do i=1,2*nf-1
            if ((aux(i+1)-aux(i)).LT.1.0d-9) then
               if(i.eq.(2*nf-1)) then
                  aux(i+1)=2.0d0*aux(i)-aux(i-1)
               else
                  aux(i+1)=0.5d0*(aux(i)+aux(i+2))
                  
                  if(indx(i+1).lt.indx(i)) then !indx(i+1) from omegu
                     numbre=indx(i+1) !indx for interpol of omegu
                     indx(i+1)=indx(i)
                     indx(i)=numbre
                  else
                     numbre=inaux(i+1)
                     inaux(i+1)=inaux(i)
                     inaux(i)=numbre
                  endif
                  
               end if
            endif
         end do
         
         
         do i=1,2*nf
            kndex(i,l)=indx(i)  !store kndex table
         end do
         do i=1,2*nf
            omegfa(i,l)=aux(i)  !store frequency vectors
         end do
         
c-----------------------
         do i=1,2*nf          !sorting while copying
            aun(i,l)=aux2(indx(i))
            afn(i,l)=aux3(inaux(i))
         end do
*----------------------------
c!        loop over ordered vector; picks out dummy fields and interpolates
         help(1)=0.0d0
         izeiger(1)=0
         izaehler=0
         do i=1,2*nf
            IF (aun(i,l).NE.dummy) THEN
               izaehler=izaehler+1
               help(2)=aun(i,l)
               izeiger(2)=i
               IF (izeiger(1).NE.(i-1)) THEN
                  mumber=i-izeiger(1)-1
                  nnn=izeiger(1)
                  DO k=1,mumber
                     kkk=nnn+k
                     if(nnn.eq.0) then
                        aun(kkk,l)=0.0d0
                     else
                        xx1=aun(i,l)-aun(nnn,l)
                        aaa=(omegfa(kkk,l)-omegfa(nnn,l))/
     &                       (omegfa(i,l)-omegfa(nnn,l))
                        aun(kkk,l)=aun(nnn,l)+xx1*aaa
                     endif
                  END DO
               END IF
               IF (izaehler.EQ.nf) THEN
                  DO kkk=i+1,2*nf
                     aun(kkk,l)=0.0d0
                  END DO
                  GOTO 111
               END IF
               help(1)=help(2)
               izeiger(1)=izeiger(2)
            END IF
         END DO
         
 111     CONTINUE
         
         
*----------------------------
c!        loop over ordered vector; picks out dummy fields and interpolates
         help(1)=0.0d0
         izeiger(1)=0
         izaehler=0
         do i=1,2*nf
            IF (afn(i,l).NE.dummy) THEN
               izaehler=izaehler+1
               help(2)=afn(i,l)
               izeiger(2)=i
               IF (izeiger(1).NE.(i-1)) THEN
                  mumber=i-izeiger(1)-1
                  nnn=izeiger(1)
                  DO k=1,mumber
                     kkk=nnn+k
                     if(nnn.eq.0) then
                        afn(kkk,l)=0.0d0
                     else
                        xx1=afn(i,l)-afn(nnn,l)
                        aaa=(omegfa(kkk,l)-omegfa(nnn,l))/
     &                       (omegfa(i,l)-omegfa(nnn,l))
                        afn(kkk,l)=afn(nnn,l)+xx1*aaa
                     endif
                  END DO
               END IF
               IF (izaehler.EQ.nf) THEN
                  DO kkk=i+1,2*nf
                     afn(kkk,l)=0.0d0
                  END DO
                  GOTO 222
               END IF
               help(1)=help(2)
               izeiger(1)=izeiger(2)
            END IF
         END DO
         
 222     CONTINUE
         
c------------------------
         DO i=1,2*nf
            if (i.eq.1) then
               depsilon(i,l)=(omegfa(2,l)-omegfa(1,l))/2.0d0
            else if (i.eq.2*nf) then
               depsilon(i,l)=(omegfa(2*nf,l)-omegfa(2*nf-1,l))/2.0d0
            else
               depsilon(i,l)=(omegfa(i+1,l)-omegfa(i-1,l))/2.0d0
            endif
            cernel(i,l)=ferm0(beta*omegfa(i,l))-
     &           ferm0(beta*(omegfa(i,l)+omegf(l)))
         END DO
c-----------------------------
         ad(l)=0.0d0
         DO i=1,2*nf
C            ad(l)=ad(l)+aun(i,l)*afn(i,l)*depsilon(i,l)
            ad(l)=ad(l)+cernel(i,l)*aun(i,l)*afn(i,l)*depsilon(i,l)
         END DO
         ad(l)=ad(l)
 1111 CONTINUE
      
      RETURN
      END
      
************************************************************************
************************************************************************
