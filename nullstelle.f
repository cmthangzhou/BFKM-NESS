      SUBROUTINE ZERORF(FCT,ABSERR,RELERR,FAB,MAXIT,IMETH,IEXTRA,
     +                  DELX,X1,X2,X3,F1,F2,F3,NUMIT,IHELP1,
     +                  IHELP2,INCL,IERR)
C
C*****************************************************************
C                                                                *
C  The SUBROUTINE ZERORF is the governing program for several    *
C  subroutines that determine zeros of a continuous real-valued  *
C  function FCT. The methods used are primarily intended for     *
C  simple zeros and zeros of odd order.                          *
C  This SUBROUTINE adapts its approach depending on whether      *
C  the starting values enclose a zero of the function FCT or not.*
C  If a zero is enclosed, an iteration sequence is constructed   *
C  using one of several iteration methods: Pegasus,              *
C  Anderson/Bjoerck, or King. These methods ensure that inclusion*
C  is preserved. All methods used work without derivatives and   *
C  each possesses a convergence order of P > 1.6.                *
C  If no two starting estimates are known that enclose a         *
C  functional zero, then two identical starting values X1, X2    *
C  may be specified.                                             *
C  The SUBROUTINE ZERORF will construct a                        *
C  second starting value X2 by adding DELX to the starting value *
C  X1. DELX has to be chosen by the user.                        *
C  If there is still no enclosure, linear extrapolation is tried,*
C  or, following the first secant step, quadratic extrapolation  *
C  using the tangent line of an interpolating parabola through   *
C  three different interpolation points is used. If enclosure    *
C  is finally achieved, the program continues with one of the    *
C  three special iterative methods as chosen by the user.        *
C  NOTE: Extensive tests have shown that the number of iteration *
C        steps required per zero for the same degree of          *
C        accuracy is the largest for the Pegasus-method          *
C        (SUBROUTINE PEG) and by far the smallest for the        *
C        combination of the Anderson/Bjoerck and King            *
C        methods (SUBROUTINE ANDBJK).                            *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FCT      : real-valued function for which a zero is to be     *
C             determined. It is declared as                      *
C                 DOUBLE PRECISION FUNCTION FCT(X)               *
C             and has to be defined as EXTERNAL within the       *
C             calling program (or as INTRINSIC if a FORTRAN      *
C             standard function is used).                        *
C  ABSERR   : ) error bounds both of which have to be >= 0.0.    *
C  RELERR   : ) Their sum has to be > 0.0. The following mixed   *
C               test is used as a break-off criterion:           *
C                   ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.         *
C               Thus if RELERR=0.0 is chosen, this tests for the *
C               absolute error, if ABSERR=0.0, this tests for the*
C               relative error.                                  *
C               The values entered for ABSERR and RELERR are     *
C               accepted unchanged by the program if they        *
C               both exceed four times the machine constant, or, *
C               if one is zero, then the other has to exceed     *
C               four times the machine constant. If this is not  *
C               the case, then both or one of the bounds is set  *
C               internally to this value.                        *
C  FAB      : break-off criterion constant for the functional    *
C             value at the last approximation of the zero. The   *
C             value for FAB can be chosen between zero and four  *
C             times the machine constant. If it is chosen        *
C             negative or larger than four times the machine     *
C             constant, it is set to that value internally.      *
C  MAXIT    : maximum number of functional evaluations.          *
C             (this equals the number of iteration steps)        *
C  IMETH    : = 1, use the Pegasus-method.                       *
C             = 2, use the King-method.                          *
C             = 3, use the Anderson/Bjoerck-method.              *
C             = 4, use the method of Anderson/Bjoerck and King.  *
C  IEXTRA   : = 0, quadratic extrapolation shall be allowed.     *
C             = 1, only linear extrapolation is acceptable (if   *
C                  e.g. a double root is expected or if setting  *
C                  IEXTRA=0 has resulted in IERR=0).             *
C             If the starting values do not enclose a zero, i.e.,*
C             if ( F1*F2 > 0.0, one may predetermine via IEXTRA, *
C             whether only linear extrapolation should be tried  *
C             or whether quadratic extrapolation is also         *
C             acceptable. If no double zero is expected, then    *
C             IEXTRA=0 should be used first.                     *
C  DELX     : is used for altering a starting value in case two  *
C             identical starting values were entered. A meaning- *
C             ful choice for DELX would be:                      *
C                     DELX=1E-6 or DELX=1E-8.                    *
C  X1,X2    : starting values for the iteration; if only one     *
C             starting value is known, setting X2=X1 is          *
C             acceptable. By adding DELX to one of the two       *
C             identical starting values, the program will create *
C             two different starting values and proceed.         *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  ABSERR   : ) error bounds actually used.                      *
C  RELERR   : )                                                  *
C  X1,X2,X3 : approximate values for the desired zero.           *
C             (see IERR).                                        *
C  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *
C  NUMIT    : number of functional evaluations performed.        *
C  INCL     : = 0, starting values with no enclosure.            *
C             = 1, starting values with known enclosure.         *
C  IERR     : = 0, zero was not found.                           *
C             = 1, zero lies between X1 and X2  (of the two      *
C                  values the one with the smaller absolute      *
C                  functional value should be chosen as the      *
C                  zero).                                        *
C                  The absolute error of the computed zero is    *
C                  smaller than or equal to ABS(X1-X2).          *
C             = 2, X2 is a zero of FCT: F2=0.0 (machine zero).   *
C             = 3, X3 is a zero with ABS(F3) < 4 * machine       *
C                  constant.                                     *
C             = 4, zero of FCT is at X2 (enclosure interval      *
C                  not definable).                               *
C             = 5, maximum number MAXIT of functional evaluations*
C                  exceeded.                                     *
C             = 6, ABSERR or RELERR are negative, or both are    *
C                  equal to zero, or MAXIT < 1.                  *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  IHELP1,IHELP2 : internal variables.                           *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: EXCHG, PEG, PEGK, ANDBJ, ANDBJK,        *
C                        SUB1, SUB2, MACHPD                      *
C                                                                *
C                                                                *
C  sources: 1. Anderson/Bjoerck, see [ANDE73].                   *
C           2. Dowell/Jarrat, see [DOWE71], [DOWE72].            *
C           3. King, see [KING73].                               *
C           4. unpublished manuscript by R. Wodicka,             *
C              RWTH Aachen.                                      *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 08.26.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  if two identical starting values are initially given (X1=X2),
C  a second starting value is internally generated.
C
      IF(X1 .EQ. X2) THEN
         X1=X1+DELX
      END IF
C
C  initializing the parameters IERR and NUMIT.
C
      IERR=1
      NUMIT=2
C
C  calculation of the machine constant FMACHP.
C
      FMACHP=1.0D0
   10 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP=2.0D0*FMACHP
C
C  testing the validity of the error bounds and MAXIT.
C
      IF(ABSERR .GE. 0.0D0 .AND. RELERR .GE. 0.0D0 .AND.
     +   ABSERR+RELERR .GT. 0.0D0 .AND. MAXIT .GE. 1) GOTO 20
      IERR=6
      RETURN
   20 DUMMY=4.0D0*FMACHP
      IF(RELERR .EQ. 0.0D0) THEN
         IF(ABSERR .LT. DUMMY) ABSERR=DUMMY
      ELSE IF(ABSERR .EQ. 0.0D0) THEN
         IF(RELERR .LT. DUMMY) RELERR=DUMMY
      ELSE
         IF(ABSERR .LT. DUMMY) ABSERR=DUMMY
         IF(RELERR .LT. DUMMY) RELERR=DUMMY
      END IF
C
C  testing the validity of the break-off parameter FAB for the
C  functional value at the approximate zero.
C
      IF(FAB .LT. 0.0D0 .OR. FAB .GT. DUMMY) FAB=DUMMY
C
C  calculating the functional values at the starting points.
C
      F1=FCT(X1)
      F2=FCT(X2)
C
C  labelling the zeros, so that ABS(F2) <= ABS(F1) holds.
C
      IF(DABS(F2) .GT. DABS(F1)) THEN
         CALL EXCHG(X1,X2,F1,F2)
      END IF
C
C  test whether X2 already is a zero of FCT,
C  in which case IERR=2.
C
      IF(F2 .EQ. 0.0D0) THEN
         IERR=2
         RETURN
      END IF
C
C  test whether the starting values X1, X2 enclose a zero of FCT.
C  If so, the internal variables are set to INCLUD=1, IHELP2=1,
C  otherwise INCLUD=0.
C  IHELP2=1 indicates that the next iteration step is to be a secant step.
C
      IF(F1*F2 .GT. 0.0D0) THEN
         INCLUD=0
         F3=F2
      ELSE
         INCLUD=1
         IHELP2=1
      END IF
      INCL=INCLUD
C
C  iteration loop to find a new approximate value X3. We account for
C  the fact here whether we have inclusion or not. The new
C  approximate value X3 is taken as the x-intercept of the straight
C  line connecting (X2, F2) and (X1, F1). Unless a secant step is
C  performed, F1 is changed: in case of enclosure according to
C  the specified method, otherwise by using quadratic extrapolation
C  (see above remarks).
C
   30 IF(INCLUD .EQ. 0) THEN
         F1DF2=F1/F2
         IF(F1DF2 .GT. 1.0D0) THEN
C
C     check whether linear extrapolation is the only acceptable
C     method or whether quadratic extrapolation is also allowed.
C
            IF(IEXTRA .EQ. 0) THEN
               IF((F1DF2-F1/F3) .GT. 1.0D0) THEN
                  G=1.0D0-F2/F3
                  F1=G*F1
               END IF
            END IF
         ELSE
C
C     now F1/F2 <= 1.0. If moreover F1/F2 < 0.0, then we have enclosure.
C
            IF(F1DF2 .LT. 0.0D0) THEN
               INCLUD=1
               IF(DABS(X1-X2) .LE. DABS(X2)*RELERR+ABSERR) THEN
                  IERR=1
                  RETURN
               ELSE
                  IHELP2=1
               END IF
            ELSE
C
C     if there is no enclosure, then a zero cannot be found
C     because of ABS(F1) <= ABS(F2).
C
               IERR=0
               RETURN
            END IF
         END IF
      END IF
C
C  calculation of the scaling factor Q for X3=X2+Q(X1-X2).
C
      Q=F2/(F2-F1)
C
C  calculation of the new approximate value.
C
      X3=X2+Q*(X1-X2)
C
C  testing whether the new approximate value X3 differs from both
C  X1 and X2. If this is not the case, an alternate value
C  X3NEW is calculated for X3. If there is no enclosure and
C  X2=X3=X3NEW, then the program is stopped with setting
C  IERR=4 (zero of FCT is at X2).
C
      IF(INCLUD .EQ. 0) THEN
         IF(X2 .EQ. X3) THEN
            X3NEW=X2+(X2-X1)/9.0D0
            IF(X2 .EQ. X3NEW) THEN
               IERR=4
               RETURN
            ELSE
               X3=X3NEW
            END IF
         ELSE
            IF(X2 .EQ. X3) THEN
               X3NEW=X2+(X1-X2)/3.0D0
               IF(X2 .EQ. X3NEW) THEN
                  IERR=1
                  RETURN
               ELSE
   40             Q=2.0D0*Q
                  X3NEW=X2+Q*(X1-X2)
                  IF(X3NEW .EQ. X2) GOTO 40
                  IF(X3NEW .EQ. X1) THEN
                     IERR=1
                     RETURN
                  ELSE
                     X3=X3NEW
                  END IF
               END IF
            ELSE
               IF(X1 .EQ. X3) THEN
                  X3NEW=X1+(X2-X1)/3.0D0
                  IF(X3NEW .EQ. X1) THEN
                     IERR=1
                     RETURN
                  ELSE
                     Q=F1/(F1-F2)
   50                Q=2.0D0*Q
                     X3NEW=X1+Q*(X2-X1)
                     IF(X3NEW .EQ. X1) GOTO 50
                     IF(X3NEW .EQ. X2) THEN
                        IERR=1
                        RETURN
                     ELSE
                        X3=X3NEW
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C  now X3 differs from both X1 and X2, and F3 has been calculated.
C
      F3=FCT(X3)
C
C  increasing the counter NUMIT.
C
      NUMIT=NUMIT+1
C
C  test whether |F3| is less than four times the
C  machine constant. In this case, X3 is a zero with IERR=3.
C
      IF(DABS(F3) .LE. FAB) THEN
         IERR=3
         RETURN
      END IF
C
C  test whether the maximum number of functional evaluations
C  allowed has been exceeded.
C
      IF(NUMIT .GE. MAXIT) THEN
         IERR=5
         RETURN
      END IF
C
C  the approximate values and their functional values
C  are relabelled by applying the SUBROUTINE EXCHG.
C  If F2*F3 > 0.0, the auxiliary variable IHELP1 is set to 0;
C  if F2*F3 < 0.0, it is set to IHELP1=1.
C
      IF(INCLUD .EQ. 0) THEN
         CALL EXCHG(X1,X2,F1,F2)
      ELSE
         IF(F2*F3 .LT. 0.0D0) THEN
            CALL EXCHG(X1,X2,F1,F2)
            IHELP1=1
         ELSE
            IHELP1=0
         END IF
      END IF
      CALL EXCHG(X2,X3,F2,F3)
C
C  if we have enclosure, we check whether the break-off criterion
C  holds for the new enclosure interval. If the break-off condition
C  is met, the iteration is stopped with IERR=1.
C  Otherwise the user specified method is used to determine F1.
C
      IF(INCLUD .EQ. 1) THEN
         IF(DABS(X1-X2) .GT. DABS(X2)*RELERR+ABSERR) THEN
            IF(IMETH .EQ. 1) THEN
               CALL PEG(IHELP1,F1,F2,F3)
            ELSE IF(IMETH .EQ. 2) THEN
               CALL PEGK(IHELP1,IHELP2,F1,F2,F3)
            ELSE IF(IMETH .EQ. 3) THEN
               CALL ANDBJ(IHELP1,F1,F2,F3)
            ELSE
               CALL ANDBJK(IHELP1,IHELP2,F1,F2,F3)
            END IF
         ELSE
            IERR=1
            RETURN
         END IF
      END IF
      GOTO 30
      END
C
C

      SUBROUTINE PEG(IHELP1,F1,F2,F3)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE uses the Pegasus-method to calculate a new    *
C  functional value F1 for the governing program ZERORF.         *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  IHELP1   : auxiliary variable, as specified by the            *
C             governing program ZERORF.                          *
C  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  F1       : new functional value at X1.                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SUB1                                    *
C                                                                *
C                                                                *
C  sources: Dowell/Jarrat, see at [DOWE71], [DOWE72].            *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 08.26.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(IHELP1 .EQ. 0) THEN
         CALL SUB1(F1,F2,F3)
      END IF
      RETURN
      END
C
C
C

      SUBROUTINE PEGK(IHELP1,IHELP2,F1,F2,F3)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE uses an improved version of the Pegasus       *
C  method by King for calculating a new functional value F1      *
C  for use in the governing program ZERORF.                      *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  IHELP1   : ) auxiliary variables that are specified by the    *
C  IHELP2   : ) governing program ZERORF.                        *
C  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETER:                                             *
C  =================                                             *
C  F1       : new functional value at X1.                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SUB1                                    *
C                                                                *
C                                                                *
C  sources: method of King, see at [KING73].                     *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 08.26.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(IHELP2 .EQ. 1) THEN
         IHELP2=0
         CALL SUB1(F1,F2,F3)
      ELSE IF(IHELP1 .EQ. 0) THEN
         CALL SUB1(F1,F2,F3)
      ELSE
         IHELP2=1
      END IF
      RETURN
      END
C
C
C

      SUBROUTINE ANDBJ(IHELP1,F1,F2,F3)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE uses the Anderson / Bjoerck method to calcu-  *
C  late a new functional value F1 for the governing program      *
C  ZERORF.                                                       *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  IHELP1   : auxiliary variable that is determined in the       *
C             governing program ZERORF.                          *
C  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETER:                                             *
C  =================                                             *
C  F1       : new functional value at X1.                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SUB2                                    *
C                                                                *
C                                                                *
C  sources: method of Anderson/Bjoerck, see at [ANDE73].        *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 08.26.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(IHELP1 .EQ. 0) THEN
           CALL SUB2(F1,F2,F3)
      END IF
      RETURN
      END
C
C

      SUBROUTINE ANDBJK(IHELP1,IHELP2,F1,F2,F3)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE uses a combination of the Anderson/Bjoerck    *
C  and King methods suggested by King to calculate a new func-   *
C  tional value F1 for use in the governing program ZERORF.      *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  IHELP1   : ) auxiliary variables set by the governing         *
C  IHELP2   : ) program ZERORF.                                  *
C  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETER:                                             *
C  =================                                             *
C  F1       : new functional value at X1.                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SUB2                                    *
C                                                                *
C                                                                *
C  sources: 1. method of Anderson/Bjoerck, see at [ANDE73].      *
C           2. method of King, see at [KING73].                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 08.26.1985                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF(IHELP2 .EQ. 1) THEN
         IHELP2=0
         CALL SUB2(F1,F2,F3)
      ELSE IF(IHELP1 .EQ. 0) THEN
         CALL SUB2(F1,F2,F3)
      ELSE
         IHELP2=1
      END IF
      RETURN
      END
C
C

      SUBROUTINE EXCHG(X,Y,FX,FY)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE exchanges X and Y, and FX and FY.             *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DUMMY=X
      X=Y
      Y=DUMMY
      DUMMY=FX
      FX=FY
      FY=DUMMY
      RETURN
      END
C
C

      SUBROUTINE SUB1(F1,F2,F3)
C
C*****************************************************************
C                                                                *
C  Auxiliary routine for subroutines PEG and PEGK.               *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      G=F3/(F2+F3)
      F1=G*F1
      RETURN
      END
C
C

      SUBROUTINE SUB2(F1,F2,F3)
C
C*****************************************************************
C                                                                *
C  Auxiliary routine for subroutines ANDBJ and ANDBJK.           *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      G=1.0D0-F2/F3
      IF(G .LE. 0.0D0) THEN
         G=0.5D0
      END IF
      F1=G*F1
      RETURN
      END

      INTEGER FUNCTION MACHPD(X)
      DOUBLE PRECISION X
      MACHPD=0
      IF (1.0D0 .LT. X) MACHPD=1
      RETURN
      END


************************************************************************
************************************************************************
*****************************LAMBDA BESTIMMUNG**************************
************************************************************************
************************************************************************
C     uses routine from selfenergies to determine update of rlambd0
      DOUBLE PRECISION FUNCTION Qbstmng(delta)
      IMPLICIT NONE
      INTEGER i,j,k,kk,l,m,n
      INTEGER ngrid,nadd,naux,new,nNeu,mnew,nnn,nvolt,nwv
      INCLUDE 'PARAMETER'
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      INTEGER iscr(ngrid)
      DOUBLE PRECISION omega(ngrid),xxx(ngrid)
      DOUBLE PRECISION delta,exom,beta,wert
      DOUBLE PRECISION ferm0
      DOUBLE PRECISION aux,wert1,wert2,freq,aux2
      DOUBLE PRECISION freq1,freq2,wertn,omex
      DOUBLE PRECISION out
      EXTERNAL ferm0
      COMMON /temp/ beta
      COMMON /whtev/ omega,xxx


      exom=delta                !exom =shift of xxx
      DO j=1,ngrid
         omex=exom+omega(j)     !frequency argument of xxx
         CALL locate(omega,ngrid,omex,m)
         iscr(j)=m              !index of epsilon on eps+exom
      END DO
      aux=0.0d0
      m=iscr(1)
      IF(m.eq.0) m=1
      IF(m.eq.ngrid) m=ngrid-1
      wert1=xxx(m)+(xxx(m+1)-xxx(m))/(omega(m+1)-omega(m))
     &     *(omega(1)+exom-omega(m)) !interpolate A_b
      DO j=1,ngrid-1                !integration loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
                                !if iscr(j)=iscr(j+1),no add. points
         l=iscr(j+1)
         IF((l.eq.0).OR.(l.eq.n)) GOTO 20
         m=iscr(j)
         kk=l-m
         wert2=xxx(l)+(xxx(l+1)-xxx(l))/(omega(l+1)-omega(l))
     &        *(omega(j+1)+exom-omega(l)) !interpolate A_b
         IF(kk.eq.0) THEN
            aux=aux+0.5*(omega(j+1)-omega(j))*
     &           (ferm0(beta*omega(j+1))*wert2+
     &           ferm0(beta*omega(j))*wert1)
c     
         ELSE
            IF((m.eq.0).OR.(m.eq.ngrid)) GOTO 20
            freq1=omega(j)
            DO k=1,kk
c     write(86,*) omega(m+k)-exom,abm(m+k)
               freq2=omega(m+k)-exom
               wertn=xxx(m+k)
               aux=aux+0.5*(freq2-freq1)*
     &              (ferm0(beta*freq2)*wertn+
     &              ferm0(beta*freq1)*wert1)
               freq1=freq2
               wert1=wertn
            END DO
            aux=aux+0.5*(omega(j+1)-freq1)* !add last interval
     &           (ferm0(beta*omega(j+1))*wert2+
     &           ferm0(beta*freq1)*wert1)
         ENDIF
 20      CONTINUE
         wert1=wert2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCnegative freq
      END DO
      Qbstmng=aux/3.1415926535898d0-1.0d0
c      Qbstmng=aux-1.0d0
      RETURN
      END

***************************************************************************
***************************************************************************
      SUBROUTINE RUPDATE(sefp,sefr,omega,beta,rlambd0,delta)
      IMPLICIT NONE
      INTEGER n,numit,ihelp1,ihelp2,incl,ierr,ngrid,naux,nadd,nNeu,new
      INTEGER i,nspin,k
      INTEGER mnew,nnn,nvolt,nwv
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535898d0)
      PARAMETER(nspin=2)
      INCLUDE 'PARAMETER'
      PARAMETER (naux=N+(N-1)*nadd,nNeu=2*naux)
      PARAMETER (nnn=nNeu+4*new)
      PARAMETER (nwv=nnn+4*nvolt)
      PARAMETER (ngrid=nwv+4*mnew)
      DOUBLE PRECISION sefp(ngrid),sefr(ngrid)
      DOUBLE PRECISION afp(ngrid),omega(ngrid)
      DOUBLE PRECISION xxx(ngrid),omega0(ngrid),beta0
      DOUBLE PRECISION beta,ABSERR,RELERR,x1,x2,x3,f1,f2,f3
      DOUBLE PRECISION ferm0,Qbstmng
      DOUBLE PRECISION rlambd0,delta
      EXTERNAL ferm0
      EXTERNAL Qbstmng
      COMMON /temp/ beta0
      COMMON /whtev/ omega0,xxx
      
      beta0=beta
      DO i=1,ngrid
         afp(i)=sefp(i)/
     &        ((omega(i)+rlambd0-sefr(i))**2+sefp(i)**2)
      END DO
     

      DO i=1,ngrid
         omega0(i)=omega(i)
         xxx(i)=nspin*afp(i)
      END DO
     
      ABSERR=1.0D-11
      RELERR=1.0d-8
      x1=-1.0d-1
      x2=1.0d-1
      CALL ZERORF(Qbstmng,ABSERR,RELERR,1.0d-99,4000,4,0,
     +                  1D-8,X1,X2,X3,F1,F2,F3,NUMIT,IHELP1,
     +                  IHELP2,INCL,IERR)

      IF(ierr.EQ.0) THEN
         PRINT*,'zero was not found!'
         CALL ROUT('nischt.dat',ngrid,omega,afp)
         STOP
      ELSEIF(ierr.EQ.1) THEN
         delta=x2
      ELSEIF(ierr.EQ.2) THEN
         delta=x2
      ELSEIF(ierr.EQ.3) THEN
         delta=x3
      ELSE
          PRINT*,'there is a problem!'
         print*,ierr
         STOP
      ENDIF
     
******************************************************************
C  IERR     : = 0, zero was not found.                           *
C             = 1, zero lies between X1 and X2  (of the two      *
C                  values the one with the smaller absolute      *
C                  functional value should be chosen as the      *
C                  zero).                                        *
C                  The absolute error of the computed zero is    *
C                  smaller than or equal to ABS(X1-X2).          *
C             = 2, X2 is a zero of FCT: F2=0.0 (machine zero).   *
C             = 3, X3 is a zero with ABS(F3) < 4 * machine       *
C                  constant.                                     *
C             = 4, zero of FCT is at X2 (enclosure interval      *
C                  not definable).                               *
C             = 5, maximum number MAXIT of functional evaluations*
C                  exceeded.                                     *
C             = 6, ABSERR or RELERR are negative, or both are    *
C                  equal to zero, or MAXIT < 1.                  *
******************************************************************
      rlambd0=rlambd0+delta
      RETURN
      END


      SUBROUTINE SHIFT(delta,omega,n,sefp,sefm,sefr,sebp,sebm,sebr)
      IMPLICIT NONE
      INTEGER n,i
      DOUBLE PRECISION delta,omega(n),sefp(n),sefm(n),sefr(n),sebp(n),
     +     sebm(n),sebr(n)
      DOUBLE PRECISION sefp0(n),sefm0(n),sebp0(n),sebm0(n)
      DOUBLE PRECISION omegas(n)

      DO i=1,n
         sefp0(i)=sefp(i)
         sefm0(i)=sefm(i)
         sebp0(i)=sebp(i)
         sebm0(i)=sebm(i)
         write(77,*) omega(i),sefp0(i)
      END DO

      DO i=1,n
         omegas(i)=omega(i)+delta
      END DO
      CALL GINDERPOL(omegas,n,omega,n,sefp0,sefp)

      CALL KRAKRO(n,omega,sefp,sefr)
      CALL KRAKRO(n,omega,sebp,sebr)
      

      RETURN
      END
