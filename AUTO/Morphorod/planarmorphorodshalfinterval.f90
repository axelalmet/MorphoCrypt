!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!           A planar morphoelastic rod upon an elastic foundation
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

 SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!---------- ---- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
 DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
 DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)
 DOUBLE PRECISION alpha, gamma, T, K, L, H, sigma, E, S, X, Y, NX, NY, Theta, M 

 ! Define the variables
 S = U(1)
 X = U(2)
 Y = U(3)
 NX = U(4)
 NY = U(5)
 Theta = U(6)
 M = U(7)

 ! Define the parameters
 T = PAR(1)
 K = PAR(2)
 L = PAR(3)
 H = PAR(4)
 ! sigma = 2.d0*((3.d0)**0.5)*0.01/0.015
 ! gamma = DEXP(DEXP(-((S - L)/sigma)**2.0)*T)
 gamma = T
 !E = 1.d0 - 0.75*DEXP(-((S - L)/sigma)**2.0)
 E = 1.d0
 alpha = 1.d0 + (NX*DCOS(Theta) + NY*DSIN(Theta))/E

   ! Define derivatives
   F(1) = L
   F(2) = alpha*gamma*L*DCOS(Theta)
   F(3) = alpha*gamma*L*DSIN(Theta)
   F(4) = K*L*(X-S)
   F(5) = K*L*(Y)
   F(6) = gamma*L*M/E
   F(7) = alpha*gamma*L*(NX*DSIN(Theta) - NY*DCOS(Theta))

 END SUBROUTINE FUNC
!----------------------------------------------------------------------

 SUBROUTINE STPNT(NDIM,U,PAR,x) 
!---------- ----- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(IN) :: x
 DOUBLE PRECISION h, w, L0, L, y0

 L0 = 0.125
 h = 0.015
 w = 0.01
 L = 2.d0*(3.d0)**0.5*L0/h
 H = 0.d0*2.d0*(3.d0)**0.5*L0/h

   PAR(1) = 1.d0 ! Set gamma to one 
   PAR(2) = 0.01  ! Set k
   PAR(3) = 0.5*L  ! Set L  
   PAR(4) = H ! Set H

 END SUBROUTINE STPNT
!----------------------------------------------------------------------

 SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!---------- ---- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
 DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
 DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
 DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
 DOUBLE PRECISION L, H

 L = PAR(3)
 H = PAR(4)

 ! These boundary conditions correspond to the rod being clamped at the boundaries:
 ! S(L) = L, X(0) = 0, X(L) = L, Y(0) = Y(L) = H, Theta(0) = Theta(L) = 0

 FB(1)=U0(1)
 FB(2)=U0(2)
 FB(3)=U1(2) - L
 FB(4)=U0(3) - H
 FB(5)=U1(5)
 FB(6)=U0(6)
 FB(7)=U1(6)

 END SUBROUTINE BCND

!----------------------------------------------------------------------

 SUBROUTINE PVLS(NDIM,U,PAR)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
 DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
 DOUBLE PRECISION, EXTERNAL :: GETP
 DOUBLE PRECISION MAX, MIN, EXTREMUM
 
 MAX = GETP('MAX', 3, U)
 MIN = GETP('MIN', 3, U)

 IF (-MIN > MAX) THEN
   EXTREMUM = -MIN
 ELSE
   EXTREMUM = MAX
 END IF
   
 PAR(5) = EXTREMUM

 END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
 SUBROUTINE ICND
 END SUBROUTINE ICND

 SUBROUTINE FOPT 
 END SUBROUTINE FOPT
!----------------------------------------------------------------------
!----------------------------------------------------------------------