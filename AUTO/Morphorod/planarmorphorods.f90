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
 DOUBLE PRECISION gamma, K, L, t, g, S, X, Y, NX, NY, Theta, M, alpha

 ! Define the variables
 S = U(1)
 X = U(2)
 Y = U(3)
 NX = U(4)
 NY = U(5)
 Theta = U(6)
 M = U(7)

 ! Define the parameters
 t = PAR(1)
 K = PAR(2)
 L = PAR(3)
 g = PAR(4)
  gamma = 1 + t
  ! alpha = 1.d0 + (15.d0)**2.0/(12.d0*(125.d0)**2.0)*NX*COS(Theta) + NY*SIN(Theta)
  ! alpha = 1.d0

   ! Define derivatives
   F(1) = L
   F(2) = gamma*L*COS(Theta)
   F(3) = gamma*L*SIN(Theta)
   F(4) = K*gamma*L*(X - S)
   F(5) = K*gamma*L*(Y)
   F(6) = gamma*L*M
   F(7) = gamma*L*(NX*DSIN(Theta) - NY*DCOS(Theta))

 END SUBROUTINE FUNC
!----------------------------------------------------------------------

 SUBROUTINE STPNT(NDIM,U,PAR,x) 
!---------- ----- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(IN) :: x
 DOUBLE PRECISION L, h, w, K, f

 h = 15.d0
 w = 10.d0
 L = 125.d0
 f = 0.01

 K = f*12.d0*(L**4.d0)/(w*(h**3.d0))
 ! K = 850.d0

   PAR(1) = 0.d0 ! Set gamma to one 
   PAR(2) = K ! Set k
   PAR(3) = 0.5  ! Set L  
   PAR(4) = 0.16 ! Set g

 END SUBROUTINE STPNT
!----------------------------------------------------------------------

 SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!---------- ---- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
 DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
 DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
 DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
 DOUBLE PRECISION L

 ! These boundary conditions correspond to the rod being clamped at the boundaries:
 ! S(L) = L, X(0) = 0, X(L) = L, Y(0) = Y(L) = H, Theta(0) = Theta(L) = 0

 ! Half-clampeds replace Y(0) = 0 with NY(0) = 0

 FB(1)=U0(1)
 FB(2)=U0(2)
 FB(3)=U1(2) - 0.5
 FB(4)=U0(5)
 FB(5)=U1(3)
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