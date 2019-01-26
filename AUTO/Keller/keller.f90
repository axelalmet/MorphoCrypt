!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!          Flaherty et al. (1972) model of a buckling, inextensible, planar ring
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

 SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!---------- ---- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
 DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
 DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,*), DFDP(NDIM,*)
 DOUBLE PRECISION P, L, Q, N, Theta, K, alpha

 ! Define the variables
 Q = U(1)
 N = U(2)
 Theta = U(3)
 K = U(4)
 
 ! Define the parameters
 P = PAR(1)
 L = PAR(2)

   ! Define derivatives
   F(1) = -L*K*N
   F(2) = K*Q + P
   F(3) = L*K
   F(4) = L*N

 END SUBROUTINE FUNC
!----------------------------------------------------------------------

 SUBROUTINE STPNT(NDIM,U,PAR,x) 
!---------- ----- 

 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NDIM
 DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
 DOUBLE PRECISION, INTENT(IN) :: x
 DOUBLE PRECISION PI, P, L

    P = 0.d0
    PI = 4.d0*DATAN(1.d0)
    L = 0.5*PI

    U(1) = P
    U(2) = 0.d0
    U(3) = L*x
    U(4) = 1.d0

    PAR(1) = P ! Set P 
    PAR(2) = L  ! Set L

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

 L = PAR(2)

 FB(1)=U0(1)
 FB(2)=U1(1)
 FB(3)=U0(3)
 FB(4)=U1(3) - L

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
   
 PAR(3) = EXTREMUM
 PAR(4) = GETP('BV0', 4, U)

 END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
 SUBROUTINE ICND
 END SUBROUTINE ICND

 SUBROUTINE FOPT 
 END SUBROUTINE FOPT
!----------------------------------------------------------------------
!----------------------------------------------------------------------