program main  
    INCLUDE 'link_fnl_shared.h'

    IMPLICIT   NONE
!                                 Declare variables
    INTEGER(kind=4) :: K
    INTEGER(kind=4),PARAMETER :: N=3
    iNTEGER(kind=4),PARAMETER :: ITMAX = 500
!
    REAL(kind=8) :: FNORM, X(N), XGUESS(N)
    REAL(kind=8),PARAMETER :: ERRREL=1E-5
    EXTERNAL :: FCN
!                                 Set values of initial guess
!                                 XGUESS = (  4.0  4.0  4.0 )
    DATA XGUESS /4.0, 4.0, 4.0/
!
!                                 Find the solution
    CALL DNEQNF(FCN,ERRREL,N,ITMAX,xguess, x, fnorm)
!                                 Output
    WRITE(*,"('The solution to the system is:',/,' X = (',3F5.1,&
        ')',/,' with FNORM = ',F5.4,//)") (X(K),K=1,N), FNORM
    
    pause
!
END PROGRAM MAIN
