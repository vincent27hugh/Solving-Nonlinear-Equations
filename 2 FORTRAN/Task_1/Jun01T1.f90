! Jun 9, 2017 in IVF & VS
! June 1,2017
! May 22-24, 2017
! 20170411 PM
! Mar18,2017
! Task #20170221
! Related to May16/Oct17,2016
! edited in Feb21,2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Jun01T1(epsilond,epsilonc,theta,alpha,exitflag)
    INCLUDE 'link_fnl_shared.h'    
    USE NEQNF_INT
    
    implicit none
    
    real(kind=8) :: epsilond
    real(kind=8) :: epsilonc
    real(kind=8) :: theta
    real(kind=8) :: alpha
    logical(kind=4) :: exitflag
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The number of equations to be solved. for NEQNF
    integer(kind=4),parameter :: N=4
    real(kind=8):: fnorm, X(N), XGUESS(N)
    ! fnorm: A scalar that has the value F(1)2 + … + F(N)2 at the point X.   (Output of NEQNF)
    ! X(N): Solution of nonlinear systems; A vector of length N.   (Output of NEQNF) 
    ! X contains the best estimate of the root found by NEQNF.
    ! XGUESS: A vector of length N.   (Input of NEQNF) 
    real(kind=8),parameter:: ERRREL = 1d-6
    ! ERRREL:Tolerance; Stopping criterion;  (Input of NEQNF)
    integer(kind=4),PARAMETER:: ITMAX = 1500
    ! The maximum allowable number of iterations.   (Input) 
    
    EXTERNAL FCN
    ! User-supplied SUBROUTINE to evaluate the system of equations to be solved. 
    ! The usage is CALL fun_solveJun01T1(X, F, N)
    
    !real ceq(N)
    !  A vector of length N. FVEC contains the functions evaluated at the point X.
    
    DATA XGUESS /-5.0,-1.0,1.0,0.5/
    
    ! Routine NEQNF is based on the MINPACK subroutine HYBRD1, 
    ! which uses a modification of M.J.D. Powell’s hybrid algorithm.
    CALL DNEQNF(FCN,ERRREL,N,ITMAX,XGUESS,X,fnorm)
    
    write(*,*) X
    
    epsilond = X(1)
    epsilonc = X(2)
    theta = X(3)
    alpha = X(4)
    
    ! Check whether the solution is valid
    if (fnorm<=ERRREL) then
        exitflag = .true.
    else
        exitflag = .false.
    end if
    
    return
end subroutine Jun01T1