program Test
    
    include 'link_fnl_shared.h'
    
    implicit none
    
    external :: F
    
    real(kind=8) :: y
    real(kind=8) :: RESULT
    real(kind=8),parameter :: ERRABS=1E-5
    real(kind=8),parameter :: ERRREL=1E-5
    real(kind=8) :: ERREST
    
    real(kind=8) :: A,B

    A=18.0
    B=25.0
    
    CALL dQDAGS(F,A,B,ERRABS,ERRREL,RESULT,ERREST)
    
    y = RESULT
    
    write(*,*) y,ERREST
    
    pause
    
end program Test
    
real(kind=8) function F(X)

    implicit none
    real(kind=8) :: X
    Real(kind=8) :: F2
    EXTERNAL F2
    
    F = 1.0-F2(X)
    
    return
end function
    
real(kind=8) function F2(X)
    implicit none
    real(kind=8) :: X
    real(kind=8) ::  sin
    intrinsic sin

    F2=sin(X)
    
    return
end function
    