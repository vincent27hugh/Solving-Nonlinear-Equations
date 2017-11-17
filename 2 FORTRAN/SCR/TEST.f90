PROGRAM TEST
    IMPLICIT NONE
    
    real(kind=8) :: y
    
    !!!!!!!!!!!
    real(kind=8),external :: fun_intF,fun_F_x
    real(kind=8) :: B1,B2
    !!!!!!!!!!!!!!!!!!!!
    B1=18.0
    B2=25.0
    
    y=fun_intF(B1,B2)   
    
    WRITE(*,*) fun_F_X(-3.0)
    
    pause
    
END PROGRAM TEST
    
real(kind=8) function fun_intF(a,b) 
    INCLUDE 'link_fnl_shared.h'
    
    implicit none
    
    external F
    
    real(kind=8) ::  a,b
    real RESULT,ERRABS,ERRREL,ERREST
    
    ERRABS=1E-3
    ERRREL=1E-3
    
    CALL QDAGS(F,real(a),real(b),ERRABS,ERRREL,RESULT,ERREST)
    
    fun_intF=RESULT
    
    return
end function

real function F(X)
    implicit none
    Real X
    Real(kind=8) :: fun_F_x
    external fun_F_x
    
    F=1-fun_F_x(X)
    
    return
end function
    
real(kind=8) function fun_F_x(x)
    
    implicit none
    
    real X
    real(kind=8) :: epsilon_u
    character(len=3) :: typen
    Real(kind=8) :: erf,sqrt,log
    intrinsic erf,sqrt,log
    
    epsilon_u = 1.0
    typen = "III"
    
    if (typen=='III') then
        fun_F_x=.5-.5*erf((log(-X+1)+.5*log(2.0))/sqrt(2.0*log(2.0)))    
    end if
    
    return
end function