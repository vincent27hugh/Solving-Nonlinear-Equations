real(kind=8) function fun_intF(a,b) 
    INCLUDE 'link_fnl_shared.h'
    
    implicit none
    
    external :: F
    
    real(kind=8) :: a,b,RESULT
    real(kind=8),parameter :: ERRABS=1E-5
    real(kind=8),parameter :: ERRREL=1E-5
    real(kind=8) :: ERREST
    
    CALL dQDAGS(F,a,b,ERRABS,ERRREL,RESULT,ERREST)
    
    fun_intF=RESULT
    
    return
end function

real(kind=8) function F(X)
    implicit none
    Real(kind=8) :: X
    Real(kind=8) :: fun_F_x
    external fun_F_x
    
    F=1.0-fun_F_x(X)
    
    return
end function