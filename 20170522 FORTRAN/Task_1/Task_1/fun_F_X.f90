real(kind=8) function fun_F_x(x)
    use point_para
    
    implicit none
    
    real(kind=8) :: X
    real(kind=8) :: epsilon_u
    character(len=3) :: typen
    Real(kind=8) ::  erf,sqrt,log
    intrinsic erf,sqrt,log
    
    epsilon_u = p_epu
    typen = p_typen
    
    if (typen=='III') then
        fun_F_x=.5-.5*erf((log(-X+1.0)+.5*log(2.0))/sqrt(2.0*log(2.0)))    
    end if
    
    return
end function