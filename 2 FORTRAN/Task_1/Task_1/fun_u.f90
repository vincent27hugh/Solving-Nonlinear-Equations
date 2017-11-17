subroutine fun_u(epsilon_d,theta,u)
    USE point_para
    
    IMPLICIT NONE
    
    REAL(kind=8),intent(in) :: epsilon_d,theta
    REAL(kind=8),intent(out) :: u
    REAL(kind=8),external :: funqtheta,fun_F_x
    ! Local variables
    REAL(kind=8) :: q_theta,F_epsilond,A,B1,B2,lambda
    
    !!!!!!!!!!!!!  
    A = p_para(1)
    B1 = p_para(14)
    B2 = p_para(15)
    lambda = p_para(9)
    
    q_theta = funqtheta(theta,A,B1,B2)
    F_epsilond = fun_F_x(epsilon_d)
    u = lambda*F_epsilond/(lambda*F_epsilond+theta*q_theta)
    
    
    return
end subroutine