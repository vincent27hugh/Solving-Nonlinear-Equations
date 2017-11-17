subroutine fun_hs(epsilon_c,epsilon_d,theta,alpha,u,h_s)
    USE point_para
    
    implicit none
    
    REAL(kind=8),intent(in) :: epsilon_c,epsilon_d,theta,alpha,u
    REAL(kind=8),intent(out) :: h_s
    REAL(kind=8),external :: funqtheta,fun_F_x
    ! LOCAL VARIABLES
    REAL(kind=8) :: q_theta,F_epsilonc,F_epsilond,temp1,temp2,A,B1,B2,lambda
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    A = p_para(1)
    B1 = p_para(14)
    B2 = p_para(15)
    lambda = p_para(9)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    q_theta = funqtheta(theta,A,B1,B2)
    F_epsilonc = fun_F_x(epsilon_c)
    F_epsilond = fun_F_x(epsilon_d)
    
    temp1 = (1-u)*(1-alpha)*lambda*(F_epsilonc-F_epsilond)
    temp2 = lambda*F_epsilonc+(1-alpha)*theta*q_theta
    
    h_s =temp1/temp2
    
    return
end subroutine
    