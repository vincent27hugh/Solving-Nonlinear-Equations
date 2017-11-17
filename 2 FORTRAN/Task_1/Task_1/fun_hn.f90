subroutine fun_hn(epsilon_c,theta,alpha,u,h_n)
    USE point_para
    IMPLICIT NONE
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(kind=8),intent(in) :: epsilon_c,theta,alpha,u
    REAL(kind=8),intent(out) :: h_n
    REAL(kind=8),external :: funqtheta,fun_F_x
    ! LOCAL VARIABLES
    REAL(kind=8) :: q_theta,F_epsilonc,temp1,temp2,A,B1,B2,lambda
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    A = p_para(1)
    B1 = p_para(14)
    B2 = p_para(15)
    lambda = p_para(9)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    q_theta = funqtheta(theta,A,B1,B2)
    F_epsilonc = fun_F_x(epsilon_c)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    temp1 = (1-u)*(alpha)*lambda*(F_epsilonc)
    temp2 = lambda*F_epsilonc+(1-alpha)*theta*q_theta
    
    h_n =temp1/temp2
    return
    end subroutine