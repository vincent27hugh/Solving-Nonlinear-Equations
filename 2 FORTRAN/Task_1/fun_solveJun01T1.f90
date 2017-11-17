subroutine FCN (X,F,N)
    use point_para
    
    implicit none
    
    integer(KIND=4) ::  N
    real(kind=8):: X(N),F(N)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=8) :: A,B1,B2,r,c_p,beta,phi,delta,sigma,lambda,pstar,b,c_f,epsilon_u,mu
    real(kind=8) :: p
    real(kind=8) :: q_theta,int_Fdu,int_Fcu
    real(kind=8) :: part11,part12,temp211,temp212,temp213,part21,temp2211,temp2212,temp221,part22,&
        part31,part32,temp321,temp322,part41,part42,temp421,temp422
    character(len=3) :: typen
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=8) :: fun_intF,fun_F_X,funqtheta
    
    external :: funqtheta,fun_intF,fun_F_X
    
    REAL       EXP, SIN
    INTRINSIC  EXP, SIN
    
    A = p_para(1)

    B1= p_para(14)
    B2 = p_para(15)
    
    r = p_para(3)
    c_p = p_para(4)
    beta = p_para(5)
    phi = p_para(6)
    delta = p_para(7)
    sigma = p_para(8)
    lambda = p_para(9)
    pstar = p_para(10)
    b = p_para(11)
    c_f = p_para(12)
    mu = p_para(13)
   
    typen = p_typen
    epsilon_u=p_epu
    
    p=pstar
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    q_theta = funqtheta(X(3),A,B1,B2)

    int_Fdu = fun_intF(X(1),epsilon_u)
    int_Fcu = fun_intF(X(2),epsilon_u)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    part11 = X(1)+lambda*int_Fdu/(r+lambda+X(3)*q_theta)

    part12 = (b-p)/sigma
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    temp211 = delta*(r+lambda-phi*(1-fun_F_x(X(2))))*(X(2)-X(1))
    temp212 = r+lambda+X(3)*q_theta
    temp213 = (lambda/(r+lambda)-lambda*delta/temp212)*int_Fcu
    part21 = X(2)-temp211/((1-phi)*temp212)+temp213

    temp2211 = ((beta+phi*(1-beta))*c_f*(1-X(4)))/((1-beta)*(1-phi))
    temp2212 = (beta*c_p*X(4))/(1-beta)
    temp221 = X(3)*(temp2211 + temp2212)
    part22 = delta*X(1)+((1-delta)*(b-p)+temp221)/delta
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    part31 = 1/q_theta

    temp321 = (1-beta)*(1-phi)/c_f
    temp322 = sigma*(epsilon_u-X(2))/(r+lambda)+&
        delta*sigma*(X(2)-X(1))/((1-phi)*(r+lambda+X(3)*q_theta))
    part32 = temp321*temp322
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    part41 = 1/q_theta

    temp421 = (1-beta)*delta*sigma*(epsilon_u-X(1))
    temp422 = c_p*(r+lambda+X(3)*q_theta)
    part42 = temp421/temp422

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    F(1) = part11-part12
    F(2) = part21-part22
    F(3) = part31-part32
    F(4) = part41-part42
    
    return
end SUBROUTINE fcn
    
