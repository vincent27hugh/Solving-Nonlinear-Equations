real(kind=8) function funqtheta(theta,A,B1,B2)
    implicit none
    
    real(kind=8) :: theta,A,B1,B2
    
    funqtheta = A*theta**(-B1/B2)
    
    return
end function