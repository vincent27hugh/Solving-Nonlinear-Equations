program Test
    implicit none
    real(kin=8) :: x,y
    integer(kind=4) :: B1,B2
    
    x=10
    B1=18
    B2=25
    
    y=x**(B1/B2)
    
    write(*,*) y
    
end program Test
    