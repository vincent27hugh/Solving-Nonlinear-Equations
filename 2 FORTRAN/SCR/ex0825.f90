program ex0825
    implicit none 
    
    interface 
    ! define the interface of func
    function random10(lbound,ubound)
    implicit none
    real :: lbound,ubound
    real :: random10(10)
    ! the return value is an array
    end function
    end interface
    
    real :: a(10)
    CALL RANDOM_SEED()
    ! subroutine in library
    ! initialize a pseudo-random number sequence
    
    a = random10(1.0,10.0)
    ! generare 10 random numbers between 1.0 and 10.0
    write(*,"(10F6.2)") a
    pause
end

! random10 would return10 random numbers
    function random10(lbound,ubound)
    implicit none
    real :: lbound,ubound
    real :: len
    real :: random10(10)
    real t
    integer i
    
    len = ubound - lbound
    ! calculate the size of scope
    do i = 1,10
        call RANDOM_NUMBER(t)
        ! t is random numbers in 0~1
        random10(i) = lbound + len*t
        ! convert t into the random numbers between lbound and ubound 
        
    end do
    
    return
    end