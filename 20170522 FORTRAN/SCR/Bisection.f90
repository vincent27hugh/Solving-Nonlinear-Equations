! Example of Bisection
module NUMERICAL
    implicit none
    real, parameter :: zero = 1d-5
    ! tolerance
    contains
    real function bisect(A,B)
        implicit none
        real A,B
        ! the guess value
        real C
        ! to get (A+B)/2
        real FA
        ! F(A)
        real FB
        ! F(B)
        real FC
        ! F(C)
        
        ! get C & F(C) at first
        C = (A+B)/2.0
        FC = func(C)
        
        ! stop loop when F(C) is smaller than ZERO
        do while (abs(FC)>ZERO)
            FA = func(A)
            FB = func(B)
            if (FA*FC<0) then
                ! F(A)*F(C)<0
                B=C
                C=(A+B)/2.0
            else
                ! F(B)*F(C)<0
                A=C
                C=(A+B)/2.0
            end if
            
            ! get new F(C)
            FC=FUNC(C)
        end do
        
        bisect = C
        return   
    end function
    
    ! the function to solve
    real function func(X)
        implicit none
        real X
        
        FUNC = (X+3)*(X-3)
        
        return
    end function
end module
    
program main
    use NUMERICAL
    implicit none
    real A,B
    ! initial values
    real ANS
    ! the answer
    
    do while (.true.)
        write(*,*) 'Input two guess values:'
        read(*,*) A,B
        
        if (FUNC(A)*FUNC(B)<0) exit
        
        write(*,*) 'No valid guess values'
    end do
    
    ANS = bisect(A,B)
    !
    write(*,"('X=',F6.3)") ANS
    
    pause 
    
    stop
end program main