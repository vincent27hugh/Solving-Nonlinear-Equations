PROGRAM ex0815
    implicit none
    real :: a=1.0
    call ShowInteger(a)
    call ShowReal(a)
    
    pause
    stop
    
    END PROGRAM ex0815
    
subroutine ShowInteger(num)
    implicit none
    integer :: num
    write(*,*) num
    return
end 
    
subroutine ShowReal(num)
    implicit none
    real :: num
    write(*,*) num
    
    return
end