module func
    integer,pointer :: p(:)  
end module
    
program ex1008
    use func
    implicit none
    integer,target :: a(8) = (/10,15,8,25,9,20,17,19/)

    p=>a(:)
    
    write(*,"('Pointers are :',8(/I2))") p
    write(*,"('Targets are:',8(/I2))") a
    
    a(2)=1
    
    call test
    
    write(*,"('Pointers are :',8(/I2))") p
    write(*,"('Targets are:',8(/I2))") a
     
    pause
    
    stop
    
    end program ex1008
    
subroutine test
    use func
    implicit none
    integer :: t
    
    write(*,"('Pointers are :',8(/I2))") p
    
    t = p(3)
    
    write(*,"('Test are :',(/I2))") t
    
    p(2)=3
    
end