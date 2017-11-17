program ex1007
    implicit none    
    integer, target :: a(8) = (/10,15,8,25,9,20,17,19/)
    integer, pointer :: p(:)
    
    interface
    function getmin(p)
    integer,pointer :: p(:)
    integer,pointer :: getmin
    end function
    end interface
    
    p=>a(1:8:2)
    
    write(*,*) getmin(p)
    
    pause
    stop
    
end program ex1007
    
function getmin(p)
    implicit none
    integer,pointer :: p(:)
    integer, pointer :: getmin
    integer i,s
    integer min
    
    s=size(p,1)
    
    min = 2**30
    ! set min as a really large number
    do i = 1,s
        if (min>p(i)) then
            min = p(i)
            getmin=>p(i)
        end if
        
    end do
    
    return
end function