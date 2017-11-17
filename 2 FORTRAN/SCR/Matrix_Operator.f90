program main
    
    include 'link_fnl_shared.h'
    
    use lin_sol_gen_int
    use rand_gen_int
    use error_option_packet
    
    implicit none
    integer,parameter :: n=32
    real(kind(1e0)),parameter :: one=1e0
    real(kind(1e0)) err
    real(kind(1e0)) A(n,n),b(n,n),x(n,n),res(n,n),y(n**2)
    
    !  generate a random matrix
    call rand_gen(y)
    A = reshape(y,(/n,n/))
    
    ! generate random right hand sides
    call rand_gen(y)
    b = reshape(y,(/n,n/))
    
    ! compute the solution matrix of Ax=b
    call lin_sol_gen(A,b,x)
    
    ! check the results
    res=b-matmul(A,x)
    err = maxval(abs(res)/sum(abs(A)+abs(b)))
    
    if (err <= sqrt(epsilon(one))) then
        write(*,*) 'Example for LIN_SOL_GEN is correct!'
    end if
    
    pause
    
end program
    