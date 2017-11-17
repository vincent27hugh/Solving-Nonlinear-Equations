! Jun 9, 2017 in IVF & VS
! June 1,2017
! May 22-24, 2017
! 20170411 PM
! Mar18,2017
! Task #20170221
! Related to May16/Oct17,2016
! edited in Feb21,2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module point_para
    ! pointer to get parameters
    real(kind=8),pointer :: p_para(:)
    character(len=:), pointer :: p_typen
    real(kind=8),pointer :: p_epu
end module
    
program main_Jun01_Task_1
    ! Include the IMSL lib
    include 'link_fnl_shared.h'
    ! use module
    use point_para

    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! declare paramaters
    character(len=7),parameter :: taskcode = "Jun01T1"
    character(len=3),parameter :: num_para = "V27"
    real(kind=4),parameter :: pi = 3.1415926
    !!!!!!!!!!!!!!!!!!!
    ! declare function
    character(len=10), external :: cell_type
    character(len=10), external :: cell_case
    
    real,intrinsic :: sqrt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! declare variable
    integer I
    ! for the loop
    
    character(len=3), target :: typen
    character(len=10):: casecode
    character(len=:),allocatable :: str_s
    character(len=:),allocatable :: str_fig
    character(len=13) :: str_para 
    
    real(kind=8),target :: epsilon_u
    
    real(kind=8),TARGET :: parameters(15)
    ! parameters=[A;B;r;c_p;beta;phi;delta;sigma;lambda;pstar;b;c_f;mu;B1;B2];
    
    real(kind=4),allocatable :: vt(:)
    real(kind=8),allocatable :: vepsilon_d(:)
    real(kind=8),allocatable :: vepsilon_c(:)
    real(kind=8),allocatable :: vtheta(:)
    real(kind=8),allocatable :: valpha(:)
    logical(kind=4),allocatable :: vexitflag(:)
    
    real(kind=8),allocatable :: vu(:)
    real(kind=8),allocatable :: vhs(:)
    real(kind=8),allocatable :: vhn(:)
    
    real(kind=8) :: pstar
    integer tt,cc
    ! tt is # of type of function F(x)
    ! cc is # of case
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p_para=>parameters
    p_typen=>typen
    P_epu=>epsilon_u
    
    str_para = 'paraJun01_'//num_para
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tt=3
    ! type of F(x)
    loop2: do cc= 8,10
    ! case #
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        typen = cell_type(tt)
        casecode = cell_case(cc)
    
        if (typen=="I") then
            epsilon_u = sqrt(3.0)
        else if (typen=="II") then
            epsilon_u=sqrt(2.0)/(sqrt(pi-2.0))
        else if (typen=="III") then
            epsilon_u=1.0
        else if (typen=="IV") then
            epsilon_u=(3.0-sqrt(3.0))/2.0
        else if (typen=="O") then
            epsilon_u=(pi/4.0)/sin(pi/4.0)
        end if
    
        ! get the parameters
        call MyParameter_Jun01V27(parameters)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! check the pointers and targets
        write(*,"('Parameters are:',15(/f8.3))") parameters
        write(*,"('Parameters are:',15(/f8.3))") P_PARA
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change parameters with different cases
        if (casecode=="pstar") then ! 1
            vt = (/ (0.01*i,i = 18,250) /) 
            str_s="pstar"
        else if (casecode=="b") then ! 2
            vt = (/ (0.01*i,i = 1,100) /)
            str_s = "b"
        else if (casecode=="phi") then ! 3
            vt = (/ (0.01*i,i = 30,95) /)
            str_s = "phi"
        else if (casecode=="sigma") then ! 4
            vt = (/ (0.01*i,i = 1,120) /)
            str_s = "sigma"
        else if (casecode=="beta") then ! 5
            vt = (/ (0.01*i,i = 40,99) /)
            str_s = "beta"
        else if (casecode=="lambda") then ! 6
            vt = (/ (0.01*i,i = 2,63) /)
            str_s = "lambda"
        else if (casecode=="c_f") then ! 7
            vt = (/ (0.01*i,i = 12,130) /)
            str_s = "c_F"
        else if (casecode=="r") then ! 8
            vt = (/ (0.001*i,i = 3,150) /)
            str_s = "r"
        else if (casecode=="c_p") then ! 9
            vt = (/ (0.001*i,i = 10,200) /)
            str_s = "c_P"
        else if (casecode=="delta") then ! 10
            vt = (/ (0.01*i,i = 21,70) /)
            str_s = "delta"
        end if
    
        str_fig=taskcode//'_'//typen//'_'//str_s
    
        open(unit=10,file=str_fig//'_'//str_para//"_ep_d"//".txt")
        open(unit=11,file=str_fig//'_'//str_para//"_ep_c"//".txt")
        open(unit=12,file=str_fig//'_'//str_para//"_theta"//".txt")
        open(unit=13,file=str_fig//'_'//str_para//"_alpha"//".txt")
        open(unit=14,file=str_fig//'_'//str_para//"_exitflag"//".txt")
    
        open(unit=15,file=str_fig//'_'//str_para//"_u"//".txt")
        open(unit=16,file=str_fig//'_'//str_para//"_hs"//".txt")
        open(unit=17,file=str_fig//'_'//str_para//"_hn"//".txt")
        !!!!!!!!!!!!!
        allocate(vepsilon_d(size(vt)))
        allocate(vepsilon_c(size(vt)))
        allocate(vtheta(size(vt)))
        allocate(valpha(size(vt)))
        allocate(vexitflag(size(vt)))
    
        allocate(vu(size(vt)))
        allocate(vhs(size(vt)))
        allocate(vhn(size(vt)))
    
        loop1: do i = 1,size(vt)
            if (casecode=="pstar") then
                parameters(10)=vt(i)
            else if (casecode=="b") then
                parameters(11)=vt(i)
            else if (casecode=="phi") then
                parameters(6)=vt(i)
            else if (casecode=="sigma") then
                parameters(8)=vt(i)
            else if (casecode=="beta") then
                parameters(5)=vt(i)
            else if (casecode=="lambda") then
                parameters(9)=vt(i)
            else if (casecode=="c_f") then
                parameters(12)=vt(i)
            else if (casecode=="r") then
                parameters(3)=vt(i)
            else if (casecode=="c_p") then
                parameters(4)=vt(i)
            else if (casecode=="delta") then
                parameters(7)=vt(i)
            end if
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(*,"('Parameters are:',15(/f8.3))") parameters
            write(*,"('Parameters are:',15(/f8.3))") P_PARA
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call Jun01T1(vepsilon_d(i),vepsilon_c(i),vtheta(i),valpha(i),vexitflag(i))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call fun_u(vepsilon_d(i),vtheta(i),vu(i))
            call fun_hs(vepsilon_c(i),vepsilon_d(i),vtheta(i),valpha(i),vu(i),vhs(i))
            call fun_hn(vepsilon_c(i),vtheta(i),valpha(i),vu(i),vhn(i))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(10,"(f16.8)") vepsilon_d(i)
            write(11,"(f16.8)") vepsilon_c(i)
            write(12,"(f16.8)") vtheta(i)
            write(13,"(f16.8)") valpha(i)
            write(14,"(L4)") vexitflag(i)
        
            write(15,"(f16.8)") vu(i)
            write(16,"(f16.8)") vhs(i)
            write(17,"(f16.8)") vhn(i)
        
        end do loop1
        !!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(vepsilon_d)
        deallocate(vepsilon_c)
        deallocate(vtheta)
        deallocate(valpha)
        deallocate(vexitflag)
    
        deallocate(vu)
        deallocate(vhs)
        deallocate(vhn)
        !!!!!!!!!!
        close(UNIT=10,STATUS='KEEP')
        close(UNIT=11,STATUS='KEEP')
        close(UNIT=12,STATUS='KEEP')
        close(UNIT=13,STATUS='KEEP')
        close(UNIT=14,STATUS='KEEP')
    
        close(UNIT=15,STATUS='KEEP')
        close(UNIT=16,STATUS='KEEP')
        close(UNIT=17,STATUS='KEEP')
    end do loop2
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.FALSE.) then 
        ! output to test
        write(*,"('Task code is: ',A10/,'Parameter code is: ',A10/&
            &,'Type is: ',A3/,'epsilon_u = ',F8.2)") taskcode,num_para,typen,epsilon_u
    
        do i = 1,size(parameters)
            write(*,"('Parameters NO.',I5,' is: ',F8.3/)") I,parameters(I)
        end do
    
        do i = 1,size(vt)
            write(*,"('vt NO.',I5,' is: ',F8.3/)") I,vt(I)
        end do
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pause
end program main_Jun01_Task_1
    
character(len=10) function cell_type(tt)
    implicit none
    integer, intent(in) :: tt
    ! tt is # of type of F(x), which 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (tt==1) then
        cell_type = "I"
    else if (tt==2) then
        cell_type = "II"
    else if (tt==3) then
        cell_type = "III"
    else if (tt==4) then
        cell_type = "IV"
    else if (tt==5) then
        cell_type = "O"
    end if
    
    return
end
    
character(len=10) function cell_case(cc)
    implicit none
    integer, intent(in) :: cc
    !
    if (cc==1) then
        cell_case = "pstar"
    elseif (cc==2) then
        cell_case = "b"
    else if (cc==3) then
        cell_case = "phi"
    else if (cc==4) then
        cell_case = "sigma"
    else if (cc==5) then
        cell_case = "beta"
    else if (cc==6) then
        cell_case = "lambda"
    else if (cc==7) then 
        cell_case = "c_f"
    else if (cc==8) then
        cell_case = "r"
    else if (cc==9) then
        cell_case = "c_p"
    else if (cc==10) then
        cell_case = "delta"
    end if
    return
end