! Jun 9, 2017 in IVF & VS
! June 1,2017
! May 22-24, 2017
! 20170411 PM
! Mar18,2017
! Task #20170221
! Related to May16/Oct17,2016
! edited in Feb21,2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MyParameter_Jun01V27(parameters)
    implicit none
    real(kind=8) :: parameters(15)
    ! internal parameters
    
    ! A
    parameters(1)=1.355
    ! B1
    parameters(14)=18.0
    ! B2
    parameters(15)=25.0

    ! B
    parameters(2)=parameters(14)/parameters(15)
    ! r
    parameters(3) = 0.012
    ! c_p
    parameters(4) = 0.06
    ! beta
    parameters(5) = 0.72
    ! phi
    parameters(6) = 0.6
    ! delta
    parameters(7) = 0.53
    ! sigma
    parameters(8) = 0.06
    ! lambda
    parameters(9) = 0.1
    ! pstar
    parameters(10) = 1.0
    ! b
    parameters(11) = 0.2
    ! c_f
    parameters(12) = 0.36
    ! mu
    parameters(13) = 0.08
    !!!!!!!!!!!!!!!!!!!!
    
end subroutine MyParameter_Jun01V27