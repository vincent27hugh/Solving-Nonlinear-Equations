program main

  INCLUDE 'link_fnl_shared.h' 
  use numerical_libraries
  use lin_sol_gen_int
  
  implicit none
  real :: A(3,3) = (/ 1,3,2,1,2,1,2,1,3 /)
  real :: B(3,1) = (/ 4,6,6 /)
  real :: X(3,1)

  call lin_sol_gen(A,B,X) ! A*X=B,Ω‚X
  write(*,"(3F5.2)") X

  stop
end program