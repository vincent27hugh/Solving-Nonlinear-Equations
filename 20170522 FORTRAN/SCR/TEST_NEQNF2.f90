!                                 User-defined subroutine
SUBROUTINE FCN (X, F, N)
    INTEGER(kind=4) :: N
    REAL(kind=8)::X(N), F(N)
!
    REAL(kind=8) ::   EXP, SIN
    INTRINSIC  EXP, SIN
!
    F(1) = X(1) + EXP(X(1)-1.0) + (X(2)+X(3))*(X(2)+X(3)) - 27.0
    F(2) = EXP(X(2)-2.0)/X(1) + X(3)*X(3) - 10.0
    F(3) = X(3) + SIN(X(2)-2.0) + X(2)*X(2) - 7.0
    
    RETURN
END SUBROUTINE FCN