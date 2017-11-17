      program main
        INCLUDE 'link_fnl_shared.h'  
        USE QAND_INT
        USE UMACH_INT
        IMPLICIT NONE
        INTEGER I, J, MAXFCN, N, NOUT
        REAL A(3), B(3), CNST, ERRABS, ERREST, ERRREL, F, RESULT
        EXTERNAL F
        ! Get output unit number
        CALL UMACH (2, NOUT)
        !
        N = 3
        MAXFCN = 100000
        ! Set error tolerances
        ERRABS = 0.0001
        ERRREL = 0.001
        !
        DO 20 I=1, 6
        CNST = I/2.0
        ! Set limits of integration
        ! As CNST approaches infinity, the
        ! answer approaches PI**1.5
        DO 10 J=1, 3
        A(J) = -CNST
        B(J) = CNST
     10 CONTINUE
        CALL QAND (F, N, A, B, RESULT, ERRABS, ERRREL, MAXFCN, ERREST)
        WRITE (NOUT,99999) CNST, RESULT, ERREST
     20 CONTINUE
  99999 FORMAT (1X, 'For CNST = ', F4.1, ', result = ', F7.3, ' with ','error estimate ', 1PE10.3)
	  END
        !
      REAL FUNCTION F (N, X)
        INTEGER N
        REAL X(N)
        REAL EXP
        INTRINSIC EXP
        F = EXP(-(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
        RETURN
      END