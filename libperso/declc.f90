       REAL(DP), INTENT(IN) :: aa,bb
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funcc) :: funk
       !
       REAL(DP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
