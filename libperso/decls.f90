       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(IN), DIMENSION(:) :: arg
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funcs) :: funk
       !
       REAL(SP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(SP), DIMENSION(2*3**(n-2)) :: x
