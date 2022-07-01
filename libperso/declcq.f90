       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funccq) :: funk
       !
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
