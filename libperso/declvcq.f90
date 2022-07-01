       REAL(QP), INTENT(IN) :: aa,bb
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT), DIMENSION(m) :: s
       INTEGER(I4B), INTENT(IN) :: m,n
       PROCEDURE(funcvcq) :: funk
       !
       REAL(QP) :: a,b,del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
