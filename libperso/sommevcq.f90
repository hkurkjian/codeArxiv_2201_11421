       if (n == 1) then
        s(:)=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg,m ),dim=1)
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it)
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s(:)=s(:)/3.0_qp+del*sum(func(x,arg,m),dim=1)
       end if
       CONTAINS
        FUNCTION func(x,arg,mm)
        INTEGER, INTENT(IN) ::  mm
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x),mm) :: func
        INTEGER is
