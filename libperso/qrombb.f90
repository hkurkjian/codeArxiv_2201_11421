MODULE qrombb
 USE nrtype; USE nrutil, ONLY :arth
 IMPLICIT NONE
       INTERFACE trapz
        MODULE PROCEDURE trapzs,trapzd,trapzq,trapzc,trapzcq
       END INTERFACE trapz
       INTERFACE qromb
        MODULE PROCEDURE qrombs,qrombd,qrombq,qrombc,qrombcq
       END INTERFACE qromb
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzs(func,a,b,ar,s,n)
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(IN), DIMENSION(:) :: ar
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(SP) :: del,fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_sp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_sp*del,del,it),ar))
       s=0.5_sp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzd(func,a,b,ar,s,n)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(DP) :: del,fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_dp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_dp*del,del,it),ar))
       s=0.5_dp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzq(func,a,b,ar,s,n)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del,fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_qp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_qp*del,del,it),ar))
       s=0.5_qp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzc(func,a,b,ar,s,n)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(DP) :: del
       COMPLEX(DPC) :: fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_dp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_dp*del,del,it),ar))
       s=0.5_dp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE trapzcq(func,a,b,ar,s,n)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
!       This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of b a f(x)dx. Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!       additional interior points. s should not be modified between sequential calls.
       REAL(QP) :: del
       COMPLEX(QPC) :: fsum
       INTEGER(I4B) :: it
       if (n == 1) then
       s=0.5_qp*(b-a)*sum(func( (/ a,b /) ,ar))
       else
       it=2**(n-2)
       del=(b-a)/it !This is the spacing of the points to be added.
       fsum=sum(func(arth(a+0.5_qp*del,del,it),ar))
       s=0.5_qp*(s+del*fsum) !This replaces s by its refined value.
       end if
       END SUBROUTINE trapzcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombs(func,a,b,ar)
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(IN), DIMENSION(:) :: ar
       REAL(SP) :: qrombs
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(SP), PARAMETER :: EPS=1.0e-7_sp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(SP), DIMENSION(JMAXP) :: h,s !These store the successive trapezoidal approximations and their relative stepsizes
       REAL(SP) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombs'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qrombs,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombs)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_sp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombd(func,a,b,ar)
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       REAL(DP) :: qrombd
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP), PARAMETER :: EPS=1.0e-6_dp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(DP), DIMENSION(JMAXP) :: h,s !These store the successive trapezoidal approximations and their relative stepsizes
       REAL(DP) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombd'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qrombd,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombd)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_dp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombq(func,a,b,ar)
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       REAL(QP) :: qrombq
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP), PARAMETER :: EPS=1.0e-6_qp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(QP), DIMENSION(JMAXP) :: h,s !These store the successive trapezoidal approximations and their relative stepsizes
       REAL(QP) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombq'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qrombq,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombq)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_qp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombc(func,a,b,ar)
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(DPC) :: qrombc
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP), PARAMETER :: EPS=1.0e-8_dp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(DP), DIMENSION(JMAXP) :: h !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(DPC), DIMENSION(JMAXP) :: s !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(DPC) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombc'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qrombc,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombc)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_dp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qrombcq(func,a,b,ar)
       INTERFACE
        FUNCTION func(x,arg)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: func
        END FUNCTION func
       END INTERFACE
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: ar
       COMPLEX(QPC) :: qrombcq
       INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP), PARAMETER :: EPS=1.0e-8_qp
!       Returns the integral of the function func from a to b. Integration is performed by Romberg’s
!       method of order 2K, where, e.g., K=2 is Simpson’s rule.
!       Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation error
!       estimate; JMAX limits the total number of steps; K is the number of points used in the
!       extrapolation.
       REAL(QP), DIMENSION(JMAXP) :: h !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(QPC), DIMENSION(JMAXP) :: s !These store the successive trapezoidal approximations and their relative stepsizes
       COMPLEX(QPC) :: dqromb 
       INTEGER(I4B) :: j
!       write(6,*)'qrombcq'
       h(1)=1.0
       do j=1,JMAX
        call trapz(func,a,b,ar,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qrombcq,dqromb)
         if (abs(dqromb) <= EPS*abs(qrombcq)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_qp*h(j) !This is a key step: The factor is 0.25 even
!        though the stepsize is decreased by only
!        0.5. This makes the extrapolation a polynomial
!        in h2 as allowed by equation (4.2.1),
!        not just a polynomial in h.
       end do
       STOP 'qromb:nbr d itérations dépassé'
       END FUNCTION qrombcq
END MODULE qrombb
