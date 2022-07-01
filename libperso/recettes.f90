!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE recettes 
        USE nrtype
        USE nrutil, ONLY : assert
        IMPLICIT NONE
        INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
        INTERFACE tri
                MODULE PROCEDURE tri_s,tri_d,tri_q
        END INTERFACE
        INTERFACE tri_pos
                MODULE PROCEDURE tri_pos_q
        END INTERFACE
        INTERFACE indexx
                MODULE PROCEDURE indexx_sp,indexx_dp,indexx_qp
        END INTERFACE
        INTERFACE laguer
                MODULE PROCEDURE laguer_d,laguer_q
        END INTERFACE
        INTERFACE zroots
                MODULE PROCEDURE zroots_d,zroots_q
        END INTERFACE
        INTERFACE rf
                 MODULE PROCEDURE rf_s,rf_d,rf_q,rf_c,rf_cq
        END INTERFACE
        INTERFACE rd
                 MODULE PROCEDURE rd_s,rd_d,rd_q,rd_c,rd_cq
        END INTERFACE
        INTERFACE rc
                 MODULE PROCEDURE rc_q
        END INTERFACE
        INTERFACE rj
                 MODULE PROCEDURE rj_q
        END INTERFACE
        INTERFACE elle
                 MODULE PROCEDURE elle_s,elle_d,elle_q,elle_cq,elle_ccq,elleim_q
        END INTERFACE
        INTERFACE ellf
                 MODULE PROCEDURE ellf_s,ellf_d,ellf_q,ellf_cq,ellf_ccq,ellfim_q
        END INTERFACE
        INTERFACE ellfmath
                 MODULE PROCEDURE ellfmath_s,ellfmath_d,ellfmath_q
        END INTERFACE
        INTERFACE ellemath
                 MODULE PROCEDURE ellemath_s,ellemath_d,ellemath_q
        END INTERFACE
        INTERFACE ellpi
                 MODULE PROCEDURE ellpi_q,ellpi_cq
        END INTERFACE
        INTERFACE argum
                 MODULE PROCEDURE argum_d,argum_q
        END INTERFACE
!        INTERFACE arth
!                MODULE PROCEDURE arth_r, arth_d, arth_q, arth_i
!        END INTERFACE
        INTERFACE rtsafe
                MODULE PROCEDURE rtsafed, rtsafeq
        END INTERFACE
        INTERFACE mnewt
                MODULE PROCEDURE mnewt_s, mnewt_d, mnewt_q
        END INTERFACE
        INTERFACE newt
                MODULE PROCEDURE newt_q
        END INTERFACE
        INTERFACE ludcmp
                MODULE PROCEDURE ludcmp_s, ludcmp_d, ludcmp_q, ludcmp_cq
        END INTERFACE
        INTERFACE lubksb
                MODULE PROCEDURE lubksb_s, lubksb_d, lubksb_q
        END INTERFACE
        INTERFACE lnsrch
                MODULE PROCEDURE lnsrch_q
        END INTERFACE
        INTERFACE fdjac
                MODULE PROCEDURE fdjac_q
        END INTERFACE
        INTERFACE fmin
                MODULE PROCEDURE fmin_s,fmin_d,fmin_q
        END INTERFACE
        ABSTRACT INTERFACE
         FUNCTION funcvs(x)
         USE nrtype
         REAL(SP), DIMENSION(:), INTENT(IN) :: x
         REAL(SP), DIMENSION(size(x)) :: funcvs
         END FUNCTION funcvs
         FUNCTION funcvd(x)
         USE nrtype
         REAL(DP), DIMENSION(:), INTENT(IN) :: x
         REAL(DP), DIMENSION(size(x)) :: funcvd
         END FUNCTION funcvd
         FUNCTION funcvq(x)
         USE nrtype
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(size(x)) :: funcvq
         END FUNCTION funcvq
        END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
        SUBROUTINE fmin_s(x,f,fm,func)
        REAL(SP), DIMENSION(:), INTENT(IN)  :: x
        REAL(SP), DIMENSION(size(x)), INTENT(OUT) :: f
        REAL(SP), INTENT(OUT) :: fm
!        Returns f = 1/2(F · F) at x. FUNCTION funcv(x) is a fixed-name, user-supplied routine that
!        returns the vector of functions at x.
        PROCEDURE(funcvs) func
        f=func(x)
        fm=0.5_sp*dot_product(f,f)
        END SUBROUTINE fmin_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE fmin_d(x,f,fm,func)
        REAL(DP), DIMENSION(:), INTENT(IN)  :: x
        REAL(DP), DIMENSION(size(x)), INTENT(OUT) :: f
        REAL(DP), INTENT(OUT) :: fm
        PROCEDURE(funcvd) func
        f=func(x)
        fm=0.5_dp*dot_product(f,f)
        END SUBROUTINE fmin_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE fmin_q(x,f,fm,func)
        REAL(QP), DIMENSION(:), INTENT(IN)  :: x
        REAL(QP), DIMENSION(size(x)), INTENT(OUT) :: f
        REAL(QP), INTENT(OUT) :: fm
        PROCEDURE(funcvq) func
        f=func(x)
        fm=0.5_qp*dot_product(f,f)
        END SUBROUTINE fmin_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE newt_q(x,tolf,check,func)
        USE nrutil, ONLY :nrerror,vabs
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: x
        REAL(QP), INTENT(IN) :: tolf
        LOGICAL, INTENT(OUT) :: check
        PROCEDURE(funcvq) func

        INTEGER, PARAMETER :: MAXITS=200
        REAL(QP), PARAMETER :: TOLMIN=1.0e-6_qp,TOLX=epsilon(x),STPMX=100.0_qp
!       Given an initial guess x for a root in N dimensions, find the root by a globally convergent
!       Newton's method. The length N vector of functions to be zeroed, called fvec in the routine
!       below, is returned by a user-supplied routine that must be called funcv and have the
!       declaration FUNCTION funcv(x). The output quantity check is false on a normal return
!       and true if the routine has converged to a local minimum of the function fmin defined
!       below. In this case try restarting from a different initial guess.
!       Parameters: MAXITS is the maximum number of iterations; TOLF sets the convergence
!       criterion on function values; TOLMIN sets the criterion for deciding whether spurious convergence
!       to a minimum of fmin has occurred; TOLX is the convergence criterion on δx;
!       STPMX is the scaled maximum step length allowed in line searches.
        INTEGER(I4B) :: its
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(QP) :: d,f,fold,stpmax
        REAL(QP), DIMENSION(size(x)) :: g,p,xold
        REAL(QP), DIMENSION(size(x)) :: fvec
        REAL(QP), DIMENSION(size(x),size(x)) :: fjac
        call fmin(x,fvec,f,func) !fvec is also computed by this call.
        if (maxval(abs(fvec(:))) < 0.01_qp*TOLF) then !Test for initial guess being a root.
         check=.false.  !Use more stringent test than simply TOLF.
         RETURN
        end if
        stpmax=STPMX*max(vabs(x(:)),real(size(x),qp)) !Calculate stpmax for line searches.
        do its=1,MAXITS !Start of iteration loop.
         call fdjac(x,fvec,fjac,func)
         !If analytic Jacobian is available, you can replace the routine fdjac below with your own routine.
         g(:)=matmul(fvec(:),fjac(:,:)) !Compute ∇f for the line search.
         xold(:)=x(:)                   !Store x,
         fold=f                         !and f.
         p(:)=-fvec(:)                  !Right-hand side for linear equations.
         call ludcmp(fjac,indx,d)       !Solve linear equations by LU decomposition.
         call lubksb(fjac,indx,p)
         call lnsrch(xold,fold,g,p,x,f,stpmax,check,fminsimp)
         !lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.
         if (maxval(abs(fvec(:))) < TOLF) then !Test for convergence on function values.
          check=.FALSE. 
          RETURN
         end if
         if (check) then !Check for gradient of f zero, i.e., spurious convergence.
          check=(maxval(abs(g(:))*max(abs(x(:)),1.0_qp)/max(f,0.5_qp*size(x))) < TOLMIN)
          RETURN !Test for convergence on δx.
         end if
         if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_qp)) < TOLX) RETURN
        end do
        call nrerror('MAXITS exceeded in newt')
        CONTAINS 
         FUNCTION fminsimp(x)
         IMPLICIT NONE
         REAL(QP) :: fminsimp
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         call fmin(x,fvec,f,func)
         fminsimp=f
         END FUNCTION fminsimp
        END SUBROUTINE newt_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lnsrch_q(xold,fold,g,p,x,f,stpmax,check,func)
        USE nrutil, ONLY :assert_eq, nrerror,vabs
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(IN) :: xold,g
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: p
        REAL(QP), INTENT(IN) :: fold,stpmax
        REAL(QP), DIMENSION(:), INTENT(OUT) :: x
        REAL(QP), INTENT(OUT) :: f
        LOGICAL, INTENT(OUT) :: check
        INTERFACE
         FUNCTION func(x)
         USE nrtype
         REAL(QP) :: func
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         END FUNCTION func
        END INTERFACE
        REAL(QP), PARAMETER :: ALF=1.0e-4_qp,TOLX=epsilon(x)
!        Given an N-dimensional point xold, the value of the function and gradient there, fold
!        and g, and a direction p, finds a new point x along the direction p from xold where the
!        function func has decreased “sufficiently.” xold, g, p, and x are all arrays of length N.
!        The new function value is returned in f. stpmax is an input quantity that limits the length
!        of the steps so that you do not try to evaluate the function in regions where it is undefined
!        or subject to overflow. p is usually the Newton direction. The output quantity check is
!        false on a normal exit. It is true when x is too close to xold. In a minimization algorithm,
!        this usually signals convergence and can be ignored. However, in a zero-finding algorithm
!        the calling program should check whether the convergence is spurious.
!        Parameters: ALF ensures sufficient decrease in function value; TOLX is the convergence
!        criterion on Δx.
        INTEGER(I4B) :: ndum
        REAL(QP) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
        ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
        check=.false.
        pabs=vabs(p(:))
        if (pabs > stpmax) p(:)=p(:)*stpmax/pabs !Scale if attempted step is too big.
        slope=dot_product(g,p)
        if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
        alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_qp)) !Compute λmin.
        alam=1.0 !Always try full Newton step first.
        do !Start of iteration loop.
         x(:)=xold(:)+alam*p(:)
         f=func(x)
         if (alam < alamin) then !Convergence on Δx. For zero finding, the calling program should verify the convergence.
          x(:)=xold(:)
          check=.true.
          RETURN
         else if (f <= fold+ALF*alam*slope) then !Sufficient function decrease.
          RETURN
         else !Backtrack.
          if (alam == 1.0) then !First time.
           tmplam=-slope/(2.0_qp*(f-fold-slope))
          else !Subsequent backtracks.
           rhs1=f-fold-alam*slope
           rhs2=f2-fold-alam2*slope
           a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
           b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
           if (a == 0.0) then
            tmplam=-slope/(2.0_qp*b)
           else
            disc=b*b-3.0_qp*a*slope
            if (disc < 0.0) then
             tmplam=0.5_qp*alam
            else if (b <= 0.0) then
             tmplam=(-b+sqrt(disc))/(3.0_qp*a)
            else
             tmplam=-slope/(b+sqrt(disc))
            end if
           end if
           if (tmplam > 0.5_qp*alam) tmplam=0.5_qp*alam !λ ≤ 0.5λ1.
          end if
         end if
         alam2=alam
         f2=f
         alam=max(tmplam,0.1_qp*alam) !λ ≥ 0.1λ1.
        end do !Try again.
        END SUBROUTINE lnsrch_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE fdjac_q(x,fvec,df,funcv)
        USE nrtype; USE nrutil, ONLY :assert_eq
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(IN) :: fvec
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: x
        REAL(QP), DIMENSION(:,:), INTENT(OUT) :: df
        PROCEDURE(funcvq) funcv

        REAL(QP), PARAMETER :: petit=1.0e-12_qp
!       Computes forward-difference approximation to Jacobian. On input, x is the point at which
!       the Jacobian is to be evaluated, and fvec is the vector of function values at the point,
!       both arrays of length N. df is the N × N output Jacobian. FUNCTION funcv(x) is a
!       fixed-name, user-supplied routine that returns the vector of functions at x.
!        Parameter: EPS is the approximate square root of the machine precision.
        INTEGER(I4B) :: j,n
        REAL(QP), DIMENSION(size(x)) :: xsav,xph,h
        n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
        xsav=x
        h=petit*abs(xsav)
        where (h == 0.0) h=petit
        xph=xsav+h !Trick to reduce finite precision error.
        h=xph-xsav
        do j=1,n
         x(j)=xph(j)
         df(:,j)=(funcv(x)-fvec(:))/h(j) !Forward difference formula.
         x(j)=xsav(j)
        end do
        END SUBROUTINE fdjac_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mnewt_s(ntrial,x,tolx,tolf,usrfun)
        INTEGER(I4B), INTENT(IN) :: ntrial
        REAL(SP), INTENT(IN) :: tolx,tolf
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
        INTERFACE
        SUBROUTINE usrfun(x,fvec,fjac)
        USE nrtype
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
        REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
        END SUBROUTINE usrfun
        END INTERFACE
!        Given an initial guess x for a root in N dimensions, take ntrial Newton-Raphson steps to
!        improve the root. Stop if the root converges in either summed absolute variable increments
!        tolx or summed absolute function values tolf.
        INTEGER(I4B) :: i
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(SP) :: d
        REAL(SP), DIMENSION(size(x)) :: fvec,p
        REAL(SP), DIMENSION(size(x),size(x)) :: fjac
        do i=1,ntrial
         call usrfun(x,fvec,fjac)
!        User subroutine supplies function values at x in fvec and Jacobian matrix in fjac.
         if (sum(abs(fvec)) <= tolf) RETURN !Check function convergence.
         p=-fvec !Right-hand side of linear equations.
         call ludcmp(fjac,indx,d) !Solve linear equations using LU decomposition.
         call lubksb(fjac,indx,p) 
         x=x+p !Update solution.
         if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
        end do
        END SUBROUTINE mnewt_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mnewt_d(ntrial,x,tolx,tolf,usrfun)
        INTEGER(I4B), INTENT(IN) :: ntrial
        REAL(DP), INTENT(IN) :: tolx,tolf
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
        INTERFACE
        SUBROUTINE usrfun(x,fvec,fjac)
        USE nrtype
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
        REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
        END SUBROUTINE usrfun
        END INTERFACE
        INTEGER(I4B) :: i
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(DP) :: d
        REAL(DP), DIMENSION(size(x)) :: fvec,p
        REAL(DP), DIMENSION(size(x),size(x)) :: fjac
        do i=1,ntrial
         call usrfun(x,fvec,fjac)
         if (sum(abs(fvec)) <= tolf) RETURN !Check function convergence.
         p=-fvec !Right-hand side of linear equations.
         call ludcmp(fjac,indx,d) !Solve linear equations using LU decomposition.
         call lubksb(fjac,indx,p) 
         x=x+p !Update solution.
         if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
        end do
        END SUBROUTINE mnewt_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mnewt_q(ntrial,x,tolx,tolf,usrfun)
        INTEGER(I4B), INTENT(IN) :: ntrial
        REAL(QP), INTENT(IN) :: tolx,tolf
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: x
        INTERFACE
        SUBROUTINE usrfun(x,fvec,fjac)
        USE nrtype
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(IN) :: x
        REAL(QP), DIMENSION(:), INTENT(OUT) :: fvec
        REAL(QP), DIMENSION(:,:), INTENT(OUT) :: fjac
        END SUBROUTINE usrfun
        END INTERFACE
        INTEGER(I4B) :: i
        INTEGER(I4B), DIMENSION(size(x)) :: indx
        REAL(QP) :: d
        REAL(QP), DIMENSION(size(x)) :: fvec,p
        REAL(QP), DIMENSION(size(x),size(x)) :: fjac
        do i=1,ntrial
         call usrfun(x,fvec,fjac)
         if (sum(abs(fvec)) <= tolf) RETURN !Check function convergence.
         p=-fvec !Right-hand side of linear equations.
         call ludcmp(fjac,indx,d) !Solve linear equations using LU decomposition.
         call lubksb(fjac,indx,p) 
         x=x+p !Update solution.
         if (sum(abs(p)) <= tolx) RETURN !Check root convergence.
        end do
        END SUBROUTINE mnewt_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE indexx_sp(arr,index)
        USE nrutil, ONLY :arth,assert_eq,nrerror,swap
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
        INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!        Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantity arr is not changed.
        REAL(SP) :: a
        INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        n=assert_eq(size(index),size(arr),'indexx_sp')
        index=arth(1,1,n)
        jstack=0
        l=1
        r=n
        do
         if (r-l < NN) then
          do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
            if (arr(index(i)) <= a) exit
            index(i+1)=index(i)
           end do
           index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
         else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
           do
            i=i+1
            if (arr(index(i)) >= a) exit
           end do
           do
            j=j-1
            if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) then
           write(6,*)'indexx_sp:NSTACK too small'
           STOP
          endif
          if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
          else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
          end if
         end if
        end do
        CONTAINS
        SUBROUTINE icomp_xchg(i,j)
        INTEGER(I4B), INTENT(INOUT) :: i,j
        INTEGER(I4B) :: swp
        if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
        end if
        END SUBROUTINE icomp_xchg
        END SUBROUTINE indexx_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE indexx_dp(arr,index)
        USE nrutil, ONLY :arth,assert_eq,nrerror,swap
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
        INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!        Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantity arr is not changed.
        REAL(DP) :: a
        INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        n=assert_eq(size(index),size(arr),'indexx_dp')
        index=arth(1,1,n)
        jstack=0
        l=1
        r=n
        do
         if (r-l < NN) then
          do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
            if (arr(index(i)) <= a) exit
            index(i+1)=index(i)
           end do
           index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
         else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
           do
            i=i+1
            if (arr(index(i)) >= a) exit
           end do
           do
            j=j-1
            if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) then
           write(6,*)'indexx_dp:NSTACK too small'
           STOP
          endif
          if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
          else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
          end if
         end if
        end do
        CONTAINS
        SUBROUTINE icomp_xchg(i,j)
        INTEGER(I4B), INTENT(INOUT) :: i,j
        INTEGER(I4B) :: swp
        if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
        end if
        END SUBROUTINE icomp_xchg
        END SUBROUTINE indexx_dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE indexx_qp(arr,index)
        USE nrutil, ONLY :arth,assert_eq,nrerror,swap
        IMPLICIT NONE
        REAL(QP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
        INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
!        Indexes an array arr, i.e., outputs the array index of length N such that arr(index(j))
!        is in ascending order for j = 1, 2, . . . ,N. The input quantity arr is not changed.
        REAL(QP) :: a
        INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        n=assert_eq(size(index),size(arr),'indexx_qp')
        index=arth(1,1,n)
        jstack=0
        l=1
        r=n
        do
         if (r-l < NN) then
          do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
            if (arr(index(i)) <= a) exit
            index(i+1)=index(i)
           end do
           index(i+1)=indext
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
         else
          k=(l+r)/2
          call swap(index(k),index(l+1))
          call icomp_xchg(index(l),index(r))
          call icomp_xchg(index(l+1),index(r))
          call icomp_xchg(index(l),index(l+1))
          i=l+1
          j=r
          indext=index(l+1)
          a=arr(indext)
          do
           do
            i=i+1
            if (arr(index(i)) >= a) exit
           end do
           do
            j=j-1
            if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
          end do
          index(l+1)=index(j)
          index(j)=indext
          jstack=jstack+2
          if (jstack > NSTACK) then
           write(6,*)'indexx_qp:NSTACK too small'
           STOP
          endif
          if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
          else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
          end if
         end if
        end do
        CONTAINS
        SUBROUTINE icomp_xchg(i,j)
        INTEGER(I4B), INTENT(INOUT) :: i,j
        INTEGER(I4B) :: swp
        if (arr(j) < arr(i)) then
        swp=i
        i=j
        j=swp
        end if
        END SUBROUTINE icomp_xchg
        END SUBROUTINE indexx_qp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE laguer_d(a,x,its)
        USE nrutil, ONLY :poly,poly_term
        INTEGER(I4B), INTENT(OUT) :: its
        COMPLEX(DPC), INTENT(INOUT) :: x
        COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
        REAL(DP), PARAMETER :: EPS=epsilon(1.0_dp)
        INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
!        Given an array of M + 1 complex coefficients a of the polynomial
!        and
!        given a complex value x, this routine improves x by Laguerre's method until it converges,
!        within the achievable roundoff limit, to a root of the given polynomial. The number of
!        iterations taken is returned as its.
!        Parameters: EPS is the estimated fractional roundoff error. We try to break (rare) limit
!        cycles with MR different fractional values, once every MT steps, for MAXIT total allowed
!        iterations.
        INTEGER(I4B) :: iter,m
        REAL(DP) :: abx,abp,abm,err
        COMPLEX(DPC) :: dx,x1,f,g,h,sq,gp,gm,g2
        COMPLEX(DPC), DIMENSION(size(a)) :: b,d
        REAL(DP), DIMENSION(MR) :: frac = &
        (/ 0.5_dp,0.25_dp,0.75_dp,0.13_dp,0.38_dp,0.62_dp,0.88_dp,1.0_dp /)
!        Fractions used to break a limit cycle.
        m=size(a)-1
        do iter=1,MAXIT !Loop over iterations up to allowed maximum.
         its=iter
         abx=abs(x)
         b(m+1:1:-1)=poly_term(a(m+1:1:-1),x) !Efficient computation of the polynomial
         d(m:1:-1)=poly_term(b(m+1:2:-1),x) ! and its first two derivatives.
         f=poly(x,d(2:m))
         err=EPS*poly(abx,abs(b(1:m+1))) !Esimate of roundoff in evaluating polynomial.
         if (abs(b(1)) <= err) RETURN !We are on the root.
         g=d(1)/b(1) !The generic case: Use Laguerre's formula.
         g2=g*g
         h=g2-2.0_dp*f/b(1)
         sq=sqrt((m-1)*(m*h-g2))
         gp=g+sq
         gm=g-sq
         abp=abs(gp)
         abm=abs(gm)
         if (abp < abm) gp=gm
         if (max(abp,abm) > 0.0) then
          dx=m/gp
         else
          dx=exp(cmplx(log(1.0_dp+abx),iter,kind=dpc))
         end if
         x1=x-dx
         if (x == x1) RETURN !Converged.
         if (mod(iter,MT) /= 0) then
          x=x1
         else !Every so often we take a fractional step, to break any limit cycle (itself a rare occurrence).
          x=x-dx*frac(iter/MT)
         end if
        end do
        write(6,*) 'Nbr d itérations dépassé dans laguer_d'
        STOP
!        Very unusual — can occur only for complex roots. Try a different starting guess for the root.
        END SUBROUTINE laguer_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE laguer_q(a,x,its)
        USE nrutil, ONLY :poly,poly_term
        INTEGER(I4B), INTENT(OUT) :: its
        COMPLEX(QPC), INTENT(INOUT) :: x
        COMPLEX(QPC), DIMENSION(:), INTENT(IN) :: a
        REAL(QP), PARAMETER :: EPS=epsilon(1.0_qp)
        INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
!        Given an array of M + 1 complex coefficients a of the polynomial
!        and
!        given a complex value x, this routine improves x by Laguerre's method until it converges,
!        within the achievable roundoff limit, to a root of the given polynomial. The number of
!        iterations taken is returned as its.
!        Parameters: EPS is the estimated fractional roundoff error. We try to break (rare) limit
!        cycles with MR different fractional values, once every MT steps, for MAXIT total allowed
!        iterations.
        INTEGER(I4B) :: iter,m
        REAL(QP) :: abx,abp,abm,err
        COMPLEX(QPC) :: dx,x1,f,g,h,sq,gp,gm,g2
        COMPLEX(QPC), DIMENSION(size(a)) :: b,d
        REAL(QP), DIMENSION(MR) :: frac = &
        (/ 0.5_qp,0.25_qp,0.75_qp,0.13_qp,0.38_qp,0.62_qp,0.88_qp,1.0_qp /)
!        Fractions used to break a limit cycle.
        m=size(a)-1
        do iter=1,MAXIT !Loop over iterations up to allowed maximum.
         its=iter
         abx=abs(x)
         b(m+1:1:-1)=poly_term(a(m+1:1:-1),x) !Efficient computation of the polynomial
         d(m:1:-1)=poly_term(b(m+1:2:-1),x) ! and its first two derivatives.
         f=poly(x,d(2:m))
         err=EPS*poly(abx,abs(b(1:m+1))) !Esimate of roundoff in evaluating polynomial.
         if (abs(b(1)) <= err) RETURN !We are on the root.
         g=d(1)/b(1) !The generic case: Use Laguerre's formula.
         g2=g*g
         h=g2-2.0_qp*f/b(1)
         sq=sqrt((m-1)*(m*h-g2))
         gp=g+sq
         gm=g-sq
         abp=abs(gp)
         abm=abs(gm)
         if (abp < abm) gp=gm
         if (max(abp,abm) > 0.0) then
          dx=m/gp
         else
          dx=exp(cmplx(log(1.0_qp+abx),iter,kind=qpc))
         end if
         x1=x-dx
         if (x == x1) RETURN !Converged.
         if (mod(iter,MT) /= 0) then
          x=x1
         else !Every so often we take a fractional step, to break any limit cycle (itself a rare occurrence).
          x=x-dx*frac(iter/MT)
         end if
        end do
        write(6,*) 'Nbr d itérations dépassé dans laguer_d'
        STOP
!        Very unusual — can occur only for complex roots. Try a different starting guess for the root.
        END SUBROUTINE laguer_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE zroots_d(a,roots,polish)
        USE nrutil, ONLY :assert_eq, poly_term
        COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(DPC), DIMENSION(:), INTENT(OUT) :: roots
        LOGICAL(LGT), INTENT(IN) :: polish
        REAL(DP), PARAMETER :: EPS=1.0e-6_dp
!        Given the array of M + 1 complex coefficients a of the polynomial
!        this
!        routine successively calls laguer and finds all M complex roots. The logical variable
!        polish should be input as .true. if polishing (also by Laguerre's method) is desired,
!       .false. if the roots will be subsequently polished by other means.
!        Parameter: EPS is a small number.
        INTEGER(I4B) :: j,its,m
        INTEGER(I4B), DIMENSION(size(roots)) :: indx
        COMPLEX(DPC) :: x
        COMPLEX(DPC), DIMENSION(size(a)) :: ad
        m=assert_eq(size(roots),size(a)-1,'zroots')
        ad(:)=a(:) !Copy of coefficients for successive deflation.
        do j=m,1,-1 !Loop over each root to be found.
         x=cmplx(0.0_dp,kind=dpc)
!        Start at zero to favor convergence to smallest remaining root.
         call laguer(ad(1:j+1),x,its) !Find the root.
         if (abs(aimag(x)) <= 2.0_dp*EPS**2*abs(real(x))) &
         x=cmplx(real(x),kind=dpc)
         roots(j)=x
         ad(j:1:-1)=poly_term(ad(j+1:2:-1),x) !Forward deflation.
        end do
        if (polish) then
         do j=1,m !Polish the roots using the undeflated coefficients.
          call laguer(a(:),roots(j),its)
         end do
        end if
        call indexx(real(roots),indx) !Sort roots by their real parts.
        roots=roots(indx)
        END SUBROUTINE zroots_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE zroots_q(a,roots,polish)
        USE nrutil, ONLY :assert_eq, poly_term
        COMPLEX(QPC), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(QPC), DIMENSION(:), INTENT(OUT) :: roots
        LOGICAL(LGT), INTENT(IN) :: polish
        REAL(QP), PARAMETER :: EPS=1.0e-14_qp
!        Given the array of M + 1 complex coefficients a of the polynomial
!        this
!        routine successively calls laguer and finds all M complex roots. The logical variable
!        polish should be input as .true. if polishing (also by Laguerre's method) is desired,
!       .false. if the roots will be subsequently polished by other means.
!        Parameter: EPS is a small number.
        INTEGER(I4B) :: j,its,m
        INTEGER(I4B), DIMENSION(size(roots)) :: indx
        COMPLEX(QPC) :: x
        COMPLEX(QPC), DIMENSION(size(a)) :: ad
        m=assert_eq(size(roots),size(a)-1,'zroots')
        ad(:)=a(:) !Copy of coefficients for successive deflation.
        do j=m,1,-1 !Loop over each root to be found.
         x=cmplx(0.0_qp,kind=qpc)
!        Start at zero to favor convergence to smallest remaining root.
         call laguer(ad(1:j+1),x,its) !Find the root.
         if (abs(aimag(x)) <= 2.0_qp*EPS**2*abs(real(x))) &
         x=cmplx(real(x),kind=qpc)
         roots(j)=x
         ad(j:1:-1)=poly_term(ad(j+1:2:-1),x) !Forward deflation.
        end do
        if (polish) then
         do j=1,m !Polish the roots using the undeflated coefficallcients.
          call laguer(a(:),roots(j),its)
         end do
        end if
        call indexx(aimag(roots),indx) !Sort roots by their real parts.
        roots=roots(indx)
        END SUBROUTINE zroots_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_s(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(SP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(SP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(SP), PARAMETER :: TINY=1.0e-20_sp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then 
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_sp/vv !Save the scaling.
        do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
        if (j /= imax) then !Do we need to interchange rows?
        call swap(a(imax,:),a(j,:)) !Yes, do so...
        d=-d !...and change the parity of d.
        vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lubksb_s(a,indx,b)
        USE nrutil, ONLY : assert_eq
        REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
!        Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!        as the original matrix A, but rather as its LU decomposition, determined by the routine
!        ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!        input as the right-hand-side vector B, also of length N, and returns with the solution vector
!        X. a and indx are not modified by this routine and can be left in place for successive calls
!        with different right-hand sides b. This routine takes into account the possibility that b will
!        begin with many zero elements, so it is efficient for use in matrix inversion.
        INTEGER(I4B) :: i,n,iii,ll
        REAL(SP) :: summ
        n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
        iii=0 !When iii is set to a positive value, it will become the index
!        of the first nonvanishing element of b. We now do
!        the forward substitution, equation (2.3.6). The only new
!        wrinkle is to unscramble the permutation as we go.
        do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (iii /= 0) then
        summ=summ-dot_product(a(i,iii:i-1),b(iii:i-1))
        else if (summ /= 0.0) then
        iii=i !A nonzero element was encountered, so from now on we will
        end if !have to do the dot product above.
        b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
        END SUBROUTINE lubksb_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_d(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(DP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(DP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(DP), PARAMETER :: TINY=1.0e-20_dp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then 
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_dp/vv !Save the scaling.
        do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
        if (j /= imax) then !Do we need to interchange rows?
        call swap(a(imax,:),a(j,:)) !Yes, do so...
        d=-d !...and change the parity of d.
        vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lubksb_d(a,indx,b)
        USE nrutil, ONLY : assert_eq
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
!        Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!        as the original matrix A, but rather as its LU decomposition, determined by the routine
!        ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!        input as the right-hand-side vector B, also of length N, and returns with the solution vector
!        X. a and indx are not modified by this routine and can be left in place for successive calls
!        with different right-hand sides b. This routine takes into account the possibility that b will
!        begin with many zero elements, so it is efficient for use in matrix inversion.
        INTEGER(I4B) :: i,n,iii,ll
        REAL(DP) :: summ
        n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
        iii=0 !When iii is set to a positive value, it will become the index
!        of the first nonvanishing element of b. We now do
!        the forward substitution, equation (2.3.6). The only new
!        wrinkle is to unscramble the permutation as we go.
        do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (iii /= 0) then
        summ=summ-dot_product(a(i,iii:i-1),b(iii:i-1))
        else if (summ /= 0.0) then
        iii=i !A nonzero element was encountered, so from now on we will
        end if !have to do the dot product above.
        b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
        END SUBROUTINE lubksb_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_q(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        REAL(QP), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(QP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(QP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(QP), PARAMETER :: TINY=1.0e-20_qp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_qp/vv !Save the scaling.
        do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
        if (j /= imax) then !Do we need to interchange rows?
        call swap(a(imax,:),a(j,:)) !Yes, do so...
        d=-d !...and change the parity of d.
        vv(imax)=vv(j) !Also interchange the scale factor.
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE lubksb_q(a,indx,b)
        USE nrutil, ONLY : assert_eq
        REAL(QP), DIMENSION(:,:), INTENT(IN) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: b
!        Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!        as the original matrix A, but rather as its LU decomposition, determined by the routine
!        ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!        input as the right-hand-side vector B, also of length N, and returns with the solution vector
!        X. a and indx are not modified by this routine and can be left in place for successive calls
!        with different right-hand sides b. This routine takes into account the possibility that b will
!        begin with many zero elements, so it is efficient for use in matrix inversion.
        INTEGER(I4B) :: i,n,iii,ll
        REAL(QP) :: summ
        n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
        iii=0 !When iii is set to a positive value, it will become the index
!        of the first nonvanishing element of b. We now do
!        the forward substitution, equation (2.3.6). The only new
!        wrinkle is to unscramble the permutation as we go.
        do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (iii /= 0) then
        summ=summ-dot_product(a(i,iii:i-1),b(iii:i-1))
        else if (summ /= 0.0) then
        iii=i !A nonzero element was encountered, so from now on we will
        end if !have to do the dot product above.
        b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
        END SUBROUTINE lubksb_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE ludcmp_cq(a,indx,d)
        USE nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
        COMPLEX(QPC), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
        REAL(QP), INTENT(OUT) :: d
!        Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!        rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!        output vector of length N that records the row permutation effected by the partial pivoting;
!        d is output as ±1 depending on whether the number of row interchanges was even or odd,
!        respectively. This routine is used in combination with lubksb to solve linear equations or
!        invert a matrix.
        REAL(QP), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(QP), PARAMETER :: TINY=1.0e-80_qp !A small number.
        INTEGER(I4B) :: j,n,imax
        n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
        d=1.0 !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling information.
        if (any(vv == 0.0))then
         write(*,*) 'singular matrix in ludcmp'
         STOP
        endif
!        There is a row of zeros.
        vv=1.0_qp/vv !Save the scaling.
        do j=1,n
         imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
         if (j /= imax) then !Do we need to interchange rows?
          call swap(a(imax,:),a(j,:)) !Yes, do so...
          d=-d !...and change the parity of d.
          vv(imax)=vv(j) !Also interchange the scale factor.
         end if
         indx(j)=imax
         if (a(j,j) == 0.0) a(j,j)=TINY
!        If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!        For some applications on singular matrices, it is desirable to substitute TINY
!        for zero.
         a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
         a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
!        Reduce remaining submatrix.
        end do
        END SUBROUTINE ludcmp_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_s(x,y,z)
        REAL(SP), INTENT(IN) :: x,y,z
        REAL(SP) :: rd_s
        REAL(SP), PARAMETER :: ERRTOL=0.0015_sp,TINY=1.0e-37_sp,BIG=4.5e37_sp,&
        C1=3.0_sp/14.0_sp,C2=1.0_sp/6.0_sp,C3=9.0_sp/22.0_sp,&
        C4=3.0_sp/26.0_sp,C5=0.25_sp*C3,C6=1.5_sp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        REAL(SP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
        'rd_s args')
        xt=x
        yt=y
        zt=z
        sum=0.0
        fac=1.0
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_sp*fac
        xt=0.25_sp*(xt+alamb)
        yt=0.25_sp*(yt+alamb)
        zt=0.25_sp*(zt+alamb)
        ave=0.2_sp*(xt+yt+3.0_sp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_sp*eb
        ee=ed+ec+ec
        rd_s=3.0_sp*sum+fac*(1.0_sp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_d(x,y,z)
        REAL(DP), INTENT(IN) :: x,y,z
        REAL(DP) :: rd_d
        REAL(DP), PARAMETER :: ERRTOL=0.0015_dp,TINY=1.0e-200_dp,BIG=4.5e200_dp,&
        C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
        C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        REAL(DP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
        'rd_d args')
        xt=x
        yt=y
        zt=z
        sum=0.0d0
        fac=1.0d0
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_dp*fac
        xt=0.25_dp*(xt+alamb)
        yt=0.25_dp*(yt+alamb)
        zt=0.25_dp*(zt+alamb)
        ave=0.2_dp*(xt+yt+3.0_dp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_dp*eb
        ee=ed+ec+ec
        rd_d=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_q(x,y,z)
        REAL(QP), INTENT(IN) :: x,y,z
        REAL(QP) :: rd_q
        REAL(QP), PARAMETER :: ERRTOL=0.0000025_qp,TINY=1.0e-200_qp,BIG=4.5e200_qp,&
        C1=3.0_qp/14.0_qp,C2=1.0_qp/6.0_qp,C3=9.0_qp/22.0_qp,&
        C4=3.0_qp/26.0_qp,C5=0.25_qp*C3,C6=1.5_qp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        REAL(QP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
        'rd_q args')
        xt=x
        yt=y
        zt=z
        sum=0.0d0
        fac=1.0d0
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_qp*fac
        xt=0.25_qp*(xt+alamb)
        yt=0.25_qp*(yt+alamb)
        zt=0.25_qp*(zt+alamb)
        ave=0.2_qp*(xt+yt+3.0_qp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_qp*eb
        ee=ed+ec+ec
        rd_q=3.0_qp*sum+fac*(1.0_qp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_c(x,y,z)
        COMPLEX*16, INTENT(IN) :: x,y,z
        COMPLEX*16 :: rd_c
        REAL(DP), PARAMETER :: ERRTOL=0.0015_dp,TINY=1.0e-200_dp,BIG=4.5e200_dp,&
        C1=3.0_dp/14.0_dp,C2=1.0_dp/6.0_dp,C3=9.0_dp/22.0_dp,&
        C4=3.0_dp/26.0_dp,C5=0.25_dp*C3,C6=1.5_dp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative (ligne de coupure ]-inf,0]), and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        COMPLEX*16 :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rd_c args')
        xt=x
        yt=y
        zt=z
        sum=dcmplx(0.0d0,0.d0)
        fac=dcmplx(1.0d0,0.d0)
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_dp*fac
        xt=0.25_dp*(xt+alamb)
        yt=0.25_dp*(yt+alamb)
        zt=0.25_dp*(zt+alamb)
        ave=0.2_dp*(xt+yt+3.0_dp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_dp*eb
        ee=ed+ec+ec
        rd_c=3.0_dp*sum+fac*(1.0_dp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rd_cq(x,y,z)
        COMPLEX(QPC), INTENT(IN) :: x,y,z
        COMPLEX(QPC) :: rd_cq
        REAL(QP), PARAMETER :: ERRTOL=0.00005_qp,TINY=1.0e-200_qp,BIG=4.5e200_qp,&
        C1=3.0_qp/14.0_qp,C2=1.0_qp/6.0_qp,C3=9.0_qp/22.0_qp,&
        C4=3.0_qp/26.0_qp,C5=0.25_qp*C3,C6=1.5_qp*C4
        !Computes Carlson's elliptic integral of the second kind, RD(x, y, z). x and y must be
        !nonnegative (ligne de coupure ]-inf,0]), and at most one can be zero. z must be positive. TINY must be at least twice
        !the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
        !times the negative 2/3 power of the machine underflow limit.
        COMPLEX(QPC) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
        ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rd_cq args')
        xt=x
        yt=y
        zt=z
        sum=cmplx(0.0,0.0,kind=QPC)
        fac=cmplx(1.0,0.0,kind=QPC)
        do
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25_qp*fac
        xt=0.25_qp*(xt+alamb)
        yt=0.25_qp*(yt+alamb)
        zt=0.25_qp*(zt+alamb)
        ave=0.2_qp*(xt+yt+3.0_qp*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        ea=delx*dely
        eb=delz*delz
        ec=ea-eb
        ed=ea-6.0_qp*eb
        ee=ed+ec+ec
        rd_cq=3.0_qp*sum+fac*(1.0_qp+ed*(-C1+C5*ed-C6*delz*ee)&
        +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
        END FUNCTION rd_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_s(x,y,z)
        REAL(SP), INTENT(IN) :: x,y,z
        REAL(SP) :: rf_s
        REAL(SP), PARAMETER :: ERRTOL=0.08_sp,TINY=1.5e-38_sp,BIG=3.0e37_sp,&
        THIRD=1.0_sp/3.0_sp,&
        C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        REAL(SP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z) >= TINY, &
        max(x,y,z) <= BIG, 'rf_s args')
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_sp*(xt+alamb)
         yt=0.25_sp*(yt+alamb)
         zt=0.25_sp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_s=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_d(x,y,z)
        REAL(DP), INTENT(IN) :: x,y,z
        REAL(DP) :: rf_d
        REAL(DP), PARAMETER :: ERRTOL=0.0025_dp,TINY=1.5e-100_dp,BIG=3.0e100_dp,&
        THIRD=1.0_dp/3.0_dp,&
        C1=1.0_dp/24.0_dp,C2=0.1_dp,C3=3.0_dp/44.0_dp,C4=1.0_dp/14.0_dp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        REAL(DP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        if(.NOT.((min(x,y,z) >= 0.d0).AND.(min(x+y,x+z,y+z)>=TINY).AND.(max(x,y,z)<=BIG)))then
          write(6,*) 'erreur ds les arguments de rf_d'
          write(6,*) 'x,y,z=',x,y,z
          STOP
        endif
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_dp*(xt+alamb)
         yt=0.25_dp*(yt+alamb)
         zt=0.25_dp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_d=(1.0_dp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_q(x,y,z)
        REAL(QP), INTENT(IN) :: x,y,z
        REAL(QP) :: rf_q
        REAL(QP), PARAMETER :: ERRTOL=0.0000025_qp,TINY=1.5e-100_qp,BIG=3.0e100_qp,&
        THIRD=1.0_qp/3.0_qp,&
        C1=1.0_qp/24.0_qp,C2=0.1_qp,C3=3.0_qp/44.0_qp,C4=1.0_qp/14.0_qp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        REAL(QP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        if(.NOT.((min(x,y,z) >= 0.0_qp).AND.(min(x+y,x+z,y+z)>=TINY).AND.(max(x,y,z)<=BIG)))then
          write(6,*) 'erreur ds les arguments de rf_q'
          write(6,*) 'x,y,z=',x,y,z
          STOP
        endif
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         zt=0.25_qp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_q=(1.0_qp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_c(x,y,z)
        COMPLEX*16, INTENT(IN) :: x,y,z
        COMPLEX*16 :: rf_c
        REAL(DP), PARAMETER :: ERRTOL=0.0025_dp,TINY=1.5e-100_dp,BIG=3.0e100_dp,&
        THIRD=1.0_dp/3.0_dp,&
        C1=1.0_dp/24.0_dp,C2=0.1_dp,C3=3.0_dp/44.0_dp,C4=1.0_dp/14.0_dp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        COMPLEX*16 :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(x)+abs(z),abs(y)+abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rf_c args')
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_dp*(xt+alamb)
         yt=0.25_dp*(yt+alamb)
         zt=0.25_dp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_c=(1.0_dp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rf_cq(x,y,z)
        COMPLEX(QPC), INTENT(IN) :: x,y,z
        COMPLEX(QPC) :: rf_cq
        REAL(QP), PARAMETER :: ERRTOL=0.000025_qp,TINY=1.5e-100_qp,BIG=3.0e100_qp,&
        THIRD=1.0_qp/3.0_qp,&
        C1=1.0_qp/24.0_qp,C2=0.1_qp,C3=3.0_qp/44.0_qp,C4=1.0_qp/14.0_qp
        !Computes Carlson's elliptic integral of the first kind, RF (x, y, z). x, y, and z must be
        !nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
        !underflow limit, BIG at most one-fifth the machine overflow limit.
        COMPLEX(QPC) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
        call assert((/.NOT.((real(x) <= 0.0).AND.(abs(imag(x))<TINY)), &
                      .NOT.((real(y) <= 0.0).AND.(abs(imag(y))<TINY)), &
                      .NOT.((real(z) <= 0.0).AND.(abs(imag(z))<TINY)), &
                      min(abs(x)+abs(y),abs(x)+abs(z),abs(y)+abs(z)) >= TINY, &
                      max(abs(x),abs(y),abs(z)) <= BIG/),'rf_cq args')
        xt=x
        yt=y
        zt=z
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         zt=0.25_qp*(zt+alamb)
         ave=THIRD*(xt+yt+zt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
        end do
        e2=delx*dely-delz**2
        e3=delx*dely*delz
        rf_cq=(1.0_qp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
        END FUNCTION rf_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rj_q(x,y,z,p)
        USE nrutil, ONLY : assert
        REAL(QP), INTENT(IN) :: x,y,z,p
        REAL(QP) :: rj_q
        REAL(QP), PARAMETER :: ERRTOL=0.000005_qp,TINY=2.5e-233_qp,BIG=9.0e241_qp,&
        C1=3.0_qp/14.0_qp,C2=1.0_qp/3.0_qp,C3=3.0_qp/22.0_qp,&
        C4=3.0_qp/26.0_qp,C5=0.75_qp*C3,C6=1.5_qp*C4,C7=0.5_qp*C2,&
        C8=C3+C3
        !Computes Carlson's elliptic integral of the third kind, RJ(x, y, z, p). x, y, and z must be
        !nonnegative, and at most one can be zero. p must be nonzero. If p < 0, the Cauchy
        !principal value is returned. TINY must be at least twice the cube root of the machine
        !underflow limit, BIG at most one-fifth the cube root of the machine overflow limit.
        REAL(QP) :: a,alamb,alpha,ave,b,bet,delp,delx,&
        dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
        sm,tau,xt,yt,zt
        call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z,abs(p)) >= TINY, &
         max(x,y,z,abs(p)) <= BIG, "rj_q args")
        sm=0.0
        fac=1.0
        if (p > 0.0) then
         xt=x
         yt=y
         zt=z
         pt=p
        else
         xt=min(x,y,z)
         zt=max(x,y,z)
         yt=x+y+z-xt-zt
         a=1.0_qp/(yt-p)
         b=a*(zt-yt)*(yt-xt)
         pt=yt+b
         rho=xt*zt/yt
         tau=p*pt/yt
        end if
        do
         sqrtx=sqrt(xt)
         sqrty=sqrt(yt)
         sqrtz=sqrt(zt)
         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
         alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
         bet=pt*(pt+alamb)**2
         sm=sm+fac*rc(alpha,bet)
         fac=0.25_qp*fac
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         zt=0.25_qp*(zt+alamb)
         pt=0.25_qp*(pt+alamb)
         ave=0.2_qp*(xt+yt+zt+pt+pt)
         delx=(ave-xt)/ave
         dely=(ave-yt)/ave
         delz=(ave-zt)/ave
         delp=(ave-pt)/ave
         if (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL) exit
        end do
        ea=delx*(dely+delz)+dely*delz
        eb=delx*dely*delz
        ec=delp**2
        ed=ea-3.0_qp*ec
        ee=eb+2.0_qp*delp*(ea-ec)
        rj_q=3.0_qp*sm+fac*(1.0_qp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
        +delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
        if (p <= 0.0) rj_q=a*(b*rj_q+3.0_qp*(rc(rho,tau)-rf(xt,yt,zt)))
        END FUNCTION rj_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rc_q(x,y)
        USE nrutil, ONLY : assert
        REAL(QP), INTENT(IN) :: x,y
        REAL(QP) :: rc_q
        REAL(QP), PARAMETER :: ERRTOL=0.0004_qp,TINY=1.69e-68_qp,&
        SQRTNY=1.3e-34_qp,BIG=3.0e37_qp,TNBG=TINY*BIG,&
        COMP1=2.236_qp/SQRTNY,COMP2=TNBG*TNBG/25.0_qp,&
        THIRD=1.0_qp/3.0_qp,&
        C1=0.3_qp,C2=1.0_qp/7.0_qp,C3=0.375_qp,C4=9.0_qp/22.0_qp
!        Computes Carlson's degenerate elliptic integral, RC(x, y). x must be nonnegative and y
!        must be nonzero. If y < 0, the Cauchy principal value is returned. TINY must be at least
!        5 times the machine underflow limit, BIG at most one-fifth the machine maximum overflow
!        limit.
        REAL(QP) :: alamb,ave,s,w,xt,yt
        call assert( (/x >= 0.0,y /= 0.0,x+abs(y) >= TINY,x+abs(y) <= BIG, &
         y >= -COMP1 .or. x <= 0.0 .or. x >= COMP2/),"rc_q")
        if (y > 0.0) then
         xt=x
         yt=y
         w=1.0
        else
         xt=x-y
         yt=-y
         w=sqrt(x)/sqrt(xt)
        end if
        do
         alamb=2.0_qp*sqrt(xt)*sqrt(yt)+yt
         xt=0.25_qp*(xt+alamb)
         yt=0.25_qp*(yt+alamb)
         ave=THIRD*(xt+yt+yt)
         s=(yt-ave)/ave
         if (abs(s) <= ERRTOL) exit
        end do
        rc_q=w*(1.0_qp+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
        END FUNCTION rc_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellemath_s(phi,m)
        REAL(SP), INTENT(IN) :: phi,m
        REAL(SP) :: ellemath_s
!        Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
!        and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
!        int_0^phi dt*sqrt(1-m*sin(t)**2)
!        if m*sin(phi)^2>1 returns the real part of the intergral
        REAL(SP) :: cc,q,s
        s=sin(phi)
        if(m*s**2>1) s=1/sqrt(abs(m)) !The integrande becomes pure imaginary beyond t=arcsin(1/sqrt(abs(m)))
        cc=1-s**2
        q=1-m*s**2
        ellemath_s=s*(rf(cc,q,1.0_sp)-m*s**2*rd(cc,q,1.0_sp)/3.0_sp)
        END FUNCTION ellemath_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellemath_d(phi,m)
        REAL(DP), INTENT(IN) :: phi,m
        REAL(DP) :: ellemath_d
        REAL(DP) :: cc,q,s
        s=sin(phi)
        if(m*s**2>1) s=1/sqrt(abs(m))
        cc=1-s**2
        q=1-m*s**2
        ellemath_d=s*(rf(cc,q,1.0_dp)-m*s**2*rd(cc,q,1.0_dp)/3.0_dp)
        END FUNCTION ellemath_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellemath_q(phi,m)
        REAL(QP), INTENT(IN) :: phi,m
        REAL(QP) :: ellemath_q
        REAL(QP) :: cc,q,s
        s=sin(phi)
        if(m*s**2>1) s=1/sqrt(abs(m))
        cc=1-s**2
        q=1-m*s**2
        ellemath_q=s*(rf(cc,q,1.0_qp)-m*s**2*rd(cc,q,1.0_qp)/3.0_qp)
        END FUNCTION ellemath_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_s(phi,ak)
        REAL(SP), INTENT(IN) :: phi,ak
        REAL(SP) :: elle_s
!        Legendre elliptic integral of the 2nd kind E(φ, k), evaluated using Carlson's functions RD
!        and RF . The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
!        int_0^phi dt*sqrt(1-ak**2*sin(t)**2)
!        if ak^2*sin(phi)^2>1 returns the real part of the intergral
        REAL(SP) :: cc,q,s
        s=sin(phi)
        if(ak**2*s**2>1) s=1/abs(ak) !The integrande becomes pure imaginary beyond t=arcsin(1/abs(ak))
        cc=1-s**2
        q=(1-s*ak)*(1+s*ak)
        elle_s=s*(rf(cc,q,1.0_sp)-((s*ak)**2)*rd(cc,q,1.0_sp)/3.0_sp)
        END FUNCTION elle_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_d(phi,ak)
        REAL(DP), INTENT(IN) :: phi,ak
        REAL(DP) :: elle_d
        REAL(DP) :: cc,q,s
        s=sin(phi)
        if(ak**2*s**2>1) s=1/abs(ak)
        cc=1-s**2
        q=(1-s*ak)*(1+s*ak)
        elle_d=s*(rf(cc,q,1.0_dp)-((s*ak)**2)*rd(cc,q,1.0_dp)/3.0_dp)
        END FUNCTION elle_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_q(phi,ak)
        REAL(QP), INTENT(IN) :: phi,ak
        REAL(QP) :: elle_q
        REAL(QP) :: cc,q,s
        s=sin(phi)
        if(ak**2*s**2>1) s=1/abs(ak)
        cc=1-s**2
        q=(1-s*ak)*(1+s*ak)
        elle_q=s*(rf(cc,q,1.0_qp)-((s*ak)**2)*rd(cc,q,1.0_qp)/3.0_qp)
        END FUNCTION elle_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elleim_q(phi,ak,complx)
        REAL(QP), INTENT(IN) :: phi,ak
        LOGICAL, INTENT(IN) :: complx !valeur sans importance, sert pour l'INTERFACE elle
        COMPLEX(QPC) :: elleim_q
        REAL(QP) :: phi1,kp,im,re
!       Integrale elliptique du 1er type au voisinage de sa ligne de coupure sur l'axe réel
!       quand ak^2*sin^2(phi)>1, l'intégrale est réelle de 0 à phi1=asin(1/abs(ak)) (partie réelle calculée par elle(phi,ak))
!       puis imaginaire pure de phi1 à phi. La partie imaginaire est calculée ci-dessous
        re=elle(phi,ak)
        if(abs(ak)>1)then 
         phi1=asin(1/abs(ak))
         kp=sqrt(ak**2/(ak**2-1))
         im=-(elle(PI/2.0_qp-phi1,kp)-elle(PI/2.0_qp-phi,kp))*sqrt(ak**2-1)
        else
         im=0.0_qp
        endif
        elleim_q=cmplx(re,im,kind=qpc)
        END FUNCTION elleim_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_cq(phi,ak)
! Calcule l'intégrale elliptique de première espèce prolongée naivement au plan complexe
! ak complexe, phi réel
        USE modsim
        IMPLICIT NONE
        REAL(QP), INTENT(IN) :: phi
        COMPLEX(QPC), INTENT(IN) :: ak
        COMPLEX(QPC) :: elle_cq
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref
        INTEGER it,nt
        
        elle_cq=qromocq(intee,0.0_qp,phi,(/bidon/),midpntcq,1.0e-18_qp)
        
        CONTAINS
         FUNCTION intee(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: intee
        
         intee=sqrt(1.0_qp-ak**2.0_qp*sin(x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION elle_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elle_ccq(phi,ak)
! Calcule l'intégrale elliptique de première espèce prolongée naivement au plan complexe
! ak complexe, phi complexe
        USE modsim
        IMPLICIT NONE
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: elle_ccq
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref
        INTEGER it,nt
        
        elle_ccq=qromocq(intee,0.0_qp,1.0_qp,(/bidon/),midpntcq,1.0e-18_qp)
        elle_ccq=phi*elle_ccq
        
        CONTAINS
         FUNCTION intee(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: intee
        
         intee=sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION elle_ccq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION elleadiab(phi,ak)
! Calcule l'intégrale elliptique de seconde espèce prolongée au plan complexe
! Déplace les lignes de coupure par relèvement adiabatique de la racine de l'intégrande
! Des points de brancherment demeurent là où l'argument de la racine carrée s'annulle
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: elleadiab
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref
        INTEGER it,nt
        
        elleadiab=cmplx(0.0_qp,0.0_qp,kind=qpc)
        nt=54000
        dt=1.0_qp/real(nt)
        int1=1.0_qp
        do it=0,nt
! Règle de Simpson
         if((it==0).OR.(it==nt))then
          pref=1.0_qp
         elseif(modulo(it,2)==1)then
          pref=4.0_qp
         elseif(modulo(it,2)==0)then
          pref=2.0_qp
         endif
         t=it*dt
         int2=int1
         int1=intee(t,(/bidon/))
         int1=int1*(-1.0_qp)**(minloc(abs((/int1+int2,int1-int2/)),DIM=1))
!Suivi adiabatique de l'intégrande. Échoue aux points de branchement où intef,intee=0
         elleadiab=elleadiab+pref*int1
        enddo
        elleadiab=phi*elleadiab*dt/3.0_qp
        
        CONTAINS
         FUNCTION intee(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC) :: intee
        
         intee=sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellfmath_s(phi,m)
        REAL(SP), INTENT(IN) :: phi,m
        REAL(SP) :: ellfmath_s
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
!        int_0^phi dt/sqrt(1-m*sin(t)**2)
        REAL(SP) :: s,cc
        s=sin(phi)
        if(m*s**2>1) s=1/sqrt(abs(m))
        cc=1-s**2
        ellfmath_s=s*rf(cc,1.0_sp-m*s**2,1.0_sp)
        END FUNCTION ellfmath_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellfmath_d(phi,m)
        REAL(DP), INTENT(IN) :: phi,m
        REAL(DP) :: ellfmath_d
        REAL(DP) :: s,cc
        s=sin(phi)
        if(m*s**2>1) s=1/sqrt(abs(m))
        cc=1-s**2
        ellfmath_d=s*rf(cc,1.0_dp-m*s**2,1.0_dp)
        END FUNCTION ellfmath_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellfmath_q(phi,m)
        REAL(QP), INTENT(IN) :: phi,m
        REAL(QP) :: ellfmath_q
        REAL(QP) :: s,cc
        s=sin(phi)
        cc=1-s**2
        ellfmath_q=s*rf(cc,1.0_qp-m*s**2,1.0_qp)
        END FUNCTION ellfmath_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_s(phi,ak)
        REAL(SP), INTENT(IN) :: phi,ak
        REAL(SP) :: ellf_s
!        Legendre elliptic integral of the 1st kind F(φ, k), evaluated using Carlson's function RF .
!        The argument ranges are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
!        int_0^phi dt/sqrt(1-ak**2*sin(t)**2)
!        if ak^2*sin(phi)^2>1 returns the real part of the intergral
        REAL(SP) :: s,cc
        s=sin(phi)
        if(ak**2*s**2>1) s=1/abs(ak)  !The integrande becomes pure imaginary beyond t=arcsin(1/abs(ak))
        cc=1-s**2
        ellf_s=s*rf(cc,(1.0_sp-s*ak)*(1.0_sp+s*ak),1.0_sp)
        END FUNCTION ellf_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_d(phi,ak)
        REAL(DP), INTENT(IN) :: phi,ak
        REAL(DP) :: ellf_d
        REAL(DP) :: s,cc
        s=sin(phi)
        if(ak**2*s**2>1) s=1/abs(ak)
        cc=1-s**2
        ellf_d=s*rf(cc,(1.0_dp-s*ak)*(1.0_dp+s*ak),1.0_dp)
        END FUNCTION ellf_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_q(phi,ak)
        REAL(QP), INTENT(IN) :: phi,ak
        REAL(QP) :: ellf_q
        REAL(QP) :: s,cc
        s=sin(phi)
        if(ak**2*s**2>1) s=1/abs(ak)
        cc=1-s**2
        ellf_q=s*rf(cc,(1.0_qp-s*ak)*(1.0_qp+s*ak),1.0_qp)
        END FUNCTION ellf_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellfim_q(phi,ak,complx)
!       Integrale elliptique du 2eme type au voisinage de sa ligne de coupure sur l'axe réel
!       quand ak^2*sin^2(phi)>1, l'intégrale est réelle de 0 à phi1=asin(1/abs(ak)) (partie réelle calculée par ellf(phi,ak))
!       puis imaginaire pure de phi1 à phi. La partie imaginaire est calculée ci-dessous
        REAL(QP), INTENT(IN) :: phi,ak
        LOGICAL, INTENT(IN) :: complx!valeur sans importance, sert pour l'INTERFACE ellf
        COMPLEX(QPC) :: ellfim_q
        REAL(QP) :: phi1,kp,im,re
        re=ellf(phi,ak)
        if(abs(ak)>1)then 
         phi1=asin(1/abs(ak))
         kp=sqrt(ak**2/(ak**2-1))
         im=(ellf(PI/2.0_qp-phi1,kp)-ellf(PI/2.0_qp-phi,kp))/sqrt(ak**2-1)
        else
         im=0.0_qp
        endif
        ellfim_q=cmplx(re,im,kind=qpc)
        END FUNCTION ellfim_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_cq(phi,ak)
! Calcule l'intégrale elliptique de seconde espèce prolongée naivement au plan complexe
! ak complexe, phi réel
        USE modsim
        IMPLICIT NONE
        REAL(QP), INTENT(IN) :: phi
        COMPLEX(QPC), INTENT(IN) :: ak
        COMPLEX(QPC) :: ellf_cq
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref
        INTEGER it,nt
        
        ellf_cq=qromocq(intef,0.0_qp,phi,(/bidon/),midpntcq,1.0e-18_qp)
        
        CONTAINS
         FUNCTION intef(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: intef
        
         intef=sqrt(1.0_qp-ak**2.0_qp*sin(x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION ellf_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellf_ccq(phi,ak)
! Calcule l'intégrale elliptique de seconde espèce prolongée naivement au plan complexe
! ak complexe, phi complexe
        USE modsim
        IMPLICIT NONE
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: ellf_ccq
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref
        INTEGER it,nt
        
        ellf_ccq=qromocq(intef,0.0_qp,1.0_qp,(/bidon/),midpntcq,1.0e-18_qp)
        ellf_ccq=phi*ellf_ccq
        
        CONTAINS
         FUNCTION intef(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), DIMENSION(:), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(size(x)) :: intef
        
         intef=sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION ellf_ccq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellfadiab(phi,ak)
! Calcule l'intégrale elliptique de seconde espèce prolongée au plan complexe
! Déplace les lignes de coupure par relèvement adiabatique de la racine de l'intégrande
! Des points de brancherment demeurent là où l'argument de la racine carrée s'annulle
        COMPLEX(QPC), INTENT(IN) :: phi,ak
        COMPLEX(QPC) :: ellfadiab
        COMPLEX(QPC) :: int1,int2
        REAL(QP) t,dt,pref
        INTEGER it,nt
        
        ellfadiab=cmplx(0.0_qp,0.0_qp,kind=qpc)
        nt=54000
        dt=1.0_qp/real(nt)
        int1=1.0_qp
        do it=0,nt
! Règle de Simpson
         if((it==0).OR.(it==nt))then
          pref=1.0_qp
         elseif(modulo(it,2)==1)then
          pref=4.0_qp
         elseif(modulo(it,2)==0)then
          pref=2.0_qp
         endif
         t=it*dt
         int2=int1
         int1=intef(t,(/bidon/))
         int1=int1*(-1.0_qp)**(minloc(abs((/int1+int2,int1-int2/)),DIM=1))
!Suivi adiabatique de l'intégrande. Échoue aux points de branchement où intef,intee=0
         ellfadiab=ellfadiab+pref*int1
        enddo
        ellfadiab=phi*ellfadiab*dt/3.0_qp
        
        CONTAINS
         FUNCTION intef(x,arg)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: x
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC) :: intef
        
         intef=sqrt(1.0_qp-ak**2.0_qp*sin(phi*x)**2.0_qp)
        
         END FUNCTION
        END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellpi_q(phi,en,ak)
        REAL(QP), INTENT(IN) :: phi,en,ak
        REAL(QP) :: ellpi_q
!        Legendre elliptic integral of the 3rd kind Π(φ, n, k), evaluated using Carlson's functions RJ
!        and RF.
!        The ranges of φ and k are 0 ≤ φ ≤ π/2, 0 ≤ k sin φ ≤ 1.
!        definition mathematica int_0^phi dt/(1-en*sin^2 t)/sqrt(1-ak*sin^2 t)
        REAL(QP) :: cc,enss,q,s
        s=sin(phi)
        enss=-en*s*s
        cc=cos(phi)**2
        q=(1.0_qp-s**2*ak)
        ellpi_q=s*(rf(cc,q,1.0_qp)-enss*rj(cc,q,1.0_qp,1.0_qp+enss)/3.0_qp)
        END FUNCTION ellpi_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION ellpi_cq(phi,en,ak)
        USE modsim
        REAL(QP), INTENT(IN) :: phi
        COMPLEX(QPC), INTENT(IN) :: en,ak
        COMPLEX(QPC) :: ellpi_cq

        ellpi_cq=qromo(fonc,0.0_qp,phi,(/bidon/),midpntcq,1.0e-12_qp)
        CONTAINS
         FUNCTION fonc(x,arg)
         REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
         COMPLEX(QPC), DIMENSION(size(x)) :: fonc
         fonc=1/(1-en*sin(x)**2)/sqrt(1-ak*sin(x)**2)
         END FUNCTION fonc
        END FUNCTION ellpi_cq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rtnewt(funcd,arg,x1,x2,xacc)
        REAL(DP), INTENT(IN) :: x1,x2,xacc
        REAL(DP), DIMENSION(:), INTENT(IN) :: arg
        REAL(DP) :: rtnewt
        INTERFACE
         SUBROUTINE funcd(x,arg,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), INTENT(OUT) :: fval,fderiv
         END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=25
!        Using the Newton-Raphson method, find the root of a function known to lie in the interval
!        [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd
!        is a user-supplied subroutine that returns both the function value and the first derivative of
!        the function.
!        Parameter: MAXIT is the maximum number of iterations.
        INTEGER(I4B) :: j
        REAL(DP) :: df,dx,f
        rtnewt=0.5_dp*(x1+x2) !Premier essai
        do j=1,MAXIT
        call funcd(rtnewt,arg,f,df)
        dx=f/df
        rtnewt=rtnewt-dx
        if ((x1-rtnewt)*(rtnewt-x2) < 0.0)then
         write (6,*) 'nrerror: rtnewt valeur hors de lintervalle'
         STOP
        endif
        if (abs(dx) < xacc) RETURN !Convergence.
        end do
        write (6,*) 'nrerreur: nombre ditérations dépassé dans rtnewt '
        STOP 
        END FUNCTION rtnewt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rtsafed(funcd,arg,x1,x2,xacc)
        REAL(DP), INTENT(IN) :: x1,x2,xacc
        REAL(DP), DIMENSION(:), INTENT(IN) :: arg
        REAL(DP) :: rtsafed
        INTERFACE
         SUBROUTINE funcd(x,arg,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(DP), INTENT(OUT) :: fval,fderiv
         END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=1800
!       Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
!       between x1 and x2. The root, returned as the function value rtsafed, will be refined until
!       its accuracy is known within ±xacc. funcd is a user-supplied subroutine that returns both
!       the function value and the first derivative of the function.
!       Parameter: MAXIT is the maximum allowed number of iterations.
        INTEGER(I4B) :: j
        REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
        call funcd(x1,arg,fl,df)
        call funcd(x2,arg,fh,df)
        if ((fl > 0.0 .and. fh > 0.0) .or. &
            (fl < 0.0 .and. fh < 0.0)) then
           write(6,*)'nerreur: racine pas encadrée dans rtsafed'
           stop
        endif
        if (fl == 0.0) then
           rtsafed=x1
           RETURN
        else if (fh == 0.0) then
           rtsafed=x2
           RETURN
        else if (fl < 0.0) then !Orient the search so that f(xl) < 0.
           xl=x1
           xh=x2
        else
           xh=x1
           xl=x2
        end if
        rtsafed=0.5_dp*(x1+x2) !Initialize the guess for root,
        dxold=abs(x2-x1)      !the “stepsize before last,”
        dx=dxold              !and the last step.
        call funcd(rtsafed,arg,f,df)
        do j=1,MAXIT !Loop over allowed iterations.
           if (((rtsafed-xh)*df-f)*((rtsafed-xl)*df-f) > 0.0 .or. &
              abs(2.0_dp*f) > abs(dxold*df) ) then !Bisect if Newton out of range, or not decreasing fast enough.
              dxold=dx
              dx=0.5_dp*(xh-xl)
              rtsafed=xl+dx
              if (xl == rtsafed) RETURN !Change in root is negligible.
           else !Newton step acceptable. Take it.
           dxold=dx
           dx=f/df
           temp=rtsafed
           rtsafed=rtsafed-dx
           if (temp == rtsafed) RETURN
           end if
           if (abs(dx) < xacc) RETURN !Convergence criterion.
           call funcd(rtsafed,arg,f,df) !One new function evaluation per iteration.
           !write(6,*)'x,f=',rtsafed,f
           if (f < 0.0) then !Maintain the bracket on the root.
              xl=rtsafed
           else
              xh=rtsafed
           end if
        end do
        write (6,*) 'nrerreur: nombre ditérations dépassé dans rtsafed '
        STOP
        END FUNCTION rtsafed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION rtsafeq(funcd,arg,x1,x2,xacc)
        REAL(QP), INTENT(IN) :: x1,x2,xacc
        REAL(QP), DIMENSION(:), INTENT(IN) :: arg
        REAL(QP) :: rtsafeq
        INTERFACE
         SUBROUTINE funcd(x,arg,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(QP), INTENT(IN) :: x !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg !x: variable sur laquelle effectuer la recherche de racine. arg: tous les autres variables muettes
         REAL(QP), INTENT(OUT) :: fval,fderiv
         END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=1800
!       Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
!       between x1 and x2. The root, returned as the function value rtsafe, will be refined until
!       its accuracy is known within ±xacc. funcd is a user-supplied subroutine that returns both
!       the function value and the first derivative of the function.
!       Parameter: MAXIT is the maximum allowed number of iterations.
        INTEGER(I4B) :: j
        REAL(QP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
        call funcd(x1,arg,fl,df)
        call funcd(x2,arg,fh,df)
        if ((fl > 0.0 .and. fh > 0.0) .or. &
            (fl < 0.0 .and. fh < 0.0)) then
           write(6,*)'nerreur: racine pas encadrée dans rtsafeq'
        endif
        if (fl == 0.0) then
           rtsafeq=x1
           RETURN
        else if (fh == 0.0) then
           rtsafeq=x2
           RETURN
        else if (fl < 0.0) then !Orient the search so that f(xl) < 0.
           xl=x1
           xh=x2
        else
           xh=x1
           xl=x2
        end if
        rtsafeq=0.5_qp*(x1+x2) !Initialize the guess for root,
        dxold=abs(x2-x1)      !the “stepsize before last,”
        dx=dxold              !and the last step.
        call funcd(rtsafeq,arg,f,df)
        do j=1,MAXIT !Loop over allowed iterations.
           if (((rtsafeq-xh)*df-f)*((rtsafeq-xl)*df-f) > 0.0 .or. &
              abs(2.0_qp*f) > abs(dxold*df) ) then !Bisect if Newton out of range, or not decreasing fast enough.
              dxold=dx
              dx=0.5_qp*(xh-xl)
              rtsafeq=xl+dx
              if (xl == rtsafeq) RETURN !Change in root is negligible.
           else !Newton step acceptable. Take it.
           dxold=dx
           dx=f/df
           temp=rtsafeq
           rtsafeq=rtsafeq-dx
           if (temp == rtsafeq) RETURN
           end if
           if (abs(dx) < xacc) RETURN !Convergence criterion.
           call funcd(rtsafeq,arg,f,df) !One new function evaluation per iteration.
           !write(6,*)'x,f=',rtsafeq,f
           if (f < 0.0) then !Maintain the bracket on the root.
              xl=rtsafeq
           else
              xh=rtsafeq
           end if
        end do
        write (6,*) 'nrerreur: nombre ditérations dépassé dans rtsafeq '
        STOP
        END FUNCTION rtsafeq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_s(arr)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER(I4B) :: i,j,n
        REAL(SP) :: a
        n=size(arr)
        do j=2,n 
          a=arr(j)
          do i=j-1,1,-1
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a
        end do
        END SUBROUTINE tri_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_d(arr)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER(I4B) :: i,j,n
        REAL(DP) :: a
        n=size(arr)
        do j=2,n 
          a=arr(j)
          do i=j-1,1,-1 
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a
        end do
        END SUBROUTINE tri_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_q(arr)
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: arr
!        Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced
!        on output by its sorted rearrangement.
        INTEGER(I4B) :: i,j,n
        REAL(QP) :: a
        n=size(arr)
        do j=2,n !Pick out each element in turn.
          a=arr(j)
          do i=j-1,1,-1 !Look for the place to insert it.
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
          end do
          arr(i+1)=a !Insert it.
        end do
        END SUBROUTINE tri_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tricol_q(arr,col,dimm)
        REAL(QP), DIMENSION(:,:), INTENT(INOUT) :: arr
        INTEGER,  INTENT(IN) :: col,dimm
        INTEGER(I4B) :: i,j,n
        REAL(QP), ALLOCATABLE, DIMENSION(:) :: a
        INTEGER m
        m=size(arr,dim=dimm)
        allocate(a(1:size(arr,dim=3-dimm)))
        do j=2,m !Pick out each element in turn.
         if(dimm==1)then
            a=arr(j,:)
            do i=j-1,1,-1 !Look for the place to insert it.
              if (arr(i,col) <= a(col)) exit
              arr(i+1,:)=arr(i,:)
            end do
            arr(i+1,:)=a!Insert it.
            do i=1,m
            enddo
         else
            a=arr(:,j)
            do i=j-1,1,-1 
              if (arr(col,i) <= a(col)) exit
              arr(:,i+1)=arr(:,i)
            end do
            arr(:,i+1)=a
         endif
        end do
        deallocate(a)
        END SUBROUTINE tricol_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE tri_pos_q(arr,pos)
        REAL(QP), DIMENSION(:), INTENT(INOUT) :: arr
        INTEGER,  DIMENSION(size(arr)), INTENT(OUT) :: pos
!        Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced
!        on output by its sorted rearrangement.
        INTEGER(I4B) :: i,j,n
        REAL(QP) :: a
        n=size(arr)
        pos(:)=1
        do j=2,n !Pick out each element in turn.
          a=arr(j)
          do i=j-1,1,-1 !Look for the place to insert it.
            if (arr(i) <= a) exit
            arr(i+1)=arr(i)
            pos(i+1)=pos(i)
          end do
          arr(i+1)=a !Insert it.
          pos(i+1)=j
        end do
        END SUBROUTINE tri_pos_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION argum_d(z)
        COMPLEX(DPC), INTENT(IN) :: z
        REAL(DP) argum_d
        argum_d=atan2(imag(z),real(z))
        END FUNCTION argum_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION argum_q(z)
        COMPLEX(QPC), INTENT(IN) :: z
        REAL(QP) argum_q
        argum_q=atan2(imag(z),real(z))
        END FUNCTION argum_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE frenel(x,s,c)
        USE nrutil, ONLY : nrerror
        REAL(QP), INTENT(IN) :: x
        REAL(QP), INTENT(OUT) :: s,c
        INTEGER(I4B), PARAMETER :: MAXIT=1000
        REAL(QP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x),BIG=huge(x)*EPS,XMIN=1.5_qp
        !Computes the Fresnel integrals S(x) and C(x) for all real x.
        !Parameters: MAXIT is the maximum number of iterations allowed; EPS is the relative error;
        !FPMIN is a number near the smallest representable floating-point number; BIG is a number
        !near the machine overflow limit; XMIN is the dividing line between using the series and
        !continued fraction.
        INTEGER(I4B) :: k,n
        REAL(QP) :: a,ax,fact,pix2,sign,sum,sumc,sums,term,test
        COMPLEX(QPC) :: b,cc,d,h,del,cs
        LOGICAL(LGT) :: odd
        ax=abs(x)
        if (ax < sqrt(FPMIN)) then !Special case: avoid failure of convergence test because of underflow.
         s= 0.0 
         c=ax
        else if (ax <= XMIN) then !Evaluate both series simultaneously.
         sum=0.0
         sums=0.0
         sumc=ax
         sign=1.0_qp
         fact=PIs2*ax*ax
         odd=.true.
         term=ax
         n=3
         do k=1,MAXIT
            term=term*fact/k
            sum=sum+sign*term/n
            test=abs(sum)*EPS
            if (odd) then
               sign=-sign
               sums=sum
               sum=sumc
            else
               sumc=sum
               sum=sums
            end if
            if (term < test) exit
            odd=.not. odd
            n=n+2
         end do
         if (k > MAXIT) call nrerror('frenel: series failed')
         s=sums
         c=sumc
        else !Evaluate continued fraction by modified Lentz's method (§5.2).
         pix2=PI*ax*ax 
         b=cmplx(1.0_qp,-pix2,kind=qpc)
         cc=BIG
         d=1.0_qp/b
         h=d
         n=-1
         do k=2,MAXIT
          n=n+2
          a=-n*(n+1)
          b=b+4.0_qp
          d=1.0_qp/(a*d+b) !Denominators cannot be zero.
          cc=b+a/cc
          del=cc*d
          h=h*del
          if (absc(del-1.0_qp) <= EPS) exit
         end do
         if (k > MAXIT) call nrerror('cf failed in frenel')
         h=h*cmplx(ax,-ax,kind=qpc)
         cs=cmplx(0.5_qp,0.5_qp,kind=qpc)*(1.0_qp-cmplx(cos(0.5_qp*pix2),sin(0.5_qp*pix2),kind=qpc)*h)
         c=real(cs)
         s=aimag(cs)
        end if
        if (x < 0.0) then !Use antisymmetry.
         c=-c
         s=-s
        end if
        CONTAINS
         FUNCTION absc(z)
         IMPLICIT NONE
         COMPLEX(QPC), INTENT(IN) :: z
         REAL(QP) :: absc
         absc=abs(real(z))+abs(aimag(z))
         END FUNCTION absc
        END SUBROUTINE frenel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        FUNCTION nlignes(fich)
        CHARACTER(len=*), INTENT(IN) :: fich
        INTEGER nlignes

        INTEGER io,ilignes

        nlignes=0
        open(1001,file=trim(fich),ACTION="READ",iostat=io)
          if(io==29)then
           nlignes=0
           return
          endif
          do ilignes=1,100000000
           read(1001,*,end=11)
           nlignes=nlignes+1
          enddo
        11 CONTINUE
        close(1001)
        END FUNCTION nlignes
END MODULE recettes
