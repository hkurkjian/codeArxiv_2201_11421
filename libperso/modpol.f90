MODULE modpol
 USE nrtype
 USE OMP_LIB
 IMPLICIT NONE
       INTERFACE iminloc
        MODULE PROCEDURE iminloc_s,iminloc_d,iminloc_q
       END INTERFACE iminloc
       INTERFACE polint 
        MODULE PROCEDURE polints,polintd,polintq,polintc,polintcq
       END INTERFACE polint
       INTERFACE polin2 
        MODULE PROCEDURE polin2s,polin2q
       END INTERFACE polin2
       INTERFACE locate
        MODULE PROCEDURE locates,located,locateq,locatei
       END INTERFACE locate
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polints(xa,ya,x,y,dy)
!       Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
!       and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) =
!       yai, i = 1, . . . ,N, then the returned value y = P(x).
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: y,dy
       INTEGER(I4B) :: m,n,ns
       REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
       INCLUDE "polint.f90"
       END SUBROUTINE polints
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polin2s(x1a,x2a,ya,x1,x2,y,dy)
       USE nrtype; USE nrutil, ONLY : assert_eq
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: x1a
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: x2a
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
       REAL(SP), INTENT(IN) :: x1,x2
       REAL(SP), INTENT(OUT) :: y,dy
!       Given arrays x1a of length M and x2a of length N of independent variables, and an M×N
!       array of function values ya, tabulated at the grid points defined by x1a and x2a, and given
!       values x1 and x2 of the independent variables, this routine returns an interpolated function
!       value y, and an accuracy indication dy (based only on the interpolation in the x1 direction,
!       however).
       INTEGER(I4B) :: j,m,ndum
       REAL(SP), DIMENSION(size(x1a))   :: ymtmp
       REAL(SP), DIMENSION(size(x2a,1)) :: yntmp,x2atmp
       m   =assert_eq(size(x1a),size(ya,2),    "polin2: m")
       ndum=assert_eq(size(x2a,1),size(ya,1),  "polin2: ndum")
       do j=1,m !Loop over rows.
        yntmp=ya(j,:) !Copy row into temporary storage.
        x2atmp=x2a(j,:)
        call polint(x2atmp,yntmp,x2,ymtmp(j),dy) !Interpolate answer into temporary storage.
       enddo
       call polint(x1a,ymtmp,x1,y,dy) !Do the final interpolation.
       END SUBROUTINE polin2s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE locates(xx,x,pos)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), INTENT(IN) :: x
       INTEGER(I4B), INTENT(OUT) :: pos
       INCLUDE "locate.f90"
       END SUBROUTINE locates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION iminloc_s(arr)
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(1) :: imin
       INTEGER(I4B) :: iminloc_s
       imin=minloc(arr(:))
       iminloc_s=imin(1)
       END FUNCTION iminloc_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintd(xa,ya,x,y,dy)
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: y,dy
       INTEGER(I4B) :: m,n,ns
       REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho
       INCLUDE "polint.f90"
       END SUBROUTINE polintd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE located(xx,x,pos)
       USE nrtype
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: xx
       REAL(DP), INTENT(IN) :: x
       INTEGER(I4B), INTENT(OUT) :: pos
       INCLUDE "locate.f90"
       END SUBROUTINE located
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION iminloc_d(arr)
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(1) :: imin
        INTEGER(I4B) :: iminloc_d
        imin=minloc(arr(:))
        iminloc_d=imin(1)
       END FUNCTION iminloc_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintq(xa,ya,x,y,dy)
       REAL(QP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(QP), INTENT(IN) :: x
       REAL(QP), INTENT(OUT) :: y,dy
       INTEGER(I4B) :: m,n,ns
       REAL(QP), DIMENSION(size(xa)) :: c,d,den,ho
       INCLUDE "polint.f90"
       END SUBROUTINE polintq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polin2q(x1a,x2a,ya,x1,x2,y,dy)
       USE nrtype; USE nrutil, ONLY : assert_eq
       IMPLICIT NONE
       REAL(QP), DIMENSION(:), INTENT(IN) :: x1a
       REAL(QP), DIMENSION(:,:), INTENT(IN) :: x2a
       REAL(QP), DIMENSION(:,:), INTENT(IN) :: ya
       REAL(QP), INTENT(IN) :: x1,x2
       REAL(QP), INTENT(OUT) :: y,dy
!       Given arrays x1a of length M and x2a of length N of independent variables, and an M×N
!       array of function values ya, tabulated at the grid points defined by x1a and x2a, and given
!       values x1 and x2 of the independent variables, this routine returns an interpolated function
!       value y, and an accuracy indication dy (based only on the interpolation in the x1 direction,
!       however).
       INTEGER(I4B) :: j,m,ndum
       REAL(QP), DIMENSION(size(x1a))   :: ymtmp
       REAL(QP), DIMENSION(size(x2a,1)) :: yntmp,x2atmp
       m   =assert_eq(size(x1a),size(ya,1),size(x2a,1),    "polin2: m")
       ndum=assert_eq(size(x2a,2),size(ya,2),              "polin2: ndum")
       do j=1,m !Loop over rows.
        yntmp=ya(:,j) !Copy row into temporary storage.
        x2atmp=x2a(:,j)
        call polint(x2atmp,yntmp,x2,ymtmp(j),dy) !Interpolate answer into temporary storage.
!        write(6,*)"j,ymtmp(j)=",j,ymtmp(j)
       enddo
       call polint(x1a,ymtmp,x1,y,dy) !Do the final interpolation.
       END SUBROUTINE polin2q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE locateq(xx,x,pos)
       USE nrtype
       IMPLICIT NONE
       REAL(QP), DIMENSION(:), INTENT(IN) :: xx
       REAL(QP), INTENT(IN) :: x
       INTEGER(I4B), INTENT(OUT) :: pos
       INCLUDE "locate.f90"
       END SUBROUTINE locateq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION iminloc_q(arr)
       REAL(QP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4B), DIMENSION(1) :: imin
       INTEGER(I4B) :: iminloc_q
       imin=minloc(arr(:))
       iminloc_q=imin(1)
       END FUNCTION iminloc_q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintc(xa,ya,x,y,dy)
       COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: ya
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa
       REAL(DP), INTENT(IN) :: x
       COMPLEX(DPC), INTENT(OUT) :: y,dy
       INTEGER(I4B) :: m,n,ns
       COMPLEX(DPC), DIMENSION(size(xa)) :: c,d,den
       REAL(DP), DIMENSION(size(xa)) :: ho
       INCLUDE "polint.f90"
       END SUBROUTINE polintc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE polintcq(xa,ya,x,y,dy)
       COMPLEX(QPC), DIMENSION(:), INTENT(IN) :: ya
       REAL(QP), DIMENSION(:), INTENT(IN) :: xa
       REAL(QP), INTENT(IN) :: x
       COMPLEX(QPC), INTENT(OUT) :: y,dy
       INTEGER(I4B) :: m,n,ns
       COMPLEX(QPC), DIMENSION(size(xa)) :: c,d,den
       REAL(QP), DIMENSION(size(xa)) :: ho
       INCLUDE "polint.f90"
       END SUBROUTINE polintcq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE locatei(xx,x,pos)
       USE nrtype
       IMPLICIT NONE
       INTEGER, DIMENSION(:), INTENT(IN) :: xx
       INTEGER, INTENT(IN) :: x
       INTEGER(I4B), INTENT(OUT) :: pos
       INCLUDE "locate.f90"
       END SUBROUTINE locatei
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROUTINE phi1 (r,r0,v)
!       
!!! PHI1 evaluates the multiquadric radial basis function.
!!  Parameters:
!!
!!    Input,  real R(N), the radial separation. 0 < R.
!!    Input,  real R0, a scale factor. Doit ếtre compris entre la distance min et max entre les points à interpoler
!!    Output, real V(N), the value of the radial basis function.
!!
!!  Reference:
!!
!!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!!    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
!!    Third Edition,
!!    Cambridge University Press, 2007,
!!    ISBN13: 978-0-521-88068-8,
!!    LC: QA297.N866.
!
!       REAL(QP), DIMENSION(:), INTENT(IN) ::  r
!       REAL(QP), INTENT(INOUT) ::  r0
!       REAL(QP), DIMENSION(size(r)), INTENT(OUT) :: v 
!       
!       v=sqrt(r**2+r0**2)
!       
!       END SUBROUTINE phi1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROUTINE phi2 (r,r0,v)
!       
!!! PHI2 evaluates the inverse multiquadric radial basis function.
!
!       REAL(QP), DIMENSION(:), INTENT(IN) ::  r
!       REAL(QP), INTENT(INOUT) ::  r0
!       REAL(QP), DIMENSION(size(r)), INTENT(OUT) :: v 
!       
!       v=1.0_qp/sqrt(r**2 + r0**2)
!       
!       END SUBROUTINE phi2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROUTINE phi3 (r,r0,v)
!       
!!! PHI3 evaluates the thin-plate spline radial basis function.
!!
!!  Discussion:
!!
!!    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
!!    it may be desirable to choose a value of R0 smaller than any possible R.
!       
!       REAL(QP), DIMENSION(:), INTENT(IN) ::  r
!       REAL(QP), INTENT(INOUT) ::  r0
!       REAL(QP), DIMENSION(size(r)), INTENT(OUT) :: v
!       
!       INTEGER it
!
!       do it = 1, size(r)
!        if(r(it).le.0.0_qp)then
!         v(it) = 0.0_qp
!        else
!         v(it) = r(it)**2*log(r(it)/r0)
!        end if
!       end do
!       
!       END SUBROUTINE phi3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROUTINE phi4 (r,r0,v)
!       
!!! PHI4 evaluates the gaussian radial basis function.
!       REAL(QP), DIMENSION(:), INTENT(IN) ::  r
!       REAL(QP), INTENT(INOUT) ::  r0
!       REAL(QP), DIMENSION(size(r)), INTENT(OUT) :: v
!       
!       v =exp(-0.5_qp*r**2/r0**2)
!       
!       END SUBROUTINE phi4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROUTINE rbf_weight(xd,r0,phi,fd,w)
!       USE recettes
!
!!! RBF_WEIGHT computes weights for radial basis function interpolation.
!!
!!  Discussion:
!!
!!    We assume that there are N (nonsingular) equations in N unknowns.
!!
!!    However, it should be clear that, if we are willing to do some kind
!!    of least squares calculation, we could allow for singularity,
!!    inconsistency, or underdetermine systems.  This could be associated
!!    with data points that are very close or repeated, a smaller number
!!    of data points than function values, or some other ill-conditioning
!!    of the system arising from a peculiarity in the point spacing.
!!
!!    Input, real ( kind = 8 ) XD(M,ND), the data points.
!!
!!    Input, real ( kind = 8 ) R0, a scale factor.  R0 should be larger than 
!!    the typical separation between points, but smaller than the maximum 
!!    separation.  The value of R0 has a significant effect on the resulting 
!!    interpolant.
!!
!!    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
!!    basis functions.
!!
!!    Input, real ( kind = 8 ) FD(ND), the function values at the data points.
!!
!!    Output, real ( kind = 8 ) W(ND), the weights.
!!
!       REAL(QP), INTENT(INOUT) ::  r0
!       INTERFACE 
!        SUBROUTINE phi(r,rr0,v)
!        USE nrtype
!        REAL(QP), DIMENSION(:), INTENT(IN) ::  r
!        REAL(QP), INTENT(INOUT) ::  rr0
!        REAL(QP), DIMENSION(size(r)), INTENT(OUT) :: v
!        END SUBROUTINE phi
!       END INTERFACE
!       REAL(QP), DIMENSION(:,:), INTENT(IN) :: xd
!       REAL(QP), DIMENSION(size(xd,dim=2)), INTENT(IN) :: fd
!       REAL(QP), DIMENSION(size(fd)), INTENT(OUT) :: w
!
!       REAL(QP), DIMENSION(size(fd),size(fd)) :: a
!       INTEGER i,j,nd
!       INTEGER(I4B), DIMENSION(size(fd)) :: indx
!       REAL(QP), DIMENSION(size(fd)) :: r,v
!       REAL(QP) d
!       
!       nd=size(fd)
!       do i = 1, nd
!       
!         do j = 1, nd
!           r(j) = sqrt(sum((xd(:,i) - xd(:,j))**2))
!         end do
!       
!         call phi(r,r0,v)
!       
!         a(i,1:nd) = v(1:nd)
!       
!       end do
!       
!       w=fd
!       call ludcmp(a,indx,d)
!       call lubksb(a,indx,w)
!       
!       END SUBROUTINE rbf_weight
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       SUBROUTINE rbf_interp(xd,r0,phi,w,xi,fi)
!       
!!! RBF_INTERP_ND evaluates a radial basis function interpolant.
!!
!!  Parameters:
!!
!!    Input, real ( kind = 8 ) XD(M,ND), the data points.
!!
!!    Input, real ( kind = 8 ) R0, a scale factor.  R0 should be larger than 
!!    the typical separation between points, but smaller than the maximum 
!!    separation.  The value of R0 has a significant effect on the resulting 
!!    interpolant.
!!
!!    Input, subroutine PHI ( N, R, R0, V ), a subroutine to evaluate the radial
!!    basis functions.
!!
!!    Input, real ( kind = 8 ) W(ND), the weights, as computed by RBF_WEIGHTS.
!!
!!    Input, real ( kind = 8 ) XI(M,NI), the interpolation points.
!!
!!    Output, real ( kind = 8 ) FI(NI), the interpolated values.
!!
!       REAL(QP), INTENT(INOUT) ::  r0
!       INTERFACE 
!        SUBROUTINE phi(r,rr0,v)
!        USE nrtype
!        REAL(QP), DIMENSION(:), INTENT(IN) ::  r
!        REAL(QP), INTENT(INOUT) ::  rr0
!        REAL(QP), DIMENSION(size(r)), INTENT(OUT) :: v
!        END SUBROUTINE phi
!       END INTERFACE
!       REAL(QP), DIMENSION(:,:), INTENT(IN) :: xd
!       REAL(QP), DIMENSION(size(xd,dim=2)), INTENT(IN) :: w
!       REAL(QP), DIMENSION(:,:), INTENT(IN) :: xi
!       REAL(QP), DIMENSION(size(xi,dim=2)), INTENT(OUT) :: fi
!       
!       INTEGER i,j
!       REAL(QP), DIMENSION(size(xd,dim=2)) :: r,v
!
!       do i=1,size(xi,dim=2) 
!       
!           do j=1,size(xd,dim=2)
!             r(j)=sqrt(sum((xi(:,i)-xd(:,j))**2))
!           end do
!       
!           call phi(r,r0,v)
!       
!           fi(i)=dot_product(v,w)
!       
!         end do
!       
!       END SUBROUTINE rbf_interp
END MODULE modpol
