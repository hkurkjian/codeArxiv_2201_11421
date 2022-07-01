MODULE modsim
 USE nrtype; USE nrutil, ONLY :arth
 USE modpol
 USE OMP_LIB
 IMPLICIT NONE
 INTEGER, PARAMETER :: mpnt=1,minf=2,msql=3,msqu=4,mexp=5,rinf=6
       INTERFACE midpnt
        MODULE PROCEDURE midpnts,midpntd,midpntq,midpntc,midpntcq,midpntvcq
       END INTERFACE midpnt
       INTERFACE midinf
        MODULE PROCEDURE midinfs,midinfd,midinfq,midinfc,midinfcq
       END INTERFACE midinf
       INTERFACE midsqu
        MODULE PROCEDURE midsqus,midsqud,midsquq,midsquc,midsqucq
       END INTERFACE midsqu
       INTERFACE midsql
        MODULE PROCEDURE midsqls,midsqld,midsqlq,midsqlc,midsqlcq
       END INTERFACE midsql
       INTERFACE midexp
        MODULE PROCEDURE midexps,midexpd,midexpq,midexpc,midexpcq
       END INTERFACE midexp
       INTERFACE racinf
        MODULE PROCEDURE racinfs,racinfd,racinfq,racinfc,racinfcq,racinfvcq
       END INTERFACE racinf
       INTERFACE qromo
        MODULE PROCEDURE qromos,qromod,qromoq,qromoc,qromocq,qromovq,qromovcq
       END INTERFACE qromo
       INTERFACE qromovfixed
        MODULE PROCEDURE qromovqfixed,qromovcqfixed
       END INTERFACE qromovfixed
       INTERFACE qromochoix
        MODULE PROCEDURE qromochoixs,qromochoixd,qromochoixq,qromochoixc,qromochoixcq,qromochoixvq,qromochoixvcq
       END INTERFACE qromochoix
       INTERFACE decoupe
        MODULE PROCEDURE decoupes,decouped,decoupeq,decoupec,decoupecq,decoupevq,decoupevcq
       END INTERFACE decoupe
       ABSTRACT INTERFACE
         FUNCTION funcs(x,arg)
          USE nrtype
          REAL(SP), DIMENSION(:), INTENT(IN) :: x
          REAL(SP), DIMENSION(:), INTENT(IN) :: arg
          REAL(SP), DIMENSION(size(x)) :: funcs
         END FUNCTION funcs
         FUNCTION funcd(x,arg)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          REAL(DP), DIMENSION(:), INTENT(IN) :: arg
          REAL(DP), DIMENSION(size(x)) :: funcd
         END FUNCTION funcd
         FUNCTION funcq(x,arg)
          USE nrtype
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          REAL(QP), DIMENSION(size(x)) :: funcq
         END FUNCTION funcq
         FUNCTION funcc(x,arg)
          USE nrtype
          REAL(DP), DIMENSION(:), INTENT(IN) :: x
          REAL(DP), DIMENSION(:), INTENT(IN) :: arg
          COMPLEX(DPC), DIMENSION(size(x)) :: funcc
         END FUNCTION funcc
         FUNCTION funccq(x,arg)
          USE nrtype
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          COMPLEX(QPC), DIMENSION(size(x)) :: funccq
         END FUNCTION funccq
         FUNCTION funcvq(x,arg,mm)
          USE nrtype
          INTEGER, INTENT(IN) ::  mm
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          REAL(QP), DIMENSION(size(x),mm) :: funcvq
         END FUNCTION funcvq
         FUNCTION funcvcq(x,arg,mm)
          USE nrtype
          INTEGER, INTENT(IN) ::  mm
          REAL(QP), DIMENSION(:), INTENT(IN) :: x
          REAL(QP), DIMENSION(:), INTENT(IN) :: arg
          COMPLEX(QPC), DIMENSION(size(x),mm) :: funcvcq
         END FUNCTION funcvcq
       END INTERFACE
 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpnts(func,a,b,arg,s,n)
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(IN), DIMENSION(:) :: arg
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funcs) :: func
!       This routine computes the nth stage of refinement of an extended midpoint rule. func is
!       input as the name of the function to be integrated between limits a and b, also input. When
!       called with n=1, the routine returns as s the crudest estimate of Subsequent
!       calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding
!       (2/3)×3n-1 additional interior points. s should not be modified between sequential calls.
       REAL(SP) :: del
       INTEGER(I4B) :: it
       REAL(SP), DIMENSION(2*3**(n-2)) :: x
       INCLUDE "somme.f90"
       END SUBROUTINE midpnts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfs(funk,aa,bb,arg,s,n)
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/x rather than in x. This allows the upper limit bb to be as large and positive
!       as the computer allows, or the lower limit aa to be as large and negative, but not both.
!       aa and bb must have the same sign.
       INCLUDE "decls.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_sp/aa !These two statements change the limits of integration accordingly
       a= 1.0_sp/bb
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=funk(1.0_sp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqus(funk,aa,bb,arg,s,n)
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in sqrt(bb-x) rather than in x. This allows to handle squareroot divergences in the upper bound
       INCLUDE "decls.f90"
       b=sqrt(bb-aa)!These two statements change the limits of integration accordingly
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)!This internal function effects the change of variable.
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqls(funk,aa,bb,arg,s,n)
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in sqrt(x-aa) rather than in x. This allows to handle squareroot divergences in the lower bound
       INCLUDE "decls.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=2.0_sp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexps(funk,aa,bb,arg,s,n)
       !This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
       !of the integral of funk from aa to bb, except that bb is assumed to be infinite (value passed
       !not actually used). It is assumed that the function funk decreases exponentially rapidly at
       !infinity.
       INCLUDE "decls.f90"
       b=exp(-aa) !These two statements change the limits of integration accordingly
       a=0.0 
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) !This internal function effects the change of variable.
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=funk(-log(x),arg)/x
        END FUNCTION func
       END SUBROUTINE midexps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfs(funk,aa,bb,arg,s,n)
!       This routine is an exact replacement for midpnt, i.e., returns as s the nth stage of refinement
!       of the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
!       points in 1/sqrt(x) rather than in x. This allows to handle 1/x**(3/2) tails at large x
       INCLUDE "decls.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_sp/sqrt(aa)
       a= 1.0_sp/sqrt(bb)
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(SP), DIMENSION(size(x)) :: func
        func=2.0_sp*funk(1.0_sp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromos(func,a,b,arg,choose,EPS)
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: arg
       REAL(SP), INTENT(IN) :: EPS
       REAL(SP) :: qromos
       PROCEDURE(funcs) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         IMPORT
         REAL(SP), INTENT(IN) :: aa,bb
         REAL(SP), DIMENSION(:), INTENT(IN) :: arg
         REAL(SP), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcs) :: funk
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
!         Romberg integration on an open interval. Returns the integral of the function func from a
!         to b, using any specified integrating subroutine choose and Romberg’s method. Normally
!         choose will be an open formula, not evaluating the function at the endpoints. It is assumed
!         that choose triples the number of steps on each call, and that its error series contains only
!         even powers of the number of steps. The routines midpnt, midinf, midsql, midsqu,
!         and midexp are possible choices for choose. The parameters have the same meaning as in qromb.
       REAL(SP), DIMENSION(JMAXP) :: h,s
       REAL(SP) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromos,dqromo)
         if ((abs(dqromo) <= EPS*abs(qromos)).AND.(abs(s(j)-s(j-1))<0.01_sp*abs(s(j)))) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_sp !This is where the assumption of step tripling and an even error series is used.
       end do
       STOP 'Nombre d itération dépassé dans qromos'
       END FUNCTION qromos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixs(func,a,b,arg,choix,EPS)
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: choix
       REAL(SP), INTENT(IN)  :: EPS
       REAL(SP) :: qromochoixs
       PROCEDURE(funcs) :: func
       if(choix==mpnt)then 
        qromochoixs=qromo(func,a,b,arg,midpnts,EPS)
       elseif(choix==minf)then
        qromochoixs=qromo(func,a,b,arg,midinfs,EPS)
       elseif(choix==msql)then
        qromochoixs=qromo(func,a,b,arg,midsqls,EPS)
       elseif(choix==msqu)then
        qromochoixs=qromo(func,a,b,arg,midsqus,EPS)
       elseif(choix==mexp)then
        qromochoixs=qromo(func,a,b,arg,midexps,EPS)
       elseif(choix==rinf)then
        qromochoixs=qromo(func,a,b,arg,racinfs,EPS)
       endif
       END FUNCTION qromochoixs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decoupes(func,bornes,arg,choix,EPS,bla) 
       !Calcule la somme des integrales de func de bornes(m) à bornes(m+1)
       !avec le changement de variable choix(m) et les arguments arg(m)
       REAL(SP), DIMENSION(:), INTENT(IN) :: bornes
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(SP), INTENT(IN) :: EPS
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funcs) :: func
       REAL(SP) decoupes
       !
       INTEGER m
       REAL(SP) I
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decoupes"
       if(size(choix)      .NE.size(bornes)-1) stop "size(choix) ne correspond pas à size(bornes)-1 dans decoupes"
       decoupes=0.0_sp
       do m=1,size(bornes)-1
        I=qromochoixs(func,bornes(m),bornes(m+1),arg(:,m),choix(m),EPS)
        decoupes=decoupes+I
        if(bla) write(6,FMT="(A1,I1,A3,1P,G20.10)") "I",m,"=  ",I
       enddo
       if(bla) write(6,FMT="(A5,1P,G20.10)") "Itot=",decoupes
       END FUNCTION decoupes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntd(func,a,b,arg,s,n)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funcd) :: func
       !
       REAL(DP) :: del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       INCLUDE "somme.f90"
       END SUBROUTINE midpntd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfd(funk,aa,bb,arg,s,n)
       INCLUDE "decld.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_dp/aa 
       a= 1.0_dp/bb
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=funk(1.0_dp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqud(funk,aa,bb,arg,s,n)
       INCLUDE "decld.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqud
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqld(funk,aa,bb,arg,s,n)
       INCLUDE "decld.f90"
       b=sqrt(bb-aa) 
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqld
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexpd(funk,aa,bb,arg,s,n)
       INCLUDE "decld.f90"
       b=exp(-aa) 
       a=0.0 
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=funk(-log(x),arg)/x
        END FUNCTION func
       END SUBROUTINE midexpd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfd(funk,aa,bb,arg,s,n)
       INCLUDE "decld.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_dp/sqrt(aa)
       a= 1.0_dp/sqrt(bb)
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(DP), DIMENSION(size(x)) :: func
        func=2.0_dp*funk(1.0_dp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromod(func,a,b,arg,choose,EPS)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), DIMENSION(:), INTENT(IN) :: arg
       REAL(DP), INTENT(IN) :: EPS
       REAL(DP) :: qromod
       PROCEDURE(funcd) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         IMPORT
         REAL(DP), INTENT(IN) :: aa,bb
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         REAL(DP), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcd) :: funk
         END SUBROUTINE choose
       END INTERFACE
       !
       INTEGER(I4B), PARAMETER :: JMAX=15,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(DP), DIMENSION(JMAXP) :: h,s
       REAL(DP) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromod,dqromo)
         if ((abs(dqromo) <= EPS*abs(qromod)).AND.(abs(s(j)-s(j-1))<0.01_dp*abs(s(j)))) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_dp
       end do
       write(6,*) 'Nombre d itération dépassé dans qromod'
       END FUNCTION qromod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixd(func,a,b,arg,choix,EPS)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: choix
       REAL(DP), INTENT(IN)  :: EPS
       REAL(DP) :: qromochoixd
       PROCEDURE(funcd) :: func
       if(choix==mpnt)then 
        qromochoixd=qromo(func,a,b,arg,midpntd,EPS)
       elseif(choix==minf)then
        qromochoixd=qromo(func,a,b,arg,midinfd,EPS)
       elseif(choix==msql)then
        qromochoixd=qromo(func,a,b,arg,midsqld,EPS)
       elseif(choix==msqu)then
        qromochoixd=qromo(func,a,b,arg,midsqud,EPS)
       elseif(choix==mexp)then
        qromochoixd=qromo(func,a,b,arg,midexpd,EPS)
       elseif(choix==rinf)then
        qromochoixd=qromo(func,a,b,arg,racinfd,EPS)
       endif
       END FUNCTION qromochoixd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decouped(func,bornes,arg,choix,EPS,bla) 
       !Calcule la somme des integrales de func de bornes(m) à bornes(m+1)
       !avec le changement de variable choix(m) et les arguments arg(m)
       REAL(DP), DIMENSION(:), INTENT(IN) :: bornes
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(DP), INTENT(IN) :: EPS
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funcd) :: func
       REAL(DP) decouped
       !
       INTEGER m
       REAL(DP) I
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decouped"
       if(size(choix)      .NE.size(bornes)-1) stop "size(choix) ne correspond pas à size(bornes)-1 dans decouped"
       decouped=0.0_dp
       do m=1,size(bornes)-1
        I=qromochoixd(func,bornes(m),bornes(m+1),arg(:,m),choix(m),EPS)
        decouped=decouped+I
        if(bla) write(6,FMT="(A1,I1,A3,1P,G25.18)") "I",m,"=  ",I
       enddo
       if(bla) write(6,FMT="(A5,1P,G25.18)") "Itot=",decouped
       END FUNCTION decouped

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntq(func,a,b,arg,s,n)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funcq) :: func
       !
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       INCLUDE "somme.f90"
       END SUBROUTINE midpntq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfq(funk,aa,bb,arg,s,n)
       INCLUDE "declq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_qp/aa
       a= 1.0_qp/bb
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=funk(1.0_qp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsquq(funk,aa,bb,arg,s,n)
       INCLUDE "declq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsquq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlq(funk,aa,bb,arg,s,n)
       INCLUDE "declq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqlq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexpq(funk,aa,bb,arg,s,n)
       INCLUDE "declq.f90"
       b=exp(-aa) 
       a=0.0 
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=funk(-log(x),arg)/x
        END FUNCTION func
       END SUBROUTINE midexpq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfq(funk,aa,bb,arg,s,n)
       INCLUDE "declq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa)
       a= 1.0_qp/sqrt(bb)
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        REAL(QP), DIMENSION(size(x)) :: func
        func=2.0_qp*funk(1.0_qp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromoq(func,a,b,arg,choose,EPS)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       REAL(QP), INTENT(IN) :: EPS
       REAL(QP) :: qromoq
       PROCEDURE(funcq) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         IMPORT
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcq) :: funk
         END SUBROUTINE choose
       END INTERFACE
       !
       INTEGER(I4B), PARAMETER :: JMAX=16,JMAXP=JMAX+1,K=5,KM=K-1
       REAL(QP), DIMENSION(JMAXP) :: h,s
       REAL(QP) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
!        write(6,*)"s(j)=",s(j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qromoq,dqromo)
!         write(6,*)"dqromo=",dqromo
         if ((abs(dqromo) <= EPS*abs(qromoq)).AND.(abs(s(j)-s(j-1))<0.01_qp*abs(s(j)))) RETURN
         if ((abs(qromoq) <= EPS*10e-10_qp)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_qp
       end do
       write(6,*) 'Nombre d itération dépassé dans qromoq'
       END FUNCTION qromoq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixq(func,a,b,arg,choix,EPS)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: choix
       REAL(QP), INTENT(IN)  :: EPS
       REAL(QP) :: qromochoixq
       PROCEDURE(funcq) :: func
       if(choix==mpnt)then 
        qromochoixq=qromo(func,a,b,arg,midpntq,EPS)
       elseif(choix==minf)then
        qromochoixq=qromo(func,a,b,arg,midinfq,EPS)
       elseif(choix==msql)then
        qromochoixq=qromo(func,a,b,arg,midsqlq,EPS)
       elseif(choix==msqu)then
        qromochoixq=qromo(func,a,b,arg,midsquq,EPS)
       elseif(choix==mexp)then
        qromochoixq=qromo(func,a,b,arg,midexpq,EPS)
       elseif(choix==rinf)then
        qromochoixq=qromo(func,a,b,arg,racinfq,EPS)
       endif
       END FUNCTION qromochoixq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decoupeq(func,bornes,arg,choix,EPS,bla) 
       !Calcule la somme des integrales de func de bornes(m) à bornes(m+1)
       !avec le changement de variable choix(m) et les arguments arg(m)
       REAL(QP), DIMENSION(:), INTENT(IN) :: bornes
       REAL(QP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(QP), INTENT(IN) :: EPS
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funcq) :: func
       REAL(QP) decoupeq
       !
       INTEGER m
       REAL(QP) I
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decoupeq"
       if(size(choix)      .NE.size(bornes)-1) stop "size(choix) ne correspond pas à size(bornes)-1 dans decoupeq"
       decoupeq=0.0_qp
       do m=1,size(bornes)-1
        I=qromochoixq(func,bornes(m),bornes(m+1),arg(:,m),choix(m),EPS)
        decoupeq=decoupeq+I
        if(bla) write(6,FMT="(A1,I1,A3,1P,G45.35)") "I",m,"=  ",I
       enddo
       if(bla) write(6,FMT="(A5,1P,G45.35)") "Itot=",decoupeq
       END FUNCTION decoupeq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntc(func,a,b,arg,s,n)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(DPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funcc) :: func
       !
       REAL(DP) :: del
       INTEGER(I4B) :: it
       REAL(DP), DIMENSION(2*3**(n-2)) :: x
       INCLUDE "somme.f90"
       END SUBROUTINE midpntc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfc(funk,aa,bb,arg,s,n)
       INCLUDE "declc.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinfc'
       b=1.0_dp/aa
       a= 1.0_dp/bb
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=funk(1.0_dp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsquc(funk,aa,bb,arg,s,n)
       INCLUDE "declc.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsquc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlc(funk,aa,bb,arg,s,n)
       INCLUDE "declc.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=2.0_dp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqlc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexpc(funk,aa,bb,arg,s,n)
       INCLUDE "declc.f90"
       b=exp(-aa) 
       a=0.0 
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=funk(-log(x),arg)/x
        END FUNCTION func
       END SUBROUTINE midexpc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfc(funk,aa,bb,arg,s,n)
       INCLUDE "declc.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_dp/sqrt(aa)
       a= 1.0_dp/sqrt(bb)
       CONTAINS
        FUNCTION func(x,arg)
        REAL(DP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(DPC), DIMENSION(size(x)) :: func
        func=2.0_dp*funk(1.0_dp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromoc(func,a,b,arg,choose,EPS)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), DIMENSION(:), INTENT(IN) :: arg
       REAL(DP), INTENT(IN)  :: EPS
       COMPLEX(DPC) :: qromoc
       PROCEDURE(funcc) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         IMPORT
         REAL(DP), INTENT(IN) :: aa,bb
         REAL(DP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(DPC), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcc) :: funk
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
       !
       REAL(DP), DIMENSION(JMAXP) :: h
       COMPLEX(DPC), DIMENSION(JMAXP) :: s
       COMPLEX(DPC) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromoc,dqromo)
         if ((abs(dqromo) <= EPS*abs(qromoc)).AND.(abs(s(j)-s(j-1))<0.01_qp*abs(s(j)))) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_dp
       end do
       STOP 'Nombre d itération dépassé dans qromoc'
       END FUNCTION qromoc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixc(func,a,b,arg,choix,EPS)
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: choix
       REAL(DP), INTENT(IN)  :: EPS
       COMPLEX(DPC) :: qromochoixc
       PROCEDURE(funcc) :: func
       if(choix==mpnt)then 
        qromochoixc=qromo(func,a,b,arg,midpntc,EPS)
       elseif(choix==minf)then
        qromochoixc=qromo(func,a,b,arg,midinfc,EPS)
       elseif(choix==msql)then
        qromochoixc=qromo(func,a,b,arg,midsqlc,EPS)
       elseif(choix==msqu)then
        qromochoixc=qromo(func,a,b,arg,midsquc,EPS)
       elseif(choix==mexp)then
        qromochoixc=qromo(func,a,b,arg,midexpc,EPS)
       elseif(choix==rinf)then
        qromochoixc=qromo(func,a,b,arg,racinfc,EPS)
       endif
       END FUNCTION qromochoixc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decoupec(func,bornes,arg,choix,EPS,bla) 
       REAL(DP), DIMENSION(:), INTENT(IN) :: bornes
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(DP), INTENT(IN) :: EPS
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funcc) :: func
       COMPLEX(DPC) decoupec
       !
       INTEGER m
       COMPLEX(DPC) I
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decoupec"
       if(size(choix)      .NE.size(bornes)-1) stop "size(choix) ne correspond pas à size(bornes)-1 dans decoupec"
       decoupec=0.0_dp
       do m=1,size(bornes)-1
        I=qromochoixc(func,bornes(m),bornes(m+1),arg(:,m),choix(m),EPS)
        decoupec=decoupec+I
        if(bla) write(6,FMT="(A1,I1,A3,1P,2G25.18)") "I",m,"=  ",I
       enddo
       if(bla) write(6,FMT="(A5,1P,2G25.18)") "Itot=",decoupec
       END FUNCTION decoupec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntcq(func,a,b,arg,s,n)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT) :: s
       INTEGER(I4B), INTENT(IN) :: n
       PROCEDURE(funccq) :: func
       !
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       INCLUDE "somme.f90"
       END SUBROUTINE midpntcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfcq(funk,aa,bb,arg,s,n)
       INCLUDE "declcq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinfcq'
       b=1.0_qp/aa
       a= 1.0_qp/bb
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=funk(1.0_qp/x,arg)/x**2
        END FUNCTION func
       END SUBROUTINE midinfcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqucq(funk,aa,bb,arg,s,n)
       INCLUDE "declcq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(bb-x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqucq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlcq(funk,aa,bb,arg,s,n)
       INCLUDE "declcq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=2.0_qp*x*funk(aa+x**2,arg)
        END FUNCTION func
       END SUBROUTINE midsqlcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexpcq(funk,aa,bb,arg,s,n)
       INCLUDE "declcq.f90"
       b=exp(-aa) 
       a=0.0 
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg) 
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=funk(-log(x),arg)/x
        END FUNCTION func
       END SUBROUTINE midexpcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfcq(funk,aa,bb,arg,s,n)
       INCLUDE "declcq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa)
       a= 1.0_qp/sqrt(bb)
       INCLUDE "somme.f90"
       CONTAINS
        FUNCTION func(x,arg)
        REAL(QP), DIMENSION(:), INTENT(IN) :: x,arg
        COMPLEX(QPC), DIMENSION(size(x)) :: func
        func=2.0_qp*funk(1.0_qp/x**2,arg)/x**3
        END FUNCTION func
       END SUBROUTINE racinfcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromocq(func,a,b,arg,choose,EPS)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       REAL(QP), INTENT(IN)  :: EPS
       COMPLEX(QPC) :: qromocq
       PROCEDURE(funccq) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,arg,s,n)
         IMPORT
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funccq) :: funk
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=17,JMAXP=JMAX+1,K=5,KM=K-1
       !
       REAL(QP), DIMENSION(JMAXP) :: h
       COMPLEX(QPC), DIMENSION(JMAXP) :: s
       COMPLEX(QPC) :: dqromo
       INTEGER(I4B) :: j
       h(1)=1.0
       do j=1,JMAX
        call choose(func,a,b,arg,s(j),j)
        if (j >= K) then
         call polint(h(j-KM:j),s(j-KM:j),0.0_qp,qromocq,dqromo)
         if ((abs(dqromo) <= EPS*abs(qromocq)).AND.(abs(s(j)-s(j-1))<0.01_qp*abs(s(j)))) RETURN
         if ((abs(qromocq) <= EPS*10e-30_qp)) RETURN
        end if
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0_qp
       end do
       STOP 'Nombre d itération dépassé dans qromocq'
       END FUNCTION qromocq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixcq(func,a,b,arg,choix,EPS)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: choix
       REAL(QP), INTENT(IN)  :: EPS
       COMPLEX(QPC) :: qromochoixcq
       PROCEDURE(funccq) :: func
       if(choix==mpnt)then 
        qromochoixcq=qromo(func,a,b,arg,midpntcq,EPS)
       elseif(choix==minf)then
        qromochoixcq=qromo(func,a,b,arg,midinfcq,EPS)
       elseif(choix==msql)then
        qromochoixcq=qromo(func,a,b,arg,midsqlcq,EPS)
       elseif(choix==msqu)then
        qromochoixcq=qromo(func,a,b,arg,midsqucq,EPS)
       elseif(choix==mexp)then
        qromochoixcq=qromo(func,a,b,arg,midexpcq,EPS)
       elseif(choix==rinf)then
        qromochoixcq=qromo(func,a,b,arg,racinfcq,EPS)
       endif
       END FUNCTION qromochoixcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decoupecq(func,bornes,arg,choix,EPS,bla) 
       !Calcule la somme des integrales de func de bornes(m) à bornes(m+1)
       !avec le changement de variable choix(m) et les arguments arg(m)
       REAL(QP), DIMENSION(:), INTENT(IN) :: bornes
       REAL(QP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(QP), INTENT(IN) :: EPS
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funccq) :: func
       COMPLEX(QPC) decoupecq
       !
       INTEGER m
       COMPLEX(QPC) I
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decoupecq"
       if(size(choix)      .NE.size(bornes)-1) stop "size(choix) ne correspond pas à size(bornes)-1 dans decoupecq"
       decoupecq=0.0_qp
       do m=1,size(bornes)-1
        I=qromochoixcq(func,bornes(m),bornes(m+1),arg(:,m),choix(m),EPS)
        decoupecq=decoupecq+I
        if(bla) write(6,FMT="(A1,I1,A3,1P,2G45.35)") "I",m,"=  ",I
       enddo
       if(bla) write(6,FMT="(A5,1P,2G45.35)") "Itot=",decoupecq
       END FUNCTION decoupecq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntvq(func,a,b,m,arg,s,n)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       REAL(QP), INTENT(INOUT), DIMENSION(m) :: s
       INTEGER(I4B), INTENT(IN) :: m,n
       PROCEDURE(funcvq) :: func
       !
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s(:)=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg,m ),dim=1)
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it)
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s(:)=s(:)/3.0_qp+del*sum(func(x,arg,m),dim=1)
       end if
       END SUBROUTINE midpntvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsquvq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "sommevq.f90"
        func=2.0_qp*funk(bb-x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)*x(is)
        enddo
        END FUNCTION func
       END SUBROUTINE midsquvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlvq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "sommevq.f90"
        func=2.0_qp*funk(aa+x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)*x(is)
        enddo
        END FUNCTION func
       END SUBROUTINE midsqlvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexpvq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvq.f90"
       b=exp(-aa) 
       a=0.0 
       INCLUDE "sommevq.f90"
        func=funk(-log(x),arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)/x(is)
        enddo
        END FUNCTION func
       END SUBROUTINE midexpvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfvq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_qp/aa
       a= 1.0_qp/bb
       INCLUDE "sommevq.f90"
        func=funk(1.0_qp/x,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)/x(is)**2
        enddo
        END FUNCTION func
       END SUBROUTINE midinfvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfvq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa)
       a= 1.0_qp/sqrt(bb)
       INCLUDE "sommevq.f90"
        func=funk(1.0_qp/x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=2.0_qp*func(is,:)/x(is)**3
        enddo
        END FUNCTION func
       END SUBROUTINE racinfvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! qromov(f,a,b,dim,arg,varchange,EPS): Romberg integration of the function f
! with values in R^dim from a to b with precision EPS. arg is a vector of 
! static parameters of the function,
! varchange is a subroutine that performs a change of variable for improper integrals:
! varchange=midpntvq when the function is integrated over a compact interval where it takes finite values
! varchange=midinfvq when b -> +oo (b can be as large as allowed by machine precision) and f decays at least as 1/x^2
! varchange=racinfvq when b -> +oo and f decays as 1/x^(3/2)
! varchange=midsquvq/midsqlvq f has a 1/sqrt(b-x) or 1/sqrt(x-a) (integrable) divergence at the upper/lower bound of the integration interval
       FUNCTION qromovq(func,a,b,m,arg,choose,EPS)
       INTEGER,  INTENT(IN) :: m !Dimension du vecteur à intégrer
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       REAL(QP), DIMENSION(m) :: qromovq
       REAL(QP), INTENT(IN)  :: EPS
       PROCEDURE(funcvq) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,mm,arg,s,n)
         IMPORT
         INTEGER, INTENT(IN) :: mm
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(mm), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcvq) :: funk
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=16,JMAXP=JMAX+1,K=5,KM=K-1
       !
       REAL(QP), DIMENSION(JMAXP,m) :: h,s
       REAL(QP) :: dqromo(m)
       INTEGER(I4B) :: j,im
       LOGICAL conv(m)
       h(1,:)=1.0
       do j=1,JMAX
        call choose(func,a,b,m,arg,s(j,:),j)
!        if(EPS>1.0e-6_qp) write(6,*)"j,s(j)=",j,s(:,1)
        if (j >= K) then
         do im=1,m
          call polint(h(j-KM:j,im),s(j-KM:j,im),0.0_qp,qromovq(im),dqromo(im))
          conv(im)=(abs(dqromo(im)) <= EPS*abs(qromovq(im))).AND.(abs(s(j,im)-s(j-1,im))<0.01_qp*abs(s(j,im)))
         enddo
!         if(EPS>1.0e-6_qp) write(6,*)"j,conv=",j,conv
!         if((j==8).AND.(EPS>1.0e-6_qp)) write(6,*)"j,conv=",j,conv,abs(dqromo(1)),EPS*abs(qromovq(1)),"intpp"
!         if((j==8).AND.(EPS<1.0e-6_qp)) write(6,*)"j,conv=",j,conv,abs(dqromo(1)),EPS*abs(qromovq(1)),"rhopp"
!         if(j==8) write(6,*)"s=",s(:,1)
         if (all(conv)) RETURN
        end if
        s(j+1,:)=s(j,:)
        h(j+1,:)=h(j,:)/9.0_qp
       end do
       write(6,*) 'Nombre d itération dépassé dans qromovq'
       END FUNCTION qromovq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Same of qromov but with romberg iterations limited to JMAX
! err=.TRUE. is return if convergence to precision EPS is not reached at JMAX
       FUNCTION qromovqfixed(func,a,b,m,arg,choose,EPS,JMAX,err)
       INTEGER,  INTENT(IN) :: m !Dimension du vecteur à intégrer
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       REAL(QP), DIMENSION(m) :: qromovqfixed
       REAL(QP), INTENT(IN)  :: EPS
       INTEGER,  INTENT(IN)  :: JMAX
       LOGICAL, INTENT(OUT)  :: err
       PROCEDURE(funcvq) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,mm,arg,s,n)
         IMPORT
         INTEGER, INTENT(IN) :: mm
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         REAL(QP), DIMENSION(mm), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcvq) :: funk
         END SUBROUTINE choose
       END INTERFACE
       !
       REAL(QP), DIMENSION(JMAX+1,m) :: h,s
       REAL(QP) :: dqromo(m)
       INTEGER(I4B) :: j,im
       LOGICAL conv(m)
       err=.FALSE.
       h(1,:)=1.0
       do j=1,JMAX
        call choose(func,a,b,m,arg,s(j,:),j)
        if (j >= JMAX) then
         do im=1,m
          call polint(h(j-(JMAX-1):j,im),s(j-(JMAX-1):j,im),0.0_qp,qromovqfixed(im),dqromo(im))
          conv(im)=(abs(dqromo(im)) <= EPS*abs(qromovqfixed(im))).AND.(abs(s(j,im)-s(j-1,im))<0.01_qp*abs(s(j,im)))
         enddo
         if (all(conv)) RETURN
        end if
        s(j+1,:)=s(j,:)
        h(j+1,:)=h(j,:)/9.0_qp
       end do
       err=.TRUE.
       END FUNCTION qromovqfixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixvq(func,a,b,m,arg,choix,EPS)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: m,choix
       REAL(QP), INTENT(IN)  :: EPS
       REAL(QP) :: qromochoixvq(1:m)
       PROCEDURE(funcvq) :: func
       if(choix==mpnt)then 
        qromochoixvq=qromovq(func,a,b,m,arg,midpntvq,EPS)
       elseif(choix==minf)then
        qromochoixvq=qromovq(func,a,b,m,arg,midinfvq,EPS)
       elseif(choix==msql)then
        qromochoixvq=qromo(func,a,b,m,arg,midsqlvq,EPS)
       elseif(choix==msqu)then
        qromochoixvq=qromo(func,a,b,m,arg,midsquvq,EPS)
       elseif(choix==mexp)then
        qromochoixvq=qromo(func,a,b,m,arg,midexpvq,EPS)
       elseif(choix==rinf)then
        qromochoixvq=qromo(func,a,b,m,arg,racinfvq,EPS)
       endif
       END FUNCTION qromochoixvq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decoupevq(func,bornes,mm,arg,choix,EPS,bla) 
       !Calcule la somme des integrales de func de bornes(m) à bornes(m+1)
       !avec le changement de variable choix(m) et les arguments arg(m)
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       REAL(QP), DIMENSION(size(choix)+1), INTENT(IN) :: bornes
       REAL(QP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(QP), INTENT(IN) :: EPS
       INTEGER, INTENT(IN) :: mm
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funcvq) :: func
       REAL(QP) decoupevq(1:mm)
       !
       INTEGER m,imm
       REAL(QP) I(1:mm)
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decoupevq"
       decoupevq(:)=0.0_qp
       do m=1,size(choix)

        if(bla)then
         write(6,*)"************************"
         write(6,*)
         write(6,*)"Intégrale de ",bornes(m)," à ",bornes(m+1)
         write(6,*)
        endif

        I=qromochoixvq(func,bornes(m),bornes(m+1),mm,arg(:,m),choix(m),EPS)
        decoupevq=decoupevq+I

        if(bla)then
         write(6,*)
         do imm=1,mm
          write(6,FMT="(A2,I1,A1,I1,A3,1P,G45.35)") "I",m,"(",imm,")=  ",I(imm)
         enddo
         write(6,*)
        endif
       enddo
       if(bla)then
         write(6,*)
         do imm=1,mm
          write(6,FMT="(A6,I1,A2,1P,G45.35)") "Itot(",imm,")=",decoupevq(imm)
         enddo
         write(6,*)
       endif

       END FUNCTION decoupevq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midpntvcq(func,a,b,m,arg,s,n)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), INTENT(IN), DIMENSION(:) :: arg
       COMPLEX(QPC), INTENT(INOUT), DIMENSION(m) :: s
       INTEGER(I4B), INTENT(IN) :: m,n
       PROCEDURE(funcvcq) :: func
       !
       REAL(QP) :: del
       INTEGER(I4B) :: it
       REAL(QP), DIMENSION(2*3**(n-2)) :: x
       if (n == 1) then
        s(:)=(b-a)*sum(func( (/0.5_qp*(a+b)/),arg,m ),dim=1)
       else
        it=3**(n-2)
        del=(b-a)/(3.0_qp*it)
        x(1:2*it-1:2)=arth(a+0.5_qp*del,3.0_qp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_qp*del
        s(:)=s(:)/3.0_qp+del*sum(func(x,arg,m),dim=1)
       end if
       END SUBROUTINE midpntvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsquvcq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvcq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "sommevcq.f90"
        func=2.0_qp*funk(bb-x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)*x(is)
        enddo
        END FUNCTION func
       END SUBROUTINE midsquvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midsqlvcq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvcq.f90"
       b=sqrt(bb-aa)
       a= 0.0
       INCLUDE "sommevcq.f90"
        func=2.0_qp*funk(aa+x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)*x(is)
        enddo
        END FUNCTION func
       END SUBROUTINE midsqlvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midexpvcq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvcq.f90"
       b=exp(-aa) 
       a=0.0 
       INCLUDE "sommevcq.f90"
        func=funk(-log(x),arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)/x(is)
        enddo
        END FUNCTION func
       END SUBROUTINE midexpvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE midinfvcq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvcq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans midinf'
       b=1.0_qp/aa
       a= 1.0_qp/bb
       INCLUDE "sommevcq.f90"
        func=funk(1.0_qp/x,arg,mm)
        do is=1,size(x)
         func(is,:)=func(is,:)/x(is)**2
        enddo
        END FUNCTION func
       END SUBROUTINE midinfvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE racinfvcq(funk,aa,bb,m,arg,s,n)
       INCLUDE "declvcq.f90"
       if(aa*bb <= 0.0) STOP 'bornes dans racinf'
       b=1.0_qp/sqrt(aa)
       a= 1.0_qp/sqrt(bb)
       INCLUDE "sommevcq.f90"
        func=funk(1.0_qp/x**2,arg,mm)
        do is=1,size(x)
         func(is,:)=2.0_qp*func(is,:)/x(is)**3
        enddo
        END FUNCTION func
       END SUBROUTINE racinfvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromovcq(func,a,b,m,arg,choose,EPS)
       INTEGER,  INTENT(IN) :: m !Dimension du vecteur à intégrer
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       COMPLEX(QPC), DIMENSION(m) :: qromovcq
       REAL(QP), INTENT(IN)  :: EPS
       PROCEDURE(funcvcq) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,mm,arg,s,n)
         IMPORT
         INTEGER, INTENT(IN) :: mm
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(mm), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcvcq) :: funk
         END SUBROUTINE choose
       END INTERFACE
       INTEGER(I4B), PARAMETER :: JMAX=16,JMAXP=JMAX+1,K=5,KM=K-1
       !
       COMPLEX(QPC), DIMENSION(JMAXP,m) :: s
       REAL(QP), DIMENSION(JMAXP,m) :: h
       COMPLEX(QPC) :: dqromo(m)
       INTEGER(I4B) :: j,im
       LOGICAL conv(m)
       h(1,:)=1.0
       do j=1,JMAX
        call choose(func,a,b,m,arg,s(j,:),j)
        if (j >= K) then
         do im=1,m
          call polint(h(j-KM:j,im),s(j-KM:j,im),0.0_qp,qromovcq(im),dqromo(im))
          conv(im)=abs(dqromo(im)) <= EPS*abs(qromovcq(im))
         enddo
!         if(omp_get_thread_num()==4) write(6,*)"j=",j
!         if(omp_get_thread_num()==4) write(6,*)"abs(dqromo)=",abs(dqromo)
!         if(omp_get_thread_num()==4) write(6,*)"qromo=",real(qromovcq(1:3))
!         if(omp_get_thread_num()==4) write(6,*)"conv=",conv
         if (all(conv)) RETURN
        end if
        s(j+1,:)=s(j,:)
        h(j+1,:)=h(j,:)/9.0_qp
       end do
       write(6,*) 'Nombre d itération dépassé dans qromovcq'
       END FUNCTION qromovcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromovcqfixed(func,a,b,m,arg,choose,EPS,JMAX,err)
       INTEGER,  INTENT(IN) :: m !Dimension du vecteur à intégrer
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       COMPLEX(QPC), DIMENSION(m) :: qromovcqfixed
       REAL(QP), INTENT(IN)  :: EPS
       INTEGER,  INTENT(IN)  :: JMAX
       LOGICAL, INTENT(OUT)  :: err
       PROCEDURE(funcvcq) :: func
       INTERFACE
         SUBROUTINE choose(funk,aa,bb,mm,arg,s,n)
         IMPORT
         INTEGER, INTENT(IN) :: mm
         REAL(QP), INTENT(IN) :: aa,bb
         REAL(QP), DIMENSION(:), INTENT(IN) :: arg
         COMPLEX(QPC), DIMENSION(mm), INTENT(INOUT) :: s
         INTEGER(I4B), INTENT(IN) :: n
         PROCEDURE(funcvcq) :: funk
         END SUBROUTINE choose
       END INTERFACE
       !
       COMPLEX(QPC), DIMENSION(JMAX+1,m) :: s
       REAL(QP), DIMENSION(JMAX+1,m) :: h
       COMPLEX(QPC) :: dqromo(m)
       INTEGER(I4B) :: j,im
       LOGICAL conv(m)
       err=.FALSE.
       h(1,:)=1.0
       do j=1,JMAX
        call choose(func,a,b,m,arg,s(j,:),j)
        if (j >= JMAX) then
         do im=1,m
          call polint(h(j-(JMAX-1):j,im),s(j-(JMAX-1):j,im),0.0_qp,qromovcqfixed(im),dqromo(im))
          conv(im)=(abs(dqromo(im)) <= EPS*abs(qromovcqfixed(im))).AND.(abs(s(j,im)-s(j-1,im))<0.01_qp*abs(s(j,im)))
         enddo
         if (all(conv)) RETURN
        end if
        s(j+1,:)=s(j,:)
        h(j+1,:)=h(j,:)/9.0_qp
       end do
       err=.TRUE.
       END FUNCTION qromovcqfixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION qromochoixvcq(func,a,b,m,arg,choix,EPS)
       REAL(QP), INTENT(IN) :: a,b
       REAL(QP), DIMENSION(:), INTENT(IN) :: arg
       INTEGER, INTENT(IN) :: m,choix
       REAL(QP), INTENT(IN)  :: EPS
       COMPLEX(QPC) :: qromochoixvcq(1:m)
       PROCEDURE(funcvcq) :: func
       if(choix==mpnt)then 
        qromochoixvcq=qromovcq(func,a,b,m,arg,midpntvcq,EPS)
       elseif(choix==minf)then
        qromochoixvcq=qromovcq(func,a,b,m,arg,midinfvcq,EPS)
       elseif(choix==msql)then
        qromochoixvcq=qromo(func,a,b,m,arg,midsqlvcq,EPS)
       elseif(choix==msqu)then
        qromochoixvcq=qromo(func,a,b,m,arg,midsquvcq,EPS)
       elseif(choix==mexp)then
        qromochoixvcq=qromo(func,a,b,m,arg,midexpvcq,EPS)
       elseif(choix==rinf)then
        qromochoixvcq=qromo(func,a,b,m,arg,racinfvcq,EPS)
       endif
       END FUNCTION qromochoixvcq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       FUNCTION decoupevcq(func,bornes,mm,arg,choix,EPS,bla) 
       !Calcule la somme des integrales de func de bornes(m) à bornes(m+1)
       !avec le changement de variable choix(m) et les arguments arg(m)
       INTEGER, DIMENSION(:), INTENT(IN) :: choix
       REAL(QP), DIMENSION(size(choix)+1), INTENT(IN) :: bornes
       REAL(QP), DIMENSION(:,:), INTENT(IN) :: arg
       REAL(QP), INTENT(IN) :: EPS
       INTEGER, INTENT(IN) :: mm
       LOGICAL, INTENT(IN) :: bla
       PROCEDURE(funcvcq) :: func
       COMPLEX(QPC) decoupevcq(1:mm)
       !
       INTEGER m,imm
       COMPLEX(QPC) I(1:mm)
       if(size(arg,  DIM=2).NE.size(bornes)-1) stop "size(arg)   ne correspond pas à size(bornes)-1 dans decoupevcq"
       decoupevcq(:)=0.0_qp
       do m=1,size(choix)

        if(bla)then
         write(6,*)"************************"
         write(6,*)
         write(6,*)"Intégrale de ",bornes(m)," à ",bornes(m+1)
         write(6,*)
        endif

        I=qromochoixvcq(func,bornes(m),bornes(m+1),mm,arg(:,m),choix(m),EPS)
        decoupevcq=decoupevcq+I

        if(bla)then
         write(6,*)
         do imm=1,mm
          write(6,FMT="(A2,I1,A1,I1,A3,1P,2G45.35)") "I",m,"(",imm,")=  ",I(imm)
         enddo
         write(6,*)
        endif
       enddo
       if(bla)then
         write(6,*)
         do imm=1,mm
          write(6,FMT="(A6,I1,A2,1P,2G45.35)") "Itot(",imm,")=",decoupevcq(imm)
         enddo
         write(6,*)
       endif

       END FUNCTION decoupevcq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       FUNCTION ecritrout(taille,config) 
       INTEGER, INTENT(IN) :: taille 
       INTEGER, DIMENSION(:), INTENT(IN) :: config 
       CHARACTER(len=90) :: ecritrout
        
       INTEGER itai 
       ecritrout=trim(ecritr(config(1))) 
       do itai=2,taille 
        ecritrout=trim(ecritrout)//"  "//trim(ecritr(config(itai))) 
       enddo 
       END FUNCTION ecritrout 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       FUNCTION ecritr(r) 
       INTEGER, INTENT(IN) :: r 
       CHARACTER(len=7) :: ecritr 
       if(r==1)then 
        ecritr="mpnt" 
       elseif(r==2)then 
        ecritr="minf" 
       elseif(r==3)then 
        ecritr="msql" 
       elseif(r==4)then 
        ecritr="msqu" 
       elseif(r==5)then 
        ecritr="mexp" 
       elseif(r==6)then 
        ecritr="rinf" 
       endif 
       END FUNCTION ecritr 
END MODULE modsim
