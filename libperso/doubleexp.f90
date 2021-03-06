MODULE doubleexp
USE nrtype
! DE-Quadrature
! Numerical Automatic Integrator for Improper Integral
!     method    : Double Exponential (DE) Transformation
!     dimension : one
!     table     : not use
! subroutines
!     intde  : integrator of f(x) over (a,b).
!     intdei : integrator of f(x) over (a,infinity), 
!                  f(x) is non oscillatory function.
!     intdeo : integrator of f(x) over (a,infinity), 
!                  f(x) is oscillatory function.
!
!
! intde
!     [description]
!         I = integral of f(x) over (a,b)
!     [declaration]
!         external f
!     [usage]
!         call intde(f, a, b, eps, i, err)
!     [parameters]
!         f         : integrand f(x) (REAL(QP) function)
!         a         : lower limit of integration (REAL(QP))
!         b         : upper limit of integration (REAL(QP))
!         eps       : relative error requested (REAL(QP))
!         i         : approximation to the integral (REAL(QP))
!         err       : estimate of the absolute error (REAL(QP))
!     [remarks]
!         function
!             f(x) needs to be analytic over (a,b).
!         relative error
!             eps is relative error requested excluding 
!             cancellation of significant digits.
!             i.e. eps means : (absolute error) / 
!                              (integral_a^b |f(x)| dx).
!             eps does not mean : (absolute error) / I.
!         error message
!             err >= 0 : normal termination.
!             err < 0  : abnormal termination (m >= mmax).
!                        i.e. convergent error is detected :
!                            1. f(x) or (d/dx)^n f(x) has 
!                               discontinuous points or sharp 
!                               peaks over (a,b).
!                               you must divide the interval 
!                               (a,b) at this points.
!                            2. relative error of f(x) is 
!                               greater than eps.
!                            3. f(x) has oscillatory factor 
!                               and frequency of the oscillation 
!                               is very high.
!
!
! intdei
!     [description]
!         I = integral of f(x) over (a,infinity), 
!             f(x) has not oscillatory factor.
!     [declaration]
!         external f
!     [usage]
!         call intdei(f, a, eps, i, err)
!     [parameters]
!         f         : integrand f(x) (REAL(QP) function)
!         a         : lower limit of integration (REAL(QP))
!         eps       : relative error requested (REAL(QP))
!         i         : approximation to the integral (REAL(QP))
!         err       : estimate of the absolute error (REAL(QP))
!     [remarks]
!         function
!             f(x) needs to be analytic over (a,infinity).
!         relative error
!             eps is relative error requested excluding 
!             cancellation of significant digits.
!             i.e. eps means : (absolute error) / 
!                              (integral_a^infinity |f(x)| dx).
!             eps does not mean : (absolute error) / I.
!         error message
!             err >= 0 : normal termination.
!             err < 0  : abnormal termination (m >= mmax).
!                        i.e. convergent error is detected :
!                            1. f(x) or (d/dx)^n f(x) has 
!                               discontinuous points or sharp 
!                               peaks over (a,infinity).
!                               you must divide the interval 
!                               (a,infinity) at this points.
!                            2. relative error of f(x) is 
!                               greater than eps.
!                            3. f(x) has oscillatory factor 
!                               and decay of f(x) is very slow 
!                               as x -> infinity.
!
!
! intdeo
!     [description]
!         I = integral of f(x) over (a,infinity), 
!             f(x) has oscillatory factor :
!             f(x) = g(x) * sin(omega * x + theta) as x -> infinity.
!     [declaration]
!         external f
!     [usage]
!         call intdeo(f, a, omega, eps, i, err)
!     [parameters]
!         f         : integrand f(x) (REAL(QP) function)
!         a         : lower limit of integration (REAL(QP))
!         omega     : frequency of oscillation (REAL(QP))
!         eps       : relative error requested (REAL(QP))
!         i         : approximation to the integral (REAL(QP))
!         err       : estimate of the absolute error (REAL(QP))
!     [remarks]
!         function
!             f(x) needs to be analytic over (a,infinity).
!         relative error
!             eps is relative error requested excluding 
!             cancellation of significant digits.
!             i.e. eps means : (absolute error) / 
!                              (integral_a^R |f(x)| dx).
!             eps does not mean : (absolute error) / I.
!         error message
!             err >= 0 : normal termination.
!             err < 0  : abnormal termination (m >= mmax).
!                        i.e. convergent error is detected :
!                            1. f(x) or (d/dx)^n f(x) has 
!                               discontinuous points or sharp 
!                               peaks over (a,infinity).
!                               you must divide the interval 
!                               (a,infinity) at this points.
!                            2. relative error of f(x) is 
!                               greater than eps.
!
!
CONTAINS
      SUBROUTINE intde(f, a, b, eps, i, err)
      REAL(QP), INTENT(IN)  :: a, b, eps
      REAL(QP), INTENT(OUT) :: i, err
      INTERFACE
       FUNCTION f(x)
       USE nrtype
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: x
       REAL(QP)             :: f
       END FUNCTION f
      END INTERFACE

      INTEGER mmax
      REAL(QP) efs, hoff
      INTEGER m
      REAL(QP) pi2, epsln, epsh, h0, ehp, ehm, epst, ba, ir, h
      REAL(QP) iback, irback, t, ep, em, xw, xa, wg, fa, fb, errt
      REAL(QP) errh, errd
! ---- adjustable parameter ----
      mmax = 256
      efs = 0.1_qp
      hoff = 8.5_qp
! ------------------------------
      pi2 = 2 * atan(1.0_qp)
      epsln = 1 - log(efs * eps)
      epsh = sqrt(efs * eps)
      h0 = hoff / epsln
      ehp = exp(h0)
      ehm = 1 / ehp
      epst = exp(-ehm * epsln)
      ba = b - a
      ir = f((a + b) * 0.5_qp) * (ba * 0.25_qp)
      i = ir * (2 * pi2)
      err = abs(i) * epst
      h = 2 * h0
      m = 1
   10 continue
          iback = i
          irback = ir
          t = h * 0.5_qp
   20     continue
              em = exp(t)
              ep = pi2 * em
              em = pi2 / em
   30         continue
                  xw = 1 / (1 + exp(ep - em))
                  xa = ba * xw
                  wg = xa * (1 - xw)
                  fa = f(a + xa) * wg
                  fb = f(b - xa) * wg
                  ir = ir + (fa + fb)
                  i = i + (fa + fb) * (ep + em)
                  errt = (abs(fa) + abs(fb)) * (ep + em)
                  if (m .eq. 1) err = err + errt * epst
                  ep = ep * ehp
                  em = em * ehm
              if (errt .gt. err .or. xw .gt. epsh) goto 30
              t = t + h
          if (t .lt. h0) goto 20
          if (m .eq. 1) then
              errh = (err / epst) * epsh * h0
              errd = 1 + 2 * errh
          else
              errd = h * (abs(i - 2 * iback) + 4 * abs(ir - 2 * irback))
          end if
          h = h * 0.5_qp
          m = m * 2
      if (errd .gt. errh .and. m .lt. mmax) goto 10
      i = i * h
      if (errd .gt. errh) then
          err = -errd * m
      else
          err = errh * epsh * m / (2 * efs)
      end if
      end
!
!
      SUBROUTINE intdei(f, a, eps, i, err)
      REAL(QP), INTENT(IN)  :: a, eps
      REAL(QP), INTENT(OUT) :: i, err
      INTERFACE
       FUNCTION f(x)
       USE nrtype
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: x
       REAL(QP)             :: f
       END FUNCTION f
      END INTERFACE

      INTEGER mmax
      REAL(QP) efs, hoff
      INTEGER m
      REAL(QP) pi4, epsln, epsh, h0, ehp, ehm, epst, ir, h, iback
      REAL(QP) irback, t, ep, em, xp, xm, fp, fm, errt, errh, errd
! ---- adjustable parameter ----
      mmax = 256
      efs = 0.1_qp
      hoff = 11.0_qp
! ------------------------------
      pi4 = atan(1.0_qp)
      epsln = 1 - log(efs * eps)
      epsh = sqrt(efs * eps)
      h0 = hoff / epsln
      ehp = exp(h0)
      ehm = 1 / ehp
      epst = exp(-ehm * epsln)
      ir = f(a + 1)
      i = ir * (2 * pi4)
      err = abs(i) * epst
      h = 2 * h0
      m = 1
   10 continue
          iback = i
          irback = ir
          t = h * 0.5_qp
   20     continue
              em = exp(t)
              ep = pi4 * em
              em = pi4 / em
   30         continue
                  xp = exp(ep - em)
                  xm = 1 / xp
                  fp = f(a + xp) * xp
                  fm = f(a + xm) * xm
                  ir = ir + (fp + fm)
                  i = i + (fp + fm) * (ep + em)
                  errt = (abs(fp) + abs(fm)) * (ep + em)
                  if (m .eq. 1) err = err + errt * epst
                  ep = ep * ehp
                  em = em * ehm
              if (errt .gt. err .or. xm .gt. epsh) goto 30
              t = t + h
          if (t .lt. h0) goto 20
          if (m .eq. 1) then
              errh = (err / epst) * epsh * h0
              errd = 1 + 2 * errh
          else
              errd = h * (abs(i - 2 * iback) + 4 * abs(ir - 2 * irback))
          end if
          h = h * 0.5_qp
          m = m * 2
      if (errd .gt. errh .and. m .lt. mmax) goto 10
      i = i * h
      if (errd .gt. errh) then
          err = -errd * m
      else
          err = errh * epsh * m / (2 * efs)
      end if
      end
!
!
      SUBROUTINE intdeo(f, a, omega, eps, i, err)
      REAL(QP), INTENT(IN)  :: a, omega, eps
      REAL(QP), INTENT(OUT) :: i, err
      INTERFACE
       FUNCTION f(x)
       USE nrtype
       IMPLICIT NONE
       REAL(QP), INTENT(IN) :: x
       REAL(QP)             :: f
       END FUNCTION f
      END INTERFACE

      INTEGER mmax, lmax
      REAL(QP) efs, enoff, pqoff, ppoff
      INTEGER n, m, l, k
      REAL(QP) pi4, epsln, epsh, frq4, per2, pp, pq, ehp, ehm, ir, h
      REAL(QP) iback, irback, t, ep, em, tk, xw, wg, xa, fp, fm, errh
      REAL(QP) tn, errd
! ---- adjustable parameter ----
      mmax = 256
      lmax = 5
      efs = 0.1_qp
      enoff = 0.40_qp
      pqoff = 2.9_qp
      ppoff = -0.72_qp
! ------------------------------
      pi4 = atan(1.0_qp)
      epsln = 1 - log(efs * eps)
      epsh = sqrt(efs * eps)
      n = int(enoff * epsln)
      frq4 = abs(omega) / (2 * pi4)
      per2 = 4 * pi4 / abs(omega)
      pq = pqoff / epsln
      pp = ppoff - log(pq * pq * frq4)
      ehp = exp(2 * pq)
      ehm = 1 / ehp
      xw = exp(pp - 2 * pi4)
      i = f(a + sqrt(xw * (per2 * 0.5_qp)))
      ir = i * xw
      i = i * (per2 * 0.5_qp)
      err = abs(i)
      h = 2
      m = 1
   10 continue
          iback = i
          irback = ir
          t = h * 0.5_qp
   20     continue
              em = exp(2 * pq * t)
              ep = pi4 * em
              em = pi4 / em
              tk = t
   30         continue
                  xw = exp(pp - ep - em)
                  wg = sqrt(frq4 * xw + tk * tk)
                  xa = xw / (tk + wg)
                  wg = (pq * xw * (ep - em) + xa) / wg
                  fm = f(a + xa)
                  fp = f(a + xa + per2 * tk)
                  ir = ir + (fp + fm) * xw
                  fm = fm * wg
                  fp = fp * (per2 - wg)
                  i = i + (fp + fm)
                  if (m .eq. 1) err = err + (abs(fp) + abs(fm))
                  ep = ep * ehp
                  em = em * ehm
                  tk = tk + 1
              if (ep .lt. epsln) goto 30
              if (m .eq. 1) then
                  errh = err * epsh
                  err = err * eps
              end if
              tn = tk
              do while (abs(fm) .gt. err)
                  xw = exp(pp - ep - em)
                  xa = xw / tk * 0.5_qp
                  wg = xa * (1 / tk + 2 * pq * (ep - em))
                  fm = f(a + xa)
                  ir = ir + fm * xw
                  fm = fm * wg
                  i = i + fm
                  ep = ep * ehp
                  em = em * ehm
                  tk = tk + 1
              end do
              fm = f(a + per2 * tn)
              em = per2 * fm
              i = i + em
              if (abs(fp) .gt. err .or. abs(em) .gt. err) then
                  l = 0
   40             continue
                      l = l + 1
                      tn = tn + n
                      em = fm
                      fm = f(a + per2 * tn)
                      xa = fm
                      ep = fm
                      em = em + fm
                      xw = 1
                      wg = 1
                      do k = 1, n - 1
                          xw = xw * (n + 1 - k) / k
                          wg = wg + xw
                          fp = f(a + per2 * (tn - k))
                          xa = xa + fp
                          ep = ep + fp * wg
                          em = em + fp * xw
                      end do
                      wg = per2 * n / (wg * n + xw)
                      em = wg * abs(em)
                      if (em .le. err .or. l .ge. lmax) goto 50
                      i = i + per2 * xa
                  goto 40
   50             continue
                  i = i + wg * ep
                  if (em .gt. err) err = em
              end if
              t = t + h
          if (t .lt. 1) goto 20
          if (m .eq. 1) then
              errd = 1 + 2 * errh
          else
              errd = h * (abs(i - 2 * iback) + pq * abs(ir - 2 * irback))
          end if
          h = h * 0.5_qp
          m = m * 2
      if (errd .gt. errh .and. m .lt. mmax) goto 10
      i = i * h
      if (errd .gt. errh) then
          err = -errd
      else
          err = err * (m * 0.5_qp)
      end if
      end
!
END MODULE doubleexp
