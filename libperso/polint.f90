       n=size(xa)
       if(size(xa).NE.size(ya)) STOP 'Taille des entrees dans polint'
       c=ya !Initialize the tableau of c’s and d’s.
       d=ya
       ho=xa-x
       ns=iminloc(abs(x-xa)) !Find index ns of closest table entry.
       y=ya(ns) !This is the initial approximation to y.
       ns=ns-1
       do m=1,n-1 !For each column of the tableau,
        den(1:n-m)=ho(1:n-m)-ho(1+m:n) !we loop over the current c’s and d’s and upidate them.
        if(any(den(1:n-m) == 0.0))then
          write(6,*)"erreur ds polint, fil:",OMP_GET_THREAD_NUM()
          STOP 'polint: entrées dégénérées'
        endif
!       This error can occur only if two input xa’s are (to within roundoff) identical.
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m) !Here the c’s and d’s are updated.
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then 
!       After each column in the tableau is completed, we decide
!       which correction, c or d, we want to add to our accumulating
!       value of y, i.e., which path to take through
!       the tableau—forking up or down. We do this in such a
!       way as to take the most “straight line” route through the
!       tableau to its apex, updating ns accordingly to keep track
!       of where we are. This route keeps the partial approximations
!       centered (insofar as possible) on the target x. The
!       last dy added is thus the error indication.
         dy=c(ns+1)
        else
         dy=d(ns)
         ns=ns-1
        end if
        y=y+dy
       enddo
