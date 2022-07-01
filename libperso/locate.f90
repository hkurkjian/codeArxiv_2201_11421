       !Given an array xx(1:N), and given a value x, returns a value j such that x is between
       !xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
       !j = N is returned to indicate that x is out of range.
       INTEGER(I4B) :: n,jl,jm,ju
       LOGICAL :: ascnd
       n=size(xx)
       ascnd = (xx(n) >= xx(1)) !True if ascending order of table, false otherwise.
       jl=0 !Initialize lower
       ju=n+1 !and upper limits.
       do
        if (ju-jl <= 1) exit!Repeat until this condition is satisfied.
        jm=(ju+jl)/2 !Compute a midpoint,
        if (ascnd .eqv. (x >= xx(jm))) then
         jl=jm !and replace either the lower limit
        else
         ju=jm !or the upper limit, as appropriate.
        end if
       end do
       if (x == xx(1)) then !Then set the output, being careful with the endpoints.
        pos=1
       else if (x == xx(n)) then
        pos=n-1
       else
        pos=jl
       end if

