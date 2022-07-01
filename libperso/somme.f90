       if (n == 1) then
        s=(b-a)*sum(func( (/0.5_sp*(a+b)/),arg ))
       else
        it=3**(n-2)
        del=(b-a)/(3.0_sp*it) !The added points alternate in spacing between del and 2*del.
        x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
        s=s/3.0_sp+del*sum(func(x,arg)) !The new sum is combined with the old integral
       end if !to give a refined integral.
