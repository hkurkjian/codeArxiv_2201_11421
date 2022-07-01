set terminal epslatex color size 12cm,10cm lw 2
set encoding utf8
set key on

tabal     ="Tabal.dat"
tabalb01  ="Tabalb01.dat"

system  "rm alpha.eps"
set output 'alpha.tex'


set key bottom left Left reverse

set multiplot
set origin 0,0
set size 1,1

set xrange [0:4.5]
set yrange [-4:10]

pi=3.14159

set arrow from first 0,0 to first 2,0 lw 1 lc rgb "black" nohead

set xlabel "\\Large $\\omega_p/\\Delta$"
set ylabel "\\Large $\\Delta\\text{Re}\\;\\alpha/\\epsilon_F$"

filtre1(x)=(x>=2.0?x:1/0)

set arrow from first 2,graph 0 to first 2,graph 1 lw 2 dt 3 lc rgb "red" nohead

set arrow from first 1.696,graph 0 to first 1.696,graph 1 lw 2 dt 3 lc rgb "black" nohead

set x2tics -1000,1000,10000
set x2tics add ("$1.696$" 1.696)

plot tabal u (filtre1($1)):(-10):(10) with filledcu fc rgb "red" notitle fs transparent solid 0.1,\
     tabal    u 1:2       w lines lw 2      lc rgb "black" title "$T=0$",\
     tabalb01 u 1:2       w lines lw 2 lc rgb "red" title "$T/\\Delta=10$",\
     tabalb01 u 1:4       w lines lw 2 dt 3 lc rgb "red" title "$T\\to T_c$",\

#     tabal u 1:(2./3./$1)              w lines dt 2 lc rgb "blue" title "$2/3\\omega_p$",\
#     tabal u 1:(-2*pi/(15*sqrt(2-$1))) w lines dt 2 lc rgb "red" title "$-2\\pi/15\\sqrt{2-\\omega_p}$",\


set origin .392,.39
set size 0.6,0.55

set xrange [2.0:6.0]
set yrange [0:2]

unset xtics
unset x2tics
unset arrow
set xtics auto
set xtics 2,1,6
set tics front

set format y "{\\small $%g$}"
set format x "{\\small $%g$}"

set xlabel "$\\omega_p/\\Delta$" offset 0,0.5
set ylabel "$\\Delta\\text{Im}\\;\\alpha/\\epsilon_F$" offset 4.5

plot tabal    u (filtre1($1)):3       w lines lw 2 lc rgb "black" notitle,\
     tabalb01 u (filtre1($1)):(-$3)   w lines lw 2 lc rgb "red" notitle,\

#plot tabal u 1:2                         w lines lw 2 lc rgb "black" title "$\\textrm{Re}\\bar\\alpha$",\
#     tabal u 1:3                         w lines lw 2 lc rgb "red"   title "$-\\textrm{Im}\\bar\\alpha$",\
#     tabal u 1:(2*pi/(15.*sqrt($1-2)))   w lines dt 2 lc rgb "red"   title "$-2\\pi/15\\sqrt{\\omega_p-2}$",\
#     tabal u 1:(6./(5.*$1))              w lines dt 2 lc rgb "black" title "$6/5\\omega_p$",\
#     tabal u 1:(16.*pi/(15.*$1**3))      w lines dt 3 lc rgb "red"   title "$16\\pi/15\\omega_p^3$",\
