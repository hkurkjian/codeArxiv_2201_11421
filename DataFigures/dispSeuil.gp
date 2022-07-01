set terminal epslatex color size 12cm,10cm lw 2
set encoding utf8
set key on

system  "rm dispSeuil.tex"
set output 'dispSeuil.tex'

om1p9sec0 ="disp19reel.dat"
om1p9sec1 ="disp19sec1.dat"
om1p9sec2 ="disp19sec2.dat"
om1p9higg ="disp19higgs.dat"

set xlabel "\\Large $q\\xi$"
set ylabel "\\Large $\\text{Re}\\omega_q/\\Delta$"

set yrange [1.5:2.5]
set xrange [0.0:3.5]

#set yrange [1.8:2.0]
#set xrange [0.45:1.0]
set arrow from first 0.2,graph 0 to first 0.2,graph 1 dt 3 lw 2 lc rgb "#000000" nohead front
set arrow from first 0.5,graph 0 to first 0.5,graph 1 dt 3 lw 2 lc rgb "#777777" nohead front
set arrow from first 1.0,graph 0 to first 1.0,graph 1 dt 3 lw 2 lc rgb "#CCCCCC" nohead front

set label at first 0.6503917,first 1.8538731 "$\\bullet$"

set label at first 1.05,1.05 "\\fbox{\\colorbox{white}{\\includegraphics[width=6.5cm]{SchemaPlasmons.pdf}}}" front

plot  om1p9sec1 u 1:(2):(2*sqrt(1+$1**2))   w filledcu        fc rgb "blue"     fs transparent solid 0.1 notitle,\
      om1p9sec1 u 1:(2*sqrt(1+$1**2)):(100) w filledcu        fc rgb "red"      fs transparent solid 0.1 notitle,\
      om1p9sec0 u 1:2                       w lines lw 2      lc rgb "black"    title "Sec I",\
      om1p9sec1 u 1:2                       w lines lw 2      lc rgb "blue"     title "Sec II",\
      om1p9higg u 1:2                       w lines lw 2 dt 2 lc rgb "blue"     title "Sec II, modulus mode",\
      om1p9sec2 u 1:2                       w lines lw 2      lc rgb "red"      title "Sec III",\
      om1p9higg u 1:(2)                     w lines lw 1 dt 3 lc rgb "black"    title "$2\\Delta$",\
      om1p9sec2 u 1:(2*sqrt(1+$1**2))       w lines lw 1 dt 3 lc rgb "black"    title "$\\omega_2$",\

system  "rm rep.tex"
system  "rm rep.eps"

set output "rep.tex"

rep1="coupe1p9q02.dat"
rep2="coupe1p9.dat"
rep3="coupe1p9q1.dat"

unset label
set xrange [1.8:3.15]
set yrange [0:5]

set xlabel "\\Large $\\omega/\\Delta$"
set ylabel "\\Large $\\chi_{\\rho\\rho}(\\omega+i0^+)$"

om2_1=2*sqrt(1+0.2**2)
om2_2=2*sqrt(1+0.5**2)
om2_3=2*sqrt(1+1.0**2)

set arrow from first 2  ,graph 0 to first 2  ,graph 1 dt 3 lc rgb "black" nohead

set arrow from first om2_1,graph 0 to first om2_1,graph 1 dt 3 lc rgb "#000000" nohead
set arrow from first om2_2,graph 0 to first om2_2,graph 1 dt 3 lc rgb "#777777" nohead
set arrow from first om2_3,graph 0 to first om2_3,graph 1 dt 3 lc rgb "#CCCCCC" nohead

set arrow from first om2_1,graph 0 to first om2_1,graph 1 dt 3 lc rgb "black" nohead
set arrow from first om2_1,graph 0 to first om2_1,graph 1 dt 3 lc rgb "red"   nohead
set arrow from first om2_1,graph 0 to first om2_1,graph 1 dt 3 lc rgb "black" nohead

set arrow from first 1.8839,graph 0 to first 1.8839,graph 1 lw 2 lc rgb "#000000" nohead front
set arrow from first 1.8572,graph 0 to first 1.8572,graph 1 lw 2 lc rgb "#777777" nohead front
set arrow from first 1.8679,graph 0 to first 1.8679,graph 1 lw 2 lc rgb "#CCCCCC" nohead front

#set label at first 3,first 28 "{\\color{green!50!black} $\\times 4$}"
#set label at first 2.1,first 140 "{\\color[HTML]{FFC500} $\\times 1/4$}"

set x2tics -1000,1000,1000

omp=1.9
q1=0.2
q2=0.5
q3=1

f1=2.*q1**2/3./omp**2
f2=2.*q2**2/3./omp**2
f3=2.*q3**2/3./omp**2

fm=2./3.


plot rep1 u 3:(f1*$6)     w l lw 2 lc rgb "#000000"  title "$q\\xi=0.2$",\
     rep2 u 3:(f2*$6)     w l lw 2 lc rgb "#777777"  title "$q\\xi=0.5$",\
     rep3 u 3:(f3*$6)     w l lw 2 lc rgb "#CCCCCC"  title "$q\\xi=1.0$",\
