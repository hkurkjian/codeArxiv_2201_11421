set terminal epslatex color size 12cm,10cm lw 2
set encoding utf8
set key on

system  "rm minplasm.tex"
set output 'minplasm.tex'

mp    ="minplasm.dat"

set xlabel  "\\Large $\\omega_p/\\Delta$"
set ylabel  "\\Large $\\omega_{q_{\\rm min}}/\\Delta$"
set y2label "\\Large ${\\color{blue} {q_{\\rm min}}\\xi}$" offset -8

set rmargin 7

set ytics nomirror

#set yrange [-0.4:0.0]
set xrange [1.6:2.9]
set yrange  [1.6:2.2]
set y2range [0.0:2.5]

set y2tics auto
set arrow from first 1.6955, graph 0 to first 1.6955,graph 1 dt 2 lc rgb "black" nohead

set xtics add ("1.70" 1.6955)

#set format y2 "$%h$}"
set format y2 "\\color{blue} %h"

set arrow from first 2.7,  first 2.05 to first 2.9, first 2.05 lw 2 lc rgb "blue"
set arrow from first 1.72, first 1.75 to first 1.6, first 1.75 lw 2 lc rgb "black"

plot mp u 1:3    w lines lw 2      lc rgb "black" notitle,\
     mp u 1:1    w lines lw 1 dt 3 lc rgb "black" notitle,\
     mp u 1:(2)  w lines lw 1 dt 3 lc rgb "black" notitle,\
     mp u 1:2   axis x1y2 w lines lw 2 lc rgb "blue" notitle,\

system  "rm minplasmq.tex"
set output 'minplasmq.tex'

set ylabel "\\Large $k_Fq_{\\rm min}/k_\\Delta^2$"

set yrange [0.0:2.5]
plot mp u 1:2            w lines lw 2 lc rgb "blue" notitle,\
