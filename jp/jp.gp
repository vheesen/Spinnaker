
reset
set key bottom left font ",20" box width -6
#set terminal epslatex color default
set terminal postscript eps enhanced 'Helvetica' 24 lw 1
set output 'jp.eps'
set tics scale 1
set xlabel '{{/Symbol n}/{/Symbol n}_{brk}}'
set ylabel '{/Helevetica-Bold Non-Thermal Intensity (I_{nt})}'

set logscale x
set logscale y
set format x "10^{%L}"
set format y "10^{%L}"

plot [1e-5:100][1e-4:1000] 'spec.dat' using ($1/20.7e9):($2*60) lt 1 lw 1 lc rgb "red" smooth csplines notitle, 'spec.dat' using ($1/20.7e9):($4*60) lt 1 lw 1 lc rgb "blue" smooth csplines notitle
