
#set terminal epslatex size 6.0, 2.5
#set output "plots.tex"
set terminal post size 6.0, 2.5
set output "awgn.eps"
set samples 10000
set autoscale			
set xtic auto   
set ytic auto   
set xlabel "$\\rho$ [dB]"
set ylabel "Bit-error-rate"
set xr [0:2]
set yr [0.00001:0.9]
set log y
set xtics 0,0.2,2
set key default
set grid 
plot \
"awgn_ru1000" using 1:3 title "RU-1000" with lines lt 1 lw 3 lc rgbcolor "blue", \
"awgn_ru10000" using 1:3 title "RU-10000" with lines lt 2 lw 3 lc rgbcolor "red", \
"awgn_ru100000" using 1:3 title "RU-100000" with lines lt 3 lw 3 lc rgbcolor "green"
