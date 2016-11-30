
set terminal epslatex size 4.5, 3
set output "myth-5-figure.tex"
#set terminal post size 6.0, 2.5
#set output "mami.eps"
set samples 10000
set autoscale	
set xtic auto   
set ytic auto   
set xlabel "$\\rho$ [dB]"
set ylabel "Bit-error-rate" offset 2
set xr [-15:-11]
set yr [0.00005:0.9]
set log y
#set xtics -15,0.2,-11
set key default
set grid 
set arrow from -13.94,0.00005 to -13.94,0.9 nohead lt 1 lw 5 lc rgbcolor "red"
plot \
"mami_ru1000" using 1:3 title "1000 bits codewords" with lines lt 1 lw 3 lc rgbcolor "blue", \
"mami_ru10000" using 1:3 title "10000 bits codewords" with lines lt 2 lw 3 lc rgbcolor "red", \
"mami_ru100000" using 1:3 title "100000 bits codewords" with lines lt 3 lw 3 lc rgbcolor "green"
