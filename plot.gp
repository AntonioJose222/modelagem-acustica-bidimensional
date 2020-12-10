cd '/home/antonio/IC/modelagem_acustica_bidimensional/data'
stat "data0.dat" name "A"
set xlabel "x"
set ylabel "z"
set grid
set xrange [A_min_x:A_max_x]                                              
set yrange [A_min_y:A_max_y]
set cbrange [-3:3]
set palette rgb 21,22,23
set output '/home/antonio/IC/modelagem_acustica_bidimensional/Plots/cerjan.gif'  
set term gif animate delay 6
do for [i=0:79] {plot 'data'.i.'.dat' with image title sprintf("t = %1.5f", i*0.00025)}
unset output
cd '/home/antonio/IC/modelagem_acustica_bidimensional/gnuplotScript'