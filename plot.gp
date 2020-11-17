cd '/home/antonio/IC/modelagem_acustica_bidimensional/data'
stat "centro1.dat" name "A"
set xlabel "x"
set ylabel "z"
set grid
set xrange [A_min_x:A_max_x]                                              
set yrange [A_min_y:A_max_y]
set cbrange [-0.5:0.5]
set palette rgb 33,13,10
set output '/home/antonio/IC/modelagem_acustica_bidimensional/Plots/cerjan.gif'  
set term gif animate delay 6
do for [i=0:39] {plot 'data'.i.'.dat' with image}
unset output