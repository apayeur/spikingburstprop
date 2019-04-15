set terminal aqua size  600, 400 font 'Helvetica, 14'

set xlabel 'Time [s]'
set ylabel 'Rate and BP'
p './test-output/burstpoisson.0.brate_seed1' u 1:3 w l lc 'blue' t 'ER',\
 '' u 1:(100*$2/$3) w l lc 'red' t 'BP',\
'' u 1:2 w l lc 'orange' t 'BR';
