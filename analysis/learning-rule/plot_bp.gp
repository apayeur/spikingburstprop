set terminal aqua font 'Helvetica, 14' size 600, 400

p './d0.1/trc.0.brate' u 1:($2/$3) w l,\
'./d0.2/trc.0.brate' u 1:($2/$3) w l,\
'./d0.3/trc.0.brate' u 1:($2/$3) w l;

