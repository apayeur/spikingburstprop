set terminal aqua font 'Helvetica, 14' size 300, 200

p './d0.1/trc.0.wsum' u 1:2 w l,\
'./d0.2/trc.0.wsum' u 1:2 w l,\
'./d0.3/trc.0.wsum' u 1:2 w l;

