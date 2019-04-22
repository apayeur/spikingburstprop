mkdir -p test-output
make
BASE_CURRENT=100.e-12
for c in 0 100.e-12 200.e-12 
do
   ./test_ficurve_burstpoisson --somacurrent $c
done

python plot_ficurves.py
open ./test-output/BPversusIdend.pdf

