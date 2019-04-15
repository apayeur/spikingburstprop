mkdir -p test-output
make
for i in 1 2 3 4 5
do
   ./test_burstpoisson --seed $i
done

python plot.py


