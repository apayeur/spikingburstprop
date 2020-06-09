# Create directories for data and results if they don't exist already
DIR="../supp-data/normalization"
mkdir -p $DIR

# Compile
make

./sim_normalization_of_bpcurves --dir $DIR --somacurrent 0. --w_pyr_to_pyr 0.02
./sim_normalization_of_bpcurves --dir $DIR --somacurrent 100.e-12 --w_pyr_to_pyr 0.02
./sim_normalization_of_bpcurves --dir $DIR --somacurrent 200.e-12 --w_pyr_to_pyr 0.02

# Plot
python ../supp-analysis/ficurve.py

