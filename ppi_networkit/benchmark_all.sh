#!/bin/bash

# setting database path
SQLITE_DB_FILE=$HOME/dev/bio/data/test_matching.sqlite


# change into build dir
cd ./build

# compile without OpenMP
echo "Compiling serial code"
cmake ../ -DBENCHMARK_OMP=OFF
make

algo=bw

echo "Running $algo serial benchmark"
./benchmark -a $algo $SQLITE_DB_FILE | tee data/benchmark_${algo}_serially.csv


# compile with OpenMP
echo "Compiling OpenMP code"
cmake ../ -DBENCHMARK_OMP=ON
make

echo "Running $algo parallel benchmark with STATIC schedule"
export OMP_SCHEDULE="static"
./benchmark -a $algo $SQLITE_DB_FILE | tee data/benchmark_${algo}_omp_static.csv

echo "Running $algo parallel benchmark with DYNAMIC schedule"
export OMP_SCHEDULE="dynamic"
./benchmark -a $algo $SQLITE_DB_FILE | tee data/benchmark_${algo}_omp_dynamic.csv

for algo in cc plp
do
    echo "Running $algo benchmark"
    ./benchmark -a $algo $SQLITE_DB_FILE | tee data/benchmark_${algo}.csv
done
