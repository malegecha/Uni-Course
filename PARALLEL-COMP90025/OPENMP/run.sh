g++ -Wall -g -fopenmp -o run utility.cpp utility.h main.cpp
export OMP_NUM_THREADS=1
time ./run
m
export OMP_NUM_THREADS=2
time ./run

export OMP_NUM_THREADS=3
time ./run

export OMP_NUM_THREADS=4
time ./run


export OMP_NUM_THREADS=8
time ./run