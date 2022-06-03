# Simon_Clustering
This repository provides the code for the paper "Improved Differential and Linear Cryptanalysis on Round-Reduced SIMON"

Computes a lower boud on the probability of a linear approximation for Simon or Simeck

1. Change the variants parameters and window size at line 31~35; then change the corresponding input output difference/mask; 

2.Change the corresponding dynamic window.

3.Compile with gcc -Wall -Wextra -O3 -march=native -fopenmp -lm
