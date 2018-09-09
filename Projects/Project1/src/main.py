# Project 1 Main - Computational Physics
# Hunter Phillips

# This program allows the user to input n, a, b, and c to solve
# Poisson's equation in a single dimension by turning it into a matrix.
# using Gaussian Jordan Elimination, Specialized Tridiagonal, and LU Decomposition
# Benchmarks and error calculations are also performed.

import sys
import numpy as nmp
from solver import *

# part b
def part_b(n, a_i, b_i, c_i):
    a = nmp.full(n+2, a_i, dtype=nmp.float64)
    b = nmp.full(n+2, b_i, dtype=nmp.float64)
    c = nmp.full(n+2, c_i, dtype=nmp.float64)
    x, v = general(f_function, a, b, c, n)
    return x, v

# part c
def part_c(n, a_i, b_i, c_i):
    a = nmp.full(n+2, a_i, dtype=nmp.float64)
    b = nmp.full(n+2, b_i, dtype=nmp.float64)
    c = nmp.full(n+2, c_i, dtype=nmp.float64)
    x, v = tridiag(f_function, a, b, c, n)
    return x, v

# part d
def part_d(n, a_i, b_i, c_i):
    a = nmp.full(n+2, a_i, dtype=nmp.float64)
    b = nmp.full(n+2, b_i, dtype=nmp.float64)
    c = nmp.full(n+2, c_i, dtype=nmp.float64)
    x, v = tridiag(f_function, a, b, c, n)
    u = u_function(x)
    eps_i = nmp.log10(abs((v[1:-1]-u[1:-1])/u[1:-1])) # strip the last and first parts of matrix
    eps = nmp.max(eps_i)
    return n, eps

# part e
def part_e(n, a_i, b_i, c_i):
    a = a_i
    b = b_i
    c = c_i
    x, v = LU(f_function, a, b, c, n)
    return  x, v

# benchmarks
def benchmarking(n, a_i, b_i, c_i):

    a = nmp.full(n+2, a_i, dtype=nmp.float64) # have to get that accuracy
    b = nmp.full(n+2, b_i, dtype=nmp.float64)
    c = nmp.full(n+2, c_i, dtype=nmp.float64)

    start_time = time.time() # in the future I would like to add more complex timing analysis
    x, v = general(f_function, a, b, c, n)
    end_time = time.time()
    general_time = end_time - start_time

    start_time = time.time()
    x, v = tridiag(f_function, a, b, c, n)
    end_time = time.time()
    tridiag_time = end_time - start_time

    a = a_i
    b = b_i
    c = c_i

    start_time = time.time()
    x, v = LU(f_function, a, b, c, n)
    end_time = time.time()
    LU_time = end_time - start_time

    return general_time, tridiag_time, LU_time


if __name__ == "__main__":

    print('This program allows the user to input n (for n x n matrix), a, b, and c to solve\nPoisson\'s equation in a single dimension by turning it into a matrix\nusing Gaussian Jordan Elimination, Specialized Tridiagonal, and LU Decomposition\nBenchmarks and error calculations are also performed.\n\n')

    # get user input
    n_i = int(raw_input('Please input the size n of your desired matrix: '))
    a_i = int(raw_input('Please input the value of a for your desired matrix: '))
    b_i = int(raw_input('Please input the value of b for your desired matrix: '))
    c_i = int(raw_input('Please input the value of c for your desired matrix: '))

    print('\nGeneral Algorithm')
    b_x, b_v = part_b(n_i, a_i, b_i, c_i)
    print(b_x)
    print(b_v)

    print('\nSpecialized Tridiagonal Algorithm')
    c_x, c_v = part_c(n_i, a_i, b_i, c_i)
    print(c_x)
    print(c_v)

    print('\nRelative Error')
    n, eps = part_d(n_i, a_i, b_i, c_i)
    print('n = ' + str(n))
    print('error = ' + str(eps))

    print('\nLinear Decomposition Algorithm')
    e_x, e_v = part_e(n_i, a_i, b_i, c_i)
    print(e_x)
    print(e_v)

    general_time, tridiag_time, LU_time = benchmarking(n_i, a_i, b_i, c_i)
    print('\nGeneral Algorithm Time to Execute')
    print(general_time)
    print('\nSpecialized Tridiagonal Algorithm Time to Execute')
    print(tridiag_time)
    print('\nLU Decomposition Time to Execute')
    print(LU_time)
    print('\n\n')
