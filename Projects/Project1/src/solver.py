# Project 1 Solvers - Computational Physics
# Hunter Phillips

import numpy as nmp
import scipy.linalg as scp
import time

def u_function(x):
    return 1 - (1 - nmp.exp(-10))*x - nmp.exp(-10*x) # analytical solution

def f_function(x):
    return 100*nmp.exp(-10*x) # source term

def general(f_function, a, b, c, n):

    # solve Using Gaussian Jordan Elimination (rref)
    x = nmp.linspace(0, 1, n+2) # set x range with points between
    h = x[1] - x[0]             # find height
    s = f_function(x)*h**2      # find s
    v = nmp.zeros(n+2)          # preallocation for speed

    # forward substitution
    for m in range(1, n+2):
        rf   =  a[m]/b[m-1]     # row factor
        b[m] -= c[m-1]*rf       # gauss magic
        s[m] -= s[m-1]*rf       # take it to the front now

    # backward substitution
    v[n+1] = s[n+1] / b[n+1]    # time to back it up

    for m in range(n, -1, -1):
        v[m] = (s[m] - c[m]*v[m+1]) / b[m]

    return x,v                  # lets make sure the variables don't get lost forever

def tridiag(f_function, a, b, c, n):

    # solve using tri diagonal method special case
    x = nmp.linspace(0, 1, n+2)
    h = x[1] - x[0]
    s = f_function(x)*h**2
    v = nmp.zeros(n+2)

    b_func = lambda i: (i+1)/i # create anonymous function
    b = nmp.zeros(n+2, dtype=nmp.float64)
    i = nmp.arange(0, n+2, dtype=nmp.float64)
    b[1:-1] = b_func(i[1:-1])

    # forward substitution
    for m in range(2, n+1):
        s[m] += (s[m-1])/(b[m-1])

    # backward substitution
    v[n] = (s[n])/(b[n])

    for m in range(n-1,0,-1):
        v[m] = (s[m]+v[m+1])/b[m]

    return x, v

def LU(f_function, a, b, c, n):

    # Solve Using LU
    A = nmp.zeros(shape = (n,n))
    A[range(1,n), range(n-1)] = a # fill matrix with a,b, and c values
    nmp.fill_diagonal(A, b)
    A[range(n-1), range(1, n)] = c


    lu, pivot = scp.lu_factor(A, overwrite_a=True, check_finite=False)
    # I would be interested in the future to try a no-pivot, partial-pivot, and full-pivot
    # lu decomp scheme to compare how it optimizes for this case.

    x = nmp.linspace(0, 1, n+2)
    h = x[1]-x[0]
    s = f_function(x)[1:-1]*h**2

    # solve the lu equation system
    v_i = scp.lu_solve((lu, pivot), s, overwrite_b=True)

    v = nmp.zeros(n+2)
    v[1:-1] = v_i[:] # strip non-important data

    return x, v
