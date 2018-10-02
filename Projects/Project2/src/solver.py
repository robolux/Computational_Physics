# Project 2 Solvers - Computational Physics
# Hunter Phillips

import numpy as nmp
import scipy.integrate as scp
import time
import numpy as np
import math
import time
import os
import sys
from operator import itemgetter

def find_max(A):
    n = A.shape[0]
    maximum = abs(A[0,1])
    max_k=0
    max_l=1
    for i in xrange(n):
        for j in xrange(i+1, n):
            if abs(A[i,j]) > maximum:
                maximum = abs(A[i,j])
                max_k = i
                max_l = j
    return maximum, max_k, max_l

def solve(A, R, tol):

    n = A.shape[0]
    iterations = 0
    maximum, k, l = find_max(A)
    while maximum > tol:
        iterations += 1
        rotate(A, R, k, l)
        maximum, k, l = find_max(A)

    return iterations, tol

def rotate(A, R, k, l):
    n = A.shape[0]
    tau = (A[l,l] - A[k,k])/(2*A[k,l])
    if tau > 0:
        t = 1./(tau + math.sqrt(1 + tau**2))
    else:
        t = 1./(tau - math.sqrt(1 + tau**2))
    c = 1 / math.sqrt(1+t**2)
    s = c*t

    a_kk = A[k,k]
    a_ll = A[l,l]

    A[k,k] = c**2*a_kk - 2*c*s*A[k,l] + s**2*a_ll
    A[l,l] = s**2*a_kk + 2*c*s*A[k,l] + c**2*a_ll
    A[k,l] = 0
    A[l,k] = 0

    for i in xrange(n):
        if i != k and i != l:
            a_ik = A[i,k]
            a_il = A[i,l]
            A[i,k] = c*a_ik - s*a_il
            A[k,i] = A[i,k]
            A[i,l] = c*a_il + s*a_ik
            A[l,i] = A[i,l]

        r_ik = R[i,k]
        r_il = R[i,l]
        R[i,k] = c*r_ik - s*r_il
        R[i,l] = c*r_il + s*r_ik
def make_matrix_case_b(n, omega, p_max):

    A = np.zeros(shape=(n,n), dtype=np.float64)

    rho_0 = 0
    rho_n = p_max
    rho = np.linspace(rho_0, rho_n, n+2)[:-1]   # quickfix
    h = rho[1]-rho[0]
    V = np.zeros(n+1)
    V[1:] = omega**2*rho[1:]**2
    d = 2/h**2 + V
    e = -1/h**2

    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = e
    A[range(n-1), range(1, n)] = e
    return A, rho

def make_matrix_case_e(n, omega, p_max):

    A = np.zeros(shape=(n,n), dtype=np.float64)

    rho_0 = 0
    rho_n = p_max
    rho = np.linspace(rho_0, rho_n, n+2)[:-1]
    h = rho[1]-rho[0]
    V = np.zeros(n+1)
    V[1:] = omega**2*rho[1:]**2 + 1/rho[1:]
    d = 2/h**2 + V
    e = -1/h**2

    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = e
    A[range(n-1), range(1, n)] = e
    return A, rho

def make_matrix_case_d(n, omega, p_max):

    A = np.zeros(shape=(n,n), dtype=np.float64)

    rho_0 = 0
    rho_n = p_max
    rho = np.linspace(rho_0, rho_n, n+2)[:-1]
    h = rho[1]-rho[0]
    V = np.zeros(n+1)
    V[1:] = omega**2*rho[1:]**2 + rho[1:]**2
    d = 2/h**2 + V
    e = -1/h**2

    A[range(n), range(n)] = d[1:]
    A[range(1, n), range(n-1)] = e
    A[range(n-1), range(1, n)] = e
    return A, rho

def norm(eigenvectors):

    for i in xrange(len(eigenvectors)):
        if sum(eigenvectors[i]) < 0:
            eigenvectors[i] = -eigenvectors[i]
    return eigenvectors

def extract(A, R, sort_R=True):
    n = A.shape[0]
    eigvals = A[range(n), range(n)]
    eigvecs = norm([R[:,i].copy() for i in xrange(n)])
    eigs = [(eigval, eigvec, i) for i, (eigval, eigvec) in enumerate(zip(eigvals, eigvecs))]
    sorted_eigs = sorted(eigs, key=itemgetter(0))
    if sort_R:
        for i in xrange(n):
            R[:,i] = sorted_eigs[i][1]
    return sorted_eigs

# n = matrix size
# interacting = 1 or 0 (basically int bool)
# omega = prescribed omega value
# p_max = rho max
def solve_low(n, interacting, omega, p_max, tol):

    if interacting == 0:
        A, rho = make_matrix_case_b(n, omega, p_max)
    elif interacting == 1:
        A, rho = make_matrix_case_e(n, omega, p_max)
    elif interacting == 2:
        A, rho = make_matrix_case_d(n, omega, p_max)
    R = nmp.eye(n)

    eig_compare = nmp.linalg.eig(A)
    eig_compare2 = nmp.linalg.eigh(A)

    # solve
    start_time = time.time()
    iterations, tol = solve(A, R, tol)
    end_time = time.time()
    general_time = end_time - start_time

    # extract eigenvalues and eigenvectors
    sorted_eigs = extract(A, R)

    return iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2
