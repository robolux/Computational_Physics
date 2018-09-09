# Project 1 Plotting - Computational Physics
# Hunter Phillips

import sys
import os
import matplotlib.pyplot as plt
from main import part_b, part_c, part_e
import numpy as nmp
from solver import u_function

save_path = '../results/plots'
n_values = ['10', '100', '1000', '10000', '100000', '1000000']
n_values_lu = ['10', '100', '1000', '10000']

# pretty plots make the world go round
rainbow = ['#FC2847', '#FF8243', '#FFCF48', '#B2EC5D', '#1CD3A2', '#1DACD6']
rainbow2 = ['#FC2847', '#FF8243', '#FFCF48', '#B2EC5D', '#1CD3A2', '#1DACD6']
rainbow3 = ['#FC2847', '#FF8243', '#FFCF48', '#B2EC5D', '#1CD3A2', '#1DACD6']

# get part b solution plots
for g, rainbow in zip(n_values, rainbow):                       # enumerate across the two lists
        b_x, b_v = part_b(int(g), -1, 2, -1)                    # pass the test case to function
        plt.plot(b_x, b_v, color = rainbow, label = str(g))     # plot with legend / color

# get exact solution
n_exact = 1000000
x_exact = nmp.linspace(0, 1, n_exact)
u_exact = u_function(x_exact)
plt.plot(x_exact, u_exact, color = '#C364C5', label = 'Exact')


plt.legend(loc='upper right')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Approximation to u by General Algorithm')
filename = os.path.join(save_path,"general.png")
plt.savefig(filename)   # save figure to .png in ../results/plots
plt.clf()               # clear figure to prepare for next algo plot

# get part c solution plots
for g, rainbow2 in zip(n_values, rainbow2):
        c_x, c_v = part_c(int(g), -1, 2, -1)
        plt.plot(c_x, c_v, color = rainbow2, label = str(g))

# plot exact solution
plt.plot(x_exact, u_exact, color = '#C364C5', label = 'Exact')


plt.legend(loc='upper right')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Approximation to u by Specialized Tridiagonal Algorithm')
filename = os.path.join(save_path,"special.png")
plt.savefig(filename)
plt.clf()

# get part e solution plots
for g, rainbow3 in zip(n_values_lu, rainbow3):
        e_x, e_v = part_c(int(g), -1, 2, -1)
        plt.plot(e_x, e_v, color = rainbow3, label = str(g))

# plot exact solution
plt.plot(x_exact, u_exact, color = '#C364C5', label = 'Exact')


plt.legend(loc='upper right')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.title('Approximation to u by LU Decomposition')
filename = os.path.join(save_path,"lu.png")
plt.savefig(filename)
