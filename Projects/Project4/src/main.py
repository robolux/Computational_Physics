# Project 2 Main - Computational Physics
# Hunter Phillips

from __future__ import print_function
import sys
import numpy as nmp
from solver import *
import matplotlib.pyplot as plt

if __name__ == "__main__":

    print('\nThis program allows the user to input:')
    print('' + '*'.__mul__(18))
    print('T (Temperature)\nL (L x L Lattice)\n# of Monte Carlo Cycles\nSpin Orientation')
    print('' + '*'.__mul__(18))
    print('\nThis is to solve the Ising Model using the')
    print('Metropolis algorithm coupled inside a Monte Carlo Simulation.\n')


    T_i = float(raw_input('What is your desired temperature?: '))
    T_f = T_i
    Tn  = 1
    bx = True

    L   = int(raw_input('Please input L for the L x L size of your desired lattice: '))
    MCc = float(raw_input('Please input the desired number of Monte Carlo Cycles (can be in scientific form ex. 1e5): '))

    bz = False
    while bz == False:
        adv = raw_input('Would you like to have random (1) or ordered (2) spin orientations? (1 or 2): ')
        adv.replace("(", "") # cleanup user input
        adv.replace(")", "")

        if (int(adv) == 1):
            rii = True
            bz = True

        elif (int(adv) == 2):
            rii = False
            bz = True

        else:
            print('User Input not recognized as option, please try again.')

    fn_i = raw_input('Please input the name of your desired output file (with .txt extension): ')
    fn_i.replace(" ", "")
    save_path = '../results/usercase/'
    file_name  = save_path + fn_i

    ins = nmp.array([L,T_i,T_f,Tn,MCc])
    test = ising_model(inscript = ins)
    test.solve(random_i=rii,store_values=True,filename=file_name,when_save=1)

    udata = nmp.genfromtxt(save_path + fn_i)
    E_bar = udata[:,0]
    M_bar_abs = udata[:,2]
    MCc = udata[:,-2]
    accepted_per = udata[:,-1]

    plt.plot(MCc, E_bar)
    plt.ylabel('Mean Energy')
    plt.xlabel('Cycles')
    plt.title('Cycles vs Mean Energy')
    plt.tight_layout()
    plt.show()
    plt.clf()
    plt.plot(MCc, M_bar_abs, 'r')
    plt.ylabel('Mean Magnetization')
    plt.xlabel('Cycles')
    plt.title('Cycles vs Mean Magnetization')
    plt.show()

    print('\n\n' + '*'.__mul__(49))
    print('All results located in: ' + save_path + fn_i + '\n\n')
