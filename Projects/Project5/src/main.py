# Project 5 Main - Computational Physics
# Hunter Phillips

import sys
import numpy as nmp
from solver import *

if __name__ == "__main__":

    print('\nThis program simulates financial transactions among')
    print('financial agents using Monte Carlo methods.\n')

    m_init      = float(raw_input('Please input the desired initial amount of money (m0): '))
    N           = int(raw_input('Please input the desired number of agents (N): '))
    tran_count  = int(raw_input('Please input the desired number of transactions: '))
    sim_count   = int(raw_input('Please input the desired number of Monte Carlo simulations: '))
    lambdaa     = float(raw_input('Please input the desired value for savings rate (lambda): '))
    alpha       = float(raw_input('Please input the desired value for neighbour interactions (alpha): '))
    gamma       = float(raw_input('Please input the desired value for probability to trade with agents\nyou have traded with before (gamma): '))

    print('\nProgram is running...\n')

    agents, tot_time = wrapper(m_init, N, tran_count, sim_count, lambdaa, alpha, gamma)

    print("Total Computational Time: " + str(tot_time) + '\n')

    nmp.savetxt('user_run.txt', agents, delimiter='\n')
