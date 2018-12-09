# Project 5 Testing - Computational Physics
# Hunter Phillips

import sys
import os.path
import numpy as nmp
from solver import *

save_path = '../results'

def unit_tests():

    for w, unit_lambda in enumerate(['0', '0.25', '0.5', '0.9'], 1):

        for w2, unit_alpha in enumerate(['0', '0.5', '1.0', '1.5', '2.0'], 1):

            for w3, unit_N in enumerate(['500', '1000'], 1):

                for w3, unit_gamma in enumerate(['0.0', '1.0', '2.0', '3.0', '4.0'], 1):

                    file_name = save_path + '/N_' + unit_N + '_L_' + unit_lambda + '_a_' + unit_alpha + '_g_' + unit_gamma + '.txt'
                    comp_time = save_path + '/timing/N_' + unit_N + '_L_' + unit_lambda + '_a_' + unit_alpha + '_g_' + unit_gamma + '_timing.txt'

                    print(file_name)

                    m_init       =  float(1000)
                    N            =  int(unit_N)
                    tran_count   =  int(1e4)
                    sim_count    =  int(1e3)
                    lambdaa      =  float(unit_lambda)
                    alpha        =  float(unit_alpha)
                    gamma        =  float(unit_gamma)

                    agents, tot_time = wrapper(m_init, N, tran_count, sim_count, lambdaa, alpha, gamma)

                    nmp.savetxt(file_name, agents, delimiter='\n')

                    with open (comp_time, 'w') as f_m:
                        f_m.write('Computational Time: ' + str(tot_time))

unit_tests()
