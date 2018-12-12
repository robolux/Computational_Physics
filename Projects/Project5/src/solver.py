# Project 5 solvers - Computational Physics
# Hunter Phillips

import numpy as nmp
import math
import random
import os
import sys
import time as rolex

# Inputs for function
# N  = Number of Agents
# tc = transaction count
# ag = agents array
# lambdaa = inputted c value (can't use just lambda because it is specified in python)
# alpha   = inputted c value
# gamma   = inputted

nmp.seterr(divide='ignore', invalid='ignore')

def trading(N, transaction_count, agents, lambdaa, alpha, gamma):

    startt = nmp.zeros((N,N), dtype=int)

    avg_tempo_old       = float(1e30000) # simulating infinity
    tempo_block_old     = float(1e30000) # ""
    cum_1               = float(0)
    cum_2               = float(0)
    block_length        = int(1e4)

    for i in range(0, transaction_count):

        ag_i = random.randint(0, N-1)
        ag_j = random.randint(0, N-1)

        m_i = agents[ag_i]
        m_j = agents[ag_j]

        rn = random.random() # random number
        previous = startt[ag_i, ag_j]

        varz1 = nmp.power(nmp.absolute(m_i - m_j), -alpha)*nmp.power(previous+1, gamma)
        while ((varz1 < rn) | (ag_i == ag_j)):

            ag_i = random.randint(0, N-1)
            ag_j = random.randint(0, N-1)

            m_i = agents[ag_i]
            m_j = agents[ag_j]

            previous = startt[ag_i, ag_j]

            rn = random.random()

        epsilon = random.random()

        delta_m = (1 - lambdaa) * (epsilon * m_j - (1 - epsilon) * m_i)
        agents[ag_i] = agents[ag_i] + delta_m
        agents[ag_j] = agents[ag_j] - delta_m

        startt[ag_i, ag_j] += 1
        startt[ag_j, ag_i] += 1

        tempo = agents
        cum_1 += tempo
        cum_2 += tempo*tempo

        if (i % block_length == 0):

            avg_tempo  = cum_1 / block_length
            avg_tempo2 = cum_2 / block_length

            tempo_block = avg_tempo2 - avg_tempo*avg_tempo

            testcond1 = nmp.absolute(avg_tempo_old - avg_tempo) / nmp.absolute(avg_tempo_old)
            testcond2 = nmp.absolute(tempo_block - tempo_block_old) / nmp.absolute(tempo_block_old)

            if ((testcond1.any() < 0.2) & (testcond2.any() < 0.5)):

                break

            else:

                avg_tempo_old = avg_tempo
                tempo_block_old = tempo_block

            cum_1  = 0
            cum_2 = 0

    return agents

# wrap up the interface for testing and UI
def wrapper(m_init, N, tran_count, sim_count, lambdaa, alpha, gamma):
    agents = nmp.zeros(N)
    totagents = nmp.zeros(N)

    random.seed(a=None) # start random set

    start = rolex.time()

    for i in range(0, sim_count):
        print('progress' + str(i))
        agents.fill(m_init)
        agents = trading(N, tran_count, agents, lambdaa, alpha, gamma)
        totagents += nmp.sort(agents)

    agents = totagents/sim_count
    end = rolex.time()
    tot_time = end - start

    return agents, tot_time
