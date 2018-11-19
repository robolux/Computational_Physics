# Project 4 solvers - Computational Physics
# Hunter Phillips

import numpy as nmp
import math
import random
import os
import sys
import time as rolex

class ising_model:
    def __init__(self,inscript=False,initial_file='noname'):

        if str(initial_file) is not 'noname':
            self.L,self.T_i,self.T_f,self.Tn,self.MCc = nmp.genfromtxt(str(initial_file))
            self.L = int(self.L)
            self.MCc = int(self.MCc)
        else:
            if inscript is not False:
                self.L,self.T_i,self.T_f,self.Tn,self.MCc = inscript
                self.L = int(self.L)
                self.MCc = int(self.MCc)
            else:
                print 'There are no parameters given\nAborting.'
                sys.exit(1)

    def init_state(self,random_i=True):

        a = nmp.array([-1.,1.])
        L = self.L
        if random_i is True:
            spin_mat = nmp.random.choice(a,size=(L+2,L+2))
        else:
            spin_mat = nmp.ones((L+2,L+2))

        # setting periodic boundary conditions
        spin_mat[:,0] = spin_mat[:,L]
        spin_mat[:,L+1] = spin_mat[:,1]
        spin_mat[0,:] = spin_mat[L,:]
        spin_mat[L+1,:] = spin_mat[1,:]

        # setting initial energy and magnetic
        for i in xrange(1,L+1):
            for j in xrange(1,L+1):
                self.E -= spin_mat[i,j]*(spin_mat[i+1,j]+spin_mat[i,j+1])
                self.M += spin_mat[i,j]
        self.spin_mat = spin_mat

    def energy(self):

        E_st = 0
        for i in xrange(1,int(self.L)+1):
            for j in xrange(1,int(self.L)+1):
                E_st -= self.spin_mat[i,j]*(self.spin_mat[i+1,j]+self.spin_mat[i,j+1])
        self.E_count.append(E_st)

    def flip(self,i,j):

        L = self.L
        self.spin_mat[i,j] *= -1.

        # update PBC's depending on edge state
        if i == L:
            self.spin_mat[0,j] = self.spin_mat[i,j]
        elif i == 1:
            self.spin_mat[L+1,j] = self.spin_mat[i,j]
        if j == L:
            self.spin_mat[i,0] = self.spin_mat[i,j]
        elif j == 1:
            self.spin_mat[i,L+1] = self.spin_mat[i,j]


    def metropolis(self,w):

        L = int(self.L)
        for iterable in xrange(L*L):
            i,j = nmp.random.randint(1,L+1,2)
            dE_t = 2.*self.spin_mat[i,j]*(self.spin_mat[i,j+1] + \
                self.spin_mat[i,j-1] + self.spin_mat[i+1,j] + self.spin_mat[i-1,j])
            if dE_t <= 0:
                self.flip(i,j)
                self.E += dE_t
                self.M += 2.*self.spin_mat[i,j]
                self.accepted += 1
            else:
                if nmp.random.random() <= w[8+int(dE_t)]:
                    self.flip(i,j)
                    self.E += dE_t
                    self.M += 2.*self.spin_mat[i,j]
                    self.accepted += 1

    def solve(self,random_i=True,store_values=False,filename='noname',when_save=1):

        MCc = self.MCc; T_i = self.T_i
        T_f = self.T_f; Tn = self.Tn
        if T_i == T_f:
            T_ = [T_i]
        else:
            T_ = nmp.arange(T_i,T_f+1E-10,Tn)
        start = rolex.time()

        for T in T_:
            self.E_count = []
            self.accepted=0
            nmp.random.RandomState()
            print "Running MC with T = %g, lattice: %gX%g\n"%(T,self.L,self.L)
            # Create initial state
            self.E = 0; self.M = 0
            self.init_state(random_i)
            avg = nmp.zeros(5)
            beta = 1./T
            w = nmp.zeros(17)
            dE = nmp.array([-8.,-4.,0,4.,8.])
            for i in range(5):
                w[int(dE[i])+8] = nmp.exp(-beta*dE[i])

            cycle = 0

            if store_values is False:
                while cycle < MCc:
                    self.metropolis(w)
                    avg[0] += self.E; avg[1] += self.E**2
                    avg[2] += self.M; avg[3] += self.M**2
                    avg[4] += abs(self.M)
                    cycle += 1
                    self.energy()
                self.average(avg,T,cycle)

            else:
                if len(T_) == 1:
                    while cycle < MCc:
                        self.metropolis(w)
                        avg[0] += self.E; avg[1] += self.E**2
                        avg[2] += self.M; avg[3] += self.M**2
                        avg[4] += abs(self.M)
                        cycle += 1
                        if cycle % when_save == 0:
                            self.average(avg,T,cycle,filename)
                        self.accepted = 0
                else:
                    while cycle < MCc:
                        self.metropolis(w)
                        avg[0] += self.E; avg[1] += self.E**2
                        avg[2] += self.M; avg[3] += self.M**2
                        avg[4] += abs(self.M)
                        cycle += 1
                    self.average(avg,T,cycle,filename)

        end = rolex.time()
        print '\nTotal time: %g s'%(end-start)
        return self.E_count,self.E_var


    def average(self,avg,T,cycle,filename='noname'):

        Mc =  float(cycle) #
        T = float(T)
        L_tot = float(self.L*self.L)
        norm = 1./Mc
        E_bar, E2_bar, M_bar, M2_bar, M_bar_abs = avg*norm
        E_var = (E2_bar - E_bar*E_bar)/L_tot # = (<E^2> - <E>^2)/numberspin
        M_abs_var = (M2_bar - M_bar_abs*M_bar_abs)/L_tot # = <M^2> - <|M|>^2 = <M>^2 - <|M|>^2
        Cv = E_var/(T*T)             # Heat capacity per spin
        Chi = M_abs_var/T            # Susceptibility per spin
        M_bar_abs *= 1./L_tot        # <|M|> per spin
        E_bar *= 1./L_tot            # <E> per spin
        self.E_var = E_var
        self.Mabsbar = M_bar_abs
        self.ebar = E_bar

        # DEFAULT OUTPUTS
        with open (filename, 'a') as f_m:
            f_m.write('%10.11f' % E_bar)
            f_m.write(' ')
            f_m.write('%10.11f' % E_var)
            f_m.write(' ')
            f_m.write('%10.11f' % M_bar_abs)
            f_m.write(' ')
            f_m.write('%10.11f' % M_abs_var)
            f_m.write(' ')
            f_m.write('%10.11f' % T)
            f_m.write(' ')
            f_m.write('%10.11f' % Mc)
            f_m.write(' ')
            f_m.write('%10.11f' % (self.accepted/L_tot))
            f_m.write('\n')

        # FOR PART B
        # with open (filename, 'a') as f_m:
        #     f_m.write('%10.11f' % E_bar)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % Cv)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % M_bar_abs)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % Chi)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % T)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % Mc)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % (self.accepted/L_tot))
        #     f_m.write('\n')

        # FOR PART E
        # with open (filename, 'a') as f_m:
        #     f_m.write('%10.11f' % E_bar)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % M_bar_abs)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % Cv)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % Chi)
        #     f_m.write(' ')
        #     f_m.write('%10.11f' % T)
        #     f_m.write('\n')
