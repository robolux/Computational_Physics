# Project 4 Testing - Computational Physics
# Hunter Phillips

import sys
import os.path
import numpy as nmp
import matplotlib.pyplot as plt
from solver import *

save_path = '../results'
run_txt_production_b   = False # computer freeze protection
run_txt_production_c_1 = False
run_txt_production_c_2 = False
run_txt_production_d   = False
run_txt_production_e   = False

def part_b():

    L = 2
    T_i = 1
    T_f = 1
    Tn = 1

    if run_txt_production_b == True:
        for w, g in enumerate(['1e3', '1e4', '1e5', '2e5', '1e6'], 1):

            file_name = save_path + '/part_b__L_2__mc__' + g + '.txt'
            with open (file_name, 'w') as f_m:
                print("")

            MCc = float(g)
            ins = nmp.array([L,T_i,T_f,Tn,MCc])
            test = ising_model(inscript = ins)
            test.solve(random_i=True,store_values=True,filename=file_name,when_save=1)

    data_b_1 = nmp.genfromtxt(save_path + '/part_b__L_2__mc__2e5.txt')

    E_bar = data_b_1[:,0]
    M_bar_abs = data_b_1[:,2]
    MCc = data_b_1[:,-2]
    accepted_per = data_b_1[:,-1]

    plt.plot(MCc, E_bar)
    plt.ylabel('Mean Energy')
    plt.xlabel('Cycles')
    plt.tight_layout()
    plt.savefig(save_path + '/part_b_Mean_Energy_Random.png', dpi = 1200)
    plt.clf()
    plt.plot(MCc, M_bar_abs, 'r')
    plt.ylabel('Mean Magnetization')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_b_Mean_Magnetization_Random.png', dpi = 1200)
    plt.clf()

def part_c_T_1():
    L = 20
    T_i = 1
    T_f = 1
    Tn = 1

    if run_txt_production_c_1 == True:
        for w, g in enumerate(['1e4'], 1):  # reduced this to 1e4, but didn't change looping
                                            # for future expandibility

            file_name  = save_path + '/part_c__random__L_2__mc__' + g + '.txt'
            file_name2 = save_path + '/part_c__ordered__L_2__mc__' + g + '.txt'
            with open (file_name, 'w') as f_m: # clear old files for appending
                print("")
            with open (file_name2, 'w') as f_m:
                print("")

            MCc = float(g)
            ins = nmp.array([L,T_i,T_f,Tn,MCc])
            test = ising_model(inscript = ins)
            test.solve(random_i=True,store_values=True,filename=file_name,when_save=1)

            test2 = ising_model(inscript = ins)
            test2.solve(random_i=False,store_values=True,filename=file_name2,when_save=1)

    data_c_1 = nmp.genfromtxt(save_path + '/part_c__random__L_2__mc__1e4.txt')
    data_c_2 = nmp.genfromtxt(save_path + '/part_c__ordered__L_2__mc__1e4.txt')

    E_bar = data_c_1[:,0]
    M_bar_abs = data_c_1[:,2]
    MCc = data_c_1[:,-2]
    accepted_per = data_c_1[:,-1]

    E_bar2 = data_c_2[:,0]
    M_bar_abs2 = data_c_2[:,2]
    MCc2 = data_c_2[:,-2]
    accepted_per2 = data_c_2[:,-1]

    # Plot Random 20x20 over 1e4
    plt.plot(MCc, E_bar)
    plt.ylabel('Mean Energy')
    plt.xlabel('Cycles')
    plt.tight_layout()
    plt.savefig(save_path + '/part_c_Mean_Energy_Random.png', dpi = 1200)
    plt.clf()
    plt.plot(MCc, M_bar_abs, 'r')
    plt.ylabel('Mean Magnetization')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_Mean_Magnetization_Random.png', dpi = 1200)
    plt.clf()

    # Plot Ordered 20x20 over 1e4
    plt.plot(MCc2, E_bar2)
    plt.ylabel('Mean Energy')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_Mean_Energy_Ordered.png', dpi = 1200)
    plt.clf()
    plt.plot(MCc2, M_bar_abs2, 'r')
    plt.ylabel('Mean Magnetization')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_Mean_Magnetization_Ordered.png', dpi = 1200)
    plt.clf()

    # plot accepted configs for random
    plt.plot(MCc, accepted_per)
    plt.ylabel('Accepted configurations')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_accepted_configs_Random.png', dpi = 1200)
    plt.clf()

    # plot accepted configs for ordered
    plt.plot(MCc2, accepted_per2)
    plt.ylabel('Accepted configurations')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_accepted_configs_Ordered.png', dpi = 1200)
    plt.clf()

def part_c_T_2():
    L = 20
    T_i = 2.4
    T_f = 2.4
    Tn = 1

    if run_txt_production_c_2 == True:
        for w, g in enumerate(['1e4'], 1):  # reduced this to xxx, but didn't change looping
                                            # for future expandibility

            file_name  = save_path + '/part_c__random__T2__L_2__mc__' + g + '.txt'
            file_name2 = save_path + '/part_c__ordered__T2__L_2__mc__' + g + '.txt'
            with open (file_name, 'w') as f_m: # clear old files for appending
                print("")
            with open (file_name2, 'w') as f_m:
                print("")

            MCc = float(g)
            ins = nmp.array([L,T_i,T_f,Tn,MCc])
            test = ising_model(inscript = ins)
            test.solve(random_i=True,store_values=True,filename=file_name,when_save=1)

            test2 = ising_model(inscript = ins)
            test2.solve(random_i=False,store_values=True,filename=file_name2,when_save=1)

    data_c_1 = nmp.genfromtxt(save_path + '/part_c__random__T2__L_2__mc__1e4.txt')
    data_c_2 = nmp.genfromtxt(save_path + '/part_c__ordered__T2__L_2__mc__1e4.txt')

    E_bar = data_c_1[:,0]
    M_bar_abs = data_c_1[:,2]
    MCc = data_c_1[:,-2]
    accepted_per = data_c_1[:,-1]

    E_bar2 = data_c_2[:,0]
    M_bar_abs2 = data_c_2[:,2]
    MCc2 = data_c_2[:,-2]
    accepted_per2 = data_c_2[:,-1]

    # Plot Random 20x20 over 1e4
    plt.plot(MCc, E_bar)
    plt.ylabel('Mean Energy')
    plt.xlabel('Cycles')
    plt.tight_layout()
    plt.savefig(save_path + '/part_c_Mean_Energy_Random_T2.png', dpi = 1200)
    plt.clf()
    plt.plot(MCc, M_bar_abs, 'r')
    plt.ylabel('Mean Magnetization')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_Mean_Magnetization_Random_T2.png', dpi = 1200)
    plt.clf()

    # Plot Ordered 20x20 over 1e4
    plt.plot(MCc2, E_bar2)
    plt.ylabel('Mean Energy')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_Mean_Energy_Ordered_T2.png', dpi = 1200)
    plt.clf()
    plt.plot(MCc2, M_bar_abs2, 'r')
    plt.ylabel('Mean Magnetization')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_Mean_Magnetization_Ordered_T2.png', dpi = 1200)
    plt.clf()

    # plot accepted configs for random
    plt.plot(MCc, accepted_per)
    plt.ylabel('Accepted configurations')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_accepted_configs_Random_T2.png', dpi = 1200)
    plt.clf()

    # plot accepted configs for ordered
    plt.plot(MCc2, accepted_per2)
    plt.ylabel('Accepted configurations')
    plt.xlabel('Cycles')
    plt.savefig(save_path + '/part_c_accepted_configs_Ordered_T2.png', dpi = 1200)
    plt.clf()



# get probability
def part_d():

    MCc = 1E4
    L = 20
    T = 1.0
    T2 = 2.4

    in_1_T1 = nmp.array([L,T,T,1,MCc])
    in_1_T2 = nmp.array([L,T2,T2,1,MCc])
    in_2_T1 = nmp.array([L,T,T,1,MCc])
    in_2_T2 = nmp.array([L,T2,T2,1,MCc])

    filename_T1_rand = save_path + '/part_d__random__L_2__mc__1e4.txt'
    filename_T2_rand = save_path + '/part_d__random__T2__L_2__mc__1e4.txt'
    filename_T1_ord =  save_path + '/part_d__ordered__L_2__mc__1e4.txt'
    filename_T2_ord =  save_path + '/part_d__ordered__T2__L_2__mc__1e4.txt'

    if run_txt_production_d == True:
        test = ising_model(inscript = in_1_T1) # Temp 1 Random
        E_count, E_var = test.solve(random_i=True)
        vals = nmp.append(E_var,E_count)
        nmp.savetxt(filename_T1_rand,nmp.c_[vals])

        test2 = ising_model(inscript = in_1_T2) # Temp 2 Random
        E_count2, E_var2 = test2.solve(random_i=True)
        vals2 = nmp.append(E_var2,E_count2)
        nmp.savetxt(filename_T2_rand,nmp.c_[vals2])

        test3 = ising_model(inscript = in_2_T1) # Temp 1 Ordered
        E_count3, E_var3 = test3.solve(random_i=False)
        vals3 = nmp.append(E_var3,E_count3)
        nmp.savetxt(filename_T1_ord,nmp.c_[vals3])

        test4 = ising_model(inscript = in_2_T2) # Temp 2 Ordered
        E_count4, E_var4 = test4.solve(random_i=False)
        vals4 = nmp.append(E_var4,E_count4)
        nmp.savetxt(filename_T2_ord,nmp.c_[vals4])

    vals1 = nmp.genfromtxt(filename_T1_rand)
    E_var1 = vals1[0]
    E_count1 = vals1[1:]

    vals2 = nmp.genfromtxt(filename_T2_rand)
    E_var2 = vals2[0]
    E_count2 = vals2[1:]

    vals3 = nmp.genfromtxt(filename_T1_ord)
    E_var3 = vals3[0]
    E_count3 = vals3[1:]

    vals4 = nmp.genfromtxt(filename_T2_ord)
    E_var4 = vals4[0]
    E_count4 = vals4[1:]

    ss = 2500 # When approximately steady state occurs

    E_count1 = nmp.array(E_count1[ss:])/400.
    plt.hist(E_count1,weights=nmp.ones_like(E_count1)/float(len(E_count1)),bins=25)
    plt.ylabel('$\mathrm{PDF}$',size=15)
    plt.xlabel('$E,\,[\mathrm{J}]$',size=15)
    print nmp.sqrt(E_var1)/20.  # E_var is given per spin, so in order to find variance, sqrt(400)/400 = 1/20
    plt.savefig(save_path + '/part_d_T1_random.png', dpi = 1200)
    plt.clf()

    E_count2 = nmp.array(E_count2[ss:])/400.
    plt.hist(E_count2,weights=nmp.ones_like(E_count2)/float(len(E_count2)),bins=25)
    plt.ylabel('$\mathrm{PDF}$',size=15)
    plt.xlabel('$E,\,[\mathrm{J}]$',size=15)
    print nmp.sqrt(E_var2)/20.  # E_var is given per spin, so in order to find variance, sqrt(400)/400 = 1/20
    plt.savefig(save_path + '/part_d_T2_random.png', dpi = 1200)
    plt.clf()

    E_count3 = nmp.array(E_count3[ss:])/400.
    plt.hist(E_count3,weights=nmp.ones_like(E_count3)/float(len(E_count3)),bins=25)
    plt.ylabel('$\mathrm{PDF}$',size=15)
    plt.xlabel('$E,\,[\mathrm{J}]$',size=15)
    print nmp.sqrt(E_var3)/20.  # E_var is given per spin, so in order to find variance, sqrt(400)/400 = 1/20
    plt.savefig(save_path + '/part_d_T1_ord.png', dpi = 1200)
    plt.clf()

    E_count4 = nmp.array(E_count4[ss:])/400.
    plt.hist(E_count4,weights=nmp.ones_like(E_count4)/float(len(E_count4)),bins=25)
    plt.ylabel('$\mathrm{PDF}$',size=15)
    plt.xlabel('$E,\,[\mathrm{J}]$',size=15)
    print nmp.sqrt(E_var4)/20.  # E_var is given per spin, so in order to find variance, sqrt(400)/400 = 1/20
    plt.savefig(save_path + '/part_d_T2_ord.png', dpi = 1200)
    plt.clf()


def part_e(): # Writes data for L=40,60,80,100 over temp steps of 0.05 in the interval [2.0,2.3]
              # when running this uncomment
              # PART E in solver.py function average

    if run_txt_production_e == True:
        T_i = 2.0; T_f = 2.3; Tn = 0.05
        MCc = 1E4
        L_ = nmp.array([40,60,80,100])
        for L in L_:
            file_name = save_path + '/part_e__L_%s.txt'%(str(L))
            ins = nmp.array([L,T_i,T_f,Tn,MCc])
            test = ising_model(inscript = ins)
            test.solve(store_values=True,filename=file_name)
            del test

    data  = nmp.genfromtxt(save_path + '/part_e__L_40.txt')
    data2 = nmp.genfromtxt(save_path + '/part_e__L_60.txt')
    data3 = nmp.genfromtxt(save_path + '/part_e__L_80.txt')
    data4 = nmp.genfromtxt(save_path + '/part_e__L_100.txt')

    T = data[:,-1]

    # Energies
    E_bar = data[:,0]
    E_bar2 = data2[:,0]
    E_bar3 = data3[:,0]
    E_bar4 = data4[:,0]

    plt.plot(T,E_bar,label='$L=40$')
    plt.plot(T,E_bar2,label='$L=60$')
    plt.plot(T,E_bar3,label='$L=80$')
    plt.plot(T,E_bar4,label='$L=100$')

    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$\\langle E \\rangle,\,[1]$',size=15)
    plt.savefig(save_path + '/part_e_Energy.png', dpi = 1200)
    plt.clf()

    # |M|
    M_bar_abs = data[:,1]
    M_bar_abs2 = data2[:,1]
    M_bar_abs3 = data3[:,1]
    M_bar_abs4 = data4[:,1]
    plt.plot(T,M_bar_abs,label='$L=40$')
    plt.plot(T,M_bar_abs2,label='$L=60$')
    plt.plot(T,M_bar_abs3,label='$L=80$')
    plt.plot(T,M_bar_abs4,label='$L=100$')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$\\langle |M| \\rangle,\,[1]$',size=15)
    plt.savefig(save_path + '/part_e_Mbar.png', dpi = 1200)
    plt.clf()

    # Specific heat capacity; Cv
    Cv = data[:,2]
    Cv2 = data2[:,2]
    Cv3 = data3[:,2]
    Cv4 = data4[:,2]
    plt.plot(T,Cv,label='$L=40$')
    plt.plot(T,Cv2,label='$L=60$')
    plt.plot(T,Cv3,label='$L=80$')
    plt.plot(T,Cv4,label='$L=100$')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$C_v,\,[1]$',size=15)
    plt.savefig(save_path + '/part_e_Cv.png', dpi = 1200)
    plt.clf()

    # Magnetic susceptibility; chi
    chi = data[:,3]
    chi2 = data2[:,3]
    chi3 = data3[:,3]
    chi4 = data4[:,3]
    plt.plot(T,chi,label='$L=40$')
    plt.plot(T,chi2,label='$L=60$')
    plt.plot(T,chi3,label='$L=80$')
    plt.plot(T,chi4,label='$L=100$')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('$T,\,[\mathrm{J}]$',size=15)
    plt.ylabel('$\chi,\,[1/J]$',size=15)
    plt.savefig(save_path + '/part_e_chi.png', dpi = 1200)
    plt.clf()

# Find Tc
def part_f():
    data  = nmp.genfromtxt(save_path + '/part_e__L_40.txt')
    data2 = nmp.genfromtxt(save_path + '/part_e__L_60.txt')
    data3 = nmp.genfromtxt(save_path + '/part_e__L_80.txt')
    data4 = nmp.genfromtxt(save_path + '/part_e__L_100.txt')

    T = data[:,-1]

    M_bar_abs = data[:,1]
    M_bar_abs2 = data2[:,1]
    M_bar_abs3 = data3[:,1]
    M_bar_abs4 = data4[:,1]

    Cv_ = data[:,2]
    Cv2 = data2[:,2]
    Cv3 = data3[:,2]
    Cv4 = data4[:,2]
    beta = 1./8

    Cv = nmp.array([Cv_,Cv2,Cv3,Cv4])
    Mbar = nmp.array([M_bar_abs,M_bar_abs2,M_bar_abs3,M_bar_abs4])
    L = nmp.array([40.0,60.0,80.0,100.0])
    Tc = nmp.zeros(len(L))
    o = 5
    a = 0
    for i in range(len(L)):
        Tc[i] = T[nmp.argmax(Cv[i,o:]) + o]
        ind = nmp.argmax(Cv[i,o:]) + o
        a += L[i]*(Mbar[i,ind]**(1./beta))
    a /= float(len(L))
    print(a)
    Tc_L_lim = 0.
    for i in range(len(L)):
        Tc_L_lim += Tc[i] - a/L[i]
    Tc_L_lim /= float(len(L))
    print('Tc = ' + str(Tc_L_lim))

if __name__ == "__main__":
    # Run All of the Project Parts
    part_b()
    part_c_T_1()
    part_c_T_2()
    part_d()
    part_e()
    part_f()
