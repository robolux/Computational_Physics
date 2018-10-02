# Project 1 Testing - Computational Physics
# Hunter Phillips

import sys
import os.path
import numpy as nmp
from main import *

save_path = '../results'

# part b
case = 1
for w, g in enumerate(['1', '5'], 1):
    for w2, g2 in enumerate(['4', '16', '64', '128'], 1):
        for w3, g3 in enumerate(['2', '6', '12'], 1):
            save_path = '../results/partb'
            filename = os.path.join(save_path, "case" + str(case) + ".txt")
            print(filename)
            case = case + 1
            with open (filename, 'w') as f_m:           # with auto closes files which is helpful here
                n_i = int(g2)
                o_i = float(g)
                p_i = int(g3)
                t_i = 1e-10
                iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2 = solve_low(n_i, 0, o_i, p_i, t_i)
                f_m.write('\n\n' + '*'.__mul__(16)+'\n')
                f_m.write('Part B'+'\n')
                f_m.write('*'.__mul__(16) + '\n'+'\n')
                f_m.write("Solution for Part B with n = %d, omega = %d, and rho max = %d took %g seconds\n" % (n_i, o_i, p_i, general_time))
                f_m.write("Reached specified tolerance of (%.00E) in %g iterations\n" % (tol, iterations))
                f_m.write("This translates to %.2f iterations/element\n" % (iterations/float(n_i**2)))
                f_m.write('\nPart B Sol    numpy .eig    numpy .eigh\n*****************************************\n')

                eig_num = [1.0]
                for m in range(0, n_i):
                    eig_num.append(sorted_eigs[m][0])
                    eig_num.append(eig_compare[0][m])
                    eig_num.append(eig_compare2[0][m])
                eig_num.pop(0) # could restart index instead, quick fix for right now

                # this checks if eig has put the values in decending order
                # and fixes the ordering if this is the case
                # if the values are genuinally different without regard to order
                # the test will still fail

                rev_back = -2
                trash_pass = [1.0]
                rev_for = 1
                dummy2 = 0
                eig_num_dummy = nmp.copy(eig_num)
                for h in range(0, nmp.size(eig_num_dummy)/3):
                    dummy2 = eig_num_dummy[nmp.size(eig_num_dummy)+rev_back]
                    trash_pass.append(dummy2)
                    rev_back = rev_back - 3
                    rev_for = rev_for + 3
                trash_pass.pop(0)

                trash_sorted = nmp.sort(trash_pass)

                rev_back = -2
                rev_for = 1
                for h in range(0, nmp.size(eig_num_dummy)/3):
                    eig_num[rev_for] = trash_sorted[h]
                    rev_back = rev_back - 3
                    rev_for = rev_for + 3

                counter = 0
                for n in eig_num:
                    if counter == 0:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('  ')
                        counter = 1
                    elif counter == 1:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('  ')
                        counter = 2
                    elif counter == 2:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('\n')
                        counter = 0

                Pass_bool = False
                for h in range(0,nmp.size(eig_num)/3):
                    if nmp.round(eig_num[0], 8) == nmp.round(eig_num[1], 8) == nmp.round(eig_num[2], 8):
                        Pass_bool = True
                    else:
                        Pass_bool = False
                    del eig_num[0:3]

                if Pass_bool == True:
                    f_m.write('The Unit Tests have PASSED')
                elif Pass_bool == False:
                    f_m.write('The Unit Tests have FAILED')

# Part D
case = 1
for w, g in enumerate(['1', '3', '7', '11'], 1):
    for w2, g2 in enumerate(['4', '16', '64', '128'], 1):
        for w3, g3 in enumerate(['2', '6', '36', '128'], 1): # make pmax higher to determine stability
            save_path = '../results/partd'
            filename = os.path.join(save_path, "case" + str(case) + ".txt")
            print(filename)
            case = case + 1
            with open (filename, 'w') as f_m:           # with auto closes files which is helpful here
                n_i = int(g2)
                o_i = float(g)
                p_i = int(g3)
                t_i = 1e-10
                iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2 = solve_low(n_i, 2, o_i, p_i, t_i)
                f_m.write('\n\n' + '*'.__mul__(16)+'\n')
                f_m.write('Part D'+'\n')
                f_m.write('*'.__mul__(16) + '\n'+'\n')
                f_m.write("Solution for Part D with n = %d, omega = %d, and rho max = %d took %g seconds\n" % (n_i, o_i, p_i, general_time))
                f_m.write("Reached specified tolerance of (%.00E) in %g iterations\n" % (tol, iterations))
                f_m.write("This translates to %.2f iterations/element\n" % (iterations/float(n_i**2)))
                f_m.write('\nPart D Sol    numpy .eig    numpy .eigh\n*****************************************\n')

                eig_num = [1.0]
                for m in range(0, n_i):
                    eig_num.append(sorted_eigs[m][0])
                    eig_num.append(eig_compare[0][m])
                    eig_num.append(eig_compare2[0][m])
                eig_num.pop(0) # could restart index instead, quick fix for right now

                rev_back = -2
                trash_pass = [1.0]
                rev_for = 1
                dummy2 = 0
                eig_num_dummy = nmp.copy(eig_num)
                for h in range(0, nmp.size(eig_num_dummy)/3):
                    dummy2 = eig_num_dummy[nmp.size(eig_num_dummy)+rev_back]
                    trash_pass.append(dummy2)
                    rev_back = rev_back - 3
                    rev_for = rev_for + 3
                trash_pass.pop(0)

                trash_sorted = nmp.sort(trash_pass)

                rev_back = -2
                rev_for = 1
                for h in range(0, nmp.size(eig_num_dummy)/3):
                    eig_num[rev_for] = trash_sorted[h]
                    rev_back = rev_back - 3
                    rev_for = rev_for + 3

                counter = 0
                for n in eig_num:
                    if counter == 0:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('  ')
                        counter = 1
                    elif counter == 1:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('  ')
                        counter = 2
                    elif counter == 2:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('\n')
                        counter = 0

                Pass_bool = False
                for h in range(0,nmp.size(eig_num)/3):
                    if nmp.round(eig_num[0], 8) == nmp.round(eig_num[1], 8) == nmp.round(eig_num[2], 8):
                        Pass_bool = True
                    else:
                        Pass_bool = False
                    del eig_num[0:3]

                if Pass_bool == True:
                    f_m.write('The Unit Tests have PASSED')
                elif Pass_bool == False:
                    f_m.write('The Unit Tests have FAILED')

# Part E
case = 1
for w, g in enumerate(['0.01', '0.5', '1.0', '5.0'], 1):
    for w2, g2 in enumerate(['4', '16', '64', '128'], 1):
        for w3, g3 in enumerate(['2', '6', '36', '128'], 1): # make pmax higher to determine stability
            save_path = '../results/parte'
            filename = os.path.join(save_path, "case" + str(case) + ".txt")
            print(filename)
            case = case + 1
            with open (filename, 'w') as f_m:           # with auto closes files which is helpful here
                n_i = int(g2)
                o_i = float(g)
                p_i = int(g3)
                t_i = 1e-10
                iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2 = solve_low(n_i, 1, o_i, p_i, t_i)
                f_m.write('\n\n' + '*'.__mul__(16)+'\n')
                f_m.write('Part E'+'\n')
                f_m.write('*'.__mul__(16) + '\n'+'\n')
                f_m.write("Solution for Part E with n = %d, omega = %d, and rho max = %d took %g seconds\n" % (n_i, o_i, p_i, general_time))
                f_m.write("Reached specified tolerance of (%.00E) in %g iterations\n" % (tol, iterations))
                f_m.write("This translates to %.2f iterations/element\n" % (iterations/float(n_i**2)))
                f_m.write('\nPart E Sol w' + str(w) + '   numpy .eig      numpy .eigh\n********************************************')

                eig_num = [1.0]
                for m in range(0, n_i):
                    eig_num.append(sorted_eigs[m][0])
                    eig_num.append(eig_compare[0][m])
                    eig_num.append(eig_compare2[0][m])
                eig_num.pop(0) # could restart index instead, quick fix for right now

                rev_back = -2
                trash_pass = [1.0]
                rev_for = 1
                dummy2 = 0
                eig_num_dummy = nmp.copy(eig_num)
                for h in range(0, nmp.size(eig_num_dummy)/3):
                    dummy2 = eig_num_dummy[nmp.size(eig_num_dummy)+rev_back]
                    trash_pass.append(dummy2)
                    rev_back = rev_back - 3
                    rev_for = rev_for + 3
                trash_pass.pop(0)

                trash_sorted = nmp.sort(trash_pass)

                rev_back = -2
                rev_for = 1
                for h in range(0, nmp.size(eig_num_dummy)/3):
                    eig_num[rev_for] = trash_sorted[h]
                    rev_back = rev_back - 3
                    rev_for = rev_for + 3

                f_m.write('\n')

                counter = 0
                for n in eig_num:
                    if counter == 0:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('    ')
                        counter = 1
                    elif counter == 1:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('    ')
                        counter = 2
                    elif counter == 2:
                        f_m.write('%.12s' % ('%.10f' % n))
                        f_m.write('\n')
                        counter = 0

                Pass_bool = False
                for h in range(0,nmp.size(eig_num)/3):
                    if nmp.round(eig_num[0], 8) == nmp.round(eig_num[1], 8) == nmp.round(eig_num[2], 8):
                        Pass_bool = True
                    else:
                        Pass_bool = False
                    del eig_num[0:3]

                if Pass_bool == True:
                    f_m.write('The Unit Tests have PASSED')
                elif Pass_bool == False:
                    f_m.write('The Unit Tests have FAILED')

                f_m.write('\n\n')
