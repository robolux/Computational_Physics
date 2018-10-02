# Project 2 Main - Computational Physics
# Hunter Phillips

from __future__ import print_function
import sys
import numpy as nmp
from solver import *

# Print iteration output
def it_output(tol, iterations, n):
   print("Reached specified tolerance of (%.00E) in %g iterations" % (tol, iterations))
   print("This translates to %.2f iterations/element" % (iterations/float(n**2)))
   return


if __name__ == "__main__":

    print('\nThis program allows the user to input n (for n x n matrix), omega, rho max, and convergance tolerance to solve Eigenvalue problems,')
    print('from the equations of a buckling beam to Schroedinger\'s equation')
    print('for two electrons in a three- dimensional harmonic oscillator well.\n')

    # get user input
    n_i = int(raw_input('Please input the size n of your desired matrix: '))
    o_i = int(raw_input('Please input the value of omega for your desired matrix: '))
    p_i = int(raw_input('Please input the value of rho max for your desired matrix: '))
    t_i = float(raw_input('Please input the desired tolerance of the solution for your matrix in scientific notation (ex. 1e-8): '))
    # it is so convienent that float can take an input with e in it

    ### Part C ###
    print('\n\n' + '*'.__mul__(16))
    print('Part C')
    print('The two unit tests required for part C are the comparison between')
    print('The solution, numpy .eig, and numpy .eigh in the other parts of the assignment')
    print('*'.__mul__(16) + '\n')

    ### Part B ###
    print('\n\n' + '*'.__mul__(16))
    print('Part B')
    print('*'.__mul__(16) + '\n')

    # non interesting case thus intersecting = 0
    iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2 = solve_low(n_i, 0, o_i, p_i, t_i)

    print("Solution for Part B with n = %d, omega = %d, and rho max = %d took %g seconds"
    % (n_i, o_i, p_i, general_time))

    it_output(tol, iterations, n_i)

    print('\nPart B Sol    numpy .eig    numpy .eigh\n*****************************************')

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
            print('%.12s' % ('%.10f' % n), end="  ")
            counter = 1
        elif counter == 1:
            print('%.12s' % ('%.10f' % n), end="  ")
            counter = 2
        elif counter == 2:
            print('%.12s' % ('%.10f' % n))
            counter = 0

    Pass_bool = False
    for h in range(0,nmp.size(eig_num)/3):
        if nmp.round(eig_num[0], 8) == nmp.round(eig_num[1], 8) == nmp.round(eig_num[2], 8):
            Pass_bool = True
        else:
            Pass_bool = False
        del eig_num[0:3]

    if Pass_bool == True:
        print('The Unit Tests have PASSED')
    elif Pass_bool == False:
        print('The Unit Tests have FAILED')


    ### Part D ###

    print('\n\n' + '*'.__mul__(16))
    print('Part D')
    print('*'.__mul__(16) + '\n')

    iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2 = solve_low(n_i, 2, o_i, p_i, t_i)

    print("Solution for Part D with n = %d, omega = %d, and rho max = %d took %g seconds"
    % (n_i, o_i, p_i, general_time))

    it_output(tol, iterations, n_i)

    print('\nPart D Sol    numpy .eig    numpy .eigh\n*****************************************')

    eig_num = [1.0]
    for m in range(0, n_i):
        eig_num.append(sorted_eigs[m][0])
        eig_num.append(eig_compare[0][m])
        eig_num.append(eig_compare2[0][m])
    eig_num.pop(0)


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
            print('%.12s' % ('%.10f' % n), end="  ")
            counter = 1
        elif counter == 1:
            print('%.12s' % ('%.10f' % n), end="  ")
            counter = 2
        elif counter == 2:
            print('%.12s' % ('%.10f' % n))
            counter = 0

    Pass_bool = False
    for h in range(0,nmp.size(eig_num)/3):
        if nmp.round(eig_num[0], 8) == nmp.round(eig_num[1], 8) == nmp.round(eig_num[2], 8):
            Pass_bool = True
        else:
            Pass_bool = False
        del eig_num[0:3]

    if Pass_bool == True:
        print('The Unit Tests have PASSED')
    elif Pass_bool == False:
        print('The Unit Tests have FAILED')

    ### Part E ###
    print('\n\n' + '*'.__mul__(16))
    print('Part E')
    print('Uses omega_r = 0.01, 0.5, 1, and 5 for the ground state only')
    print('*'.__mul__(16) + '\n')


    omega_array = ['0.01', '0.5', '1', '5']
    for w, g in enumerate(omega_array, 1):

        iterations, tol, sorted_eigs, eig_compare, A, general_time, eig_compare2 = solve_low(n_i, 1, float(g), p_i, t_i)

        print('Solution for Part E Sol w' + str(w) + ' with n = %d, omega = %.2f, and rho max = %d took %g seconds'
        % (n_i, float(g), p_i, general_time))

        it_output(tol, iterations, n_i)

        print('\nPart E Sol w' + str(w) + '   numpy .eig      numpy .eigh\n********************************************')

        eig_num = [1.0]
        for m in range(0, n_i):
            eig_num.append(sorted_eigs[m][0])
            eig_num.append(eig_compare[0][m])
            eig_num.append(eig_compare2[0][m])
        eig_num.pop(0)

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
                print('%.12s' % ('%.10f' % n), end="    ")
                counter = 1
            elif counter == 1:
                print('%.12s' % ('%.10f' % n), end="    ")
                counter = 2
            elif counter == 2:
                print('%.12s' % ('%.10f' % n))
                counter = 0

        Pass_bool = False
        for h in range(0,nmp.size(eig_num)/3):
            if nmp.round(eig_num[0], 8) == nmp.round(eig_num[1], 8) == nmp.round(eig_num[2], 8):
                Pass_bool = True
            else:
                Pass_bool = False
            del eig_num[0:3]

        if Pass_bool == True:
            print('The Unit Tests have PASSED')
        elif Pass_bool == False:
            print('The Unit Tests have FAILED')

        print('\n')
