# Project 1 Testing - Computational Physics
# Hunter Phillips

import sys
import os.path
import numpy as nmp
from main import part_b, part_c, part_d, part_e, benchmarking

save_path = '../results'

# part b
for w, g in enumerate(['10', '100', '1000', '10000', '100000', '1000000'], 1):
    filename = os.path.join(save_path, "general_" + str(g) + ".txt")
    print(filename)
    with open (filename, 'w') as f_m:           # with auto closes files which is helpful here
        b_x, b_v = part_b(int(g), -1, 2, -1)
        s_b_x = ', '.join(str(h) for h in b_x)  # convert arrays to a string for output
        f_m.write(s_b_x)                        # write to a file the data
        f_m.write("\n\n")
        s_b_v = ', '.join(str(h) for h in b_v)
        f_m.write(s_b_v)

# part c
for w, g in enumerate(['10', '100', '1000', '10000', '100000', '1000000'], 1):
    filename = os.path.join(save_path, "special_" + str(g) + ".txt")
    print(filename)
    with open (filename, 'w') as f_m:
        c_x, c_v = part_c(int(g), -1, 2, -1)
        s_c_x = ', '.join(str(h) for h in c_x)
        f_m.write(s_c_x)
        f_m.write("\n\n")
        s_c_v = ', '.join(str(h) for h in c_v)
        f_m.write(s_c_v)

# part d
for w, g in enumerate(['10', '100', '1000', '10000', '100000', '1000000'], 1):
    filename = os.path.join(save_path, "error_" + str(g) + ".txt")
    print(filename)
    with open (filename, 'w') as f_m:
        n, eps = part_d(int(g), -1, 2, -1)
        f_m.write("%f" % n)
        f_m.write("\n\n")
        f_m.write("%f" % eps)

# part e
for w, g in enumerate(['10', '100', '1000', '10000'], 1):
    filename = os.path.join(save_path, "lu_" + str(g) + ".txt")
    print(filename)
    with open (filename, 'w') as f_m:
        e_x, e_v = part_e(int(g), -1, 2, -1)
        s_e_x = ', '.join(str(h) for h in e_x)
        f_m.write(s_e_x)
        f_m.write("\n\n")
        s_e_v = ', '.join(str(h) for h in e_v)
        f_m.write(s_e_v)

# benchmarking
for w, g in enumerate(['10', '100', '1000', '10000', '100000', '1000000'], 1):
    filename = os.path.join(save_path, "bench_" + str(g) + ".txt")
    print(filename)
    with open (filename, 'w') as f_m:
        g_bench, t_bench, lu_bench = benchmarking(int(g), -1, 2, -1)
        f_m.write("%f" % g_bench)
        f_m.write("\n\n")
        f_m.write("%f" % t_bench)
        f_m.write("\n\n")
        f_m.write("%f" % lu_bench)
