# Project 5 Plotting - Computational Physics
# Hunter Phillips

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from matplotlib import rc

# Plot 1 - Lambda 0, 0.25, 0.5, 0.9
plot1_1 = np.loadtxt("../results/N_500_L_0_a_0_g_0.0_1.txt")
plot1_2 = np.loadtxt("../results/N_500_L_0.25_a_0_g_0.0_1.txt")
plot1_3 = np.loadtxt("../results/N_500_L_0.5_a_0_g_0.0_1.txt")
plot1_4 = np.loadtxt("../results/N_500_L_0.9_a_0_g_0.0_1.txt")
length_1 = len(plot1_1)

sz = 0.125
selection_1=int(max(plot1_1)/sz)
selection_2=int(max(plot1_2)/sz)
selection_3=int(max(plot1_3)/sz)
selection_4=int(max(plot1_4)/sz)

dh1, edge = np.histogram(plot1_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.plot(cen, dh1/(float(length_1)))

dh2, edge = np.histogram(plot1_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.plot(cen, dh2/(float(length_1)))

dh3, edge = np.histogram(plot1_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.plot(cen, dh3/(float(length_1)))

dh4, edge = np.histogram(plot1_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.plot(cen, dh4/(float(length_1)))

plt.xlim([0, 3])
plt.show()
plt.clf()

# Plot 2 - log(wm) as function of m
plot2_1 = np.loadtxt("../results/N_500_L_0_a_0_g_0.0_1.txt")
plot2_2 = np.loadtxt("../results/N_500_L_0.25_a_0_g_0.0_1.txt")
plot2_3 = np.loadtxt("../results/N_500_L_0.5_a_0_g_0.0_1.txt")
plot2_4 = np.loadtxt("../results/N_500_L_0.9_a_0_g_0.0_1.txt")

beta_1_2 = 1/np.mean(plot2_1)
beta_2_2 = 1/np.mean(plot2_2)
beta_3_2 = 1/np.mean(plot2_3)
beta_4_2 = 1/np.mean(plot2_4)
omega_1_2 = beta_1_2*np.exp(-beta_1_2 * plot2_1)
omega_2_2 = beta_2_2*np.exp(-beta_2_2 * plot2_2)
omega_3_2 = beta_3_2*np.exp(-beta_3_2 * plot2_3)
omega_4_2 = beta_4_2*np.exp(-beta_4_2 * plot2_4)

plt.semilogy(plot2_1,omega_1_2)
plt.semilogy(plot2_2,omega_2_2)
plt.semilogy(plot2_3,omega_3_2)
plt.semilogy(plot2_4,omega_4_2)
plt.show()
plt.clf()

# Plot 3 - Parameterization
plot3_1 = np.loadtxt("../results/N_500_L_0.9_a_0_g_0.0_1.txt")
length_3 = len(plot3_1)

lambda_3 = 0.9
x = plot3_1 / np.mean(plot3_1)
n = 1 + (3*lambda_3)/(1-lambda_3)
an = (n**n) / gamma(n)
P = an*(x**(n-1))*np.exp(-n*x)

plt.loglog(x,P,'r')

sz = 0.1
selection_1=int(max(plot3_1)/sz)
dh4, edge = np.histogram(plot3_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(2*selection_1)))
plt.show()
plt.clf()

# Plot 4 - Parameterization
plot4_1 = np.loadtxt("../results/N_500_L_0_a_0_g_0.0_1.txt")
length_4 = len(plot4_1)

lambda_3 = 0.0
x = plot4_1 / np.mean(plot4_1)
n = 1 + (3*lambda_3)/(1-lambda_3)
an = (n**n) / gamma(n)
P = an*(x**(n-1))*np.exp(-n*x)

plt.loglog(x,P,'b')

sz = 0.5
selection_1=int(max(plot4_1)/sz)
print(selection_1)
dh4, edge = np.histogram(plot4_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_4)),'r')
plt.show()
plt.clf()

# Plot 5 - Pareto Power Law Tails
plot2_1 = np.loadtxt("../results/N_500_L_0_a_0_g_0.0_1.txt")
m = np.linspace(min(plot2_1),max(plot2_1),len(plot2_1))

powerr = lambda alpha, m: m**-(1+alpha)

plt.loglog(m,powerr(2.0, m))
plt.loglog(m,powerr(1.8, m))
length_1 = len(plot2_1)
sz = 1
selection_1=int(max(plot2_1)/sz)
dh4, edge = np.histogram(plot2_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.plot(cen, dh4/(float(length_1)))
plt.xlim([0.5, 7])
plt.ylim([0.002, 15])
plt.show()
plt.clf()

# Plot 6 - Varying alpha with no savings N = 500
plot6_1 = np.loadtxt("../results/N_500_L_0_a_0_g_0.0_2.txt")
plot6_2 = np.loadtxt("../results/N_500_L_0_a_0.5_g_0.0_2.txt")
plot6_3 = np.loadtxt("../results/N_500_L_0_a_1.0_g_0.0_2.txt")
plot6_4 = np.loadtxt("../results/N_500_L_0_a_1.5_g_0.0_2.txt")
plot6_5 = np.loadtxt("../results/N_500_L_0_a_2.0_g_0.0_2.txt")
length_6 = len(plot6_1)

sz = 0.2
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot6_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_6)))

dh2, edge = np.histogram(plot6_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_6)))

dh3, edge = np.histogram(plot6_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_6)))

dh4, edge = np.histogram(plot6_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_6)))

dh5, edge = np.histogram(plot6_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_6)))

plt.show()
plt.clf()

# Plot 7 - Varying alpha with savings N = 500
plot7_1 = np.loadtxt("../results/N_500_L_0.9_a_0_g_0.0_2.txt")
plot7_2 = np.loadtxt("../results/N_500_L_0.9_a_0.5_g_0.0_2.txt")
plot7_3 = np.loadtxt("../results/N_500_L_0.9_a_1.0_g_0.0_2.txt")
plot7_4 = np.loadtxt("../results/N_500_L_0.9_a_1.5_g_0.0_2.txt")
plot7_5 = np.loadtxt("../results/N_500_L_0.9_a_2.0_g_0.0_2.txt")
length_7 = len(plot7_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot7_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_7)))

dh2, edge = np.histogram(plot7_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_7)))

dh3, edge = np.histogram(plot7_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_7)))

dh4, edge = np.histogram(plot7_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_7)))

dh5, edge = np.histogram(plot7_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_7)))

plt.show()
plt.clf()

# Plot 8 - Varying alpha without savings N = 1000
plot8_1 = np.loadtxt("../results/N_1000_L_0_a_0_g_0.0_2.txt")
plot8_2 = np.loadtxt("../results/N_1000_L_0_a_0.5_g_0.0_2.txt")
plot8_3 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_0.0_2.txt")
plot8_4 = np.loadtxt("../results/N_1000_L_0_a_1.5_g_0.0_2.txt")
plot8_5 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_0.0_2.txt")
length_8 = len(plot8_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot8_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_8)))

dh2, edge = np.histogram(plot8_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_8)))

dh3, edge = np.histogram(plot8_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_8)))

dh4, edge = np.histogram(plot8_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_8)))

dh5, edge = np.histogram(plot8_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_8)))

plt.show()
plt.clf()

# Plot 9 - Varying alpha with savings N = 1000
plot9_1 = np.loadtxt("../results/N_1000_L_0.9_a_0_g_0.0_2.txt")
plot9_2 = np.loadtxt("../results/N_1000_L_0.9_a_0.5_g_0.0_2.txt")
plot9_3 = np.loadtxt("../results/N_1000_L_0.9_a_1.0_g_0.0_2.txt")
plot9_4 = np.loadtxt("../results/N_1000_L_0.9_a_1.5_g_0.0_2.txt")
plot9_5 = np.loadtxt("../results/N_1000_L_0.9_a_2.0_g_0.0_2.txt")
length_9 = len(plot9_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot9_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_9)))

dh2, edge = np.histogram(plot9_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_9)))

dh3, edge = np.histogram(plot9_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_9)))

dh4, edge = np.histogram(plot9_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_9)))

dh5, edge = np.histogram(plot9_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_9)))

plt.show()
plt.clf()

# Plot 10 - Pareto Power Law Tails Varying Alpga
plot10_1 = np.loadtxt("../results/N_1000_L_0_a_0_g_0.0_2.txt")
plot10_2 = np.loadtxt("../results/N_1000_L_0_a_0.5_g_0.0_2.txt")
plot10_3 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_0.0_2.txt")
plot10_4 = np.loadtxt("../results/N_1000_L_0_a_1.5_g_0.0_2.txt")
plot10_5 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_0.0_2.txt")
length_10 = len(plot10_1)
m10_1 = np.linspace(min(plot10_2),max(plot10_2),len(plot10_2))
m10_2 = np.linspace(min(plot10_3),max(plot10_3),len(plot10_3))
m10_3 = np.linspace(min(plot10_4),max(plot10_4),len(plot10_4))
m10_4 = np.linspace(min(plot10_5),max(plot10_5),len(plot10_5))

powerr = lambda alpha, m: m**-(1+alpha)

plt.loglog(m10_1,powerr(0.5, m10_1), color='#C0C0C0')
plt.loglog(m10_2,powerr(1.0, m10_2), color='#808080')
plt.loglog(m10_3,powerr(1.5, m10_3), color='#696969')
plt.loglog(m10_4,powerr(2.0, m10_4), color='#000000')

sz = 2
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot10_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_10)))

dh2, edge = np.histogram(plot10_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_10)))

dh3, edge = np.histogram(plot10_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_10)))

dh4, edge = np.histogram(plot10_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_10)))

dh5, edge = np.histogram(plot10_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_10)))
plt.xlim([1, 10])
plt.ylim([0.001, 2])
plt.show()
plt.clf()

# Plot 11 - Varying gamma without savings N = 1000
plot11_1 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_0.0_3.txt")
plot11_2 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_1.0_3.txt")
plot11_3 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_2.0_3.txt")
plot11_4 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_3.0_3.txt")
plot11_5 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_4.0_3.txt")
length_11 = len(plot11_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot11_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_11)))

dh2, edge = np.histogram(plot11_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_11)))

dh3, edge = np.histogram(plot11_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_11)))

dh4, edge = np.histogram(plot11_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_11)))

dh5, edge = np.histogram(plot11_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_11)))

plt.show()
plt.clf()


# Plot 12 - Varying gamma with savings N = 1000
plot12_1 = np.loadtxt("../results/N_1000_L_0.9_a_1.0_g_0.0_3.txt")
plot12_2 = np.loadtxt("../results/N_1000_L_0.9_a_1.0_g_1.0_3.txt")
plot12_3 = np.loadtxt("../results/N_1000_L_0.9_a_1.0_g_2.0_3.txt")
plot12_4 = np.loadtxt("../results/N_1000_L_0.9_a_1.0_g_3.0_3.txt")
plot12_5 = np.loadtxt("../results/N_1000_L_0.9_a_1.0_g_4.0_3.txt")
length_12 = len(plot12_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot12_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_12)))

dh2, edge = np.histogram(plot12_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_12)))

dh3, edge = np.histogram(plot12_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_12)))

dh4, edge = np.histogram(plot12_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_12)))

dh5, edge = np.histogram(plot12_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_12)))

plt.show()
plt.clf()


# Plot 13 - Varying gamma without savings N = 1000
plot13_1 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_0.0_3.txt")
plot13_2 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_1.0_3.txt")
plot13_3 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_2.0_3.txt")
plot13_4 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_3.0_3.txt")
plot13_5 = np.loadtxt("../results/N_1000_L_0_a_2.0_g_4.0_3.txt")
length_13 = len(plot13_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot13_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_13)))

dh2, edge = np.histogram(plot13_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_12)))

dh3, edge = np.histogram(plot13_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_13)))

dh4, edge = np.histogram(plot13_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_13)))

dh5, edge = np.histogram(plot13_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_13)))

plt.show()
plt.clf()


# Plot 14 - Varying gamma with savings N = 1000
plot14_1 = np.loadtxt("../results/N_1000_L_0.9_a_2.0_g_0.0_3.txt")
plot14_2 = np.loadtxt("../results/N_1000_L_0.9_a_2.0_g_1.0_3.txt")
plot14_3 = np.loadtxt("../results/N_1000_L_0.9_a_2.0_g_2.0_3.txt")
plot14_4 = np.loadtxt("../results/N_1000_L_0.9_a_2.0_g_3.0_3.txt")
plot14_5 = np.loadtxt("../results/N_1000_L_0.9_a_2.0_g_4.0_3.txt")
length_14 = len(plot14_1)

sz = 0.5
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)
selection_3=int(max(plot6_3)/sz)
selection_4=int(max(plot6_4)/sz)
selection_5=int(max(plot6_5)/sz)

dh1, edge = np.histogram(plot14_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_14)))

dh2, edge = np.histogram(plot14_2,bins=selection_2)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh2/(float(length_14)))

dh3, edge = np.histogram(plot14_3,bins=selection_3)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh3/(float(length_14)))

dh4, edge = np.histogram(plot14_4,bins=selection_4)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh4/(float(length_14)))

dh5, edge = np.histogram(plot14_5,bins=selection_5)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh5/(float(length_14)))

plt.show()
plt.clf()


# Plot 15 - Pareto Power Law Tails Varying Alpga
plot15_1 = np.loadtxt("../results/N_1000_L_0_a_1.0_g_2.0_3.txt")
length_15 = len(plot15_1)
m15_1 = np.linspace(min(plot15_1),max(plot15_1),len(plot15_1))

powerr = lambda alpha, m: m**(-(1+alpha))

plt.loglog(m15_1,powerr(1.0, m15_1))

sz = 1
selection_1=int(max(plot6_1)/sz)
selection_2=int(max(plot6_2)/sz)

dh1, edge = np.histogram(plot15_1,bins=selection_1)
cen = 0.5*(edge[1:]+edge[:-1])
delta_b = cen[1]-cen[0]
plt.loglog(cen, dh1/(float(length_15)))
plt.xlim([1, 10])
plt.ylim([0.001, 2])

plt.show()
