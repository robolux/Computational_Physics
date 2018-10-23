# Project 3 Testing - Computational Physics
# Hunter Phillips

from __future__ import print_function
import sys
import os.path
import time
import numpy as nmp
from solver import solar_system
import matplotlib.pyplot as plt
import numpy as np
import warnings
from operator import itemgetter
import math

# suppress matplotlib axis version error
# going to be corrected in next mpl release according to stackexchange
warnings.filterwarnings("ignore", module="matplotlib") # suppress it

# constants
G = 4*nmp.pi**2

if __name__ == "__main__":

    # the earth-sun system
    solar_system_a = solar_system()
    solar_system_a.create_planetary_object(0,0, 0,0, 1) #Sun
    solar_system_a.create_planetary_object(1,0, 0,2*nmp.pi, 3.003e-6) #Earth

    p, v = solar_system_a.fill_sol(int(1e6), 50,integrator=solar_system_a.forward_euler)
    plt.axes(aspect="equal")
    plt.plot(p[:,0,0], p[:,0,1], "yo")
    plt.plot(p[:,1,0], p[:,1,1], "b-")
    plt.plot(p[0,1,0], p[0,1,1], "bx")
    plt.plot(p[-1,1,0], p[-1,1,1], "bo")
    plt.axis([-1.5,1.5,-1.5,1.5])
    plt.title("Earth Sun system over 50 years using Forward Euler Method")
    plt.legend(["Sun","Earth"])
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.savefig("../results/forward_euler.png", dpi = 1200)
    plt.clf()

    p, v = solar_system_a.fill_sol(int(1e6), 50,integrator=solar_system_a.velocity_verlet)
    plt.axes(aspect="equal")
    plt.plot(p[:,0,0], p[:,0,1], "yo")
    plt.plot(p[:,1,0], p[:,1,1], "b-")
    plt.plot(p[0,1,0], p[0,1,1], "bx")
    plt.plot(p[-1,1,0], p[-1,1,1], "bo")
    plt.axis([-1.5,1.5,-1.5,1.5])
    plt.title("Earth Sun system over 50 years using Velocity Verlet Method")
    plt.legend(["Sun","Earth"])
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.savefig("../results/velocity_verlet.png", dpi = 1200)
    plt.clf()


    # time step tests
    time_solar_system = solar_system()
    time_solar_system.create_planetary_object(0, 0, 0, 0, 1)
    time_solar_system.create_planetary_object(1, 0, 0, 6.281989, 3.003e-6)

    euler_200 = time_solar_system.fill_sol(250, 1, integrator=time_solar_system.forward_euler)[0]
    ver_200 = time_solar_system.fill_sol(250, 1)[0]
    euler_50 = time_solar_system.fill_sol(50, 1, integrator=time_solar_system.forward_euler)[0]
    ver_50 = time_solar_system.fill_sol(50, 1)[0]
    euler_25 = time_solar_system.fill_sol(25, 1, integrator=time_solar_system.forward_euler)[0]
    ver_25 = time_solar_system.fill_sol(25, 1)[0]
    plt.axes(aspect = 'equal')
    plt.plot(euler_200[:,1,0], euler_200[:,1,1], "-", color = '#00cccc')
    plt.plot(euler_50[:,1,0],  euler_50[:,1,1],  "-", color = '#69a6ea')
    plt.plot(euler_25[:,1,0],  euler_25[:,1,1],  "-", color = '#736ace')
    plt.plot(ver_200[:,1,0],   ver_200[:,1,1],   "-", color = '#30c441')
    plt.plot(ver_50[:,1,0],    ver_50[:,1,1],    "-", color = '#000000')
    plt.plot(ver_25[:,1,0],    ver_25[:,1,1],    "-", color = '#956815')
    plt.plot(0,0,"o", color = "#FDB813")
    plt.axis([-7,1.5,-1.5,2])
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.title("Differing Timesteps with Velocity Verlet\nand Forward Euler Integrators")
    plt.legend(["dt=1/250 year, FE","dt=1/50   year, FE", "dt=1/25   year, FE","dt=1/250 year, VV","dt=1/50   year, VV","dt=1/25   year, VV"], loc='upper left')
    plt.savefig("../results/time_step.png", dpi = 1200)
    plt.clf()
    
    # benchmarking forward euler vs velocity verlet with the earth-sun System
    solar_system_bench = solar_system()
    solar_system_bench.create_planetary_object(0,0, 0,0, 1)                 # Sun
    solar_system_bench.create_planetary_object(1,0, 0,2*nmp.pi, 3.003e-6)   # Earth

    pathobj = "../results/bench_earth_sun.txt"
    with open (pathobj, 'w') as f_m:
        f_m.write("Benchmarking the Earth-Sun System with varying dt and solvers\n\n")

    for w, g in enumerate(['1e3', '1e4', '1e5', '1e6'], 1):
        for w2, g2 in enumerate(['Forward Euler', 'Velocity Verlet']):

            dt_1 = 100.0/float(g)

            if g2 == 'Forward Euler':
                pre = time.clock()
                p, v = solar_system_bench.fill_sol(int(float(g)), 100,integrator=solar_system_bench.forward_euler)
                post = time.clock()
                b_time = post - pre

            if g2 == 'Velocity Verlet':
                pre = time.clock()
                p, v = solar_system_bench.fill_sol(int(float(g)), 100,integrator=solar_system_bench.velocity_verlet)
                post = time.clock()
                b_time = post - pre

            pathobj = "../results/bench_earth_sun.txt"
            with open (pathobj, 'a') as f_m:
                f_m.write(g2 + " with dt = " + str(dt_1))
                f_m.write("\nComputational Time: " + str(b_time) + "\n\n")


    # energy conservation
    solar_system_Energy = solar_system()
    solar_system_Energy.create_planetary_object(0, 0, 0, 0, 1)
    solar_system_Energy.create_planetary_object(1, 0, 0, 2*nmp.pi, 3.003e-6)

    P, V = solar_system_Energy.fill_sol(int(1e6), 10)
    t = nmp.linspace(0,10,int(1e6)+1)
    KineticEnergyEarth = 0.5*solar_system_Energy.planetary_masses[1] * (V[:,1,0]**2 + V[:,1,1]**2) #SolarMasses*AU**2/yr**2
    KineticEnergySun = 0.5*solar_system_Energy.planetary_masses[0] * (V[:,0,0]**2 + V[:,0,1]**2)
    KineticEnergy = KineticEnergySun + KineticEnergyEarth

    CenterOfMass = P[:,0,:]*solar_system_Energy.planetary_masses[0] + P[:,1,:]*solar_system_Energy.planetary_masses[1]
    distance = nmp.sqrt( (P[:,1,0] - P[:,0,0])**2 + (P[:,0,1] - P[:,1,1])**2 )
    PotentialEnergy = -G*solar_system_Energy.planetary_masses[0]*solar_system_Energy.planetary_masses[1]/distance

    plt.plot(t, KineticEnergy, 'b')
    plt.plot(t, PotentialEnergy, 'g')
    plt.plot(t, KineticEnergy+PotentialEnergy, 'c')
    plt.axis([0,10,-1.5e-4,1.5e-4])
    plt.title("Energy of Circular Planet-Sun System over 10 years")
    plt.xlabel("Time (years)")
    plt.ylabel(r"Energy ($SolarMasses*\frac{{AU}^{2}}{{Year}^{2}}$)")
    plt.legend(["Kinetic Energy","Potential Energy","Total Energy"])
    plt.savefig("../results/energy.png", dpi = 1200)
    plt.clf()

    print("Center of mass at beginning of simulation, and after 10 years:")
    print(CenterOfMass[0], "\n", CenterOfMass[-1])

    AngularMomentum = solar_system_Energy.planetary_masses[1]*nmp.cross(P[:,1], V[:,1])
    plt.plot(t, AngularMomentum)
    plt.title("Angular Momentum of Planet over 10 years")
    plt.xlabel("Time (years)")
    plt.ylabel("Angular Momentum (L)")
    print("Relative error in angular momentum over 10 years: %e" % (( nmp.min(AngularMomentum) - nmp.max(AngularMomentum) ) / nmp.min(AngularMomentum)))
    plt.axis([0,10,0,4e-5])
    plt.savefig("../results/angular_momentum.png", dpi = 1200)
    plt.clf()


    # escape velocity
    for m in range(4,10):
        solar_system_b = solar_system()
        print(m)
        solar_system_b.create_planetary_object(0,0, 0,0, 1) #Sun
        solar_system_b.create_planetary_object(1,0, 0,(nmp.sqrt(m))*nmp.pi, 3.003e-6)
        p, v = solar_system_b.fill_sol(int(1e3), 50)
        plt.axis([-60,10,-20,40])
        plt.axes(aspect='equal')
        plt.plot(p[:,1,0], p[:,1,1], "-")
    plt.title("Testing velocities in Earth Sun System")
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.plot(p[:,0,0], p[:,0,1], "yo")
    plt.legend(["Planet1, v = $\sqrt{4}$$\pi$","Planet2, v = $\sqrt{5}$$\pi$","Planet3, v = $\sqrt{6}$$\pi$", "Planet4, v = $\sqrt{7}$$\pi$", "Planet5, v = $\sqrt{8}$$\pi$", "Planet6, v = $\sqrt{9}$$\pi$","Sun"])
    plt.savefig("../results/escape_vel.png", dpi = 1200)
    plt.clf()

    # escape velocity with acc equation modified by beta in set [2,3] testing end condition
    for m in range(4,10):
        solar_system_b = solar_system()
        print(m)
        solar_system_b.create_planetary_object(0,0, 0,0, 1) #Sun
        solar_system_b.create_planetary_object(1,0, 0,(nmp.sqrt(m))*nmp.pi, 3.003e-6)
        p, v = solar_system_b.fill_sol(int(1e3), 50, acc_method=solar_system_b.acc_beta)
        plt.axis([-60,10,-20,40])
        plt.axes(aspect='equal')
        plt.plot(p[:,1,0], p[:,1,1], "-")
    plt.title("Testing velocities in Earth Sun System")
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.plot(p[:,0,0], p[:,0,1], "yo")
    plt.legend(["Planet1, v = $\sqrt{4}$$\pi$","Planet2, v = $\sqrt{5}$$\pi$","Planet3, v = $\sqrt{6}$$\pi$", "Planet4, v = $\sqrt{7}$$\pi$", "Planet5, v = $\sqrt{8}$$\pi$", "Planet6, v = $\sqrt{9}$$\pi$","Sun"])
    plt.savefig("../results/escape_vel_beta.png", dpi = 1200)
    plt.clf()

    # three body problem with changing jupiter mass
    for w, g in enumerate(['1', '10', '1000'], 1):
        solar_system_c = solar_system()
        solar_system_c.create_planetary_object(0, 0, 0, 0, 1) #Sun
        solar_system_c.create_planetary_object(1, 0, 0, 29.8*0.210805, 3.003e-6, adjust_sun=False) #Earth
        solar_system_c.create_planetary_object(5.20, 0, 0, 13.1*0.210805, 954.7e-6*float(g), adjust_sun=False) #Jupiter

        p, v = solar_system_c.fill_sol(int(1e4), 25)

        plt.plot(p[:,0,0], p[:,0,1], "y-")
        plt.plot(p[:,1,0], p[:,1,1], "b-")
        plt.plot(p[:,2,0], p[:,2,1], "r-")
        plt.plot(p[0,0,0], p[0,0,1], "yx")
        plt.plot(p[-1,0,0], p[-1,0,1], "yo")
        plt.plot(p[0,2,0], p[0,2,1], "rx")
        plt.plot(p[-1,2,0], p[-1,2,1], "ro")
        plt.plot(p[0,1,0], p[0,1,1], "bx")
        plt.plot(p[-1,1,0], p[-1,1,1], "bo")
        if g == '1':
            plt.axis('equal')
            plt.axis([-10,10,-10,10])
        elif g == '10':
            plt.axis('equal')
            plt.axis([-10,10,-10,10])
        elif g == '1000':
            plt.axis('equal')
            plt.axis([-10,30,-5,50])
        plt.xlabel("AU")
        plt.ylabel("AU")
        plt.title("Sun-Earth-Jupiter system over 25 years\n with Jupiter mass multiplier of " + str(g))
        plt.legend(["Sun", "Earth", "Jupiter"])
        plt.savefig("../results/three_body_" + str(g) + ".png", dpi = 1200)
        plt.clf()

    masses = { 'mercury' : 0.166e-6, 'venus' : 2.081e-6, 'earth' : 3.003e-6, 'mars' : 0.323e-6, 'jupiter' : 954.7e-6, 'saturn' : 285.8e-6, 'uranus' : 43.6e-6, 'neptune' : 51.5e-6}
    distances = { 'mercury' : 0.39, 'venus' : 0.72, 'earth' : 1, 'mars' : 1.52, 'jupiter' : 5.20, 'saturn' : 9.58, 'uranus' : 19.23, 'neptune' : 30.10}
    speeds_kms = { 'mercury' : 47.4, 'venus' : 35.0, 'earth' : 29.8, 'mars' : 24.1, 'jupiter' : 13.1, 'saturn' : 9.7, 'uranus' : 6.8, 'neptune' : 5.4}
    speeds = {key: speeds_kms[key]*0.210805 for key in speeds_kms}

    solar_system = solar_system()
    solar_system.create_planetary_object(0, 0, 0, 0, 1)
    bodies = ['sun']

    # sort planets after distance
    body_dist = [(key, distances[key]) for key in distances]
    sort_obj = [item[0] for item in sorted(body_dist, key=itemgetter(1))]

    for key in sort_obj:
        solar_system.create_planetary_object(distances[key], 0, 0, speeds[key], masses[key])
        bodies.append(key)

    p, v = solar_system.fill_sol(int(1e4), 250, integrator = solar_system.velocity_verlet, acc_method = solar_system.Acc)

    for i in range(len(masses)+1):
        plt.plot(p[:,i,0], p[:,i,1])
    plt.axis([-40,60,-35,35])

    plt.title("Solar System Simulated over 250 years")
    plt.legend(bodies, loc='best')
    plt.savefig("../results/complete_solar_system.png", dpi = 1200)
    plt.clf()


    # perihelion precession of Mercury
    def find_perihelion_alt(minima_x, minima_y):
        # div = nmp.divide(minima_y, minima_x)
        print(minima_y)
        print(minima_x)
        angles = nmp.arctan(minima_y / minima_x)
        return angles

    n = int(1e8)     # total time steps
    years = 100.0    # earth years
    perihelion_minima = 416  # max number of minima to store
    run_perihelion = False   # be careful young jedi

    dt = years/float(n)
    print("dt = %g" % dt)
    radians_to_arcseconds = 206265.0

    mercurymass = 0.1652e-6
    mercuryspeed = 12.44
    mercury_perihelion_distance = 0.3075

    if run_perihelion == True:

        solar_system_d = solar_system()
        solar_system_d.create_planetary_object(0, 0, 0, 0, 1)
        solar_system_d.create_planetary_object(mercury_perihelion_distance, 0, 0, mercuryspeed, mercurymass)

        pre = time.clock()
        p_newton, v_newton, minima_x_n, minima_y_n, minima_time_n = solar_system_d.fill_sol_perihelion(n, years, perihelion_minima=perihelion_minima)
        post = time.clock()
        elapsed_time_n = post - pre

        angles_newton = find_perihelion_alt(minima_x_n, minima_y_n)

        pre = time.clock()
        p_einstein, v_einstein, minima_x_e, minima_y_e, minima_time_e = solar_system_d.fill_sol_perihelion(n, years, acc_method=solar_system_d.acc_rel, perihelion_minima=perihelion_minima)
        post = time.clock()
        elapsed_time_e = post - pre

        angles_einstein = find_perihelion_alt(minima_x_e, minima_y_e)

        time_newton   = minima_time_n
        time_einstein = minima_time_e

        pathobj = "../results/perihelion.txt"
        print(pathobj)
        with open (pathobj, 'w') as f_m:

            f_m.write('Newton Angles\n')
            nmp.savetxt(f_m, angles_newton, delimiter='\n')

            f_m.write('\nNewton Times\n')
            nmp.savetxt(f_m, time_newton, delimiter='\n')

            f_m.write('\nEinstein Angles\n')
            nmp.savetxt(f_m, angles_einstein, delimiter='\n')

            f_m.write('\nEinstein Times\n')
            nmp.savetxt(f_m, time_einstein, delimiter='\n')

            f_m.write('\nNewton Computational Time: ')
            f_m.write(str(elapsed_time_n))
            f_m.write('\n\nEinstein Computational Time: ')
            f_m.write(str(elapsed_time_e))

    else:

        pathobj = "../results/perihelion.txt"
        with open (pathobj, 'r') as f_m:
            blob = f_m.read().splitlines()

        angle_newton_2_float = [float(i) for i in blob[1:417]]  # when you get the chance one-liners
        angles_newton   = nmp.asarray(angle_newton_2_float)     # in python are so beautiful

        time_newton_2_float = [float(i) for i in blob[419:835]]
        time_newton     = nmp.asarray(time_newton_2_float)

        angle_einstein_2_float = [float(i) for i in blob[837:1253]]
        angles_einstein = nmp.asarray(angle_einstein_2_float)

        time_einstein_2_float = [float(i) for i in blob[1255:1671]]
        time_einstein   = nmp.asarray(time_einstein_2_float)

        plt.plot(time_newton, angles_newton*radians_to_arcseconds, '-')
        plt.plot(time_einstein, angles_einstein*radians_to_arcseconds, '-')

        plt.title("Change in Perihelion Angle over a Century\n$\Delta t$ = %g" % dt)
        plt.legend(["Classical Mechanics Case", "Relativistic Case"], loc="best")
        plt.xlabel("Time (years)")
        plt.ylabel("Angle (seconds of arc)")
        plt.savefig("../results/perihelion_angle.png", dpi = 1200)
        plt.clf()

        precession_newton = angles_newton*radians_to_arcseconds
        precession_einstein = angles_einstein*radians_to_arcseconds
        print("Precession per 100 years (Newton):   %g" % (precession_newton[-1]))
        print("Precession per 100 years (Einstein): %g" % (precession_einstein[-1]))
        print("Observed perihelions (Newton):   %d" % len(angles_newton))
        print("Observed perihelions (Einstein): %d" % len(angles_einstein))
