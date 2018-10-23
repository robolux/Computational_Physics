# Project 3 Solvers - Computational Physics
# Hunter Phillips

import math
import numpy as nmp
import matplotlib.pyplot as plt
import sys
import time
import os

# init constants
G = 4*nmp.pi**2
c = 63197.8

# start of object oriented code
class solar_system:
    def __init__(self):
        self.n_planets = 0

    def create_planetary_object(self, x0, y0, vx0, vy0, mass, adjust_sun=True):
        if self.n_planets == 0:
            self.planetary_positions = nmp.array( [[x0, y0]] , dtype=nmp.float64)
            self.planetary_velocities = nmp.array( [[vx0, vy0]] , dtype=nmp.float64)
            self.planetary_masses = nmp.array( mass , dtype=nmp.float64)
        else:
            self.planetary_positions = nmp.append( self.planetary_positions, [[x0, y0]], axis=0 )
            self.planetary_velocities = nmp.append( self.planetary_velocities, [[vx0, vy0]], axis=0 )
            self.planetary_masses = nmp.append( self.planetary_masses, mass )
            if adjust_sun:
                self.repop_sun()
        self.n_planets += 1

    # adjust the sun to the coordinate frame [0, 0]
    def repop_sun(self):
        solar_mass = self.planetary_masses[0]
        new_planetary_object_mass = self.planetary_masses[-1]
        mass_ratio = new_planetary_object_mass / solar_mass
        self.planetary_positions[0,:] -= self.planetary_positions[-1,:] * mass_ratio
        self.planetary_velocities[0,:] -= self.planetary_velocities[-1,:] * mass_ratio

    # Classical Acceleration
    def Acc(self, Positions, Velocity, target, Masses ):
        x_acc = 0
        y_acc = 0
        global distance
        global x_distance
        global y_distance
        for i in xrange(self.n_planets):
            if i != target:
                x_distance = Positions[target,0] - Positions[i,0]
                y_distance = Positions[target,1] - Positions[i,1]
                distance = math.sqrt( x_distance**2 + y_distance**2 )
                x_acc -= G*Masses[i]*x_distance/distance**3
                y_acc -= G*Masses[i]*y_distance/distance**3
        return nmp.array([x_acc, y_acc])

    # classical acceleration with modified beta in set [2,3] testing end condition {3}
    def acc_beta(self, Positions, Velocity, target, Masses ):
        x_acc = 0
        y_acc = 0
        for i in xrange(self.n_planets):
            if i != target:
                x_distance = Positions[target,0] - Positions[i,0]
                y_distance = Positions[target,1] - Positions[i,1]
                distance = math.sqrt( x_distance**2 + y_distance**2 )
                x_acc -= G*Masses[i]*x_distance/distance**4
                y_acc -= G*Masses[i]*y_distance/distance**4
        return nmp.array([x_acc, y_acc])

    # relativistic acceleration
    def acc_rel(self, Positions, Velocity, target, Masses):
        x_acc = 0
        y_acc = 0
        global distance
        for i in xrange(self.n_planets):
            if i != target:

                x_distance = Positions[target,0] - Positions[i,0]
                y_distance = Positions[target,1] - Positions[i,1]
                distance = math.sqrt( x_distance**2 + y_distance**2 )

                l = Positions[target,0]*Velocity[1] \
                    - Positions[target,1]*Velocity[0]
                rel_fac = 1 + ( (3*l**2) / (distance**2*c**2) )

                x_acc -= G*Masses[i]*x_distance/distance**3*rel_fac
                y_acc -= G*Masses[i]*y_distance/distance**3*rel_fac
        return nmp.array([x_acc, y_acc])

    # forward euler solver
    def forward_euler(self, P, V, p_new, V_new, dt, acc_method):
        length = len(P)
        for n in xrange(length):
            p_new[n] = P[n] + V[n]*dt
            V_new[n] = V[n] + acc_method(P, V[n], n, self.planetary_masses)*dt

    # velocity verlet solver
    def velocity_verlet(self, P, V, p_new, V_new, dt, acc_method):
        length = len(P)
        acc_p = nmp.zeros((length, 2), dtype = nmp.float64)
        for n in xrange(length):
            acc_p[n] = acc_method(P, V[n], n, self.planetary_masses)
            p_new[n] = P[n] + V[n]*dt + 0.5*acc_p[n]*dt**2
        for n in xrange(length):
            acc_p_new = acc_method(p_new, V[n], n, self.planetary_masses)
            V_new[n] = V[n] + 0.5*(acc_p[n] + acc_p_new)*dt

    def fill_sol(self, steps, years, integrator = None, acc_method = None):
        if integrator == None:
            integrator = self.velocity_verlet
        if acc_method == None:
            acc_method = self.Acc
        num_objects = self.n_planets
        dt = float(years)/steps
        print "dt = %g" % dt
        p = nmp.zeros(shape = (steps+1, num_objects, 2))
        v = nmp.zeros(shape = (steps+1, num_objects, 2))
        p[0] = self.planetary_positions[:][:]
        v[0] = self.planetary_velocities[:][:]

        for i in xrange( steps ):
            integrator(p[i], v[i], p[i+1], v[i+1], dt, acc_method)

        return p, v

    def fill_sol_perihelion(self, steps, years, integrator = None, acc_method = None, perihelion_minima=None):
        if integrator == None:
            integrator = self.velocity_verlet
        if acc_method == None:
            acc_method = self.Acc
        num_objects = self.n_planets
        dt = float(years)/steps
        print "dt = %g" % dt
        p = nmp.zeros( shape = ( steps+1, num_objects, 2 ) )
        v = nmp.zeros( shape = ( steps+1, num_objects, 2 ) )
        p[0] = self.planetary_positions[:][:]
        v[0] = self.planetary_velocities[:][:]

        minima_x = nmp.zeros((perihelion_minima, 1), dtype=nmp.float64)
        minima_y = nmp.zeros((perihelion_minima, 1), dtype=nmp.float64)
        minima_time = nmp.zeros((perihelion_minima, 1), dtype=nmp.float64)

        if perihelion_minima > 1:
            curr_distance = 1.0
            prev_distance = curr_distance
            prev_prev_distance = curr_distance

            curr = 1.0
            prev = curr
            prev_prev = curr
            countt = 0

        for i in xrange(steps):

            integrator(p[i], v[i], p[i+1], v[i+1], dt, acc_method)
            prev_prev = prev
            prev = curr
            curr = distance

            if ((prev_prev > prev) & (prev < curr)):
                print(p[i][1][0]) # this is x
                print(p[i][1][1]) # this is y
                minima_x[countt] = p[i][1][0]
                minima_y[countt] = p[i][1][1]
                minima_time[countt] = countt*dt
                countt = countt + 1
                print(countt)

        return p, v, minima_x, minima_y, minima_time
