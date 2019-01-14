# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 15:54:31 2018
Class to solve 1D Finite Volume Unsteady Heat Conduction problems.
Currenlty only setup to use the explicit method.

@author: David
"""

import numpy as np
import time


class FVM_UnsteadyConduction_Explicit:

    def __init__(self, geom, material, mesh, boundary_temp,
                 initial_temp):
        '''
        geom is a class which contains the following functions:
            length: returns the length of the geometry in m.
            section_area(x): returns the cross sectional area at the point x in
            m^2.
            section_volume(x1, x2): returns the volume bounded by x1 and x2 in
            m^3.
        material is a class which defines the following variables:
            material.rho - the density in kg/m^3.
            material.c  - the specific heat capacity in J/kg.K
            material.k   - the conductivity in W/m.K
        mesh is a class which defines the following variables:
            mesh.nodes   - the x-cordinates of the nodes in m
            mesh.dxE     - dx value in the +ve x direction
            mesh.dxW     - dx value in the -ve x direction
            mesh.AE      - section area in the +ve x direction from the nodes
            mesh.AW      - section area in the -ve x direction from the nodes
            mesh.dV      - volume of control volumes at each node
        boundary_temp is a function which returns the temperature in Celcius in
            format T_left, T_right.
        initial_temp is a function which returns the temperature in Celcius
            based on the x co-ordinate
        '''
        self._geom = geom
        self._mat = material
        self._mesh = mesh
        # Generate or copy the array of initial temperatures
        self._T0 = [initial_temp(x) for x in mesh.nodes]
        self._bound = boundary_temp     # Function varying with time
        # Array which will store the solution in form [time, [Tp]]
        # Time is in seconds and Tp is temperatures of nodes in celcius
        self.solution = [[0, self._T0]]
        self.dQ_RMS = [0]
        return

    def solve(self, t_end, t_step='auto'):
        '''
        Performs the calculations to solve the Unsteady Heat Conduction problem
        using the explicit, finite volume method. Solution is found based on
        given geometry, mesh, boundary and initial conditions.

        t_end is the time to end the simulation in seconds.
        t_step is the dt between each step, either a number in seconds or auto
            to let the solver choose the best time step.
        '''
        # Copy variables from _mesh object
        dxE = self._mesh.dxE
        dxW = self._mesh.dxW
        AE = self._mesh.AE
        AW = self._mesh.AW
        dV = self._mesh.dV
        # Calculate values to generate matrix
        alpha = self._mat.k/(self._mat.rho*self._mat.c)
        # Check the time step is appropriate
        dt_max = 1/(alpha*((AE/dxE)+(AW/dxW))/dV)
        if t_step == 'auto':
            t_step = np.min(dt_max)*0.1
            print('Time step set to {:.2E} s'.format(t_step))
        elif t_step > np.min(dt_max):
            s = input('Given time stop of {:.2E} is greater than recommended '
                      'maximum of {:.2E} for stability.\n'
                      'Continue anyway?(y/n) '.format(t_step, np.min(dt_max)))
            if s.lower() == 'y':
                print('Continuing run...')
            else:
                print('Cancelling run...')
                return
        C_p1 = t_step*alpha*AE/(dV*dxE)
        C_p3 = t_step*alpha*AW/(dV*dxW)
        C_p2 = 1 - C_p1 - C_p3
        # Generate matrix to solve the problem
        M = self.generate_matrix(C_p1, C_p2, C_p3)
        # Generate the solution
        TA, TB = self._bound(0)
        T0 = np.concatenate(([TA], self._T0, [TB]))
        T_new = T0
        start = time.time()
        # Runs till time is greater than t_end
        for n in range(0, int(np.ceil(t_end/t_step))):
            T_old = T_new
            T_new = np.matmul(M, T_old)
            T_new[0], T_new[-1] = self._bound(t_step*(n+1))
            # Find the dQ_RMS value
            dQ_RMS, dQ = self.net_heat_transfer(T_old, T_new, t_step)
            self.dQ_RMS.append(dQ_RMS)
            # Add the new solution to the array
            self.solution.append([(n+1)*t_step, T_new[1:-1]])
            if time.time()-start > 5:
                print('Iteration {}   simulation time={:.2E} s   dQ_RMS={:.2E}'
                      .format(n, (n+1)*t_step, self.dQ_RMS[-1]))
                start = time.time()
        return

    def solve_steady(self, t_step='auto', target=1E-5, max_iterations=None):
        '''
        Performs the calculations to solve the Unsteady Heat Conduction problem
        using the explicit, finite volume method. Solution is found based on
        given geometry, mesh, boundary and initial conditions. S

        t_step is the dt between each step, either a number in seconds or auto
            to let the solver choose the best time step.
        target is the dQ_RMS value at which steady is determined to have been
            reached
        max_iterations is the maximum number of steps the solver will complete
        '''
        # Copy variables from _mesh object
        dxE = self._mesh.dxE
        dxW = self._mesh.dxW
        AE = self._mesh.AE
        AW = self._mesh.AW
        dV = self._mesh.dV
        # Calculate values to generate matrix
        alpha = self._mat.k/(self._mat.rho*self._mat.c)
        # Check the time step is appropriate
        dt_max = 1/(alpha*((AE/dxE)+(AW/dxW))/dV)
        if t_step == 'auto':
            t_step = np.min(dt_max)*0.1
            print('Time step set to {:.2E} s'.format(t_step))
        elif t_step > np.min(dt_max):
            s = input('Given time stop of {:.2E} is greater than recommended '
                      'maximum of {:.2E} for stability.\n'
                      'Continue anyway?(y/n) '.format(t_step, np.min(dt_max)))
            if s.lower() == 'y':
                print('Continuing run...')
            else:
                print('Cancelling run...')
                return
        C_p1 = t_step*alpha*AE/(dV*dxE)
        C_p3 = t_step*alpha*AW/(dV*dxW)
        C_p2 = 1 - C_p1 - C_p3
        # Generate matrix to solve the problem
        M = self.generate_matrix(C_p1, C_p2, C_p3)
        # Generate the solution
        TA, TB = self._bound(0)
        T0 = np.concatenate(([TA], self._T0, [TB]))
        T_new = T0
        # Runs till time is greater than t_end
        n = 0   # Tracks the number of iterations
        self.dQ_RMS[0] = target+1
        start = time.time()
        while self.dQ_RMS[-1] > target:
            T_old = T_new
            T_new = np.matmul(M, T_old)
            T_new[0], T_new[-1] = self._bound(t_step*(n+1))
            # Find the dQ_RMS value
            dQ_RMS, dQ = self.net_heat_transfer(T_old, T_new, t_step)
            self.dQ_RMS.append(dQ_RMS)
            # Add the new solution to the array
            self.solution.append([(n+1)*t_step, T_new[1:-1]])
            n = n+1
            if time.time()-start > 5:
                print('Iteration {}   simulation time={:.2E} s   dQ_RMS={:.2E}'
                      .format(n, (n+1)*t_step, self.dQ_RMS[-1]))
                start = time.time()
            if max_iterations is not None and n >= max_iterations:
                print('Note: Solver reached iteration limit before converging')
                break
        return

    def generate_matrix(self, C_p1, C_p2, C_p3):
        '''
        Generates the matrix for solving each time step of the from:
            [1    0    0    0    ....]
            [C_01 C_02 C_03 0    ....]
            [0    C_11 C_12 C_13 ....]
            [.... .... .... .... ....]
        '''
        # First generate an empy matrix
        M = np.zeros([len(C_p1)+2, len(C_p1)+2])
        M[0][0] = 1
        M[-1][-1] = 1
        # Now populate the matrix
        for i in range(1, len(C_p1)+1):
            M[i][i-1] = C_p1[i-1]
            M[i][i] = C_p2[i-1]
            M[i][i+1] = C_p3[i-1]
        return M

    def net_heat_transfer(self, t_old, t_new, t_step):
        '''
        Calculates the overall heat transfer per unit volume into each control
        volume.
        returns dQ_RMS, dQ
        Where dQ is an array of values for each node and dQ_RMS is the root
        mean square.
        '''
        dQ = self._mat.rho*self._mat.c*(t_new[1:-1]-t_old[1:-1])/t_step
        dQ_RMS = np.sqrt(np.mean(dQ**2))
        return dQ_RMS, dQ

    def temp_at_point(self, x):
        '''
        Returns the temperature at a point (x) on the geometry based on data
        from the solution. Temperature is  found using linear interpolation
        from the two nearest points.
        '''
        if len(self.solution[-1][1]) != len(self._mesh.nodes):
            raise ValueError('Cannot give temperature before solver is run')
        return np.interp(x, self._mesh.nodes, self.solution[-1][1])
