# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 16:59:32 2018

@author: David
"""
from geom import TruncatedCone
from material import Material
from mesh import FVM_1D_Mesh
from FVM_UnsteadyConduction import FVM_UnsteadyConduction_Explicit as Solver
import matplotlib.pyplot as plt
import time


fig = plt.figure()
ax = plt.axes(xlim=(0, 0.2), ylim=(0, 200))


def boundary(t):
    '''
    Returns the temperature at time t seconds, in celcius. Return format is:
    left_temperature, right_temperature
    '''
    TA = 200
    TB = 100
    return TA, TB


def initial_cond(x):
    '''
    Returns the temperature in cecius at a point x metres when t=0
    '''
    return 25


if __name__ == '__main__':
    geom = TruncatedCone(left_radius=0.04, right_radius=0.08, length=0.2)
    material = Material(rho=7800, cp=490, k=54)
    # Mesh the geometry
    mesh = FVM_1D_Mesh(geom)
    mesh.generate_linspace_mesh(300)
    # Setup the solver
    solver = Solver(geom, material, mesh, boundary, initial_cond)
    # Solve for 10 seconds with auto time-step
    start = time.time()
    solver.solve(t_end=10, t_step='auto')
    end = time.time()
    print('Total run time: {:3} s'.format(end-start))
    #solver.solve_steady(target=1E-5)
    # Plot the end solution
    ax.plot(solver._mesh.nodes, solver.solution[-1][1])
    plt.show()
    print(solver.solution[-1][1][0])
