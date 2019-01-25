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
from matplotlib import animation
import numpy as np
import time

# Setup the geometry
geom = TruncatedCone(left_radius=0.04, right_radius=0.08, length=0.2)

# Setup the materisl
material = Material(rho=7800, c=490, k=54)

# Mesh the geometry
n_nodes = 200
mesh = FVM_1D_Mesh(geom)
mesh.generate_linspace_mesh(n_nodes)


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


# Setup solver conditions
end_time = 50
end_condition = 1E-2    # Only needed if end_time is set to 'steady'
t_step = 'auto'         # Let the solver decide for us

# Variables for setting up the plot
n_cont = 71
cont_min = 25
cont_max = 200

# Variables for making the animation (if desired)
make_animation = False
n_frames = 100

# Generates colours for contour plot
colours = []
for i in np.linspace(0, 1, n_cont):
    if i < 0.25:
        colours.append([0, i*4, 1])
    elif i < 0.5:
        i = i-0.25
        colours.append([0, 1, 1-4*i])
    elif i < 0.75:
        i = i-0.5
        colours.append([4*i, 1, 0])
    else:
        i = i-0.75
        colours.append([1, 1-4*i, 0])

fig = plt.figure()
ax = plt.axes()


def animate(i, fargs):
    solver = fargs
    x = np.array([[p, p] for p in solver._mesh.nodes])
    y = np.array([[solver._geom.radius(p), -solver._geom.radius(p)]
                 for p in solver._mesh.nodes])
    z = np.array([[p, p] for p in solver.solution[int(i)][1]])
    ax.cla()
    solver._geom.plot_shape(ax)
    ax.contourf(x, y, z,
                colors=colours,
                levels=np.linspace(cont_min, cont_max, n_cont))
    ax.text(0.01, 0.07, 't={:.1F} s'.format(solver.solution[int(i)][0]))
    ax.set_title("Temperature over time")
    ax.set_xlabel("x-coordinate (m)")
    ax.set_ylabel("y-coordinate (m)")


if __name__ == '__main__':
    print('Running...')
    # Solve the problem
    solver = Solver(geom, material, mesh, boundary, initial_cond)
    start = time.time()
    if isinstance(end_time, int):
        solver.solve(t_end=end_time, t_step=t_step)
    elif isinstance(end_time, str) and end_time == 'steady':
        solver.solve_steady(target=end_condition)
    else:
        raise ValueError("end_time must be an intiger or the string 'steady'")
    end = time.time()
    print('Total run time: {:.3F} s'.format(end-start))
    # Plot the end solution

    if make_animation is True:
        # Animate the solution (500 frames, 10 seconds)
        print('Making animation of solution...')
        n = len(solver.solution)
        x = np.array([[p, p] for p in mesh.nodes])
        y = np.array([[geom.radius(p), -geom.radius(p)]
                     for p in mesh.nodes])
        z = np.array([[p, p] for p in solver.solution[0][1]])
        geom.plot_shape(ax)
        c = ax.contourf(x, y, z,
                        colors=colours,
                        levels=np.linspace(cont_min, cont_max, n_cont))
        fig.colorbar(c, label="Temperature, Celcius")
        anim = animation.FuncAnimation(fig, animate,
                                       frames=np.linspace(0, n-1, n_frames),
                                       interval=20, blit=False, fargs=[solver])
        anim.save('animation.mp4', fps=50)

    # Uncomment to make contour plot of temperature distribution

    geom.plot_shape(ax)
    x = np.array([[p, p] for p in mesh.nodes])
    y = np.array([[geom.radius(p), -geom.radius(p)]
                     for p in mesh.nodes])
    z = np.array([[p, p] for p in solver.solution[-1][1]])
    c = ax.contourf(x, y, z, colors=colours,
                    levels=np.linspace(25, 200, n_cont))
    fig.colorbar(c, label="Temperature, Celcius")
    ax.text(0.01, 0.07, 't={:.2F} s'.format(solver.solution[-1][0]))
    ax.set_title("Temperature after {:.2F} seconds".format(
            solver.solution[-1][0]))
    ax.set_xlabel("x-coordinate (m)")
    ax.set_ylabel("y-coordinate (m)")

    # Uncomment to plot dQ_RMS against time
    '''
    t = [s[0] for s in solver.solution]
    ax.semilogy(t, solver.dQ_RMS, 'k-')
    ax.set_title("Average Net Heat Transfer over Time")
    ax.set_xlabel("time, s")
    ax.set_ylabel("dQ_RMS, W/m^3.s")
    fig.savefig('LineGraph.png')
    '''
    # Uncomment to plot line graph of temperature against x-coordinate
    '''
    ax.plot(mesh.nodes, solver.solution[-1][1], 'k-')
    ax.text(0.02, 175, 't={:.2F} s'.format(solver.solution[int(-1)][0]))
    ax.set_title("Temperature after 10 seconds")
    ax.set_xlabel("x-coordinate (m)")
    ax.set_ylabel("Temperature, Celcius")
    fig.savefig('LineGraph.png')
    '''
