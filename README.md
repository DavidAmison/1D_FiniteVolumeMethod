# 1D_FiniteVolumeMethod
Python code for solving transient heat conduction problems for simple geometries (i.e. Cuboid, Cylinder, Cone etc.)
Currenlty only the explicit method is implemented.

## How it Works

### Geometry

The geom.py file contains classes which define different, simple, 3D geometries. Hopefully this will contain more as time goes on and feel free to contribute your own geometry classes if you come up with a new one.

Each class must define (at a minimum) the following functions:
  1) length(): Which returns the overall length of the shape in the x direction.
  2) section_area(x): Which takes an x co-ordinate and returns the cross-sectional area of the shape at that point.
  3) section_volume(x1, x2): Which takes two x co-ordinates and returns the volume between them.

### Mesh

The mesh.py file contains the FVM_1D_Mesh class which can be used to generate a one-dimensional, finite volume mesh for the geometry. On initialisation it must be passed the geometry which is being analysed.
You can then either have the mesh automatically generated based on the desired number of control volumes/nodes by calling generate_linspace_mesh(n) where n is the number of nodes; or you can generate a mesh based on your own set of points, stored in an array and passed to the generate_mesh_from_array(nodes) function.

### Material Properties

The material.py file defines the Material class which is used to store density, heat capacity and thermal conductivity. Currenlty they are set as constant but this could be edited to make the values dependant on x-coordinate or time (however this would also need to be implemented in the solver class)

### Boundary and Initial Conditions

Since we are working with only 1-Dimensional problems for now the solver requires onlly 2 boundary conditions at the extremes of the geometry (i.e. left (x=0) and right (x=L)). This must be defined in a function which accepts one input (time) and returns the temperatures. The following is an example:

def boundary_condtion(t):
  Ta = f(t)
  Tb = g(t)
  return Ta, Tb

Initial conditions are likewise defined by a function that accepts one input (x coordinate this time) and returns the temperature at that position. i.e.

def initial_temp(x):
  return f(x)

### Setting up the solver

To actually run your problem you need to first set-up the solver with all your inputs on whichever solver you want (geometry, material, mesh, boundary conditions and initial conditions) then run the required solver method (either transient based on an end time and time step or steady state based on an end requirement for root-mean-square value of heat flow into the control volumes).

Note: if time_step is set to auto the solver will calculate and choose an appropriate value based on the governing equations.

#### Current Solver Classes

FVM_UnsteadyCondiction_Explicit:
  Uses the explicit method (relatively simple method and works great as long the mesh isn't too fine; for fine meshes the maximum time step required for stability becomes very small leading to long run times, currently runs at about (n)^3.8 efficiency)


  
