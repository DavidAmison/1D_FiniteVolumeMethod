# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 18:10:37 2018

@author: David
"""

import numpy as np
from geom import TruncatedCone

class FVM_1D_Mesh():

    def __init__(self, geom):
        '''
        geom is a file defining the geometry to be meshed and contains the
        functions:
            length: returns the length of the geometry in m.
            section_area(x): returns the cross sectional area at the point x in
            m^2.
            section_volume(x1, x2): returns the volume bounded by x1 and x2 in
            m^3.
        '''
        self._geom = self.check_geometry(geom)  # Defines the geometry for mesh
        self.nodes = None   # For defining node co-ordinates
        self._mesh_edges = None     # For defining mesh edge co-ordinates
        self.dxE = None     # Defines dx values to east of nodes
        self.dxW = None     # Defines dx values to west of nodes
        self.AE = None      # Defines area of surface east of control volumes
        self.AW = None      # Defines area of surface west of control volumes
        self.dV = None      # Defines volume of control volumes
        return

    def check_geometry(self, geom):
        '''
        Checks that geometry has all the required functions. Raises exception
        if function doesn't exist, returns geom if geometry object is good.
        '''
        # Check for length function
        if hasattr(geom, 'length'):
            if not callable(geom.length):
                raise ValueError('geom object does not contain length'
                                 ' function')
        else:
            raise ValueError('geom object does not contain length'
                             ' function')
        # Check for section_area function
        if hasattr(geom, 'section_area'):
            if not callable(geom.section_area):
                raise ValueError('geom object does not contain section_area'
                                 ' function')
        else:
            raise ValueError('geom object does not contain section_area'
                             ' function')
        # Check for section_volume function
        if hasattr(geom, 'section_volume'):
            if not callable(geom.section_volume):
                raise ValueError('geom object does not contain section_volume'
                                 ' function')
        else:
            raise ValueError('geom object does not contain section_volume'
                             ' function')
        return geom

    def generate_linspace_mesh(self, n):
        '''
        Generates a mesh of equally spaces nodes/volumes where n is the number
        of nodes/volumes in the mesh.
        '''
        # Check that the value passed is an intiger
        if type(n) is int:
            length = self._geom.length()
            offset = length/(2*n)
            # Generate array defining x co-ordinate of nodes
            self.nodes = np.linspace(offset, length-offset, n)
            # Generate mesh metrics (edges, volumes, areas and dx values)
            self.generate_mesh_metrics()
            return
        else:
            raise ValueError('Cannot generate mesh, n must be an intiger')

    def generate_mesh_from_array(self, nodes):
        '''
        Takes an array as input which defines the points of nodes on the mesh.
        Note: nodes do not have to be sorted
        '''
        # First convert to a numpy array and sort in ascending order
        p = np.array(nodes).sort()
        # Check that the largest value in p is less than length of the geom
        if p[-1] > self._geom.length():
            raise ValueError('Mesh nodes defined outsie of geometry length')
        if p[0] < 0:
            raise ValueError('Mesh nodes defined outside of geometry length')
        # Store node values
        self.nodes = p
        # Generate other mesh values
        self.generate_mesh_metrics()
        return

    def generate_mesh_metrics(self):
        '''
        Uses current mesh nodes to generate mesh_edges, dxE, dxW, AE, AW and V
        '''
        self.generate_mesh_edges()
        self.area_east()
        self.area_west()
        self.dx_east()
        self.dx_west()
        self.section_volumes()
        return

    def generate_mesh_edges(self):
        '''
        Uses the mesh nodes to find the x co-ordinated of the edges
        '''
        # Caclulate dx between nodes
        self._mesh_edges = np.array([(self.nodes[n+1]+self.nodes[n])/2
                                     for n in range(0, len(self.nodes)-1)])
        # Insert dx value between left of geometry and first node
        self._mesh_edges = np.insert(self._mesh_edges, 0, 0)
        # Insert co-ordinate of right geometry
        self._mesh_edges = np.append(self._mesh_edges, self._geom.length())
        return

    def section_volumes(self):
        '''
        Finds the volume of every element in the mesh
        '''
        dV = [self._geom.section_volume(self._mesh_edges[n], self._mesh_edges[n+1])
              for n in range(0, len(self._mesh_edges)-1)]
        self.dV = np.array(dV)
        return

    def area_east(self):
        '''
        generates an array of all the areas at to the east (i.e. +ve x) of the
        nodes
        '''
        Ae = [self._geom.section_area(self._mesh_edges[n+1])
              for n in range(0, len(self._mesh_edges)-1)]
        self.AE = np.array(Ae)
        return

    def area_west(self):
        '''
        generates an array of all the areas to the west (i.e. -ve x) of the
        nodes
        '''
        Aw = [self._geom.section_area(self._mesh_edges[n])
              for n in range(0, len(self._mesh_edges)-1)]
        self.AW = np.array(Aw)
        return

    def dx_east(self):
        '''
        dx values to the east (i.e. +ve x) of the nodes
        '''
        dxE = [self.nodes[n+1]-self.nodes[n]
               for n in range(0, len(self.nodes)-1)]
        dxE.append(self.nodes[0])
        self.dxE = np.array(dxE)
        return

    def dx_west(self):
        '''
        dx values to the west (i.e. -ve x) of the nodes
        '''
        dxW = [self.nodes[n]-self.nodes[n-1]
               for n in range(1, len(self.nodes))]
        dxW.insert(0, self.nodes[0])
        self.dxW = np.array(dxW)
        return
