# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 16:43:09 2018

Defines the material properties required for heat conduction:
    density in kg/m^3
    heat capacity in J/kg.K
    conductivity in W/m.K

@author: David
"""


class Material:

    def __init__(self, rho, c, k):
        '''
        rho is density in kg/m^3
        c is heat capacity in J/kg.K
        k is conductivity in W/m.K
        '''
        self.rho = rho
        self.c = c
        self.k = k

    def set_density(self, rho):
        self.rho = rho

    def set_heat_capacity(self, c):
        self.c = c

    def set_conductivity(self, k):
        self.k = k
