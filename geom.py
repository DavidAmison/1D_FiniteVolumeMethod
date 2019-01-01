# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 15:27:09 2018
Class which defines geometry for the Process Modelling Individual Assignment
Shape is a truncated cone defined by a base radius, top radius and height.

Class defines functions for calculating radius, section areas and sectio volume
given x co-ordiantes. Base is assumed to be at x=0.

@author: David Amison
"""

import numpy as np
import unittest


class TruncatedCone:

    def __init__(self, left_radius, right_radius, length):
        self._ra = left_radius
        self._rb = right_radius
        self._L = length

    def length(self):
        return self._L

    def radius(self, x):
        if x < 0:
            raise ValueError('Cannot find radius at co-ordinate less than 0')
        elif x > self._L:
            raise ValueError('Cannot find radius at co-ordinate greater than'
                             ' cone length')
        else:
            return self._ra - x*(self._ra-self._rb)/self._L

    def section_area(self, x):
        return np.pi*self.radius(x)**2

    def section_volume(self, x1, x2):
        r1 = self.radius(x1)
        r2 = self.radius(x2)
        return np.pi*(r1**2+r1*r2+r2**2)*abs(x2-x1)/3


class TestTruncatedCone(unittest.TestCase):

    def test_length(self):
        cone = TruncatedCone(8, 4, 20)
        self.assertEqual(cone.length(), 20)

    def test_base_radius(self):
        cone = TruncatedCone(8, 4, 20)
        self.assertEqual(cone.radius(0), 8)

    def test_top_radius(self):
        cone = TruncatedCone(8, 4, 20)
        self.assertEqual(cone.radius(20), 4)

    def test_radius(self):
        cone = TruncatedCone(4, 8, 20)
        self.assertEqual(cone.radius(15), 7)

    def test_section_area(self):
        cone = TruncatedCone(8, 4 ,20)
        self.assertEqual(cone.section_area(10), np.pi*6**2)

    def test_section_volume(self):
        cone = TruncatedCone(8, 4, 20)
        self.assertAlmostEqual(cone.section_volume(5, 10), 664.9704, places=4)


if __name__ == '__main__':
    # Perform unit tests if desired
    unit_tests = True
    if unit_tests:
        unittest.main()
