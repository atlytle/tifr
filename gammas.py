# Define gamma matrices.  Need to CHECK overlap & hisq conventions.

import sys
import numpy as np

# Kentucky format gamma matrices.
# The overlap propagators are written in this convention.

gx = np.array([[0, 0, 0, -1j],
               [0, 0, -1j, 0],
               [0, 1j, 0, 0],
               [1j, 0, 0, 0]])

gy = -np.array([[0, 0, 0, 1],
               [0, 0, -1, 0],
               [0, -1, 0, 0],
               [1, 0, 0, 0]])

gz = np.array([[0, 0, -1j, 0],
               [0, 0, 0, 1j],
               [1j, 0, 0, 0],
               [0, -1j, 0, 0]])

gt = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, -1, 0],
               [0, 0, 0, -1]])
               
g5 = np.array([[0, 0, 1, 0],
               [0, 0, 0, 1],
               [1, 0, 0, 0],
               [0, 1, 0, 0]])

id4 = np.identity(4, dtype=complex)
id3 = np.identity(3, dtype=complex)


# Transformation from Kentucky to "Wilson"
T = np.array([[0, -1, 0, 1],
              [1, 0, -1, 0],
              [0, 1, 0, 1],
              [-1, 0, -1, 0]])/np.sqrt(2)


# MILC matrices.
gxM = np.array([[0,0,0,1j],
                [0,0,1j,0],
                [0,-1j,0,0],
                [-1j,0,0,0]])

gyM = np.array([[0,0,0,-1],
                [0,0,1,0],
                [0,1,0,0],
                [-1,0,0,0]])

gzM = np.array([[0,0,1j,0],
                [0,0,0,-1j],
                [-1j,0,0,0],
                [0,1j,0,0]])

gtM = np.array([[0,0,1,0],
                [0,0,0,1],
                [1,0,0,0],
                [0,1,0,0]])

g5M = np.array([[1,0,0,0],
                [0,1,0,0],
                [0,0,-1,0],
                [0,0,0,-1]])