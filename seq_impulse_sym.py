#!/usr/bin/env python
# //==============================================================================
# /*
#     Software License Agreement (BSD License)
#     Copyright (c) 2022

#     All rights reserved.

#     Redistribution and use in source and binary forms, with or without
#     modification, are permitted provided that the following conditions
#     are met:

#     * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.

#     * Neither the name of authors nor the names of its contributors may
#     be used to endorse or promote products derived from this software
#     without specific prior written permission.

#     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#     FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
#     COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#     INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#     BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#     LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#     LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#     ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#     POSSIBILITY OF SUCH DAMAGE.

#     \author    <amunawa2@jh.edu>
#     \author    Adnan Munawar
# */
# //==============================================================================

# Notes:
# https://allenchou.net/2013/12/game-physics-resolution-contact-constraints/
# https://box2d.org/files/ErinCatto_SequentialImpulses_GDC2006.pdf
# http://www.mft-spirit.nl/files/articles/ImpulseSolverBrief.pdf
# https://danielchappuis.ch/download/ConstraintsDerivationRigidBody3D.pdf
# https://allenchou.net/2013/12/game-physics-constraints-sequential-impulse/

# General Equations
# P -> Impulse (or Sequence Impulse)
# J -> Constraint Jacobian
# M -> Mass matrix (Accumulated Mass Matrix)
# v -> Constraint velocity
# dt -> Time step
# bias -> Bias factor
# K -> Constraint Mass Matrix
# s(t) -> Combined state vector of two bodies between which there is a constraint i.e. [x_1 q_1 x_2 q_t]'
# C(s) -> Constraint, which is a function of the state vector.

# For 3D. E is identity Matrix (3x3 for 3D), m1 and m2 masses, I1 and I2 Inertias
# M = diag(m1 * E, I1, m2 * E, I2)


# P = transpose(J) * lambda_prime
# lambda = -inverse(K) * ( J * v_prime + bias)
# lambda_prime = lambda / dt
# K = J * inverse(M) * transpose(J)
# v_prime = v + inverse(M) * F_external * dt (For constraint, no need to put external force F_external)
# bias = b * C / dt
# b -> [0, 1]

# For non-penetrating contact
# n is the contact normal and p_1 and p_2 are contact points on body 1 and 2 in world coordinates. r_1 and r_2 are
# the vectors from the COM of body 1 and 2 to their respective contact points.
# e = p2 - p1
# C = dot(e, n) = dot((x_2 + r_2 - x_1 - r_1), n)

# dC/dt = (-n', -cross(r_1, n)', n_1', cross(r_2, n)') * (v_1, w_1, v_2, w_2)'
# Here a' indicates transpose of a vector a.

# The first part of dC/dt is the Jacobian for non penetration contact
# J = (-n', -cross(r_1, n)', n_1', cross(r_2, n)')

# From above
# K = (-n', -cross(r_1, n)', n_1', cross(r_2, n)') * inverse(M) * (-n', -cross(r_1, n)', n_1', cross(r_2, n)')'


import sympy as sp
from sympy import pprint


m1, ix1, iy1, iz1, m2, ix2, iy2, iz2 = sp.symbols('m1 ix1, iy1, iz1, m2, ix1, iy2, iz2')
M = sp.diag(m1 * sp.eye(3), sp.diag(ix1, iy1, iz1), m2 * sp.eye(3), sp.diag(ix2, iy2, iz2))
M_inv = M.inv()
# pprint(M)
# print('----')
# pprint(M_inv)

rx1, ry1, rz1, rx2, ry2, rz2, nx, ny, nz = sp.symbols('rx1, ry1, rz1, rx2, ry2, rz2, nx, ny, nz')
r1 = sp.Matrix([rx1, ry1, rz1])
r2 = sp.Matrix([rx2, ry2, rz2])
n = sp.Matrix([nx, ny, nz])
# pprint(n)

vx1, vy1, vz1, wx1, wy1, wz1, vx2, vy2, vz2, wx2, wy2, wz2 = \
    sp.symbols('vx1 vy1 vz1 wx1, wy1 wz1 vx2 vy2 vz2 wx2 wy2 wz2')
V = sp.Matrix([vx1, vy1, vz1, wx1, wy1, wz1, vx2, vy2, vz2, wx2, wy2, wz2])
# pprint(V)

x1, y1, z1, wx1, wy1, wz1, x2, y2, z2, wx2, wy2, wz2 = \
    sp.symbols('x1, y1, z1, wx1, wy1, wz1, x2, y2, z2, wx2, wy2, wz2')

P1 = sp.Matrix([x1, y1, z1])
P2 = sp.Matrix([x2, y2, z2])

C = (-P1 - r1 + P2 + r2).dot(n)
# pprint(C)

J = sp.Matrix([[-n.T, -r1.cross(n).T, n.T, r2.cross(n).T]])
# pprint(J)

K = J * M_inv * J.T
# pprint(K)

inv_K = K.inv()

dt, b = sp.symbols('dt b')

bias = -b / dt * C
# pprint(bias)

lam = -inv_K * (sum(J * V, bias))
pprint(lam.shape)

P = J.T * lam / dt
pprint(P)