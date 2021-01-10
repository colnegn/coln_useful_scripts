#!/usr/bin/env python

import sys
import numpy as np
from numpy import linalg as LA


if len(sys.argv) != 5:
  sys.exit('usage: '
           '1.first AO-to-MO coef matrix; '
           '2.second AO-to-MO coef matrix; '
           '3.AO overlap matrix; '
           '4.number of occupied spatial orbitals (1/2 the # of electrons)')


C1   = np.genfromtxt(sys.argv[1])
C2   = np.genfromtxt(sys.argv[2])
S    = np.genfromtxt(sys.argv[3])
nocc = int(sys.argv[4])


## only want to use occupied columns of coef matrices:
oC1 = C1.T[:nocc].T
oC2 = C2.T[:nocc].T


#!# ignore first few AOs:
nAOs = 5
#print(oC1.T)
#print(oC2.T)
for moi in range(len(oC1.T)):
  for aoi in range(nAOs):
    oC1[aoi][moi] = 0.0
  oC1[:,moi] /= pow(np.matmul(np.matmul(oC1[:,moi].T, S), oC1[:,moi]), 0.5)
for moi in range(len(oC2.T)):
  for aoi in range(nAOs):
    oC2[aoi][moi] = 0.0
  oC2[:,moi] /= pow(np.matmul(np.matmul(oC2[:,moi].T, S), oC2[:,moi]), 0.5)
#print(oC1.T)
#print(oC2.T)
#!##################

## < \psi_1 | \psi_2 > = det| C_1^{\dag} S C_2 |
wf_overlap = pow(LA.det(np.matmul(np.matmul(oC1.T, S), oC2)), 0.5)

print(abs(wf_overlap))

