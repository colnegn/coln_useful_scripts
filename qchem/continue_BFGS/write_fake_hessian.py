#!/usr/bin/env python

import sys
import numpy as np

if len(sys.argv) != 2:
  sys.exit('usage: 1.path to HESS file from scratch to convert to binary 132.0 file')

with open(sys.argv[1],'r') as f:
  file_list = [li.strip().split() for li in f.readlines()]

## dimension printed in HESS file:
dim = int(file_list[1][1])

## flatten hessian:
hess_list = []
for rw in file_list[2:-1]:
  for cl in rw:
    hess_list.append(float(cl))

## empty Hessian array to fill:
hess = np.zeros((dim,dim))
## index for hess_list:
i = 0
## only upper triangle of Hessian stored:
for rwi in range(dim):
  for cli in range(rwi+1):
    ## diagonal or off-diagonal:
    hess[rwi][cli] += hess_list[i]
    ## if off-diagonal, need to symmetrize:
    if rwi != cli:
      hess[cli][rwi] += hess_list[i]

    ## increment index for hess_list:
    i += 1

#print(hess)
## write fake hessian file (132.0):
out_hess_file = './fake_hess_132.0'
hess.tofile(out_hess_file)

