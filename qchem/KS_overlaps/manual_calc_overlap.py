#!/usr/bin/env python

import sys
import numpy as np

###########################
def calc_overlap(p_i, p_j):
  ### NOTE: integral of {(4*pi) * (r^2) exp(-b_i * (r^2)) dr} from 0 to inf = (pi / b)^(3/2)

  ## value of integral:
  S = 0.0

  ## want to compute < \psi_i | \psi_j >, so double-loop:
  for b_c_i in p_i:
    for b_c_j in p_j:

      ## overlap between i and j:
      s_ij = (b_c_i[1] * b_c_j[1]) * pow(np.pi / (b_c_i[0] + b_c_j[0]), 1.5)

      S += s_ij

  return S
###########################




## parameters for contracted orbital i, indexed by gaussian: [0] ~ exponent (b_i); [1] ~ coefficient (c_i)
params_i = np.genfromtxt(sys.argv[1])
params_j = np.genfromtxt(sys.argv[2])

## make sure that parameter arrays are 2-index:
try:
  dummy = params_i.shape[1]
except IndexError:
  params_i = np.array([params_i])

try:
  dummy = params_j.shape[1]
except IndexError:
  params_j = np.array([params_j])



## calculate normalization
N_i = 1.0 / pow(calc_overlap(params_i, params_i), 0.5)
N_j = 1.0 / pow(calc_overlap(params_j, params_j), 0.5)

## overlap between i and j:
S_ij = calc_overlap(params_i, params_j) * N_i * N_j


print S_ij


