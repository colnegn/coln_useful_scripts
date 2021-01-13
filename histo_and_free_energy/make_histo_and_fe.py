#!/usr/bin/env python

import sys
import math
import numpy as np
from copy import deepcopy as copy


##################################################
###  only keep bins with enough sampling  ########

def enough_sampling(h, binedges):
  '''
  maybe make another option to return smoothed few-samples
  bins in separate array to make a dotted line (secondary) plot thing
  '''

  #!# HARD-CODED minimum number of samples (going, coming back once is not enough, so want more than 2 samples):
  min_num_samples = 3
  ## list of bin indices
  fewsamples_inds = [xi for xi, x in enumerate(h) if x < min_num_samples]

  #!# if enough sampling:
  if len(fewsamples_inds) == 0:
    return h, binedges

  #!# find islands of too few samples:
  #!# HARD-CODED minimum separation between few-samples bins:
  min_few_sep = 5  #!#
  ## list of few-samples island bouds (+1 for upper bound):
  fewsamples_island_bounds = [[fewsamples_inds[0], fewsamples_inds[0]+1]]
  for fxi, fx in enumerate(fewsamples_inds[:-1]):
    ## see if fewsamples bin inds are separated
    if abs(fx - fewsamples_inds[fxi+1]) < min_few_sep:
      ## extend outer-bound of most recent island:
      fewsamples_island_bounds[-1][1] = fewsamples_inds[fxi+1] + 1
    ## else, next fx starts a new island:
    else:
      new_island_bounds = [fewsamples_inds[fxi+1], fewsamples_inds[fxi+1]+1]
      fewsamples_island_bounds.append(copy(new_island_bounds))


  ## now make histogram only for bins with enough sampling:
  to_use_inds  = []
  ## indices lower than lowest-indexed island (if exists):
  to_use_inds += list(range(0, fewsamples_island_bounds[0][0]))
  ## loop over island bounds:
  for ibi, ib in enumerate(fewsamples_island_bounds[:-1]):
    to_use_inds += list(range(ib[1]+1, fewsamples_island_bounds[ibi+1][0]))
  ## indices higher than highest-indexed island (if exists):
  to_use_inds += list(range(fewsamples_island_bounds[-1][1]+1, nbins))


  ## new histogram with enough sampling:
  out_h = h[to_use_inds]
  out_bindeges = binedges[to_use_inds]

  return out_h, out_bindeges

#================================================#
#================================================#



#print('\n\nTAKING ABSOLUTE VALUE OF DATA\n\n', file=sys.stderr)
print('\n\nALSO NORMALIZING DATA\n\n', file=sys.stderr)

if len(sys.argv) != 7:
  sys.exit('usage: 1.data file; 2.column index to bin (indices start with 1); 3.number of bins; 4.output histogram file (to overwrite); 5.output FE file; 6.bool to symmetrize axis (1~do symmetrization, 0~NO symmetrization)')

d = np.genfromtxt(sys.argv[1])
icol = int(sys.argv[2]) - 1
nbins = int(sys.argv[3])
outhistf = sys.argv[4]
outfef = sys.argv[5]
sym_bool = bool(int(sys.argv[6]))


## only need one column:
d = d[:,icol]
nsamples = float(len(d))

if sym_bool:
  d = np.concatenate((d, -1. * d))


## make histogram:
h, tmpedges = np.histogram(d, nbins)

binedges = []
for i in range(len(tmpedges) - 1):
  binedges.append((tmpedges[i] + tmpedges[i+1]) / 2.0)
binedges = np.asarray(binedges)

## make sure our histogram has enough sampling:
h, binedges = enough_sampling(h, binedges)


h_normd = [float(x) / nsamples for x in h] #!# normalization

with open(outhistf, 'w') as f:
  for di in range(len(h)):
    f.write('{0:>12.6f}{1:>12.6f}\n'.format(binedges[di], h_normd[di]))



#!# also do FE:
kelvin_to_kcalmol = 0.00198719
#!# assume RT:
kT = 298.15 * kelvin_to_kcalmol

FE = np.zeros((len(h_normd),2))
FE[:,0] += binedges
FE[:,1] += np.asarray([-1.0 * kT * math.log(x) for x in h_normd])
## zero:
FE[:,1] -= min(FE[:,1])

## save FE:
np.savetxt(outfef, FE)


