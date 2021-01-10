#!/usr/bin/env python

import sys
import gxyz_io as xyz_io
import numpy as np
from copy import deepcopy as copy

if len(sys.argv) != 6:
  sys.exit('usage: 1.xyz file; 2.normal mode disp file (corresponding to xyz); 3.output xyz file to overwrite; 4.number of output configs to generate; 5.maximum displacement to generate')


## read xyz file into 'xyz' object:
in_xyz = xyz_io.read_xyz(sys.argv[1])[0]
orig_xyz = copy(in_xyz)

## read dvec as a 2-index array:
d = np.genfromtxt(sys.argv[2]).reshape((-1,3))
N = int(sys.argv[4])

## maximum displacement (and velocity of distortion):
max_disp = float(sys.argv[5])


#!# HARD-CODED do we want sine-displacements?
sine_bool = True
#sine_bool = False

with open(sys.argv[3], 'w') as f:

  if sine_bool:
    for i in list(np.linspace(0, (2*np.pi), N)):
      c = max_disp * np.sin(i)
      in_xyz.crds = orig_xyz.crds + (c * d)
      in_xyz.write_xyz_file(out_file=f, new_comment='{:0.3f}'.format(c), append_comment=False)

  else:
    for i in list(np.linspace(-1, 1, N)):
      c = max_disp * i
      in_xyz.crds = orig_xyz.crds + (c * d)
      in_xyz.write_xyz_file(out_file=f, new_comment='{:0.3f}'.format(c), append_comment=False)

