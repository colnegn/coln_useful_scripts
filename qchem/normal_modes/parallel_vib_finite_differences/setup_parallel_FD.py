#!/usr/bin/env python

import sys
import numpy as np
import periodictable as pt
from copy import deepcopy as cp

import gxyz_io
import COMPAX

## seems to work without compax coords:
## go to: /home/ckegan/9-junction/junction_CuPc_in_Ag_4/7_17_2019_read_orbs_NEWER/f9.8/f9.8_optd_read_orbs/prepare_xyzs/z-example
## ... and run: ./do_nms_from_FD.py FD/for_FD_0_parallel_FD.xyz 4.5.6 FD/dlist.dlist; head -4615 ref_freq/output | tail -12
want_compax_bool = False

## step size used by Q-Chem:
ss = 0.001

if len(sys.argv) != 4:
  sys.exit('usage: 1.input xyz file; 2.string of (period-delimited) atom indices for normal mode calculation (indices starting with 1); 3.prefix for output files to write (suffix: "_parallel_FD.xyz")')

## read input xyz file into 'xyz' objects:
cfs = gxyz_io.read_xyz(sys.argv[1])

## atoms to displace in FD calc:
#FD_atoms = list(map(int, sys.argv[2].strip().split('.')))
FD_atoms = sorted(list(map(lambda x: int(x) - 1, sys.argv[2].strip().split('.'))))

## output name prefix:
out_pref = sys.argv[3]


## loop over configuration indices:
for cfi in range(len(cfs)):


  #!# do we want to put configs into COMPAX coordinates?
  if want_compax_bool:

    ## make temporary xyz for cfi with FD_atoms:
    tmp_xyz = cfs[cfi].fragment_xyz(permd_frag_inds=FD_atoms)

    ## first translate cfi into COM coordinates:
    com_crds = COMPAX.COM_shift(tmp_xyz.crds, tmp_xyz.masses)
    ## get COM from difference:
    com_cfi = tmp_xyz.crds[0] - com_crds[0]

    #!# get rotation matrix into principal axis coordinates:
    tmp_xyz.crds, rotmat_cfi = COMPAX.init_coord_frame(com_crds, tmp_xyz.masses, return_rotmat_bool=True)
    #tmp_xyz.crds, rotmat_cfi = COMPAX.init_coord_frame(tmp_xyz.crds, tmp_xyz.masses, return_rotmat_bool=True)

    #!# translate entire system (including frozens) to COMPAX:
    cfs[cfi].crds -= com_cfi
    cfs[cfi].crds = np.matmul(cfs[cfi].crds, rotmat_cfi)

  #!####################################################

  ## string to keep displaced xyzs for cfi:
  out_str = ''

  ## get xyz for minimum energy structure of cfi:
  out_str += cfs[cfi].write_xyz_file(return_flag=True)
  out_str += '\n'

  ## loop over atoms to displace for FD:
  for at in FD_atoms:
    ## loop over x,y,z dimensions for at:
    for xy in range(3):
      ## loop over positive and negative displacements:
      for posneg in [1.0, -1.0]:
        dvec = np.zeros(cfs[cfi].crds.shape)
        dvec[at][xy] += posneg * ss
        tmp_cf = cp(cfs[cfi])
        tmp_cf.crds += dvec

        ## get xyz file string:
        out_str += tmp_cf.write_xyz_file(return_flag=True)
        out_str += '\n'

  ## write xyz file for cfi:
  with open(out_pref + '_' + str(cfi) + '_parallel_FD.xyz', 'w') as f:
    f.write(out_str)


