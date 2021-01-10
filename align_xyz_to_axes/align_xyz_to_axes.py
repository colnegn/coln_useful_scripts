#!/usr/bin/env python

import sys
import numpy as np
from numpy import linalg as LA
from copy import deepcopy as copy

##################################################
###  array to xyz file  ##########################

def array_to_xyz(xyz, symbols, comment):

  print(len(symbols))
  print(str(comment))

  ## loop over atoms:
  for i, at in enumerate(xyz):
    print(symbols[i] + ' '*4 + \
     '   '.join(map(lambda x: '{:>14s}'.format(str('{:0.8f}'.format(x))), at)))

  return

#================================================#
#================================================#


##################################################
###  inputs n usage  #############################

max_nargs = 4

## need at least 2 arguments, but fewer than max_nargs:
if len(sys.argv) > (max_nargs + 1) or len(sys.argv) < (2 + 1):
  sys.exit('\n\n\tusage:\n\t\t'
           '1.input xyz file;\n\t\t'
           '2.index of atom to center at origin (INDICES START WITH 1 !!!!);\n\t\t'
           '3.(optional) index of atom to align to x-axis\n\t\t'
           '4.(optional) index of atom to align to z-axis\n\n')

in_xyz  = sys.argv[1]
orig_at = int(sys.argv[2]) - 1

## assume optional arguments are true:
do_x_align = True
do_z_align = True
## optional arguments:
try:
  x_ax_at = int(sys.argv[3]) - 1
except IndexError:
  do_x_align = False
try:
  z_ax_at = int(sys.argv[4]) - 1
except IndexError:
  do_z_align = False

#================================================#
#================================================#


##################################################
###  initialize xyz  #############################

#with open(in_xyz, 'r') as f:
#  xyz = np.array(map(lambda x: [float(at) for at in x.strip().split()[1:]], f.readlines()[2:]))
#
#with open(in_xyz, 'r') as f:
#  symbols = map(lambda x: str(x.strip().split()[0]), f.readlines()[2:])

with open(in_xyz, 'r') as f:
  in_xyz_file = f.readlines()

## number of atoms from top of xyz file:
natoms = int(in_xyz_file[0].strip())
## number of configs (assuming same number of atoms in each config):
nconfigs = int(len(in_xyz_file) / (natoms + 2))

## list to hold all xyz coordinates, atom symbols, and comment lines:
all_xyzs = []
all_symbols = []
all_comments = []
for cf in range(nconfigs):
  ## line-index for top of cf:
  topline = cf * (natoms + 2)
  ## lower- and upper-indices for xyz coordinates:
  lw = topline + 2
  up = lw + natoms
  all_xyzs.append(np.array(list(map(lambda x: [float(at) for at in x.strip().split()[1:]], in_xyz_file[lw:up]))))
  all_symbols.append(list(map(lambda x: str(x.strip().split()[0]), in_xyz_file[lw:up])))
  all_comments.append(in_xyz_file[topline + 1].strip())

#xyz     = np.array(map(lambda x: [float(at) for at in x.strip().split()[1:]], in_xyz_file[2:]))
#symbols = map(lambda x: str(x.strip().split()[0]), in_xyz_file[2:])
#comment = in_xyz_file[1].strip()

# test:
#print(xyz)
#print(symbols)
#print(comment)
#exit()

#================================================#
#================================================#


#!# now loop over configurations:
for cfi, xyz in enumerate(all_xyzs):
  symbols = all_symbols[cfi]
  comment = all_comments[cfi]

  ##################################################
  ###  translate to origin  ########################

  ## translate coordinates to center orig_at at origin:
  xyz -= xyz[orig_at]

  #================================================#
  #================================================#


  ## do we stop after translating to origin?
  if not do_x_align:
    ## write output xyz file to stdout:
    array_to_xyz(xyz, symbols, comment)
    continue


  ##################################################
  ###  align to x-axis  ############################

  ## want to align x_ax_at to x-axis:
  x_axis = np.array([1.,0.,0.])


  ## find vector normal to plane containing x-axis, orig_at, and x_ax_at:
  nvec = np.cross( (xyz[x_ax_at] / LA.norm(xyz[x_ax_at])), x_axis )


  ## angle (in radians) between a-axis and x_ax_at:
  ang = np.arcsin(LA.norm(nvec))


  ### prepare for Euler-Rodriguez rotation:
  '''
  see:
   https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula#Definition

  in particular:
   https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula#Vector_formulation
   https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula#Rotation_angle_and_rotation_axis
  '''

  if np.sign(xyz[x_ax_at][0]) < 0:
    ang += np.pi

  ## "a" factor:
  a = np.cos(ang / 2)

  ## if 0-rotation:
  if LA.norm(nvec) == 0.0:
    skip_rotate_x_bool = True
  else:
    ## "omega" vector:
    omega = np.sign(xyz[x_ax_at][0]) * np.sin(ang / 2) * (nvec / LA.norm(nvec))
    skip_rotate_x_bool = False


  ### now do rotation:

  ## 0-rotation (already on x-axis):
  if skip_rotate_x_bool:
    nxyz = copy(xyz)

  ## do rotation (not already on x-axis):
  else:
    ## declare new xyz array:
    nxyz = np.zeros(xyz.shape)

    ## loop over atoms, & rotate each via cross-product formula for Euler-Rodriguez:
    for at in range(len(xyz)):
      t1 = 2 * a * np.cross(omega, xyz[at])
      t2 = 2     * np.cross(omega, np.cross(omega, xyz[at]))

      nxyz[at] = xyz[at] + t1 + t2

  #================================================#
  #================================================#


  ## do we stop after aligning x-axis?
  if not do_z_align:
    ## write output xyz file to stdout:
    array_to_xyz(nxyz, symbols, comment)
    continue


  ##################################################
  ###  align to z-axis  ############################

  #ang2 = np.arctan((xyz[z_ax_at][1] - xyz[x_ax_at][1]) / (xyz[z_ax_at][2] - xyz[x_ax_at][2]))
  ang2 = np.arccos(nxyz[z_ax_at][2] / LA.norm(np.array([nxyz[z_ax_at][1], nxyz[z_ax_at][2]]))) # - xyz[x_ax_at]))
  #if np.sign(nxyz[z_ax_at][1]) < 0:
  #  ang2 += np.pi

  #print(cfi, 'UNDER ERROR')
  #sys.stdout.flush()


  a2 = np.cos(ang2 / 2)
  omega2 = np.sign(nxyz[z_ax_at][1]) * np.sin(ang2 / 2) * np.array([1.0, 0.0, 0.0])

  nnxyz = np.zeros(nxyz.shape)
  for at in range(len(nxyz)):
    t12 = 2 * a2 * np.cross(omega2, nxyz[at])
    t22 = 2      * np.cross(omega2, np.cross(omega2, nxyz[at]))
    nnxyz[at] = nxyz[at] + t12 + t22

  #================================================#
  #================================================#


  ## write output xyz file to stdout:
  array_to_xyz(nnxyz, symbols, comment)


