#!/usr/bin/env python

import sys
import math
import numpy as np
from numpy import linalg as LA


plot_bool = False


####################
def distance(a1,a2):
  dist = LA.norm(a1 - a2)
  return dist
####################
####################

####################
def angle(a1,a2,a3):
  b1  = a1 - a2
  b3  = a3 - a2
  ang = math.acos( np.dot(b1,b3) / ( LA.norm(b1) * LA.norm(b3) ) )

  return (180. / np.pi) * ang
  #return ang
####################
####################

####################
def dihedral(a1,a2,a3,a4):

  b1 = a2 - a1
  b2 = a3 - a2
  b3 = a4 - a3

  n1 = np.cross(b1, b2) / LA.norm(np.cross(b1, b2))
  n2 = np.cross(b2, b3) / LA.norm(np.cross(b2, b3))

  m1 = np.cross(n1, (b2 / LA.norm(b2)))

  x = np.dot(n1, n2)
  y = np.dot(m1, n2)


  dihed = np.arctan2(y, x)

  return (180. / np.pi) * dihed
  #return dihed
####################
####################

####################
def pt_coord(o1, o2, h):
  return LA.norm(o1 - h) - LA.norm(o2 - h)
####################
####################

###
def bisector_dih(c1,f11,f12,c2,f21,f22):

  c1_f11 = f11 - c1
  c1_f12 = f12 - c1

  c1_f11 /= LA.norm(c1_f11)
  c1_f12 /= LA.norm(c1_f12)

  b1 = (c1_f11 + c1_f12) / 2.0
  b1 += c1

  c2_f21 = f21 - c2
  c2_f22 = f22 - c2
  c2_f21 /= LA.norm(c2_f21)
  c2_f22 /= LA.norm(c2_f22)

  b2 = (c2_f21 + c2_f22) / 2.0
  b2 += c2


  return dihedral(b1,c1,c2,b2)

###


####################
def lin_combo_f(a1, v1, a2, v2):
  return (a1 * v1) + (a2 * v2)
####################
####################


########################################
def array_to_xyz(xyz, symbols, comment):

  print(len(symbols))
  print(str(comment))

  ## loop over atoms:
  for i, at in enumerate(xyz):
    print(symbols[i] + ' '*4 + \
     '   '.join(list(map(lambda x: '{:>10s}'.format(str('{:0.9f}'.format(x))), at))))

  return
########################################


if len(sys.argv) != 3:
  sys.exit("usage: 1.input xyz file; 2.input coord-parameter file (atom indices start with 1; "
           "'L' ~ distance; 'A' ~ angle; 'D' ~ dihedral; 'P' ~ PT-coord (o1, o2, h); 'B' ~ bisector dihedral (center1, flanking11, f12, c2, f21, f22))")

xyz_file_path  = sys.argv[1]
internals_path = sys.argv[2]

## convert xyz file into list (rows) of lists (columns):
with open(xyz_file_path) as f:
  in_xyz = list(map(lambda x: x.strip().split(), f.readlines()))

## convert internals file into list of lists:
with open(internals_path) as f:
  internals = [li.strip().split() for li in f.readlines()
                                  if '#' not in li and len(li.strip()) != 0]
#                 ignore comment lines ^^^          ... and empty lines ^^^


#################
### write header:

## header string:
#out_head = 'index;' + ';'.join(['.'.join(i)
#                      for ii, i in enumerate(internals)])

head_max_len = 50
line_len = 0
out_str = '#! COLUMNS:  '
### loop over fields in header:
#for fii, fi in enumerate([i for i, ch in enumerate(out_head)
#                          if ch == ';']):
for fii, fi in enumerate(['cf'] + ['.'.join(i) for i in internals]):
  ## if about to exceed maximum header-line length:
  if (line_len + len(fi)) > head_max_len:
    ## print current line:
    print(out_str.strip())
    ## adjust max_len:
    line_len = 0
    ## reset out_str:
    out_str = '#!           '

  ## add field, fi, to out_str:
  out_str += str(fii + 1) + ') ' + fi + ';  '
  ## add length for field, fi, to line_len:
  line_len += len(fi) + len(str(fii + 1)) + 5

## print the last line of header to stdout:
print(out_str.strip())

#################



## keep needed atom indices in a set:
at_inds_needed = set()
## convert atom indices in internals into ints:
for i in range(len(internals)):
  for j in range(1, len(internals[i])):
    internals[i][j] = int(internals[i][j]) - 1  #!# forcing indices to start with 1
    at_inds_needed.add(internals[i][j])


## dictionary of coord-functions:
cfunctions = {'L': distance, 'A': angle, 'D': dihedral, 'P': pt_coord, 'B': bisector_dih}


### loop through all configurations from xyz file:

## dict to hold np.arrays of coordinates for needed atoms:
at_dict = {str(at): None for at in at_inds_needed}
## keep internal coordinates in list of lists:
crd_list = []
## line (row) counter for while-loop:
li_count = 0
## configuration counter:
conf_count = 0
while li_count < len(in_xyz):

  ## new crd_list sublist for each configuration:
  crd_list.append([conf_count])

  ## new out_str for each configuration (first column is configuration index):
  out_str = '{:<5s}'.format(str(conf_count))

  ## get number of atoms for current config:
  n_atoms = int(in_xyz[li_count][0])


  ## get needed coordinates:
  for at in at_dict:
    #print at, np.array(map(float, in_xyz[li_count + 2 + int(at)][1:]))
    at_dict[at] = np.array(list(map(float, in_xyz[li_count + 2 + int(at)][1:])))

  ## now loop over internal coordinates:
  for i in internals:
    atoms_for_i = []
    for at in i[1:]:
      atoms_for_i.append(at_dict[str(at)])

    crd_val = cfunctions[i[0]](*atoms_for_i)


    crd_list[-1].append(crd_val)
    #out_str += '   ' + '{:0.4f}'.format(crd_val)
    out_str += '{:>16s}'.format('{:0.9f}'.format(crd_val))
    #out_str += '   ' + '{:0.4f}'.format(cfunctions[i[0]](*atoms_for_i))

  print(out_str)

#  ###################################################################
#  #!!# print xyz for configs with both O-H lengths smaller than 0.98:
#  if all(oh_l < 0.98 for oh_l in crd_list[-1][1:]):
#    print '\n'.join(['   '.join(li) for li in in_xyz[li_count:(li_count + n_atoms + 2)]])
#  ###################################################################


  ## increment li_count:
  li_count += (n_atoms + 2)
  ## increment conf_count:
  conf_count += 1


