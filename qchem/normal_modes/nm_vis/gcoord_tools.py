#!/usr/bin/python

from __future__ import division
import sys
import numpy as np
from numpy import linalg as LA
import math
from copy import deepcopy as copy


########################################
def array_to_xyz(xyz, symbols, comment):

  print (len(symbols))
  print (str(comment))

  ## loop over atoms:
  for i, at in enumerate(xyz):
    print (symbols[i] + ' '*4 + \
     '   '.join(list(map(lambda x: '{:>10s}'.format(str('{:0.5f}'.format(x))), at))))

  return
########################################

####################
def distance(a1,a2):
  dist = LA.norm(a1 - a2)
  return dist
####################
####################

####################
def pt_coord(don,acc,htr):
  d_t = LA.norm(don - htr)
  a_t = LA.norm(acc - htr)
  return d_t - a_t
####################
####################

####################
def angle(a1,a2,a3):
  b1  = a1 - a2
  b3  = a3 - a2
  ang = math.acos( np.dot(b1,b3) / ( LA.norm(b1) * LA.norm(b3) ) )
  return ang
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

  return dihed
####################
####################

####################
def eval_hbond(R,b):
  '''
  see J. Chem Phys 126 204107 (2007) p6
  '''
  ## convert b (rad.) to b_deg (deg.):
  b_deg = 180 * (b / 3.1415926536897932383)

  ## evaluate hb criteria parabola at beta-angle, b:
  Rb_hb_parabola = 3.3 - 0.00044*(pow(b_deg,2))

  return bool(R < Rb_hb_parabola)
####################
####################




#!# previously from vector_tools.py:

#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
######################################
## AXIS-ANGLE ROTATION  ##############

def axis_angle(xyz,axis,angle):

  """ Uses Euler-Rodriguez matrix to rotate the xyz
     around axis, 'axis', by angle, 'angle'

     WARNING: angle in units of: radians/2pi (i.e. fraction of unit circle)

     see: https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula#Rotation_angle_and_rotation_axis
  """


  ## unit-circle circumference:
  tau = 6.2831853071795864769253

  angle = angle*tau

  ##########################
  ## setup Euler-Rodriguez:

  ## first normalize axis vector:
  axis = ( axis / LA.norm(axis) )

  ## need for cross-product formula:
  a = np.cos(angle/2)
  omega = copy(axis) * np.sin(angle/2)


  ##################################################
  ## now use cross-product form of Euler-Rodriguez:

  ## initialize numpy array for newly rotated xyz:
  nxyz = np.zeros(xyz.shape)

  ## loop over atom vectors in (molecule) xyz:
  for at in range(len(xyz)):
    tmp_term2 = 2*a*(np.cross(omega,xyz[at]))
    tmp_term3 = 2*(np.cross( omega,(np.cross(omega,xyz[at])) ))

    nxyz[at] = xyz[at] + tmp_term2 + tmp_term3

  ## return rotated xyz as a np array:
  return np.array(nxyz)


#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
######################################
## TRANSLATE COORDINATES  ############

def translate_xyz(xyz,vect):

  '''
  arguments:
  0.xyz == [[coords] atoms]; <np.array(2x2), floats>
  1.vect == displacememnt vector (for translation)
  '''

  nxyz = np.zeros((len(xyz),3))
  ## other approaches:
  #nxyz = [[[] for co in range(3)] for at in range(len(xyz))]
  #nxyz = [[] for at in range(len(xyz))]


  for at in range(len(xyz)):
    for co in range(len(xyz[at])):
      nxyz[at][co] = xyz[at][co] + vect[co]
      ## other approach:
      #nxyz[at].append((xyz[at][co] + vect[co]))

# test:
  #print "from translate"
  #print xyz
##########################

  return nxyz


#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
######################################
## INCREMENTAL ROTATION  #############

def inc_rot_mat(ax,angle):

  '''
  arguments:
  0.ax    == axis of rotation
  1.angle == angle of (single) rotation; i.e. an incremental angle...
           ... ang IN UNITS OF FRACTION OF UNIT CIRCLE !!!!!!
  '''

  ## unit-circle circumference:
  tau = 6.2831853071795864769253

  ## rotation cosine, & sine:
  r_c = np.cos(tau*angle)
  r_s = np.sin(tau*angle)


  #####################
  ## rotation matrices:

  ## (pos) z-axis rotation:
  if ax == 'z':
    rot_pz = np.array([[ r_c, r_s, 0.0],
                       [-r_s, r_c, 0.0],
                       [ 0.0, 0.0, 1.0]])
    return rot_pz
#---------------------------------------------

  ## (pos) y-axis rotation:
  elif ax == 'y':
    rot_py = np.array([[ r_c, 0.0,-r_s],
                       [ 0.0, 1.0, 0.0],
                       [ r_s, 0.0, r_c]])
    return rot_py
#---------------------------------------------

  ## (pos) x-axis rotation:
  elif ax == 'x':
    rot_px = np.array([[ 1.0, 0.0, 0.0],
                       [ 0.0, r_c, r_s],
                       [ 0.0,-r_s, r_c]])
    return rot_px
#---------------------------------------------

  ## input axis is a vector:
  elif type(ax) is np.array \
    and len(ax) == 3:

    ## first normalize axis vector:
    ax = ( ax / LA.norm(ax) )

    ## setup Euler-Rodriguez:
    a = np.cos(angle/2)
    b = ax[0]*(np.sin(angle/2))  # note: axis[:] = [k_x, k_y, k_z]
    c = ax[1]*(np.sin(angle/2))
    d = ax[2]*(np.sin(angle/2))

    a2 = a*a; b2 = b*b; c2 = c*c; d2 = d*d
    ab = a*b; ac = a*c; ad = a*d
    bc = b*c; bd = b*d
    cd = c*d

    ## matrix form of Euler-Rodriguez:
    rot_v = np.array([[a2+b2-c2-d2,   2*(bc+ad),   2*(bd-ac)],
                      [  2*(bc-ad), a2+c2-b2-d2,   2*(cd+ab)],
                      [  2*(bd+ac),   2*(cd-ab), a2+d2-b2-c2]])

    return rot_v

  else:
    print ("axis wrong... ")
    exit()


