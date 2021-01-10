## imports
from __future__ import division

from copy import deepcopy as copy
import math
import numpy as np
from numpy import linalg as LA



#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
######################################
## COM SHIFT  ########################

def COM_shift(xyz,masses):

  ## mass-weighted vector sum (of coordinates):
  mw_vsum = np.zeros(3)
  ## denominator just the sum of masses (a scalar)
  denom = 0.0

  ## loop over atoms:
  for at in range(len(xyz)):
    denom   += masses[at]
    mw_vsum += masses[at] * xyz[at]

  ## center of mass vector:
  COM = (mw_vsum / denom)

  ## translate xyz coords to COM-coordinate frame
  xyz_COM = xyz - COM


# test:
  #print "\nxyz type: ",type(xyz), "\nxyz:\n",xyz
  #print "\nxyz_COM type: ",type(xyz_COM), "\nxyz_COM:\n",xyz_COM
#######

  return xyz_COM

#===================================#
#===================================#


#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
######################################
## ROTATE TO PRINCIPAL AXIS FRAME  ###

def P_AX_rotate(crds,masses,return_rotmat_bool=False):


# test:
  #print "crds before P_AX:"
  #for at in range(len(crds)):
    #print crds[at][0]
    #print crds[at][1]
    #print crds[at][2]
    #print "\n"

  #print "masses:"
  #for at in range(len(masses)):
    #print masses[at]
#######


  xx,yy,zz = 0.0,0.0,0.0
  xy,xz,yz = 0.0,0.0,0.0

  for at in range(len(crds)):
    xx += math.pow(crds[at][0],2)*masses[at]
    yy += math.pow(crds[at][1],2)*masses[at]
    zz += math.pow(crds[at][2],2)*masses[at]

    xy += crds[at][0]*crds[at][1]*masses[at]
    xz += crds[at][0]*crds[at][2]*masses[at]
    yz += crds[at][1]*crds[at][2]*masses[at]

  mi_tensor = np.array([[(yy + zz), -xy, -xz],
                        [-xy, (xx + zz), -yz],
                        [-xz, -yz, (xx + yy)]])

# test:
  #print "\n\nmi_tensor:"
  #print mi_tensor
#######


  ## eigh lists evals in ascending order; i.e. principal axis ~ z-axis
  #Ivals, Ivects = LA.eigh(mi_tensor)
  Ivals, Ivects = LA.eig(mi_tensor)

  crds_p_ax = np.matmul(crds,Ivects)
#  crds_p_ax2 =np.matmul(crds_p_ax,Ivects)


# test:
  #print "\n\nIvects:\n", Ivects
  #print "Ivals:\n", Ivals
  #print "crds\n", crds
  #print "crds_p_ax:\n", crds_p_ax, "\n"
#######


  ## do we want to return rotation matrix?
  if return_rotmat_bool:
    return crds_p_ax, Ivects

  else:
    return crds_p_ax

#===================================#
#===================================#



#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
######################################
## INITIALIZE MONOMER COORDINATES  ###

def init_coord_frame(xyz,masses,return_rotmat_bool=False):
  '''
  arguments:
    + xyz:    Nx3 np.array of coordinates
    + masses: list of masses in same order as xyz

  returns coordinates of system in center of mass/principal
  axes (COMPAX) coordinate frame

  also returns PAX rotation matrix if
  return_rotmat_bool == True
  '''

  ## translate molecule to center of mass (COM) coords:
  xyz = COM_shift(copy(xyz),masses)

  ## return COMPAX coords, and rotmat (if bool == True):
  return P_AX_rotate(xyz,masses,return_rotmat_bool)

#===================================#
#===================================#



#######################################
### test code:

#f_xyz = np.array([[2.4,2.2,1.8],
#                  [5.1,2.5,0.0],
#                  [1.1,3.2,1.6]])
#
#f_xyz = np.array([[ 1.4, 0.0, 0.0],
#                  [-0.6, 0.0, 0.0]])
#
#print f_xyz, '\n\n'
#
#f_masses = [1.007, 1.007] #, 2.014]
#
#print init_coord_frame(f_xyz,f_masses)
#
#print '\n\n', f_xyz


