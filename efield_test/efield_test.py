#!/usr/bin/env python

import sys
import numpy as np
from numpy import linalg as LA



##################################################
###  init unit conversions  ######################

def init_unit_conversions(SI_bool=True):

  if SI_bool:
    ## permittivity of free space:
    eps0 = 8.8541878128e-12  #  ~ F/m = (s^2 C^2) / (kg m^3)
#        = 8987551792.26118

    ## Coulomb constant:
    k = 1. / (4. * np.pi * eps0)

    ## Coulomb/e
    charge_conv = 1.602176634e-19

    ## meter/Angstrom
    dist_conv = 1.0e-10

    ## polarizability conversion
    alpha_conv = pow(dist_conv, 3) / k  #  ~ (F m^2) / Ang^3


  ## else use internal units:
  else:
    k = 1.
    charge_conv = 1.  #  ~ e
    dist_conv = 1.    #  ~ Ang.
    alpha_conv = 1.   #  ~ Ang.^3


  return k, charge_conv, dist_conv, alpha_conv

#================================================#
#================================================#




##################################################
###  Get Permanent Field  ########################

def GetE0(crds, charges, k, debug_bool=False):

  N   = len(crds)
  dim = len(crds[0])
  ## initialize E0 vector:
  E0 = np.zeros(dim * N)

  ## loop over pairs of sites:
  for s1i, s1 in enumerate(crds):
    for s2i, s2 in enumerate(crds):
      if s1i == s2i:
        continue

      ## take difference vector:
      d12 = s1 - s2
      ## distance:
      R12 = LA.norm(d12)

      ## add contribution of E0 at s1 due to charge at s2:
      E0[(s1i * dim):(s1i * dim)+3] += charges[s2i] * d12 / pow(R12, 3)


  ## 1/(4 pi eps0)
  E0 *= k

  if debug_bool:
    print('\n\nE0:\n', E0)


  return E0

#================================================#
#================================================#


##################################################
###  Get 2nd Rank Interaction Tensor, T2  ########

def GetT2(crds, k, debug_bool=False):

  N   = len(crds)
  dim = len(crds[0])
  ## initialize T2:
  T2 = np.zeros((dim * N, dim * N))

  ## keep (dim*N)-dimensional intermediate matrices for debugging:
  RaRb = np.zeros((dim * N, dim * N))
  diag = np.zeros((dim * N, dim * N))

  ## loop over pairs of sites:
  for s1i, s1 in enumerate(crds):
    ## lower indices for (dim*N)-dimensional tensors:
    lwri = s1i * dim
    for s2i, s2 in enumerate(crds):
      lwci = s2i * dim
      if s1i == s2i:
        continue

      ## take difference vector:
      d12 = s1 - s2
      ## distance:
      R12 = LA.norm(d12)

      RaRb[lwri:lwri+3, lwci:lwci+3] += 3. * np.outer(d12,d12)
      diag[lwri:lwri+3, lwci:lwci+3] += np.diag(np.repeat(pow(R12,2), dim))

      t12 = pow(R12, -5) * (RaRb[lwri:lwri+3, lwci:lwci+3] - diag[lwri:lwri+3, lwci:lwci+3])

      #RaRb = 3. * np.outer(d12,d12)
      #diag = np.diag(np.repeat(pow(R12,2), dim))
      #t12 = pow(R12, -5) * (RaRb - diag)

      ## add contribution of T2 at s1 due to site s2:
      T2[lwri:lwri+3, lwci:lwci+3] += t12


  ## 1/(4 pi eps0)
  T2 *= k


  if debug_bool:
    print('\n\nT2:\n', T2)
    print('\n\nRaRb:\n', RaRb)
    print('\n\ndiag:\n', diag)


  return T2

#================================================#
#================================================#




##################################################
###  MAIN  #######################################

### MAIN:
def main(incrds, incharges, inalphas, SI_bool=True, debug_bool=False):

  k, charge_conv, dist_conv, alpha_conv = init_unit_conversions(SI_bool)



  charges = np.genfromtxt(incharges)
  charges *= charge_conv

  ## cartesian coordinates of sites:
  crds = np.genfromtxt(incrds)
  crds *= dist_conv


  if debug_bool:
    print('\n\ncharges:\n', charges)
    print('\n\ncoordinates:\n', crds)


  N   = len(crds)
  dim = len(crds[0])



  isotropic_bool = True

  ### make polarizability matrix for isotropic sites:
  isotropicalphas = np.genfromtxt(inalphas)

  ## (dim * N)-dimensional, square matrices for polarizabilities:
  alphas_dN = np.zeros((dim * N, dim * N))
  invalphas_dN = np.zeros((dim * N, dim * N))
  for Ni in range(N):
    for di in range(dim):
      rwi = (Ni * dim) + di
      if isotropic_bool:
        ## isotropic => diagonal:
        cli = rwi
        alphas_dN[rwi][cli] += alpha_conv * isotropicalphas[Ni]
        ## inverse of diagonal matrix:
        invalphas_dN[rwi][cli] += 1. / (alpha_conv * isotropicalphas[Ni])
      else:
        sys.exit('ERROR: only isotropic polarizabilities')


  if debug_bool:
    print('\n\nalpha matrix:\n', alphas_dN)
    print('\n\ninverse alpha matrix:\n', invalphas_dN)



  ## get (dim * N)-dimensional vector of electric field due to permanent point charges at each site:
  E0_dN = GetE0(crds, charges, k, debug_bool)

  ## get 2nd rank interaction tensor:
  T2 = GetT2(crds, k, debug_bool)


###############
#!# BACKGROUND:
#
# + Mu is (dim * N)-dimensional vector of dipole components:
#     Mu = [mu1_x, mu1_y, mu1_z, mu2_x, ... muN_y, muN_z]^T
#
# + E0 is (dim * N)-dimensional vector of total permanent electric field components at each site:
#     E0 = [E0_tot_x(r1), E0_tot_y(r1), E0_tot_z(r1), E0_tot_x(r2), ... E0_tot_z(rN)]^T
#
# + alphas is (dim * N) x (dim * N) block-diagonal matrix of dipole-polarizability tensors:
#     alphas = [[{a1}, {0}, {0}, ... {0}],
#               [{0}, {a2}, {0}, ... {0}],     with: {a1} = [[a1_xx, a1_xy, a1_xz],
#               [ .         .           ],                   [a1_yx, a1_yy, a1_yz],
#               [ .            .        ],                   [a1_zx, a1_zy, a1_zz]]
#               [ .               .     ],
#               [{0}, {0}, ...      {aN}]]     and with {0} = 3x3 zero matrix
#
# + T2 is (dim * N) x (dim * N) 2nd-rank interaction tensor:
#     T2 = [[{0}, {t12}, {t13}, ... {t1N}],
#           [{t21}, {0}, {t23}, ... {t2N}],     with: {t12} = 3x3 tensor between sites 1 & 2
#           [ .           .              ],     such that:
#           [ .                .         ],         {t12} [mu2_x, mu2_y, mu2_z]^T
#           [ .                     .    ],       = [Eind_x(r1), Eind_y(r1), Eind_z(r1)]^T
#           [{tN1}, {tN2}, ...        {0}]]
#
#
#!# DERIVATION:
#
#   Mu = alphas (E0 + (T2 Mu))
#   (alphas^-1 Mu) = E0 + (T2 Mu)
#   (alphas^-1 Mu) - (T2 Mu) = E0
#   Mu = (alphas^-1 - T2)^-1 E0
#
###############

  ## build matrix to invert:
  A = invalphas_dN - T2

  ## invert A:
  invA = LA.inv(A)

  if debug_bool:
    ## test matrix inversion:
    print('\n\ntesting AA^-1:\n', np.dot(A,invA))


  ## get (dim * N)-dimensional vector of induced dipoles at each site:
  Mu_dN = np.dot(invA, E0_dN)

  if debug_bool:
    print('\n\ninduced dipole:\n', Mu_dN)


  ## get (dim * N)-dimensional vector of electric field due to induced dipoles at each site:
  Eind_dN = np.dot(T2, Mu_dN)



  ## finally add permanent and induced electric fields:
  Etot_dN = E0_dN + Eind_dN


  ## print(in N x 3 format:)
  print('\n\npermanent E-field:\n', E0_dN.reshape(N, dim))
  print('\n\ninduced E-field:\n', Eind_dN.reshape(N, dim))
  print('\n\ntotal E-field:\n', Etot_dN.reshape(N, dim))


  return E0_dN.reshape(N, dim), Eind_dN.reshape(N, dim), Etot_dN.reshape(N, dim)

#================================================#
#================================================#


##################################################
###  usage  ######################################

def usage():

  if len(sys.argv) != 4:
    sys.exit('usage: 1.file containing coordinates (no atom labels); 2.file containing charges; 3.file containing isotropic polarizabilities (1 value per site)')

  incrds = sys.argv[1]
  incharges = sys.argv[2]
  inalphas = sys.argv[3]


  return incrds, incharges, inalphas

#================================================#
#================================================#





if __name__ == '__main__':

  ## parse input arguments:
  incrds, incharges, inalphas = usage()


  #!# HARD-CODED:
  #debug_bool = True
  debug_bool = False
  SI_bool = False
  #SI_bool = True


  E0_dN, Eind_dN, Etot_dN = main(incrds, incharges, inalphas, SI_bool, debug_bool)


  if debug_bool:
    E0_dN_SI, Eind_dN_SI, Etot_dN_SI = main(incrds, incharges, inalphas, True, False)
    E0_dN_int, Eind_dN_int, Etot_dN_int = main(incrds, incharges, inalphas, False, False)

    print('\n\nE0 SI/internal:\n', E0_dN_SI[:,0] / E0_dN_int[:,0])
    print('\n\nEind SI/internal:\n', Eind_dN_SI[:,0] / Eind_dN_int[:,0])




