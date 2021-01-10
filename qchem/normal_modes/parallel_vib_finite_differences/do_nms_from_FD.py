#!/usr/bin/env python

import sys
import numpy as np
from numpy import linalg as LA
import periodictable as pt
from copy import deepcopy as copy

import gxyz_io
import COMPAX

np.set_printoptions(7, linewidth=120, suppress=True)

### speed of light in natural units (~ hbar/E_h):
#c = 137.035999
### (~ Angstroms / Bohr):
bohr2ang = 0.52917721067

#!# from QCElemental (NIST 2014):
## Avogadro's constant (~ 1/mol):
N_A = 6.022140857e23  # 'avogadro constant'
## (~ Joule / Hartree):
hartree2joule = 4.359744650e-18  # 'hartree-joule relationship'
## speed of light in Hz (?):
c =  299792458  # 'inverse meter-hertz relationship'
## (~ meter / Bohr):
bohr2meter = 0.52917721067e-10 # 'bohr radius'

#!# from psi4 (see: psi4/psi4/driver/qcdb/vib.py)
## unit conversion for Hessian eigenvalues in a.u.: E_h/((a_0)^2) --> cm^-1:
uconv_cm_1 = (np.sqrt(N_A * hartree2joule * 1.0e19) / (2 * np.pi * c * bohr2ang))  # psi4
#            1.0e19 is (Ang / a_0)^2 * 1.e3  ^^^

##################################################
###  PREPARATION  ################################

def prepare_xyz(arg1, arg2):

  ## atoms to displace in FD calc:
  FD_atom_inds = sorted(list(map(lambda x: int(x) - 1, arg2.strip().split('.'))))

  ## read input xyz file into 'xyz' object (including frozens):
  full_FD_cfs = gxyz_io.read_xyz(arg1)
  ## get fragment of optd config without frozen atoms:
  optd_cf = full_FD_cfs[0].fragment_xyz(permd_frag_inds=FD_atom_inds)
  ## get fragments of FD (distortion) configs without frozen atoms:
  FD_cfs = [cf.fragment_xyz(permd_frag_inds=FD_atom_inds) for cf in full_FD_cfs[1:]]

  ## get FD displacement vectors in same order as FD_cfs:
  FD_disps = [cf.crds.flatten() - optd_cf.crds.flatten() for cf in FD_cfs]
  #!# make sure all displacements same magnitude:
  for dii, di in enumerate(FD_disps):
    if abs(LA.norm(di) - LA.norm(FD_disps[0])) > pow(10,-8):  #!#
      print(dii, LA.norm(di), '!=', LA.norm(di[0]))
      raise Exception('non-equivalent displacements')
  ## set displacement (distortion) step size IN ANGSTROMS:
  ss = LA.norm(FD_disps[0])
  ##!# convert ss into natural units (Bohr):
  #ss *= 1.88973
  #print(ss)
  #exit()

  ## get indices of displaced (distorted) coordinate for each FD_disp (same order as FD_disp):
  coi_FD_disps = [np.argsort(abs(di))[-1] for di in FD_disps]
  #*# get order of displacements in order of xyz atom-coordinate-order:
  pos_FD_disp_order = []
  neg_FD_disp_order = []
  #!# indices????????????????
  for dii in np.argsort(coi_FD_disps):
    if FD_disps[dii][coi_FD_disps[dii]] > 0.:
      #!# NEED +1 because we don't want to include optd_cf:
      pos_FD_disp_order.append(dii + 1)
    elif FD_disps[dii][coi_FD_disps[dii]] < 0.:
      #!# NEED +1 because we don't want to include optd_cf:
      neg_FD_disp_order.append(dii + 1)
  FD_disp_order = {'pos':pos_FD_disp_order, 'neg':neg_FD_disp_order}
  #old:FD_disp_order = np.argsort(coi_FD_disps)


  #!# MATCH Q-CHEM ????!!!!!????!!!!
  optd_cf.masses = np.array([15.99491, 1.00783, 1.00783])


  #!# translate cf into COM coordinates:
  optd_cf.crds = COMPAX.COM_shift(optd_cf.crds, optd_cf.masses)

  #!# get eigenvectors of moment of inertia tensor, X:
  compax_cf_crds, X = COMPAX.init_coord_frame(optd_cf.crds, optd_cf.masses, return_rotmat_bool=True)
# ^^^ don't need compax cf


  return FD_atom_inds, optd_cf, X, ss, FD_disp_order

#================================================#
#================================================#


##################################################
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
###  TRANSLATIONS, ROTATIONS FROM PSI4  ##########
'''
NOTE: THIS WAS TAKEN DIRECTLY FROM PSI4 vib.py
look here: psi4/psi4/driver/qcdb/vib.py
'''
def _get_TR_space(m, geom, space='TR', tol=None, verbose=1):
    """Form the idealized translation and rotation dof from geometry `geom` and masses `m`.
    Remove any linear dependencies and return an array of shape (3, 3) for atoms, (5, 3 * nat) for linear `geom`,
    or (6, 3 * nat) otherwise. To handle noisy linear geometries, pass `tol` on the order of max deviation.

    m1 = np.asarray([1.])
    m2 = np.asarray([1., 1.])
    m3 = np.asarray([1., 1., 1.])
    m4 = np.asarray([1., 1., 1., 1.])
    g4 = np.asarray([[ 1.,  1., 0.],
          [-1.,  1., 0.],
          [-1., -1., 0.],
          [ 1., -1., 0.]])
    g2 = np.asarray([[ 1.,  1., 0.],
          [-1., -1., 0.]])
    g3 = np.asarray([[3., 3., 3.],
          [4., 4., 4.,],
          [5., 5., 5.]])
    g3noisy = np.asarray([[3., 3.001, 3.],
          [4., 4.001, 4.,],
          [5., 5., 5.01]])
    g33 = np.asarray([[0., 0., 0.],
                     [1., 0., 0.],
                     [-1., 0., 0.]])
    g1 = np.asarray([[0., 0., 0.]])
    g11 = np.asarray([[1., 2., 3.]])
    noise = np.random.normal(0, 1, 9).reshape(3, 3)
    noise = np.divide(noise, np.max(noise))

    assert(_get_TR_space(m4, g4).shape == (6, 12))
    assert(_get_TR_space(m2, g2).shape == (5, 6))
    assert(_get_TR_space(m3, g3).shape == (5, 9))
    assert(_get_TR_space(m3, g33).shape == (5, 9))
    assert(_get_TR_space(m1, g1).shape == (3, 3))
    assert(_get_TR_space(m1, g11).shape == (3, 3))
    assert(_get_TR_space(m3, g3noisy, tol=1.e-2).shape == (5, 9))
    for ns in range(2, 6):
        tol = 10. ** -ns
        gnoisy = g3 + tol * noise
        assert(_get_TR_space(m3, gnoisy, tol=10*tol).shape == (5, 9))

    """
    sqrtmmm = np.repeat(np.sqrt(m), 3)
    xxx = np.repeat(geom[:, 0], 3)
    yyy = np.repeat(geom[:, 1], 3)
    zzz = np.repeat(geom[:, 2], 3)

    z = np.zeros_like(m)
    i = np.ones_like(m)
    ux = np.ravel([i, z, z], order='F')
    uy = np.ravel([z, i, z], order='F')
    uz = np.ravel([z, z, i], order='F')

    # form translation and rotation unit vectors
    T1 = sqrtmmm * ux
    T2 = sqrtmmm * uy
    T3 = sqrtmmm * uz
    R4 = sqrtmmm * (yyy * uz - zzz * uy)
    R5 = sqrtmmm * (zzz * ux - xxx * uz)
    R6 = sqrtmmm * (xxx * uy - yyy * ux)

    TRspace = []
    if 'T' in space:
        TRspace.append([T1, T2, T3])
    if 'R' in space:
        TRspace.append([R4, R5, R6])
    if not TRspace:
        # not sure about this, but it runs
        ZZ = np.zeros_like(T1)
        TRspace.append([ZZ])

    TRspace = np.vstack(TRspace)

    def orth(A, tol=tol):
        u, s, vh = np.linalg.svd(A, full_matrices=False)
        if verbose >= 2:
            print(s)
        M, N = A.shape
        eps = np.finfo(float).eps
        if tol is None:
            tol = max(M, N) * np.amax(s) * eps
        num = np.sum(s > tol, dtype=int)
        Q = u[:, :num]
        return Q

    TRindep = orth(TRspace.T)
    TRindep = TRindep.T

    if verbose >= 2:
        print(TRindep.shape, '<--', TRspace.shape)
        print(np.linalg.norm(TRindep, axis=1))
        print('-' * 80)

    return TRindep

#================================================#
#================================================#


##################################################
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
###  REPHASE EIGENVECTORS FROM PSI4 (???)  #######

def _phase_cols_to_max_element(arr, tol=1.e-2, verbose=1):
    """Returns copy of 2D `arr` scaled such that, within cols, max(fabs)
    element is positive. If max(fabs) is pos/neg pair, scales so first
    element (within `tol`) is positive.

    """
    arr2 = np.copy(arr)

    rephasing = []
    for v in range(arr.shape[1]):
        vextreme = 0.0
        iextreme = None

        # find most extreme value
        for varr in arr[:, v]:
            vextreme = max(np.absolute(varr), vextreme)

        # find the first index whose fabs equals that value, w/i tolerance
        for iarr, varr in enumerate(arr[:, v]):
            if (vextreme - np.absolute(varr)) < tol:
                iextreme = iarr
                break

        sign = np.sign(arr[iextreme, v])
        if sign == -1.:
            rephasing.append(str(v))
        arr2[:, v] *= sign

    if rephasing and verbose >= 2:
        print('Negative modes rephased:', ', '.join(rephasing))

    return arr2

#================================================#
#================================================#



if len(sys.argv) != 4:
  sys.exit('usage:\n'
      '1.FD xyz file with optd first, and then each FD (distortion) step;\n'
      '2.string of (period-delimited) atom indices for normal mode calculation (indices starting with 1);\n'
      '3.directory list of directories containing GRAD scratch files (INCLUDING OPTD FIRST) IN SAME ORDER AS input-xyz')


#!# get: 1) indices of non-frozen atoms;
#!#      2) optd_cf in COM coordinates;
#!#      3) moment of inertia eigenvectors;
#!#      4) FD step size
#!#      5) coordinate order of FD displacements relative to cf input order:
FD_atom_inds, optd_cf, X, ss, FD_disp_order = prepare_xyz(sys.argv[1], sys.argv[2])

## convert FD_atom_inds to FD_coord_inds:
FD_coord_inds = []
for ati in FD_atom_inds:
  for coi in range(3):
    FD_coord_inds.append((ati * 3) + coi)


oldgrad = 1
newgrad = 0

## read in gradient vector files:
unsorted_FD_grads = []
## loop over directories containing gradient vector files:
with open(sys.argv[3],'r') as f:
  for d in f.readlines():
    if oldgrad:
      ## read in each gradient vector file, throwing away frozen-atom components:
      unsorted_FD_grads.append(np.genfromtxt(d.strip() + '/gradient.gvec')[FD_coord_inds])
      ##!#
      #unsorted_FD_grads.append(np.genfromtxt(d.strip() + '/gradient.gvec'))
      #print(d)
      #print(np.genfromtxt(d.strip() + '/gradient.gvec')[FD_coord_inds])
      #print('\n')


    elif newgrad:
      tmp_grad = []
      ## open GRAD file in directory, d:
      with open(d.strip() + '/GRAD', 'r') as gf:
        ingrad = gf.readlines()

        found_grad = False
        for li in ingrad:
          if '$gradient' in li.strip():
            found_grad = True
            continue
          if found_grad:
            if '$end' in li.strip():
              break
            tmp_grad += list(map(float, li.strip().split()))

        unsorted_FD_grads.append(np.array(tmp_grad)[FD_coord_inds])
        ##!#
        #unsorted_FD_grads.append(np.array(tmp_grad))
        #print(d)
        #print(np.array(tmp_grad)[FD_coord_inds])
        #print('\n')

#exit()


## convert to np.array for easy indexing/permuting:
unsorted_FD_grads = np.array(unsorted_FD_grads)

#print(unsorted_FD_grads)
#exit()


## get pos and neg gradients in correct order:
FD_grads = {'pos':unsorted_FD_grads[FD_disp_order['pos']],
            'neg':unsorted_FD_grads[FD_disp_order['neg']]}
#  includes +1 to skip optd_cf (in prepare_xyz) ^^^

#print(3 * optd_cf.n_atoms)
#print(len(FD_grads['pos']))
#print(len(FD_grads['pos'][0]))
#print(len(FD_grads['neg']))
#print(len(FD_grads['neg'][0]))
#exit()

## do finite differences:
FD_appx_Hess = np.zeros((3 * optd_cf.n_atoms, 3 * optd_cf.n_atoms))
## 2 * step size IN BOHR:
N = (2. * ss) / bohr2ang
for coi in range(len(FD_appx_Hess)):
  FD_appx_Hess[coi] += (FD_grads['pos'][coi] - FD_grads['neg'][coi]) / N

## make sure Hessian, f_CART, is symmetric:
f_CART = 0.5 * (FD_appx_Hess + FD_appx_Hess.T)
#!f_CART = FD_appx_Hess

#print(f_CART)
#exit()

## transform Cartesian Hessian into mass-weighted coordinates:
f_MWC = np.zeros(f_CART.shape)
for at1 in range(len(optd_cf.masses)):
  m1 = optd_cf.masses[at1]
  for at2 in range(len(optd_cf.masses)):
    m2 = optd_cf.masses[at2]
    for co1 in range(3):
      for co2 in range(3):
        rw = (at1 * 3) + co1
        cl = (at2 * 3) + co2
        f_MWC[rw][cl] = pow(m1*m2, -0.5) * f_CART[rw][cl]

#print(f_MWC)
#exit()

TRspace = _get_TR_space(optd_cf.masses, optd_cf.crds)

#!# also 'borrowed' source code/ideas from psi4 here...
#!# see: psi4/psi4/driver/qcdb/vib.py
P = np.identity(3 * optd_cf.n_atoms)  # psi4
for irt in TRspace:                   # psi4
  P -= np.outer(irt, irt)             # psi4




## project-out translations, rotations from f_MWC:
#!#f_INT = np.dot(P.T, f_MWC).dot(P)  # psi4
f_INT = f_MWC



## solve for eigenvectors, eigenvalues:
f_INT_evals, f_INT_evects = LA.eigh(f_INT)
#f_INT_evals, f_INT_evects = LA.eig(f_INT)

frequency_cm_1 = np.lib.scimath.sqrt(f_INT_evals) * uconv_cm_1  # psi4

#print(f_INT_evals)
#print(f_INT_evects)
#exit()

#f_INT_evects = _phase_cols_to_max_element(f_INT_evects)  # psi4

#print(f_INT_evals)
#print(f_INT_evects)
#print('\n')
#exit()


## apparently this works...
M_mat = np.diag(np.repeat(pow(optd_cf.masses, -0.5), 3))
CART_evects = np.matmul(M_mat, f_INT_evects)

#print(M_mat)
#print('\n')


#print(CART_evects)
#print('\n')
#exit()

CART_evects = CART_evects.T

#print(CART_evects)
#exit()


##!# from psi4:
#reduced_mass_u = np.divide(1.0, np.linalg.norm(CART_evects, axis=0)**2)  # psi4
#
#print(reduced_mass_u)
#exit()
#
#nms = np.sqrt(reduced_mass_u) * CART_evects  # psi4

#print(nms)
#exit()



## reduced mass of each nm is square of normalization factor:
red_mass_nm = []
for nm in CART_evects:
  N = pow(sum([pow(x,2) for x in nm]), -0.5)
  nm *= N
  red_mass_nm.append(pow(N,2))

np.set_printoptions(6, linewidth=120, suppress=True)

#print(optd_cf.masses)

num_nm_skip = 0
#print(f_INT_evals)
for i in range(len(frequency_cm_1) - num_nm_skip):
  print(frequency_cm_1[i + num_nm_skip], CART_evects[i + num_nm_skip])
  #print(frequency_cm_1[i + num_nm_skip], nms[i + num_nm_skip])



