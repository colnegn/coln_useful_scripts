#!/usr/bin/env python

import cmath
import re
import sys
import math
import numpy as np
from numpy import linalg as LA
from copy import deepcopy as copy

import gcoord_tools as ctools

##################################################
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
###  xyz class  ##################################

class xyz(object):


##################################################
###  initialize  #################################

  def __init__(self, in_xyz_list=None):
    '''
    in_xyz_list is a list (rows) of lists (columns) from
    an xyz file, for a *single* configuration.
    '''

    ## empty 'xyz' object (need for fragment_xyz):
    if in_xyz_list is None:
      return

    ## initialize attributes from in_xyz_list:
    else:
      ## number of atoms:
      self.n_atoms = int(in_xyz_list[0][0])
      ## comment line:
      self.comment = ' '.join(in_xyz_list[1])
      ## atom symbols:
      self.symbols = np.array([c[0] for c in in_xyz_list[2:]])
      ## coordinates:
      self.crds = np.array([list(map(float, c[1:])) for c in in_xyz_list[2:]])

      ## dictionary of lists of atom-indices for each element:
      self.at_inds = self.get_at_inds()


      ## get masses from symbols:
      try:
        self.masses = self.get_masses()
      except:
        self.masses = None


      ## dictionary for energies:
      self.energy = None

    return

#================================================#
#================================================#


##################################################
###  get masses  #################################

  def get_masses(self, symbol_convert_dict=None):

    ################################
    ## convert symbols for xyz file:
    if symbol_convert_dict is None:
      symbol_convert_dict = {sy: sy for sy in self.symbols}
    ################################


    out_masses = []

    try:
      import periodictable as pt
      for sy in self.symbols:
        #!# get mass of the most abundant isotope IN DALTONS (1.6605390666E-27 kg):
        #out_masses.append(getattr(pt, symbol_convert_dict[sy].capitalize()).mass)
        out_masses.append(sorted([[iso.abundance, iso.mass] for iso in getattr(pt,symbol_convert_dict[sy].capitalize())], key=lambda x: x[0], reverse = True)[0][1])


    except:
      import datmass
      mass_dict = datmass.get_masses()

      out_masses = []
      for sy in self.symbols:
        out_masses.append(mass_dict[symbol_convert_dict[sy]]['mass'])

    return np.array(out_masses)

#================================================#
#================================================#


##################################################
###  fragment xyz  ###############################

  def fragment_xyz(self, permd_frag_inds, new_comment=None, append_comment=True):
    '''
    returns a new 'xyz' object with fragment configuration

    =========
    arguments

    + permd_frag_inds :  <list>  indices (permuted) of fragment atoms from self
    + new_comment     :  <str>   to for comment line of xyz-file
    + append_comment  :  <bool>  append new_comment to self.comment in xyz-file
    '''

    ##################
    ## handle comment:
    if new_comment is None:
      frag_comment = copy(self.comment)
    else:
      if append_comment:
        frag_comment = copy(self.comment) + ' | ' + str(new_comment)
      else:
        frag_comment = str(new_comment)
    ##################

    #################################
    ## get fragment 'xyz' attributes:
    frag_symbols = copy(self.symbols[permd_frag_inds])
    frag_crds    = copy(self.crds[permd_frag_inds])
    #################################


    ## now initialize (empty) 'xyz' fragment object:
    frag_xyz = xyz()
    ## set attributes:
    frag_xyz.set_members(frag_comment,frag_symbols,frag_crds)


    ## return 'xyz' object for fragment xyz:
    return frag_xyz

#================================================#
#================================================#


##################################################
###  get atom index dictionary  ##################

  def get_at_inds(self):

    ## initialize dictionary of empty lists for each element in self:
    out_dict = {sy: [] for sy in set(self.symbols)}

    ## loop over atoms (symbols):
    for i, sy in enumerate(self.symbols):
      out_dict[sy].append(i)

    return out_dict

#================================================#
#================================================#


##################################################
###  set members separately  #####################

  def set_members(self, in_comment, in_symbols, in_crds):
    self.n_atoms = len(in_symbols)
    self.comment = str(in_comment)
    self.symbols = np.array(in_symbols)
    self.crds    = np.array(in_crds)

    self.at_inds = self.get_at_inds()

    try:
      self.masses = self.get_masses()
    except:
      self.masses = None

    return

#================================================#
#================================================#


##################################################
###  write xyz file  #############################

  def write_xyz_file(self, return_flag=False, plist_in=None, out_file=None,
                     new_comment=None, append_comment=True,
                     symbol_convert_dict=None, new_permdsymbols=None):
    '''
    writes an xyz-file from 'xyz' object (self).

    =========
    arguments

    + return_flag         : <bool> returns 'str' instead of writing to file
    + plist               : <list> permuted indices to permute atoms in xyz-file
    + out_file            : <file> to write xyz-file to (must be open)
    + new_comment         : <str>  to for comment line of xyz-file
    + append_comment      : <bool> append new_comment to self.comment in xyz-file
    + symbol_convert_dict : <dict> keys are symbol strings to convert (to elements)
    + new_permdsymbols    : <list> symbols in permuted order to replace self.symbols (not compatible with convert_dict)
    '''

    ####################
    ## atom-permutation:
    if plist_in is None:
      ## 'identity' permutation list:
      plist = list(range(self.n_atoms))
    else:
      ## in case we get a string of indices:
      plist = [int(i) for i in plist_in]
    ####################

    ################################
    ## convert symbols for xyz file:
    if symbol_convert_dict is None:
      symbol_convert_dict = {sy: sy for sy in self.symbols}
    ################################

    ########################
    ## new permuted symbols:
    if new_permdsymbols is None:
      out_symbols = [symbol_convert_dict[sy] for sy in self.symbols[plist]]
    else:
      out_symbols = np.array(new_permdsymbols)
    ########################

    #######################
    ## handle comment line:
    if new_comment is None:
      out_comment = self.comment
    else:
      if append_comment:
        out_comment = self.comment + ' || ' + str(new_comment)
      else:
        out_comment = str(new_comment)
    #######################


    ## prepare xyz file out-string:
    out_str = str(len(plist)) + '\n' + out_comment

    ####################################
    ## loop over (permuted) coordinates:
    for i, at in enumerate(self.crds[plist]):
    #for i, at in enumerate(self.crds): #!#
      out_str += '\n' + out_symbols[i] + ' '*4 + \
       '   '.join(list(map(lambda x: '{:>10s}'.format(str('{:0.6f}'.format(x))), at)))
    ####################################

    #################################################
    ## only writes or prints if return_flag == False:
    if return_flag:
      return out_str
    #################################################

    ##################
    ## write xyz file:
    if out_file is None:
      print(out_str)
    else:
      out_file.write(out_str + '\n')
    ##################

    return

#================================================#
#================================================#


#================================================#
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#================================================#




##################################################
###  read xyz-file  ##############################
def read_xyz(xyz_file_path):

  # test:
  #debug_flag = True
  debug_flag = False


  ## convert xyz file into list (rows) of lists (columns):
  with open(xyz_file_path) as f:
    in_xyz = list(map(lambda x: x.strip().split(), f.readlines()))


  ### gather all configurations from xyz file:

  ## list of configurations:
  in_configs = []
  ## line (row) counter for while-loop:
  li_count = 0
  while li_count < len(in_xyz):

    ## get number of atoms for current config:
    n_atoms = int(in_xyz[li_count][0])
    ## append sub-list of config to in_configs:
    in_configs.append(xyz(in_xyz[li_count:(li_count + n_atoms + 2)]))

        #configuration(xyz(in_xyz[li_count:(li_count + n_atoms + 2)]), 'xyz object')
        #             )

    ## increment li_count:
    li_count += (n_atoms + 2)


  #######
  # test:
  if debug_flag:
    for cf in in_configs:
      print (cf.symbols)
      print (cf.crds)
      print ('\n====\n')
  #######


  return in_configs
#================================================#
#================================================#


##################################################
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
###  MAIN  #######################################

if __name__ == '__main__':

  sys.exit('please call me from some separate \'main\' script')

#================================================#
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#================================================#





