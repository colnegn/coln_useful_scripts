#!/usr/bin/env python

import sys
import numpy as np


##################################################
###  readxyz  ####################################

## taken from: https://stackoverflow.com/questions/14645789/numpy-reading-file-with-filtering-lines-on-the-fly
def filter_lines(f, natoms):
  ## loop over lines in file, f:
  for lii, li in enumerate(f):
    ## == False for atom number & comment lines:
    if lii % (natoms + 2) > 1:
      ## generator:
      yield li


def readxyz(inxyz):
  ## first get number of atoms from first line:
  with open(inxyz) as f:
    natoms = int(f.readline())
    ## throw away comment line:
    f.readline()
    ## next get atom labels:
    atom_labels = []
    for ati in range(natoms):
      atom_labels.append(f.readline().strip().split()[0])

  ## close previous file instance, and re-open:
  with open(inxyz) as f:
    ## read xyz coordinates into np.array:
    crds = np.genfromtxt(filter_lines(f, natoms),
                         dtype='f',
                         usecols=(1, 2, 3))


  ## reshape so indices are: [frame][atom][coordinate]
  nframes = int(len(crds) / natoms)
  crds = crds.reshape((nframes, natoms, 3))


  return natoms, atom_labels, crds

#================================================#
#================================================#






if __name__ == '__main__':

  if len(sys.argv) != 2:
    sys.exit('usage: 1.xyz file')

  inxyz = sys.argv[1]

  natoms, atom_labels, crds = readxyz(inxyz)

  print(crds)


