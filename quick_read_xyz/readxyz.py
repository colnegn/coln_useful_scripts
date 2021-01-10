#!/usr/bin/env python

import sys
import numpy as np


## taken from: https://stackoverflow.com/questions/14645789/numpy-reading-file-with-filtering-lines-on-the-fly
def filter_lines(f, natoms):
  ## loop over lines in file, f:
  for lii, li in enumerate(f):
    ## == False for atom number & comment lines:
    if lii % (natoms + 2) > 1:
      ## generator:
      yield li


def readxyz(inxyz, natoms):
  with open(inxyz) as f:
    ## read xyz coordinates into np.array:
    crds = np.genfromtxt(filter_lines(f, natoms),
                         dtype='f',
                         usecols=(1, 2, 3))


  nframes = int(len(crds) / natoms)
  crds = crds.reshape((nframes, natoms, 3))

  return crds




if __name__ == '__main__':

  if len(sys.argv) != 3:
    sys.exit('usage: 1.xyz file; 2.number of atoms')

  inxyz = sys.argv[1]
  natoms = int(sys.argv[2])

  crds = readxyz(inxyz, natoms)

  print(crds)


