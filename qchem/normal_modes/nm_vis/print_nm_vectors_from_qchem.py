#!/usr/bin/env python

import sys
import numpy as np


if len(sys.argv) != 4:
  sys.exit('usage: 1.qchem output file for freq calculation; 2.number of atoms; 3.prefix for nm vector files')

natoms = int(sys.argv[2])
pref   = sys.argv[3]

################################
### get qchem output file:

with open(sys.argv[1],'r') as f:
  ## open file as list of strings:
  out_file = f.readlines()

  find_str = 'VIBRATIONAL ANALYSIS'
  ## find line number containing "VIBRATIONAL ANALYSIS" in out_file:
  top_li = [i for i, li in enumerate(out_file) if find_str in li][0]

  ## only need lines below top_li:
  out_file = out_file[top_li:]
################################



find_str = '               X      Y      Z'
## find all blocks of normal mode vectors:
top_blocks = [i for i, li in enumerate(out_file) if find_str in li]


## list of normal modes:
nms = []
## list of vibrational frequencies:
freqs = []

## loop over blocks:
for bl in top_blocks:
  ## split up strings (lines) containing normal mode vectors:
  nm_bl = np.asarray([[float(lii) for lii in li.strip().split()[1:]]
                      for li in out_file[(bl+1):(bl+natoms+1)]]).T
#                transpose so nm vector components are together ^^

  ## also try to get frequencies for block:
  freqs += [float(frq) for frq in out_file[bl-6].strip().split()[1:]]


  #!# make sure nm block has correct number of columns:
  if len(nm_bl) % 3 != 0:
    print('problem?')
    print(nm_bl, '\n\n')
    raise Exception('normal mode block wrong size')


  ## separate the normal modes in nm_bl:
  for nm in range(int(float(len(nm_bl)) / 3.0)):
#        3-dimensional space ^^^

    ## 3*natoms-dimensional vector for nm:
    nms.append(nm_bl[(nm*3):((nm+1)*3)].T)
    #nms.append(nm_bl[(nm*3):((nm+1)*3)].flatten())

nms = np.array(nms)


## format for frequencies:
frq_fmt = '{:0.2f}'

for i in range(len(nms)):
  with open(pref + '_' + str(i) + '_' + frq_fmt.format(freqs[i]) + 'cm-1.nmvec', 'w') as f:
    for at in nms[i]:
      #f.write('  '.join(map(str, at)) + '\n')
      f.write(' '.join(['{:>12.6f}'.format(xy) for xy in at]) + '\n')


