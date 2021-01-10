#!/usr/bin/env python



class main_functions(object):

##################################################
###  initialize  #################################

  def __init__(self):

    ## list of functions:
    self.function_list = {f: getattr(self,f)
                           for f in dir(self)
                             if callable(getattr(self,f))
                             and f[:2] != '__'}

    self.loaded_modules_bool = False

    return

#================================================#
#================================================#

  def load_modules(self):

    if not self.loaded_modules_bool:

      global          np
      import numpy as np
      global scipy
      import scipy
      global              itts
      import itertools as itts

      global                       copy
      from copy import deepcopy as copy
      global                      LA
      from numpy import linalg as LA
      global            signal
      from scipy import signal

      global            xyz_io
      import gxyz_io as xyz_io
      global                 ctools
      import gcoord_tools as ctools
      global COMPAX
      import COMPAX


      self.loaded_modules_bool = True

    return

##################################################
###  help  #######################################

  def Help(self, method_name=''):

    help_dict = {}

#    help_dict['Initialize_to_COMPAX'] = ' + Initialize_to_COMPAX:\n' +\
#        '    1.in_xyz_path;

    help_str = '\n\n"Help" followed by method name give usage\n\n========\nMETHODS:\n\n'
    help_str += ' + Initialize_to_COMPAX\n'
    help_str += ' + Rotate_molecule\n'
    help_str += ' + Rotate_scan\n'
    help_str += ' + Translate_molecule\n'
    help_str += ' + Extract_fragments\n'
    help_str += ' + Calculate_RDF\n'
    help_str += ' + Calculate_RADF\n'
    help_str += ' + Calculate_Z_C2_RADF\n'
    help_str += ' + Calculate_Z_ACN_RRDF\n'
    help_str += ' + Calculate_kinetic_energy\n'
    help_str += ' + Bin_zundel_dihedral\n'
    help_str += ' + ACN_OCF\n'
    help_str += ' + ACN_ACF\n'
    #!#
    help_dict[''] = help_str
    help_str = ''


    help_str += ' + Initialize_to_COMPAX:\n'
    help_str += ' 1.in_xyz_path;\n 2.out_COMPAX_xyz_path;\n 3.bool to print out PAX-rotation (1~True, [def:]0~False)\n\n'
    #!#
    help_dict['Initialize_to_COMPAX'] = help_str
    help_str = ''


    help_str += ' + Rotate_molecule:\n'
    help_str += ' 1.in_xyz_path;\n 2.axis_path;\n 3.angle (units of radians/2pi);\n 4.path to list of integers(starting with 0) of atoms to rotate (other atoms remain fixed)\n\n'
    #!#
    help_dict['Rotate_molecule'] = help_str
    help_str = ''


    help_str += ' + Rotate_scan:\n'
    help_str += ' 1.in_xyz_path;\n 2.axis_path;\n 3.incremental angle (units of radians/2pi);\n 4.number of incremental rotations;\n 5.path to list of integers(starting with 0) of atoms to rotate (other atoms remain fixed)\n 6.inclusive_bool, 0 or 1 (defaults to 1)\n\n'
    #!#
    help_dict['Rotate_scan'] = help_str
    help_str = ''


    help_str += ' + Translate_molecule:\n'
    help_str += ' 1.in_xyz_path;\n 2.increment_disp_3N_vector_path;\n 3.nsteps;\n 4.inclusive_bool, 0 or 1 (defaults to 1)\n\n'
    #!#
    help_dict['Translate_molecule'] = help_str
    help_str = ''


    help_str += ' + Extract_fragments:\n'
    help_str += ' 1.in_xyz_path;\n 2.natoms_per_mol_list_path;\n 3.nmols_per_frag (int);\n'
    help_str +=    ' 4.inner_cufoff_distance (float);\n 5.outer_cutoff_distance (float);\n'
    help_str +=    ' 6.lbox (float) assumes cubic box (set to 0.0 if constant pressure);\n 7.ref_mol_index_path (0 index) (defaults to all);\n 8.const_pressure_bool (1~True, [def:]0~False) ASSUMING I-PI FORMAT (defaults to False)\n'
    help_str +=    ' 9.bool for all pairs (1~include molecules before and after reference molecule, [def:]0~only get molecules after reference) careful with double counting...;\n 10.bool for only get closest fragment (1~True, [def:]0~False);\n'
    help_str +=    ' 11.path to list of: comma-delimited lists of one or more atom indices for each molecule to calculate average distances (0 is the first atom of any molecule), defaults to first atom of each molecule (default is same as list with: 0 0 0 0...)\n\n'
    #!#
    help_dict['Extract_fragments'] = help_str
    help_str = ''


    help_str += ' + Calculate_RDF:\n'
    help_str += ' 1.xyz file (i-pi format comment, CUBIC BOX ONLY);\n 2.at1_list path (indices start with 1);\n 3.at2_list path;\n 4.minimum distance;\n 5.maximum distance;\n 6.number of bins;\n 7.output path (to overwrite);\n 8.bool for same atom type (1~true, [def:]0~false)\n\n'
    #!#
    help_dict['Calculate_RDF'] = help_str
    help_str = ''


    help_str += ' + Calculate_RADF:\n'
    help_str += ' 1.xyz file (i-pi format comment, CUBIC BOX ONLY);\n 2.at1_list path (indices start with 1);\n 3.at23_list_path with two columns (atom 2 and atom 3, paired);\n 4.number of angular bins;\n 5.minimum distance;\n 6.maximum distance;\n 7.number of radial bins;\n 8.output path (to overwrite);\n 9.bool for same 1-2 atom type (1~true, [def:]0~false);\n 10.bool for same 1-3 atom type;\n 11.bool for same 2-3 atom type\n\n'
    #!#
    help_dict['Calculate_RADF'] = help_str
    help_str = ''


    help_str += ' + Calculate_Z_C2_RADF:\n'
    help_str += ' 1.xyz file (i-pi format comment, CUBIC BOX ONLY);\n 2.at1_list path (indices start with 1);\n 3.zundel inds list path (assuming order: O1,O2,H1f,H1f,Hs,H2f,H2f);\n 4.number of angular bins;\n 5.minimum distance;\n 6.maximum distance;\n 7.number of radial bins;\n 8.output path (to overwrite);\n9.prefix for residence time outputs;\n10.initial atom1 molecule index for this split (inds start with 1)\n\n'
    #!#
    help_dict['Calculate_Z_C2_RADF'] = help_str
    help_str = ''


    help_str += ' + Calculate_kinetic_energy:\n'
    help_str += ' 1.input velocity xyz-format trajectory file;\n 2.mass conversion factor (to multiply by AMU);\n 3.velocity conversion factor (to multiply by not-squared velocity vector from trajectory;\n 4.bool to print average ([def:]1~print avg, 0~print list)\n\n'
    #!#
    help_dict['Calculate_kinetic_energy'] = help_str
    help_str = ''


    help_str += ' + Bin_zundel_dihedral:\n'
    help_str += ' 1.input xyz path (i-pi formatted comment lines);\n 2.path to dihedral atom index path INDICES START WITH 0 (2 lines: H11 H12 O1, H21 H22 O2);\n 3.number of bins;\n 4.prefix for output xyzs\n\n'
    #!#
    help_dict['Bin_zundel_dihedral'] = help_str
    help_str = ''


    help_str += ' + ACN_OCF:\n'
    help_str += ' 1.input xyz path (i-pi formatted comment lines);\n 2.path to atom1 (ref) index list (INDICES START WITH 1);\n 3.path to atom2 (distance) index list;\n 4.path to atom3 (angle) index list;\n 5.minimum distance between atom1 and atom2;\n 6.maximum distance;\n 7.box length (only constant V at the moment);\n 8.minimum number of steps for OCF segments to average;\n 9.timestep between configurations (in ps);\n 10.molecule index for first atoms 2, 3 in this split (indices start with 1; use 1 if xyz not split)\n\n'
    #!#
    help_dict['ACN_OCF'] = help_str
    help_str = ''


    help_str += ' + ACN_ACF:\n'
    help_str += ' 1.input xyz path (i-pi formatted comment lines);\n 2.path to atom1 (ref) index list (INDICES START WITH 1);\n 3.path to atom2 (distance) index list;\n 4.path to atom3 (angle) index list;\n 5.minimum distance between atom1 and atom2;\n 6.maximum distance;\n 7.box length (only constant V at the moment);\n 8.minimum number of steps for ACF segments to average;\n 9.timestep between configurations\n\n'
    #!#
    help_dict['ACN_ACF'] = help_str
    help_str = ''


    ## some possible error? (e.g. wrong method name)
    try:
      print(help_dict[method_name])
    except:
      print(help_dict[''])

#    help_str += '\n'
#    print(help_str)


    return

#================================================#
#================================================#


##################################################
###  calculate kinetic energy from velocities  ###

  def Calculate_kinetic_energy(self,
                               in_velxyz_path,
                               mass_conversion_factor,
                               velocity_conversion_factor,
                               average_bool=True):


    self.load_modules()

    mass_conversion_factor = float(mass_conversion_factor)
    velocity_conversion_factor = float(velocity_conversion_factor)


    vel_traj = xyz_io.read_xyz(in_velxyz_path)


    ## list to store total kinetic energies of each trajectory frame:
    tot_KEs = []

    ## loop over frames of trajectory:
    for fr in vel_traj:
      ## want total kinetic energy:
      tot_KEs.append(0.0)
      ## loop over atom indices and velocities:
      for ati, atvel in enumerate(fr.crds):
        ## converted mass of atom:
        convd_mass_ati = fr.masses[ati] * mass_conversion_factor
        ## convert velocity vector for atom, ati:
        convd_vel_ati = atvel * velocity_conversion_factor
        ## calculate CONVERTED v^2:
        vv = np.dot(convd_vel_ati, convd_vel_ati)

        ## add kinetic energy of atom, ati, to total KE:
        tot_KEs[-1] += 0.5 * convd_mass_ati * vv


    if bool(int(average_bool)):
      print(np.mean(tot_KEs))
    else:
      ## print output to stdout:
      print('\n'.join([str(fri) + ' ' + str(ke) for fri, ke in enumerate(tot_KEs)]))

    return

#================================================#
#================================================#


##################################################
###  initialize coordinates to COMPAX  ###########

  def Initialize_to_COMPAX(self,
                           in_xyz_path,
                           out_COMPAX_xyz_path,
                           print_rotmats_bool=False):

    self.load_modules()
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    print_rotmats_bool = bool(int(print_rotmats_bool))

    ## list to keep rotation matrices:
    #PAX_rotmats = []

    ## open output file:
    with open(out_COMPAX_xyz_path, 'w') as outf:
      ## orient each configuration into COMPAX frame:
      for xy in xyz_obj_list:
        xy.crds, rotmat = COMPAX.init_coord_frame(xy.crds, xy.masses, return_rotmat_bool=True)

        #PAX_rotmats.append(rotmat)

        ## write xyz file with configurations in COMPAX:
        xy.write_xyz_file(out_file=outf)

        ## do we want to print rotation matrix to stdout?
        if print_rotmats_bool:
          print(rotmat)

    return

#================================================#
#================================================#


#################################################
###  PBC reimage coordinates  ####################

  def PBC_reimage_coordinates(self,
                              reference_coords,
                              coords_to_reimage,
                              lbox):

    self.load_modules()
    #!# ASSUMES CUBIC BOX !!!!!!!!!

    #debug_bool = True
    debug_bool = False

    reference_coords  = np.asarray(reference_coords)

    ## in case 1-index array:
    coords_to_reimage = np.asarray(coords_to_reimage)
    if coords_to_reimage.ndim == 1:
      ## make it a 2-index array:
      coords_to_reimage = np.array([coords_to_reimage])
    ## ???
    if coords_to_reimage.ndim > 2:
      sys.exit('ERROR: too many indices in input to PBC_reimage_coordinates')

    ## array for reimaged coordinates:
    out_crds = np.zeros(coords_to_reimage.shape)

    ## loop over atoms to be reimaged:
    for ati, at in enumerate(coords_to_reimage):
      diffvec = at - reference_coords
      ## reimaging vector from diffvec:
      reimvec = np.round(diffvec / lbox) * lbox

      #######
      # test:
      if debug_bool:
        print('ref:', reference_coords)
        print('at: ', at)
        print('non-imaged dist:', LA.norm(at - reference_coords))
        print('lbox / 2:       ', lbox / 2.0)
        print('diffvec:  ', diffvec)
        print('reimcells:', np.round(diffvec / lbox))
        print('reimvec:  ', reimvec)
        print('imaged-coords:', at - reimvec)
        print('imaged dist:', LA.norm(at - reimvec - reference_coords))
        print('\n')
      #######

      out_crds[ati] += at - reimvec

    return out_crds

#================================================#
#================================================#


#################################################
###  Get PBC reimage vector  ####################

  def Get_PBC_reimage_vector(self,
                             reference_coords,
                             single_at_coords_to_reimage,
                             lbox):

    self.load_modules()
    #!# ASSUMES CUBIC BOX !!!!!!!!!

    #debug_bool = True
    debug_bool = False

    reference_coords = np.asarray(reference_coords)
    single_at_coords_to_reimage = np.asarray(single_at_coords_to_reimage)


    ## array for reimaged coordinates:
    #out_crds = np.zeros(coords_to_reimage.shape)

    diffvec = single_at_coords_to_reimage - reference_coords
    ## reimaging vector from diffvec:
    reimvec = np.round(diffvec / lbox) * lbox

    #######
    # test:
    if debug_bool:
      print('ref:', reference_coords)
      print('at: ', single_at_coords_to_reimage)
      print('non-imaged dist:', LA.norm(single_at_coords_to_reimage - reference_coords))
      print('lbox / 2:       ', lbox / 2.0)
      print('diffvec:  ', diffvec)
      print('reimcells:', np.round(diffvec / lbox))
      print('reimvec:  ', reimvec)
      print('imaged-coords:', single_at_coords_to_reimage - reimvec)
      print('imaged dist:', LA.norm(single_at_coords_to_reimage - reimvec - reference_coords))
      print('\n')
    #######

    #out_crds[ati] = single_at_coords_to_reimage - reimvec

    return reimvec

#================================================#
#================================================#



##################################################
###  reimage molecule, unbreaking bonds  #########

  def Reimage_molecule_PBC(self,
                           reference_coords,
                           mol_reference_coords,
                           mol_other_coords,
                           lbox):

    self.load_modules()
    ## in case 1-index array:
    mol_other_coords = np.asarray(mol_other_coords)
    if mol_other_coords.ndim == 1:
      ## make it a 2-index array:
      mol_other_coords = np.array([mol_other_coords])


    ## total number of atoms in molecule:
    num_ats_mol = len(mol_other_coords) + 1
    ## dimension of space:
    ndim = len(reference_coords)

    ## first we reimage rest of molecule relative to mol_ref:
    reim_mol_crds = np.zeros((num_ats_mol, ndim))
    reim_mol_crds[0] += np.asarray(mol_reference_coords)
    reim_mol_crds[1:] += self.PBC_reimage_coordinates(mol_reference_coords, mol_other_coords, lbox)

    ## now get the reimaging vector for mol_ref relative to ref:
    reimvec = self.Get_PBC_reimage_vector(reference_coords, mol_reference_coords, lbox)

    ## return reimaged molecule coordinates with mol_ref as first atom:
    return reim_mol_crds - reimvec

#================================================#
#================================================#



##################################################
###  calculate RDF  ##############################

  def Calculate_RDF(self,
                    in_xyz_path,
                    at1_list_path,
                    at2_list_path,
                    min_dist,
                    max_dist,
                    num_bins,
                    output_path,
                    same_at_type_bool=False):


    self.load_modules()
    ## convert trajectory into xyz objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    ## get box length for each configuration from comment line:
    lbox_list = [float(cf.comment.split(' ')[2]) for cf in xyz_obj_list]

    at1_list = np.genfromtxt(at1_list_path, dtype=int) - 1
    at2_list = np.genfromtxt(at2_list_path, dtype=int) - 1

    ## list for number of atoms for each (for normalization):
    nat_eachtype_list = []
    ## check if single element lists:
    try:
      nat_eachtype_list.append(len(at1_list))
    except TypeError:
      at1_list = np.array([at1_list])
      nat_eachtype_list.append(len(at1_list))

    try:
      nat_eachtype_list.append(len(at2_list))
    except TypeError:
      at2_list = np.array([at2_list])
      nat_eachtype_list.append(len(at2_list))


    ## normalization:
    if bool(int(same_at_type_bool)):
      norm = 1.0 / float(nat_eachtype_list[0] * (nat_eachtype_list[1] - 1))
    else:
      norm = 1.0 / float(nat_eachtype_list[0] * nat_eachtype_list[1])

    norm /= float(len(xyz_obj_list))


    min_dist = float(min_dist)
    max_dist = float(max_dist)
    num_bins = int(num_bins)
    dR       = (max_dist - min_dist) / float(num_bins)

    RDF = np.array([[min_dist + ((i + 0.5) * dR), 0.0] for i in range(num_bins)])


    for cfi, cf in enumerate(xyz_obj_list):
      lbox = lbox_list[cfi]

      for at1 in at1_list:
        for at2 in at2_list:
          if at1 == at2:
            continue

          dvec = cf.crds[at1] - cf.crds[at2]
          ## PBC:
          dist12 = LA.norm(dvec - (np.array([round(xy / lbox) for xy in dvec]) * lbox))

          ## skip if out of range:
          if dist12 > max_dist or dist12 < min_dist:
            continue

          ## find bin:
          binindex = int(np.floor((dist12 - min_dist) / dR))

          RDF[binindex][1] += 1.0



    dR2 = dR / 2.0
    for bn in RDF:
      #!# ASSUMES CUBIC BOX !!!!!!
      bn[1] *= (norm * pow(lbox_list[0], 3) * 3.0) / (4.0 * np.pi * (pow(bn[0] + dR2, 3) - pow(bn[0] - dR2, 3)))


    ## write RDF to file:
    #np.savetxt(output_path, RDF, fmt='%.10f')
    np.savetxt(output_path, RDF)

    return

#================================================#
#================================================#


##################################################
###  radial/angular distribution function  #######

  def Calculate_RADF(self,
                     in_xyz_path,
                     at1_list_path,
                     at23_list_path,
                     num_ang_bins,
                     min_dist,
                     max_dist,
                     num_rad_bins,
                     output_path,
                     same_12_at_type_bool=False,
                     same_13_at_type_bool=False,
                     same_23_at_type_bool=False):



    self.load_modules()

    ## convert trajectory into xyz objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    ## get box length for each configuration from comment line:
    lbox_list = [float(cf.comment.split(' ')[2]) for cf in xyz_obj_list]

    at1_list = np.genfromtxt(at1_list_path, dtype=int) - 1
    at23_list = np.genfromtxt(at23_list_path, dtype=int) - 1

    ## list for number of atoms for each (for normalization):
    nat_eachtype_list = []
    ## check if single element lists:
    try:
      nat_eachtype_list.append(len(at1_list))
    except TypeError:
      at1_list = np.array([at1_list])
      nat_eachtype_list.append(len(at1_list))

    #!# expecting at least two values, so make sure to convert to 2-index array:
    if len(at23_list.shape) == 1:
      at23_list = np.asarray([at23_list])

    nat_eachtype_list.append(len(at23_list))



    ## normalization (remember: at2, at3 are paired, so same number of atoms):
    if   bool(int(same_12_at_type_bool)) and bool(int(same_13_at_type_bool)):
      norm = 1.0 / float(nat_eachtype_list[0] * (nat_eachtype_list[1] - 1))
    elif bool(int(same_12_at_type_bool)) or bool(int(same_13_at_type_bool)):
      norm = 1.0 / float(nat_eachtype_list[0] * (nat_eachtype_list[1] - 1))
    else:
      norm = 1.0 / float(nat_eachtype_list[0] * nat_eachtype_list[1])

    norm /= float(len(xyz_obj_list))


    num_ang_bins = int(num_ang_bins)
    dA           = 180.0 / float(num_ang_bins)

    min_dist     = float(min_dist)
    max_dist     = float(max_dist)
    num_rad_bins = int(num_rad_bins)
    dR           = (max_dist - min_dist) / float(num_rad_bins)


    RADF = np.array([[min_dist + (float(ri) + 0.5) * dR, (float(ai) + 0.5) * dA, 0.0] for ri in range(num_rad_bins) for ai in range(num_ang_bins)])


    for cfi, cf in enumerate(xyz_obj_list):
      lbox = lbox_list[cfi]

      ## start with at2 since we reimage 1 and 3 relative to 2:
      for ats23 in at23_list:
        at2 = ats23[0]
        at3 = ats23[1]
        crds2 = cf.crds[at2]


        for at1 in at1_list:
          if at1 == at2:
            continue

          ## reimage at1:
          reimvec = self.Get_PBC_reimage_vector(crds2, cf.crds[at1], lbox)
          reim_crds1 = cf.crds[at1] - reimvec

          ## 12 distance:
          dist12 = ctools.distance(crds2, reim_crds1)

          ## skip if out of range:
          if dist12 > max_dist or dist12 < min_dist:
            continue


          ## reimage at3:
          reimvec = self.Get_PBC_reimage_vector(crds2, cf.crds[at3], lbox)
          reim_crds3 = cf.crds[at3] - reimvec

          ## calculate angle:
          ang123 = (180.0 / np.pi) * ctools.angle(reim_crds1, crds2, reim_crds3)

          ## find bin:
          radbinindex = int(np.floor((dist12 - min_dist) / dR))
          angbinindex = int(np.floor(ang123 / dA))

          binindex = (radbinindex * num_ang_bins) + angbinindex

          RADF[binindex][2] += 1.0


    ## divide by volume and normalize:
    for ri in range(num_rad_bins):
      r_lw = min_dist + (ri * dR)
      r_up = r_lw + dR

      #!# ASSUMES CUBIC BOX !!!!!!

      ## volume of entire radial shell:
      V_r = (4.0 * np.pi / 3.0) * (pow(r_up, 3) - pow(r_lw, 3))

      for ai in range(num_ang_bins):

        ## fraction of radial shell for angular bin, ai (from area of sherical cap):
        frac_V = 0.5 * (np.cos(ai * dA * np.pi / 180.) - np.cos((ai+1) * dA * np.pi / 180.))

        Vol_bin = V_r * frac_V

        ## to first order (from Jacobian determinant, after integrating out az-angle):
        #Vol_bin = 360.0 * pow(r_lw + (dR / 2.0), 2) * dR * np.sin((ai + 0.5) * (dA * np.pi / 180.0)) * dA
#                 ^^^^^ 360 from integrating out azimuthal angle (in degrees)

        #!# ASSUMES FIXED VOLUME:
        coef_ri_ai = norm * pow(lbox_list[0], 3) / Vol_bin


        binindex = (ri * num_ang_bins) + ai
        RADF[binindex][2] *= coef_ri_ai


    ## write RADF to file:
    #np.savetxt(output_path, RADF, fmt='%.10f')
    np.savetxt(output_path, RADF)

    return

#================================================#
#================================================#



##################################################
###  Z C2 radial/angular distribution function  ##

  def Calculate_Z_C2_RADF(self,
                          in_xyz_path,
                          at1_list_path,
                          z_inds_path,
                          num_ang_bins,
                          min_dist,
                          max_dist,
                          num_rad_bins,
                          output_path,
                          respref,
                          initial_at1mi):

    self.load_modules()

    ## to give correct molecule index in residence times outputs:
    initial_at1mi = int(initial_at1mi)

    ## convert trajectory into xyz objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    ## get box length for each configuration from comment line:
    lbox_list = [float(cf.comment.split(' ')[2]) for cf in xyz_obj_list]

    at1_list = np.genfromtxt(at1_list_path, dtype=int) - 1
    z_inds_list = np.genfromtxt(z_inds_path, dtype=int) - 1

    ## list for number of atoms for each (for normalization):
    nat_eachtype_list = []
    ## check if single element lists:
    try:
      nat_eachtype_list.append(len(at1_list))
    except TypeError:
      at1_list = np.array([at1_list])
      nat_eachtype_list.append(len(at1_list))



    norm = 1.0 / float(nat_eachtype_list[0])
    norm /= float(len(xyz_obj_list))


    num_ang_bins = int(num_ang_bins)
    dA           = 180.0 / float(num_ang_bins)

    min_dist     = float(min_dist)
    max_dist     = float(max_dist)
    num_rad_bins = int(num_rad_bins)
    dR           = (max_dist - min_dist) / float(num_rad_bins)


    RADF = np.array([[min_dist + (float(ri) + 0.5) * dR, (float(ai) + 0.5) * dA, 0.0] for ri in range(num_rad_bins) for ai in range(num_ang_bins)])



    ## get zundel indices assuming order: O1,O2,H1f,H1f,Hs,H2f,H2f
    oinds = copy(z_inds_list[:2])
    hfinds1 = copy(z_inds_list[2:4])
    hfinds2 = copy(z_inds_list[5:])
    hfinds = np.asarray([hfinds1, hfinds2])
    hsind = z_inds_list[4]

    zinds = np.concatenate((oinds, hfinds1, np.asarray([hsind]), hfinds2))



    #!# keep track of symmetric-Hf swapping:
    list_sym_Hf_swap_cfis = []
    ## permutation list to get correct symmetry for basis coefficients:
    sym_plist = np.array([[0,1],[0,1]])


    #!# for residence times (need one entry for each at1):
    within_range_bools = [False for at1i in range(len(at1_list))]
    num_prev_within    = [0 for at1i in range(len(at1_list))]


    for cfi, cf in enumerate(xyz_obj_list):
      lbox = lbox_list[cfi]

      #reimcrds = self.PBC_reimage_coordinates(cf.crds[at2], cf.crds, lbox)

      for at1i, at1 in enumerate(at1_list):
        #
        at1mi = at1i + initial_at1mi

        ## reimage z coordinates relative to at1 maintaining bonds:
        reim_z_vec = self.Get_PBC_reimage_vector(cf.crds[at1], cf.crds[hsind], lbox)
        reim_z = cf.crds[zinds] - reim_z_vec

        ## test:
        #print(zinds)
        #print(cf.crds[zinds], '\n')
        #print(reim_z, '\n\n\n')
        #exit()


        #!#! get approximate C2 axis:

        ## get O-Hf (position) direction vectors:
        oh_directions = np.zeros((2,2,3))
        ## loop over oxygens:
        for oii, oi in enumerate(oinds):
          ## loop over flanking-hydrogens:
          for hii, hi in enumerate(hfinds[oii]):
            #?oh_directions[oi][hii] = (reim_z[hi] - reim_z[oi]) / LA.norm(reim_z[hi] - reim_z[oi])
            oh_directions[oii][hii] += (reim_z[hi] - reim_z[oi]) / LA.norm(reim_z[hi] - reim_z[oi])


        #!#################################
        #!#! (9/18/20) let's try to find symmetry of Hf (assuming C2 minimum-like geom):


        ## project-out components of oh_directions parallel to oo:
        oo  = reim_z[oinds[1]] - reim_z[oinds[0]]
        oo /= LA.norm(oo)

        oh_dirs_perp_oo = copy(oh_directions)
        for oi in range(2):
          for hi in range(2):
            oh_dirs_perp_oo[oi][hi] -= oo * np.dot(oh_dirs_perp_oo[oi][hi], oo)
            oh_dirs_perp_oo[oi][hi] /= LA.norm(oh_dirs_perp_oo[oi][hi])


        ## symmetry is correct if cross products (each within same o) point in opposite directions:
        cross1 = np.cross(oh_dirs_perp_oo[0][sym_plist[0][0]], oh_dirs_perp_oo[0][sym_plist[0][1]])
        cross2 = np.cross(oh_dirs_perp_oo[1][sym_plist[1][0]], oh_dirs_perp_oo[1][sym_plist[1][1]])

        cross1 /= LA.norm(cross1)
        cross2 /= LA.norm(cross2)

        ## if we need to permute:
        if np.sign(np.dot(cross1, cross2)) > 0.:
          ## add config index to list to keep a record:
          list_sym_Hf_swap_cfis.append(cfi)
          ## permute list for o1:
          sym_plist[0] = sym_plist[0][[1,0]]

          # test:
          #if cfi > 0:
          #  print(cfi, 'sym_plist:\n', sym_plist)
          #  print(cfi, '\ncross products:', cross1, cross2)
          #  print('dot product:', np.dot(cross1, cross2))
          #  #print('prev crds:\n', zcrds[cfi-1])
          #  #print('curr crds:\n', zcrds[cfi], '\n\n=====\n')

          #  #exit()


          #!# TEST:
          cross1 = np.cross(oh_dirs_perp_oo[0][sym_plist[0][0]], oh_dirs_perp_oo[0][sym_plist[0][1]])
          cross2 = np.cross(oh_dirs_perp_oo[1][sym_plist[1][0]], oh_dirs_perp_oo[1][sym_plist[1][1]])
          if np.sign(np.dot(cross1, cross2)) > 0.:
            sys.exit('permutation did not work???')



        #!#! find approximate C2 axis (for shared-H decomposition):
        appx_C2_1 = np.zeros(3)
        appx_C2_2 = np.zeros(3)
        ## vector sums of (approximately) symmetric Hf OH-directions (perp to oo):
        appx_C2_1 += oh_dirs_perp_oo[0][sym_plist[0][0]] + oh_dirs_perp_oo[1][sym_plist[1][0]]
        appx_C2_1 /= LA.norm(appx_C2_1)

        appx_C2_2 += oh_dirs_perp_oo[0][sym_plist[0][1]] + oh_dirs_perp_oo[1][sym_plist[1][1]]
        appx_C2_2 /= LA.norm(appx_C2_2)

        ## constructive sum:
        C2_sum_sign = np.sign(np.dot(appx_C2_1, appx_C2_2))
        appx_C2  = appx_C2_1 + (C2_sum_sign * appx_C2_2)
        appx_C2 /= LA.norm(appx_C2)

        #!# also choose sign with makes C2 point in direction most similar to Hf:
        #C2_direction_dot_prods = []
        C2_direction = 0.0
        for oii in range(2):
          for hii in range(2):
            C2_direction += np.dot(oh_dirs_perp_oo[oi][hi], appx_C2)
            #C2_direction_dot_prods.append(np.dot(oh_dirs_perp_oo[oi][hi], appx_C2))

        #!# C2_direction should be positive if C2 points with Hf, negative if points away:
        appx_C2 *= np.sign(C2_direction)


        ## vector pointing from Hs to at1:
        vec_hs_at1 = cf.crds[at1] - reim_z[hsind]
        ## 12 distance:
        dist12 = LA.norm(vec_hs_at1)


        ## project-out component parallel to oo:
        vec_hs_at1_perp_oo  = copy(vec_hs_at1)
        vec_hs_at1_perp_oo -= oo * np.dot(vec_hs_at1, oo)
        ## normalize:
        vec_hs_at1_perp_oo /= LA.norm(vec_hs_at1_perp_oo)


        ## calculate angle:
        ang123 = (180.0 / np.pi) * ctools.angle(appx_C2, np.zeros(3), vec_hs_at1_perp_oo)


        #!# for residence times (5.5 and 120 based on RADFs):
        if dist12 < 5.5 and ang123 > 120.:
          within_range_bools[at1i] = True
          num_prev_within[at1i]   += 1
        else:
          within_range_bools[at1i] = False

        ## keep track of molecule indices within range:
        if num_prev_within[at1i] > 0 and not within_range_bools[at1i]:

          #!# list of residence time within bounds:
          #with open(respref + str(at2i + ind_first_ats23_mol) + '.dat', 'a') as f:
          with open(respref + '_at1_' + str(at1mi) + '.dat', 'a') as f:
            ## inclusive initial, final indices within bounds:
            indi = cfi - num_prev_within[at1i]
            ## remember at2 was within bounds up to (not including) current cfi:
            indf = cfi - 1
            ## note: '{{' is escaped to a single '{' in an f-string:
            f.write(f'{{{indi}..{indf}}}\n')

          ## reset num_prev:
          num_prev_within[at1i] = 0



        ## skip if out of range:
        if dist12 > max_dist or dist12 < min_dist:
          continue


        ## find bin:
        radbinindex = int(np.floor((dist12 - min_dist) / dR))
        angbinindex = int(np.floor(ang123 / dA))

        binindex = (radbinindex * num_ang_bins) + angbinindex

        RADF[binindex][2] += 1.0


    ## divide by volume and normalize:
    for ri in range(num_rad_bins):
      r_lw = min_dist + (ri * dR)
      r_up = r_lw + dR

      #!# ASSUMES CUBIC BOX !!!!!!

      ## volume of entire radial shell:
      V_r = (4.0 * np.pi / 3.0) * (pow(r_up, 3) - pow(r_lw, 3))

      for ai in range(num_ang_bins):

        #!# need this if polar-like angle (not azimuthal-like)
        #!### fraction of radial shell for angular bin, ai (from area of sherical cap):
        #!#frac_V = 0.5 * (np.cos(ai * dA * np.pi / 180.) - np.cos((ai+1) * dA * np.pi / 180.))

        ## azimuthal-like angle:
        frac_V = 1. / num_ang_bins  ## only need to do this once

        Vol_bin = V_r * frac_V

        ## to first order (from Jacobian determinant, after integrating out az-angle):
        #Vol_bin = 360.0 * pow(r_lw + (dR / 2.0), 2) * dR * np.sin((ai + 0.5) * (dA * np.pi / 180.0)) * dA
#                 ^^^^^ 360 from integrating out azimuthal angle (in degrees)

        #!# ASSUMES FIXED VOLUME:
        coef_ri_ai = norm * pow(lbox_list[0], 3) / Vol_bin


        binindex = (ri * num_ang_bins) + ai
        RADF[binindex][2] *= coef_ri_ai



    ## write RADF to file:
    #np.savetxt(output_path, RADF, fmt='%.10f')
    np.savetxt(output_path, RADF)

    return

#================================================#
#================================================#


##################################################
###  bin Zundel dihedral  ########################

  def Bin_zundel_dihedral(self,
                          in_xyz_path,
                          dihedral_ind_path,
                          num_bins,
                          outxyz_prefix):

    self.load_modules()
    ## convert trajectory into xyz objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    ## get box length for each configuration from comment line:
    lbox_list = [float(cf.comment.split(' ')[2]) for cf in xyz_obj_list]

    ## get atom indices for dihedral:
    at_inds_dih = np.genfromtxt(dihedral_ind_path, dtype=int)

    ## get angle increments:
    num_bins = int(num_bins)
    dA       = 180.0 / float(num_bins)
    ## bins for molecule indices:
    dih_mol_inds = [[] for n in range(num_bins)]


    ## indices for dihedral:
    h11i = at_inds_dih[0][0]
    h12i = at_inds_dih[0][1]
    o1i  = at_inds_dih[0][2]
    h21i = at_inds_dih[1][0]
    h22i = at_inds_dih[1][1]
    o2i  = at_inds_dih[1][2]
    ## loop over configurations:
    for cfi, cf in enumerate(xyz_obj_list):

      ## need to reimage needed atoms relative to any atom:
      reimcrds = self.PBC_reimage_coordinates(cf.crds[h11i], cf.crds[[h11i, h12i, o1i, h21i, h22i, o2i]], lbox_list[cfi])

      h11 = reimcrds[0]
      h12 = reimcrds[1]
      o1  = reimcrds[2]

      h21 = reimcrds[3]
      h22 = reimcrds[4]
      o2  = reimcrds[5]

      ## find bisectors:
      oh1_1 = h11 - o1
      oh2_1 = h12 - o1
      oh1_1 /= LA.norm(oh1_1)
      oh2_1 /= LA.norm(oh2_1)
      b1 = (oh1_1 + oh2_1) / 2.0
      b1 += o1

      oh1_2 = h21 - o2
      oh2_2 = h22 - o2
      oh1_2 /= LA.norm(oh1_2)
      oh2_2 /= LA.norm(oh2_2)
      b2 = (oh1_2 + oh2_2) / 2.0
      b2 += o2

      dih = ctools.dihedral(b1, o1, o2, b2) * (180.0 / np.pi)


      ## find bin index:
      binindex = int(np.floor(dih / dA))
      dih_mol_inds[binindex].append(cfi)



    ## loop over bins:
    for bni, bn in enumerate(dih_mol_inds):
      A_str = '{:0.5f}'.format(bni * dA)
      outxyz = outxyz_prefix + '_' + A_str + '.xyz'

      with open(outxyz, 'w') as f:
        ## loop over configuration indices in bin:
        for cfi in bn:
          xyz_obj_list[cfi].write_xyz_file(out_file=f)

#================================================#
#================================================#


##################################################
###  calculate C_2(t) slowly  ####################

  def Calculate_C2_slowly(self,
                          vector_list):

    out_C2 = np.zeros(len(vector_list))
    ## loop over vectors:
    for v1i, v1 in enumerate(vector_list):
      ## loop over time delays relative to v:
      for tau, v2 in enumerate(vector_list[v1i:]):
        #!# ASSUMES VECTORS ARE NORMALIZED !!!!!
        out_C2[tau] += 0.5 * ((3. * pow(np.dot(v1,v2),2)) - 1)

    return out_C2

#================================================#
#================================================#


##################################################
###  (test) calculate C_1(t) slowly  #############

  def Calculate_C1_slowly(self,
                          vector_list):

    out_C1 = np.zeros(len(vector_list))
    ## loop over vectors:
    for v1i, v1 in enumerate(vector_list):
      ## loop over time delays relative to v:
      for tau, v2 in enumerate(vector_list[v1i:]):
        #!# ASSUMES VECTORS ARE NORMALIZED !!!!!
        out_C1[tau] += np.dot(v1,v2)

    return out_C1

#================================================#
#================================================#


##################################################
###  ACN OCF  ####################################

  def ACN_OCF(self,
              in_xyz_path,
              at1_list_path,
              at2_list_path,
              at3_list_path,
              min_dist,
              max_dist,
              lbox,
              minimum_length_ocf,
              timestep,
              ind_first_ats23_mol):


    self.load_modules()

    min_dist = float(min_dist)
    max_dist = float(max_dist)

    ## read xyz-file into list of 'xyz' objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    at1_list = np.genfromtxt(at1_list_path, dtype=int) - 1
    at2_list = np.genfromtxt(at2_list_path, dtype=int) - 1
    at3_list = np.genfromtxt(at3_list_path, dtype=int) - 1

    if len(at1_list.shape) == 0:
      at1_list = np.asarray([at1_list])
    if len(at2_list.shape) == 0:
      at1_lis2 = np.asarray([at2_list])
    if len(at3_list.shape) == 0:
      at3_list = np.asarray([at3_list])


    minimum_length_ocf = int(minimum_length_ocf)
    timestep = float(timestep)

    #!# keep list of residence time within distance bounds for each molecule:
    ind_first_ats23_mol = int(ind_first_ats23_mol)
    respref = 'residence_time_from_ocf_'

    #!# constant pressure:
    lbox_list = [float(lbox) for cf in xyz_obj_list]


    #!# do we want fixed coordinate frame?
    #COORD_FRAME_BOOL = True
    COORD_FRAME_BOOL = False


    ## keep list of all OCF segments:
    all_OCF_segments = []

    for at1i, at1 in enumerate(at1_list):
      #!# also need index of other at1 (ASSUMES ONLY 2)
#!#Hs      other_at1 = at1_list[(at1i+1)%2]

      for at2i, at2 in enumerate(at2_list):
        if at1 == at2:
          continue

        #!# list of residence time within bounds:
        with open(respref + 'at1_' + str(at1i) + '_' + str(at2i + ind_first_ats23_mol) + '.dat', 'a') as f:
          f.write(f'#!#!#!#! START (timestep: {timestep} ps; inds start with 0; inds inclusive) #!#!#!#!\n')

        ## index of at3 corresponding to at2:
        at3 = at3_list[at2i]

        ## assume at2 not in 1st solvation shell of at1:
        ss1_bool = False
        ## list of at3 direction (unit) vectors in chosen coordinate system:
        at3_direction_traj = []

        for cfi, cf in enumerate(xyz_obj_list):
          lbox = lbox_list[cfi]

          ## minimum image convention to determine if in 1st solvation shell:
          dvec = cf.crds[at1] - cf.crds[at2]
          ## PBC:
          dist12 = LA.norm(dvec - (np.array([round(xy / lbox) for xy in dvec]) * lbox))


          #!# make sure not 1st solvation shell to other at1:

          ## reimage at2 and other_at1 relative to at1
          reim_at2     = self.PBC_reimage_coordinates(cf.crds[at1], cf.crds[at2], lbox)[0]
#!#Hs          reim_other_at1 = self.PBC_reimage_coordinates(cf.crds[at1], cf.crds[other_at1], lbox)[0]
#!#Hs          ## angle:
#!#Hs          angle_other_at1_at2_at1 = ctools.angle(reim_other_at1, cf.crds[at1], reim_at2)
#!#Hs          ## degrees:
#!#Hs          angle_other_at1_at2_at1 *= 180. / np.pi


          #!# HARD-CODED minimum angle (to see if in 1st ss of other at1):
          min_ang_at1_ss = 64.0

          ## skip if out of range:
          if dist12 > max_dist or dist12 < min_dist:
#!#Hs          if (dist12 > max_dist) or (dist12 < min_dist) or (angle_other_at1_at2_at1 < min_ang_at1_ss):

            ## if at2 in 1st solvation shell in previous steps in trajectory:
            if ss1_bool:

              #!# list of residence time within bounds:
              #with open(respref + str(at2i + ind_first_ats23_mol) + '.dat', 'a') as f:
              with open(respref + 'at1_' + str(at1i) + '_' + str(at2i + ind_first_ats23_mol) + '.dat', 'a') as f:
                ## inclusive initial, final indices within bounds:
                indi = cfi - len(at3_direction_traj)
                ## remember at2 was within bounds up to (not including) current cfi:
                indf = cfi - 1
                ## note: '{{' is escaped to a single '{' in an f-string:
                f.write(f'{{{indi}..{indf}}}\n')


              ## enforce minimum length of OCF:
              if len(at3_direction_traj) >= minimum_length_ocf:

                at3_direction_traj = np.asarray(at3_direction_traj)

                ## calculate OCF for this segment:
                tmp_OCF = self.Calculate_C2_slowly(at3_direction_traj)
                #tmp_OCF_array = signal.fftconvolve(at3_direction_traj, at3_direction_traj[::-1], mode='full', axes=0)
                ### equivalent to doing dot product inside CF calculation:
                #tmp_OCF = np.sum(tmp_OCF_array, axis=1)[len(at3_direction_traj)-1:]
#                        only want positive time delays ^^^^^^^^^^^^^^^^^^^^^^^^^^


                #!# maybe don't normalize here? we need to add these to others...
                #norm_OCF = np.array([float(len(at3_direction_traj) - T) for T in range(len(at3_direction_traj))])
                ## add to list:
                all_OCF_segments.append(tmp_OCF)
                #all_OCF_segments.append(tmp_OCF / norm_OCF)


              ## reset in case at2 returns to 1st solvation shell w/r/t at1:
              ss1_bool = False
              at3_direction_traj = []

            continue

          else:
            ss1_bool = True


          ## do we want to keep fixed coordinate frame?
          if COORD_FRAME_BOOL:
            ## reimage at2 and at3 wrt at1, unbreaking at2-at3 bond:
            reim_crds = np.zeros((4, len(cf.crds[at1])))
            ## at1 unchanged:
            reim_crds[0] += cf.crds[at1]
            #!# also need to reimage "other at1":
            reim_crds[1] += self.PBC_reimage_coordinates(cf.crds[at1], cf.crds[other_at1], lbox)[0]

            ##!# test:
            #print('\n\n\n',reim_crds,'\n')
            #print(reim_crds[2:],'\n')
            #print(cf.crds[at1],'\n')
            #print(cf.crds[at2],'\n')
            #print(cf.crds[at3],'\n')
            #print(self.Reimage_molecule_PBC(cf.crds[at1], cf.crds[at2], cf.crds[at3], lbox))


            reim_crds[2:] += self.Reimage_molecule_PBC(cf.crds[at1], cf.crds[at2], cf.crds[at3], lbox)

            ## translate so at1 at origin & rotate so other_at1 on x-axis, and at2 aligned to z-axis:
            reim_crds = self.Align_to_axes(reim_crds, 0, 1, 2)
            ## now translate so at2 at origin:
            reim_crds -= reim_crds[2]

            ## save direction (unit) vector of at3 in this coordinate system:
            at3_direction_traj.append(reim_crds[3] / LA.norm(reim_crds[3]))


          ## else don't want fixed coordinate frame:
          else:

            ## reimage at3 relative to at2:
            reim_at3_crds = self.PBC_reimage_coordinates(cf.crds[at2], cf.crds[at3], lbox)[0]
            ## take difference vector:
            diff32 = reim_at3_crds - cf.crds[at2]
            ## save direction (unit) vector:
            at3_direction_traj.append(diff32 / LA.norm(diff32))




          #!# if end of trajectory and in 1st solvation shell, calculate OCF:
          if cfi == len(xyz_obj_list)+1:
            if ss1_bool:

              ## enforce minimum length of OCF:
              if len(at3_direction_traj) >= minimum_length_ocf:

                at3_direction_traj = np.asarray(at3_direction_traj)

                ## calculate OCF for this segment:
                tmp_OCF = self.Calculate_C2_slowly(at3_direction_traj)
                #tmp_OCF_array = signal.fftconvolve(at3_direction_traj, at3_direction_traj[::-1], mode='full', axes=0)
                ### equivalent to doing dot product inside CF calculation:
                #tmp_OCF = np.sum(tmp_OCF_array, axis=1)[len(at3_direction_traj)-1:]
#                        only want positive time delays ^^^^^^^^^^^^^^^^^^^^^^^^^^

                #!# maybe don't normalize here? we need to add these to others...
                #norm_OCF = np.array([float(len(at3_direction_traj) - T) for T in range(len(at3_direction_traj))])
                ## add to list:
                all_OCF_segments.append(tmp_OCF)
                #all_OCF_segments.append(tmp_OCF / norm_OCF)



    ### finally, average OCF segments:

    try:
      ## find largest segment length:
      max_seglength = max([len(ocfseg) for ocfseg in all_OCF_segments])
    except ValueError:
      exit()

    sum_ocfsegs = np.zeros(max_seglength)
    for ocfseg in all_OCF_segments:
      #print(len(ocfseg), type(ocfseg), ocfseg[:4])
      for ii, ocfval in enumerate(ocfseg):
        sum_ocfsegs[ii] += ocfval
      #sum_ocfsegs += copy(ocfseg).resize(max_seglength)

    ## number of segments:
    nsegs = float(len(all_OCF_segments))

    #print('numberofsegments:' + str(nsegs))

    ## write to stdout:
    for ti, ocfval in enumerate(sum_ocfsegs):
      t   = ti * timestep
      val = ocfval #/ nsegs
      print(f'{t}  {val}')


    return


#================================================#
#================================================#



##################################################
###  ACN ACF  ####################################

  def ACN_ACF(self,
              in_xyz_path,
              at1_list_path,
              at2_list_path,
              at3_list_path,
              min_dist,
              max_dist,
              lbox,
              minimum_length_acf,
              timestep):
              #timestep,
              #const_pressure_bool=False):


    self.load_modules()
    min_dist = float(min_dist)
    max_dist = float(max_dist)

    ## read xyz-file into list of 'xyz' objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    at1_list = np.genfromtxt(at1_list_path, dtype=int) - 1
    at2_list = np.genfromtxt(at2_list_path, dtype=int) - 1
    at3_list = np.genfromtxt(at3_list_path, dtype=int) - 1

    minimum_length_acf = int(minimum_length_acf)
    timestep = float(timestep)


    ##!# constant pressure:
    #if bool(int(const_pressure_bool)):
    #  sys.exit('ERROR: not implemented... exiting...')
    #!#
    lbox_list = [float(lbox) for cf in xyz_obj_list]




    ## keep list of all ACF segments:
    all_ACF_segments = []

    for at1 in at1_list:
      for at2i, at2 in enumerate(at2_list):
        if at1 == at2:
          continue

        ## assume at2 not in 1st solvation shell of at1:
        ss1_bool = False
        ## new angle trajectory list:
        ang_traj = []

        for cfi, cf in enumerate(xyz_obj_list):
          lbox = lbox_list[cfi]


          dvec = cf.crds[at1] - cf.crds[at2]
          ## PBC:
          #dist12 = LA.norm(dvec - (round(dvec / lbox) * lbox))
          #!#
          dist12 = LA.norm(dvec - (np.array([round(xy / lbox) for xy in dvec]) * lbox))

          ## skip if out of range:
          if dist12 > max_dist or dist12 < min_dist:

            ## if at2 in 1st solvation shell in previous steps in trajectory:
            if ss1_bool:

              ## enforce minimum length of ACF:
              if len(ang_traj) >= minimum_length_acf:

                ## calculate ACF for this segment:
                tmp_ACF = signal.fftconvolve(ang_traj, ang_traj[::-1], mode='full')[len(ang_traj)-1:]
#                                                    only want positive time delays ^^^^^^^^^^^^^^^^

                norm_ACF = np.array([float(len(ang_traj) - T) for T in range(len(ang_traj))])
                ## add to list:
                all_ACF_segments.append(tmp_ACF / norm_ACF)


              ## reset in case at2 returns to 1st solvation shell w/r/t at1:
              ss1_bool = False
              ang_traj = []

            continue

          else:
            ss1_bool = True



          ## reimage at3 within same (ACN) molecule relative to at2:
          reim1_at3 = self.PBC_reimage_coordinates(cf.crds[at2], np.array([cf.crds[at3_list[at2i]]]), lbox)[0]

          ## get reimaging vector:
          reim_vec = self.Get_PBC_reimage_vector(cf.crds[at1], cf.crds[at2], lbox)

          ### now reimage relative to at1:
          #reim_acn = self.PBC_reimage_coordinates(cf.crds[at2], np.array([reim_acn]), lbox)[0]

          reim_acn = np.zeros((2,3))
          reim_acn[0] = cf.crds[at2] - reim_vec
          reim_acn[1] = reim1_at3 - reim_vec

          ## calculate cosine of angle:
          cos_ang = np.cos(ctools.angle(cf.crds[at1], reim_acn[0], reim_acn[1]))
          ## add to trajectory and continue to next step:
          ang_traj.append(cos_ang)


          #!# if end of trajectory and in 1st solvation shell, calculate ACF:
          if cfi == len(xyz_obj_list)+1:
            if ss1_bool:

              ## enforce minimum length of ACF:
              if len(ang_traj) >= minimum_length_acf:

                ## calculate ACF for this segment:
                tmp_ACF = signal.fftconvolve(ang_traj, ang_traj[::-1], mode='full')[len(ang_traj)-1:]
#                                                    only want positive time delays ^^^^^^^^^^^^^^^^

                norm_ACF = np.array([float(len(ang_traj) - T) for T in range(len(ang_traj))])
                ## add to list:
                all_ACF_segments.append(tmp_ACF / norm_ACF)


    ### finally, average ACF segments:

    ## find largest segment length:
    max_seglength = max([len(acfseg) for acfseg in all_ACF_segments])

    sum_acfsegs = np.zeros(max_seglength)
    for acfseg in all_ACF_segments:
      #print(len(acfseg), type(acfseg), acfseg[:4])
      for ii, acfval in enumerate(acfseg):
        sum_acfsegs[ii] += acfval
      #sum_acfsegs += copy(acfseg).resize(max_seglength)

    ## number of segments:
    nsegs = float(len(all_ACF_segments))

    ## write to stdout:
    for ti, acfval in enumerate(sum_acfsegs):
      t   = ti * timestep
      val = acfval / nsegs
      print(f'{t}  {val}')


    return


#================================================#
#================================================#


##################################################
###  extract fragments  ##########################

  def Extract_fragments(self,
                        in_xyz_path,
                        natoms_per_mol_list_path,
                        nmols_per_frag,
                        inner_cutoff_distance,
                        outer_cutoff_distance,
                        lbox,
                        ref_mol_index_path=None,
                        const_pressure_bool=False,
                        all_pairs_bool=False,
                        CLOSEST_MOLS_BOOL=False,
                        atom_ind_dist_list_path=None): #,


    self.load_modules()
    #!# ASSUMES CUBIC BOX !!!!!!!!!

    nmols_per_frag        = int(nmols_per_frag)
    inner_cutoff_distance = float(inner_cutoff_distance)
    outer_cutoff_distance = float(outer_cutoff_distance)

    ## read xyz-file into list of 'xyz' objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)


    #!# constant pressure:
    if bool(int(const_pressure_bool)):
      ## get box length for each configuration from comment line:
      lbox = [float(cf.comment.split(' ')[2]) for cf in xyz_obj_list]
    #!# constant volume, so same lbox for each configuration:
    else:
      lbox = [float(lbox) for cf in xyz_obj_list]


    ## read file with number of atoms per molecule:
    natoms_per_mol_list = np.genfromtxt(natoms_per_mol_list_path, dtype=int)
    try:
      nmols = len(natoms_per_mol_list)
    #!# genfromtxt yields float when 1 entry in file:
    except TypeError:
      nmols = int(xyz_obj_list[0].n_atoms / natoms_per_mol_list)
      natoms_per_mol_list = np.array([natoms_per_mol_list for mo in range(nmols)])

    ## do we want distances between first atoms of each molecule?
    if atom_ind_dist_list_path is None:
      atom_ind_dist_list = np.array([0 for moi in range(nmols)])
    else:
      ### else get list from file

      ## assume comma-delimited, and make list of lists:
      at_inds_for_dist = []
      with open(atom_ind_dist_list_path, 'r') as f:
        for moi, li in enumerate(f.readlines()):
          ## atom index for last atom in previous molecule:
          prev_ati = int(sum(natoms_per_mol_list[:moi]))
          at_inds_for_dist.append([prev_ati + int(ati) for ati in li.strip().split(',')])



    ## read 'reference' molecule index list:
    if ref_mol_index_path is None or ref_mol_index_path == 'all':
      ref_mol_list = list(range(nmols-1))
    else:
      ref_mol_list = sorted(np.genfromtxt(ref_mol_index_path, dtype=int).flatten())


    ### now loop over configs in xyz:
    for xyi, xy in enumerate(xyz_obj_list):

      ## first reimage each molecule relative to distance atom in that molecule:
      nobrokenbonds_crds = copy(xy.crds)
      for moi in range(nmols):
        ## upper and lower bounds for atoms in molecule, moi:
        lw_ati = sum(natoms_per_mol_list[:moi])
        up_ati = lw_ati + natoms_per_mol_list[moi]


        ##!# debug:
        #print(lw_ati, up_ati)
        #print(nobrokenbonds_crds[lw_ati:up_ati])

        #!# use first atom in at_inds_for_dist inner-list for reimage:
        nobrokenbonds_crds[lw_ati:up_ati] = self.PBC_reimage_coordinates(nobrokenbonds_crds[at_inds_for_dist[moi][0]], nobrokenbonds_crds[lw_ati:up_ati], lbox[xyi])

        ##!# debug:
        #print(nobrokenbonds_crds[lw_ati:up_ati])


      ## loop over 'reference' molecule indices (center of reimage):
      for ref_moi in ref_mol_list:

        ## index of first atom in molecule, ref_moi:
        #!# use first atom in at_inds_for_dist inner-list for reimage:
        ref_ati = at_inds_for_dist[ref_moi][0]
        ## copy xyz object to fill with reimaged coordinates:
        reim_xy = copy(xy)


        #!# THIS WILL NOT BREAK BONDS:
        for moi in range(nmols):  # NOTE: also want to do this for ref (if bonds broken)
          #!# use first atom in at_inds_for_dist inner-list for reimage:
          reimvec = self.Get_PBC_reimage_vector(nobrokenbonds_crds[ref_ati], nobrokenbonds_crds[at_inds_for_dist[moi][0]], lbox[xyi])

          ## lower atom index in moi:
          lw_ati = sum(natoms_per_mol_list[:moi])
          up_ati = lw_ati + natoms_per_mol_list[moi]

          ### loop over atoms in molecule moi:
          reim_xy.crds[lw_ati:up_ati] = nobrokenbonds_crds[lw_ati:up_ati] - reimvec

        reim_crds = reim_xy.crds


        ## list of indices of molecules within threshold:
        mol_inds_within_cutoff = []
        ## also keep list of average distances for sorting:
        avg_dist_for_thresh_list = []
        #!# if False, then only get molecules indexed AFTER reference, otherwise get ALL (careful w/ double-counting):
        if bool(int(all_pairs_bool)):
          ### now find distances between cf.crds[ref_ati] and 1st atom of EACH molecule, not including reference:

          #!# use average distance for threshold:
          for moi, ati_list in enumerate(at_inds_for_dist):
            ## get list of all needed distances to average:
            avg_dist_to_check = np.mean([ctools.distance(reim_crds[ati1], reim_crds[ati2])
                                         for ati1 in at_inds_for_dist[ref_moi]
                                         for ati2 in ati_list])

            ## check threshold:
            if avg_dist_to_check < outer_cutoff_distance and avg_dist_to_check > inner_cutoff_distance:
              mol_inds_within_cutoff.append(moi)
              ## add avg distance to list:
              avg_dist_for_thresh_list.append(avg_dist_to_check)


        else:
          ### now find distances between cf.crds[ref_ati] and 1st atom of each molecule AFTER (if all_pairs_bool=False) reference:

          #!# use average distance for threshold:
          for moi, ati_list in enumerate(at_inds_for_dist[(ref_moi + 1):]):
            ## get list of all needed distances to average:
            avg_dist_to_check = np.mean([ctools.distance(reim_crds[ati1], reim_crds[ati2])
                                         for ati1 in at_inds_for_dist[ref_moi]
                                         for ati2 in ati_list])

            ## check threshold:
            if avg_dist_to_check < outer_cutoff_distance and avg_dist_to_check > inner_cutoff_distance:
              mol_inds_within_cutoff.append(moi+ref_moi+1)
              ## add avg distance to list:
              avg_dist_for_thresh_list.append(avg_dist_to_check)


        #!# now sort mol_inds_within_cutoff by distance:
        ## only saving avg dists for mols within thresh, so use argsort as permutation list:
        mol_inds_within_cutoff = np.array(mol_inds_within_cutoff)[np.argsort(avg_dist_for_thresh_list)]


        ## list of atom indices for reference molecule ## each fragment (starting with atoms in first molecule):
        ref_at_inds = list(range(int(sum(natoms_per_mol_list[:ref_moi])), int(sum(natoms_per_mol_list[:ref_moi+1]))))


        #!# track skipped frames because too few (non-ref) molecules within cutoff:
        if len(mol_inds_within_cutoff) < (nmols_per_frag - 1):
          with open('skipped_frames.log', 'a') as skipped_log:
            skipped_log.write('config: ' + str(xyi) + ';  ref mol: ' + str(ref_moi) + '\n')
          #!#
          continue


        ## now make nmols_per_frag-combinations of molecule indices (within cutoff):
        for combo in itts.combinations(mol_inds_within_cutoff, nmols_per_frag-1):
          other_at_inds = []
          ## atom indices for other molecules:
          for moi in combo:
            other_at_inds += list(range(int(sum(natoms_per_mol_list[:moi])), int(sum(natoms_per_mol_list[:moi+1]))))


          ## comment with timestep number and molecule combo:
          comment = 'timestep: ' + str(xyi) + ';  molecules: ' + str([ref_moi] + list(combo))
          ## write xyz file for this combo using the corresponding atom indices:
          reim_xy.write_xyz_file(plist_in = ref_at_inds + other_at_inds, new_comment = comment, append_comment=True)

          #!# if we only want CLOSEST MOLECULES, then skip other combos:
          if bool(int(CLOSEST_MOLS_BOOL)):
            ## break out of 'combo' loop:
            break

    return

#================================================#
#================================================#


##################################################
###  translate molecule  #########################

  def Translate_molecule(self,
                         in_xyz_path,
                         increment_disp_3N_vector_path,
                         nsteps,
                         inclusive_bool=True):

    self.load_modules()
    ## read xyz-file into list of 'xyz' objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    #!# read axis (space-delimited):
    inc_dvec = np.genfromtxt(increment_disp_3N_vector_path)

    #!# make sure axis has same shape as (FIRST) xyz.crds:
    inc_dvec = inc_dvec.reshape(xyz_obj_list[0].crds.shape)

    ## loop over configurations in xyz_obj_list:
    for xy in xyz_obj_list:
      ## make copy of 'xyz' object:
      tmp_xyz = copy(xy)

      ## do we want to print 0th step?
      if bool(int(inclusive_bool)):
        tmp_xyz.write_xyz_file()

      for s in range(int(nsteps)):
        tmp_xyz.crds += inc_dvec
        tmp_xyz.write_xyz_file()


#================================================#
#================================================#


##################################################
###  align to axes  ##############################

  def Align_to_axes(self,
                    coords_to_align,
                    origin_ati,
                    x_ati=None,
                    z_ati=None):

    self.load_modules()
    ## make a copy of input coordinates:
    out_coords = np.array(copy(coords_to_align))


    ## translate to origin:
    out_coords -= out_coords[origin_ati]
    ## are we done?
    if x_ati == None:
      return out_coords


    ### align to x-axis:
    x_axis = np.array([1.,0.,0.])

    ## find vector normal to plane containing x-axis, orig_at, and x_ax_at:
    nvec = np.cross((out_coords[x_ati] / LA.norm(out_coords[x_ati])), x_axis)

    ## don't rotate if already on x-axis:
    if LA.norm(nvec) != 0.0:
      ## angle (in radians) between a-axis and x_ax_at:
      ang = np.arcsin(LA.norm(nvec))

      #!# make sure to rotate into positive x-axis??
      if np.sign(out_coords[x_ati][0]) < 0:
        ang += np.pi

      ## "a" factor:
      a = np.cos(ang / 2)
      ## "omega" vector:
      omega = np.sign(out_coords[x_ati][0]) * np.sin(ang / 2) * (nvec / LA.norm(nvec))

      ## loop over atoms, & rotate each via cross-product formula for Euler-Rodriguez:
      for at in out_coords:
        t1 = 2 * a * np.cross(omega, at)
        t2 = 2     * np.cross(omega, np.cross(omega, at))

        at += t1 + t2

    ## are we done?
    if z_ati == None:
      return out_coords


    ### align to z-axis:

    ang = np.arccos(out_coords[z_ati][2] / LA.norm(np.array([out_coords[z_ati][1], out_coords[z_ati][2]])))

    a = np.cos(ang / 2)
    omega = np.sign(out_coords[z_ati][1]) * np.sin(ang / 2) * x_axis

    ## loop over atoms, & rotate each via cross-product formula for Euler-Rodriguez:
    for at in out_coords:
      t1 = 2 * a * np.cross(omega, at)
      t2 = 2     * np.cross(omega, np.cross(omega, at))

      at += t1 + t2


    return out_coords

#================================================#
#================================================#


##################################################
###  rotate molecule  ############################

  def Rotate_molecule(self,
                      in_xyz_path,
                      axis_path,
                      angle,
                      atom_indices_to_rotate_path=None):

    self.load_modules()

    #!#  angle MUST BE IN UNITS: radians/2pi
    #!#
    #!#  i.e. fraction of unit circle !!!!!!!!!

    ## read xyz-file into list of 'xyz' objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    #!# read axis (space-delimited):
    with open(axis_path) as f:
      axis = np.array(list(map(float, f.readline().strip().split())))


    ## convert angle into float:
    angle = float(angle)

    if atom_indices_to_rotate_path == None:
      atom_indices_to_rotate = list(range(xyz_obj_list[0].n_atoms))
    else:
      ## read file with atom indices (ints starting with 0) to rotate:
      atom_indices_to_rotate = np.genfromtxt(atom_indices_to_rotate_path, dtype=int)

      #!# if only want to rotate single atom:
      try:
        fakevar = atom_indices_to_rotate[1]
      except:
        atom_indices_to_rotate = np.array([atom_indices_to_rotate])


    ## loop over 'xyz' objects:
    for xy in xyz_obj_list:
      ## now overwrite xyz_obj.crds with rotated coordinates:
      xy.crds[atom_indices_to_rotate,:] = ctools.axis_angle(xy.crds[atom_indices_to_rotate,:], axis, angle)
      ## now write rotated xyz to stdout:
      xy.write_xyz_file()

    return

#================================================#
#================================================#


##################################################
###  rotate scan  ################################

  def Rotate_scan(self,
                  in_xyz_path,
                  axis_path,
                  incremental_angle,
                  nsteps,
                  atom_indices_to_rotate_path=None,
                  inclusive_bool=True):

    self.load_modules()

    #!#  angle MUST BE IN UNITS: radians/2pi
    #!#
    #!#  i.e. fraction of unit circle !!!!!!!!!

    ## read xyz-file into list of 'xyz' objects:
    xyz_obj_list = xyz_io.read_xyz(in_xyz_path)

    #!# read axis (space-delimited):
    with open(axis_path) as f:
      axis = np.array(list(map(float, f.readline().strip().split())))


    ## convert angle into float:
    incremental_angle = float(incremental_angle)

    if atom_indices_to_rotate_path == None:
      atom_indices_to_rotate = list(range(xyz_obj_list[0].n_atoms))
    else:
      ## read file with atom indices (ints starting with 0) to rotate:
      atom_indices_to_rotate = np.genfromtxt(atom_indices_to_rotate_path, dtype=int)

      #!# if only want to rotate single atom:
      try:
        fakevar = atom_indices_to_rotate[1]
      except:
        atom_indices_to_rotate = np.array([atom_indices_to_rotate])


    ## loop over 'xyz' objects:
    for xy in xyz_obj_list:

      ## make copy of 'xyz' object:
      tmp_xyz = copy(xy)

      ## do we want to print 0th step?
      if bool(int(inclusive_bool)):
        tmp_xyz.write_xyz_file()

      for s in range(int(nsteps)):
        angle = (s+1) * incremental_angle

        ## now overwrite xyz_obj.crds with rotated coordinates:
        tmp_xyz.crds[atom_indices_to_rotate,:] = ctools.axis_angle(xy.crds[atom_indices_to_rotate,:], axis, angle)
        ## now write rotated xyz to stdout:
        tmp_xyz.write_xyz_file()

    return

#================================================#
#================================================#





##################################################
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
###  MAIN  #######################################

import sys

if __name__ == '__main__':

  ## read input arguments:
  if len(sys.argv) < 2:
    sys.exit('1st input arg specifies which main function to run, '
             'all others have inputs for that function.  '
             '\nfor list of main functions, run: \'main.py Help\'')


  ## initialize main_functions object:
  main = main_functions()


  #DEBUG_BOOL = True
  DEBUG_BOOL = False
  if DEBUG_BOOL:
    main.function_list[sys.argv[1]](*sys.argv[2:])
    exit()

  try:
    #!# call function from 1st input arg, and use others as input for function:
    main.function_list[sys.argv[1]](*sys.argv[2:])
  #!# if wrong number of arguments from sys:
  except TypeError:
    ## included correct(?) method name, but wrong arguments?
    main.function_list['Help'](sys.argv[1])
  #except KeyError:
  except:
    print('\nFirst argument is name of function to use. List of functions:')
    main.function_list['Help']()


#================================================#
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#================================================#


