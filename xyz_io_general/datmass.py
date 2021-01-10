
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
## GET MASSES  #######################

def get_masses():

# make these named tuples:

  H_dict = dict()
  H_dict['symbol']  = 'H'
  H_dict['mass']  =  1.007825

  O_dict = dict()
  O_dict['symbol']  = 'O'
  O_dict['mass']  = 15.994914

  D_dict = dict()
  D_dict['symbol']  = 'D'
  D_dict['mass']  =  2.014101

  C_dict = dict()
  C_dict['symbol']  = 'C'
  C_dict['mass']  = 12.000000

  N_dict = dict()
  N_dict['symbol']  = 'N'
  N_dict['mass']  = 14.003074

  Li_dict = dict()
  Li_dict['symbol'] = 'Li'
  Li_dict['mass'] =  6.015123

  Be_dict = dict()
  Be_dict['symbol'] = 'Be'
  Be_dict['mass'] =  9.012183

  F_dict = dict()
  F_dict['symbol']  = 'F'
  F_dict['mass']  = 18.998403

  Na_dict = dict()
  Na_dict['symbol'] = 'Na'
  Na_dict['mass'] = 22.989769

  Mg_dict = dict()
  Mg_dict['symbol'] = 'Mg'
  Mg_dict['mass'] = 23.985041

  Si_dict = dict()
  Si_dict['symbol'] = 'Si'
  Si_dict['mass'] = 27.976927

  Cl_dict = dict()
  Cl_dict['symbol'] = 'Cl'
  Cl_dict['mass'] = 34.968853

  K_dict = dict()
  K_dict['symbol']  = 'K'
  K_dict['mass']  = 38.963706

  Ca_dict = dict()
  Ca_dict['symbol'] = 'Ca'
  Ca_dict['mass'] = 39.962591

  Co_dict = dict()
  Co_dict['symbol'] = 'Co'
  Co_dict['mass'] = 58.933194

  Br_dict = dict()
  Br_dict['symbol'] = 'Br'
  Br_dict['mass'] = 78.918338

  Rb_dict = dict()
  Rb_dict['symbol'] = 'Rb'
  Rb_dict['mass'] = 84.911790

  Sr_dict = dict()
  Sr_dict['symbol'] = 'Sr'
  Sr_dict['mass'] = 87.905613

  I_dict = dict()
  I_dict['symbol']  = 'I'
  I_dict['mass']  = 126.904472

  Cs_dict = dict()
  Cs_dict['symbol'] = 'Cs'
  Cs_dict['mass'] = 132.905451

  Ba_dict = dict()
  Ba_dict['symbol'] = 'Ba'
  Ba_dict['mass'] = 137.905247

  At_dict = dict()
  At_dict['symbol'] = 'At'
  At_dict['mass'] = 209.987148

  Fr_dict = dict()
  Fr_dict['symbol'] = 'Fr'
  Fr_dict['mass'] = 223.019736

  Ra_dict = dict()
  Ra_dict['symbol'] = 'Ra'
  Ra_dict['mass'] = 224.020212


  mass_dict = dict()
  mass_dict['H']  = H_dict
  mass_dict['O']  = O_dict
  mass_dict['D']  = D_dict
  mass_dict['C']  = C_dict
  mass_dict['N']  = N_dict
  mass_dict['Li'] = Li_dict
  mass_dict['Be'] = Be_dict
  mass_dict['F']  = F_dict
  mass_dict['Na'] = Na_dict
  mass_dict['Mg'] = Mg_dict
  mass_dict['Si'] = Si_dict
  mass_dict['Cl'] = Cl_dict
  mass_dict['K']  = K_dict
  mass_dict['Ca'] = Ca_dict
  mass_dict['Co'] = Co_dict
  mass_dict['Br'] = Br_dict
  mass_dict['Rb'] = Rb_dict
  mass_dict['Sr'] = Sr_dict
  mass_dict['I']  = I_dict
  mass_dict['Cs'] = Cs_dict
  mass_dict['Ba'] = Ba_dict
  mass_dict['At'] = At_dict
  mass_dict['Fr'] = Fr_dict
  mass_dict['Ra'] = Ra_dict

  ## all lowercase keys:
  mass_dict['h']  = H_dict
  mass_dict['o']  = O_dict
  mass_dict['d']  = D_dict
  mass_dict['c']  = C_dict
  mass_dict['n']  = N_dict
  mass_dict['li'] = Li_dict
  mass_dict['be'] = Be_dict
  mass_dict['f']  = F_dict
  mass_dict['na'] = Na_dict
  mass_dict['mg'] = Mg_dict
  mass_dict['si'] = Si_dict
  mass_dict['cl'] = Cl_dict
  mass_dict['k']  = K_dict
  mass_dict['ca'] = Ca_dict
  mass_dict['co'] = Co_dict
  mass_dict['br'] = Br_dict
  mass_dict['rb'] = Rb_dict
  mass_dict['sr'] = Sr_dict
  mass_dict['i']  = I_dict
  mass_dict['cs'] = Cs_dict
  mass_dict['ba'] = Ba_dict
  mass_dict['at'] = At_dict
  mass_dict['fr'] = Fr_dict
  mass_dict['ra'] = Ra_dict

  ## all uppercase keys:
  mass_dict['LI'] = Li_dict
  mass_dict['BE'] = Be_dict
  mass_dict['NA'] = Na_dict
  mass_dict['MG'] = Mg_dict
  mass_dict['CL'] = Cl_dict
  mass_dict['SI'] = Si_dict
  mass_dict['CA'] = Ca_dict
  mass_dict['CO'] = Co_dict
  mass_dict['BR'] = Br_dict
  mass_dict['RB'] = Rb_dict
  mass_dict['SR'] = Sr_dict
  mass_dict['CS'] = Cs_dict
  mass_dict['BA'] = Ba_dict
  mass_dict['AT'] = At_dict
  mass_dict['FR'] = Fr_dict
  mass_dict['RA'] = Ra_dict


  return mass_dict




