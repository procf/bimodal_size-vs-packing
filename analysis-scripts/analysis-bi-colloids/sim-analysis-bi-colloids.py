## Analyze the results of a BD colloid simulation
## NOTE: requires matching Fortran module and solvopt module
## NOTE: this code assumes 2 colloid types (typeid=0, typeid=1)
## NOTE: to select specific analyses, scroll to the bottom
##       of this file and comment out unwanted analyses in the
##       RUN CHECKS ON A SIMULATION section
"""
## This code performs the following analyses: 
##    - extracts temperature and pressure (the negative of stress) from all frames of a GSD file 
##   - calculates:
##        * colloid coordination number:
##           - coordination number distribution (Z counts) for each frame
##           - average coordination number (<Z>) for each frame
##        * mean squared displacement (MSD) of colloids for each frame
##        * void size distribution for the final frame OR a selection of frames, 
##          using two methods: 
##           - Torquato’s Pore Size Distribution
##           - Gubbins’s Pore Size Distribution 
##           (requires the solvopt algorithm, included as a separate f90 module)
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import gsd.hoomd
import math
import fortranmod as module
from fortranmod import void_size_calculation
import os
import re
import glob
import sys

##########################
""" INPUT PARAMETERS """
##########################

# manual simulation parameters
filepath = '../Gelation.gsd' 
data_outpath = 'data'
period = 10000
kT = 0.1
eta0 = 0.3
L_X = 70
Lbox_shortest = L_X # the shortest side of the simulation box
R_C1 = 1
R_C2 = 2
colloid1_typeid = 1
colloid2_typeid = 2
kappa = 60

# set min and max particle size
R_min = min(R_C1, R_C2)
R_max = max(R_C1, R_C2)

###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)

#######
## extract thermodynamic properties from GSD for all frames
"""
# get temperature and pressure (AKA negative of stress), including all virial components
"""
def extract_properties_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)
  log = traj[-1].log

  # Check if any log exists
  if not log:
      print('WARNING: No log found, skipping extract_properties_py')
      return

  # Check required keys
  required_keys = [
      'md/compute/ThermodynamicQuantities/pressure_tensor',
      'md/compute/ThermodynamicQuantities/virial_ind_tensor',
      'md/compute/ThermodynamicQuantities/potential_energy',
      'md/compute/ThermodynamicQuantities/kinetic_temperature',
      'Simulation/tps'
  ]

  missing = [k for k in required_keys if k not in log]

  if missing:
      print('WARNING: Missing log keys:', missing)
      print('Skipping extract_properties_py')
      return

  else:
    # create a file to save the thermodynamic properties
    f=open(data_outpath+'/gsd-properties.txt', 'w')
    f.write('simframe Virial-Pressure Vr_CONS Vr_DISS Vr_RAND Vr_SQUE Vr_CONT PE kT tps\n')

    # for each frame
    for i in range(0, nframes):
      simframe = i
      # extract the total "force virial contribution" of the pressure tensor
      Pr=float(traj[i].log['md/compute/ThermodynamicQuantities/pressure_tensor'][1])
      # extract the decomposed virial compoenents of the pressure tensor
      Vr=np.asarray(traj[i].log['md/compute/ThermodynamicQuantities/virial_ind_tensor'],dtype=float)
      # extract the potential energy and scale to kilo-units (1/1000)
      Pe=float(traj[i].log['md/compute/ThermodynamicQuantities/potential_energy'][0])*0.001
      # extract the kinetic temperature (kT)
      kT=float(traj[i].log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0])
      # extract the transactios per secont (tps) for assessing speed/efficiency of the simulation
      tps=int(traj[i].log['Simulation/tps'][0])

      # write these values to the output file:
      # raw values
      f.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format(simframe, Pr, 
        Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, kT, tps))
		
      # rounded values
      #f.write('%f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %d %0.2f %d\n'%((i+1)*t1, Pr, 
      #  Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, kT, tps))
    print("Data extracted from GSD file")
#######


#######
## coordination number AKA contact number 
"""
# calculate the coordination number distribution (Z counts) and 
# the average coordination number (<Z>) for each frame
#
# for a given particle, Z is the number of other particles touching it
# (we define "contact" from the "attraction range" set with kappa; 
#   for a gel the average Z should plateau as the network is formed, 
#    and the distribution is usually centered around an average of Z=6)
"""

# use kappa to calculate attraction range (cut-off distance)
cut_off = round(3/kappa,2)

def coordination_number_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filepath,'r')
  # get the number of frames
  nframes = len(traj)
  # use the last frame to get simulation box size [L_X, L_Y, L_Z]
  Lbox = traj[-1].configuration.box[:3]
  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  # use last frame to get the radius of each colloid
  R_C = 0.5*traj[-1].particles.diameter[colloids]
  # use last frame to get the typeid of each colloid
  typeid = traj[-1].particles.typeid[colloids]

  # generate list of all types of interaction pairs from the data
  ntypes = int(len(np.unique(typeid)))
  pairs = []
  pair_names = []
  for colloid_type_a in range(ntypes):
    for colloid_type_b in range(ntypes):
      pairs.append([colloid_type_a+1, colloid_type_b+1])
  for pair in pairs:
    pair_names.append('Z_'+str(pair[0])+str(pair[1]))
  npairs = len(pairs)
  # and the number of total combos, Z_any, Z_11, Z_12..., Z_1any, ...
  n_allcombos = int(1+npairs+np.sqrt(npairs))

  ## 1. CALCULATE Z DISTRIBUTION

  # gather data for the whole simulation
  allpos_allframe = np.zeros((nframes,ncolloids,3))
  m_xys = np.zeros(nframes)
  for frame in range(nframes):
    # get all particle positions
    allpos_allframe[frame] = traj[frame].particles.position[colloids]
    # get the xy tilt factor (square=0.0, sheared-right=0.45)
    m_xys[frame] = traj[frame].configuration.box[3]

  print('Calculating Z-distribution data for '+str(len(pairs)+1)+' types of colloid-colloid pairs in '+str(nframes)+' times...')

  # "pair type" labels start at 1 (i.e. Z_11, Z_12, etc), but GSD typeids may not; account for this in Fortran
  if 0 in np.unique(typeid):
    fortran_pairs = []
    for pair in pairs:
      new_pair = []
      for i in pair:
        BD_i = i-1
        new_pair.append(BD_i)
      fortran_pairs.append(new_pair)
  else:
    fortran_pairs = pairs.copy()

  # run the fortran module
  Zs_array = module.coordination_number(nframes,Lbox,ncolloids,typeid,R_C,m_xys,allpos_allframe,cut_off,npairs,n_allcombos,fortran_pairs)

  # vectorize
  nframes, ncolloids, _ = Zs_array.shape

  z_labels = ['Z_any'] + pair_names + [f'Z_{t+1}any' for t in range(ntypes)]

  flat = Zs_array.reshape(nframes * ncolloids, -1)

  total_Zsdf = pd.DataFrame(flat, columns=z_labels)

  total_Zsdf.insert(0, 'colloidID', np.tile(np.arange(ncolloids), nframes))
  total_Zsdf.insert(0, 'frame', np.repeat(np.arange(nframes), ncolloids))
  total_Zsdf.insert(2, 'typeid', np.tile(typeid, nframes))

  total_Zsdf.to_csv(data_outpath + '/Z-counts.csv', index=False)

  print("...Coordination number calculated for all colloids in "+str(nframes)+" frames")

  ## 2. CALCULATE THE AVERAGE Z

  print("Calculating the average coordination number for "+str(len(pairs)+1)+" types of colloid-colloid pairs in "+str(nframes)+" frames...")

  # get the tags and number of colloids for each colloid type
  colloids_subpops = []
  ncolloids_subpops = []
  for colloid_type in np.unique(typeid):
    tags = np.where(traj[-1].particles.typeid == [colloid_type])[0]
    colloids_subpops.append(tags)
    number = len(tags)
    ncolloids_subpops.append(number)

  if sum(ncolloids_subpops) != ncolloids:
    print("ERROR: total number of colloids in <Z> calculation does not equal sum of subpopulations; <Z> was NOT calculated")
    return

  # create an array to save the averages
  Zavgs_array = np.zeros((nframes,n_allcombos))
  #Zavgs_array = np.zeros((nframes,len(pairs)+1))

  #caluclate the average coordination numbers
  for frame in range(nframes):
    any_sum = sum(Zs_array[frame,:,0])
    any_nc = ncolloids
    Zavgs_array[frame][0] = any_sum / any_nc

    for pair_index in range(len(pairs)):
      sub_sum = sum(Zs_array[frame,:,pair_index+1])
      primary_colloid_typeid = pairs[pair_index][0]
      sub_nc = ncolloids_subpops[primary_colloid_typeid-1]
      Zavgs_array[frame][pair_index+1] = sub_sum/sub_nc

    # and the Z_1any, Z_2any, etc...
    for type_index in range(ntypes):
      subany_sum = sum(Zs_array[frame,:,1+npairs+type_index])
      subany_nc = ncolloids_subpops[type_index]
      Zavgs_array[frame][1+npairs+type_index] = subany_sum/subany_nc 

  # convert array to data frame for easy saving
  Zavgs_df = pd.DataFrame(Zavgs_array)

  z_labels = (
      ['Z_any'] +
      pair_names +
      [f'Z_{t+1}any' for t in range(ntypes)]
  )
  Zavgs_df = pd.DataFrame(Zavgs_array, columns=z_labels)
  Zavgs_df.insert(0, 'simframe', np.arange(nframes))


  Zavgs_df.to_csv(data_outpath+'/Zavg.csv', index=False)

  print("...Average coordination number calculated for "+str(nframes)+" frames")
#######


#######
## mean squared displacement (MSD)
"""
# calculate the mean squared displacement (MSD) for colloids in all frames,
# as well as the sample standard deviation of MSD. 
#
# MSD measures how far a particle has moved from it's original position over a period of time 
# (if the data forms a diagonal line with a fixed slope, this indicates the particle is 
#   moving steadily (fluid); but if the data forms a flat/horizontal line, the particle 
#   has stopped moving (solid) -- MSD tells you about particle dynamics) 
"""

# calculate the time it takes the smallest (fastest) particle to diffuse half the shortest box length
D = kT/(6*math.pi*eta0*R_min) # diffusion coefficient = r^2 / tau
tau_to_half = (Lbox_shortest/2)**2 / D # diffusion time to half-box (L/2)

def msd_py(filename):

  if period >= tau_to_half:
    print('ERROR: the GSD file\'s recording timestep is too large for accurate MSD calculations.\n\n' +
          'To create data that you can use to accurately calculate the MSD, you should set the'
          ' period/trigger to a small enough value that a particle will not passively move'
          ' a distance equal to 1/2 the shortest box length in between frames.\n\n' + 
          'Rerun the simulation with a smaller trigger/period before calculating MSD.')
  else:
    # open the simulation GSD file as "traj" (trajectory)
    traj = gsd.hoomd.open(filename, 'r')
    # get the simulation box size [L_X, L_Y, L_Z] from the last frame 
    Lbox = traj[-1].configuration.box[:3]
    inv_Lbox = 1.0/Lbox
    # get the number of frames
    nframes = len(traj)	
    # get all the colloid particles in the last frame, and use this to count the total number of colloids
    colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
    ncolloids = len(colloids)
    # do the same for all the colloid subpopulations
    colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
    ncolloid1 = len(colloid1)
    colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
    ncolloid2 = len(colloid2)
    if (ncolloid1 + ncolloid2) != ncolloids:
      print("ERROR: ncolloid1 + ncolloid2 != ncolloids in MSD calc; MSD was NOT calculated")
      return

    # create an empty array for xyz positon of all colloids in all frames    
    allcollpos = np.zeros((nframes,ncolloids,3))
    # get the initial colloid positions from the first frame
    allcollpos[0,:,:] = traj[0].particles.position[colloids]
    # correct the change in position for colloids crossing a box boundary 
    for i in range(1,nframes):
      # calculate the change in position since the previous frame
      delpos = traj[i].particles.position[colloids] - traj[i-1].particles.position[colloids]
      # if it is more than one box length (i.e. the particle crossed a 
      # boundary) correct the position to wrap inside the sim box
      # (matching our periodic boundaries)
      delpos -= Lbox*np.rint(delpos*inv_Lbox)
      # update and record the position of the colloids in this frame
      allcollpos[i,:,:] = allcollpos[i-1,:,:] + delpos

    # create an empty array for xyz positon of all type 1 colloids in all frames    
    allcollpos1 = np.zeros((nframes,ncolloid1,3))
    # get the initial type 1 colloid positions from the first frame
    allcollpos1[0,:,:] = traj[0].particles.position[colloid1]
    # correct the change in position for type 1 colloids crossing a box boundary 
    for i in range(1,nframes):
      # calculate the change in position since the previous frame
      delpos1 = traj[i].particles.position[colloid1] - traj[i-1].particles.position[colloid1]
      # if it is more than one box length (i.e. the particle crossed a 
      # boundary) correct value the position to be inside the sim box
      delpos1 -= Lbox*np.rint(delpos1*inv_Lbox)
      # update and record the position of the colloids in this frame
      allcollpos1[i,:,:] = allcollpos1[i-1,:,:] + delpos1

    # create an empty array for xyz positon of all type 2 colloids in all frames    
    allcollpos2 = np.zeros((nframes,ncolloid2,3))
    # get the initial type 2 colloid positions from the first frame
    allcollpos2[0,:,:] = traj[0].particles.position[colloid2]
    # correct the change in position for colloids crossing a box boundary 
    for i in range(1,nframes):
      # calculate the change in position since the previous frame
      delpos2 = traj[i].particles.position[colloid2] - traj[i-1].particles.position[colloid2]
      # if it is more than one box length (i.e. the particle crossed a 
      # boundary) correct value the position to be inside the sim box
      delpos2 -= Lbox*np.rint(delpos2*inv_Lbox)
      # update and record the position of the colloids in this frame
      allcollpos2[i,:,:] = allcollpos2[i-1,:,:] + delpos2

    # run the module to calculate MSD
    # (output file "msd.txt" is created in the Fortran module)      
    module.msd_calculation(data_outpath,nframes,ncolloids,ncolloid1,ncolloid2,allcollpos,allcollpos1,allcollpos2)

    print("MSD calculation complete")
#######


#######
## void size calculation
"""
# calculates a distribution of the void space (approximated as spheres) in between particle clusters
# 
# There are two methods that are used:
#   - Torquato’s Pore Size Distribution (volume where the center of a particle can fit in between clusters)
#   - Gubbins’s Pore Size Distribution (volume occupied by a whole particlein in between clusters) 
#       (Gubbin's PSD requires solvopt, the Solver For Local Nonlinear Optimization Problems, 
#        incuded as a secondary fortran module) 
# 
# These methods are described in the Section 4.2 and the Appendix of
# Sorichetti, Hugouvieux, and Kob 2020, DOI: 10.1021/acs.macromol.9b02166
#
# This code assumes you are analyzing a porous medium made of uniform particles (like a colloidal gel),
#  and uses particle trajectories in gsd format; It takes particle and probe size as inputs.
#  Then it uses Linked-list method to compute minimum distance of a point in the void space from 
#  nearby porous medium particles. And then uses solvopt non-linear optimization code to compute Gubbin's void size.
#
# (we usually use Gubbin's PSD and expect a plot of probability vs void diameter to peak at the size of the most common voids)
"""

# calculate voidsize for the last frame only (or select a different frame)
framechoice_ps = [-1]

def void_size_calc_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)

  # convert negative framechoice values into specific frames 
  for i in range(len(framechoice_ps)):
    if framechoice_ps[i] < 0:
      framechoice_ps[i] = nframes - abs(framechoice_ps[i])
  # replace the total nframes with the desired nframes to analyze
  nframes = len(framechoice_ps)

  # get all the colloid particles in the last frame, and use this to count the total number of colloids
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  # find the radii of every colloid
  radii = 0.5*traj[-1].particles.diameter[colloids]

  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in gofr calc; gofr was NOT calculated")
    return
  radii_c1 = 0.5*traj[-1].particles.diameter[colloid1]
  radii_c2 = 0.5*traj[-1].particles.diameter[colloid2]


  ## 1. CALCULATE THE VOID SIZE DISTRIBUTION

  print('Calculating void size distribution data for '+str(len(framechoice_ps))+' frame(s)...')

  # set the number of random points used to explore void size
  nprobe = 10000 # can test quickly at 1
  # set the size of the cells used in the linked list
  dcell_init = 5.0

  # create empty arrays for holding data for all frames
  box_length = np.zeros((nframes,3))
  rxi = np.zeros((nframes,ncolloids))
  ryi = np.zeros((nframes,ncolloids))
  rzi = np.zeros((nframes,ncolloids))

  # fill arrays with data from each frame
  for i in range(nframes):
    frame = framechoice_ps[i]
    box_length[i,:]=traj[frame].configuration.box[:3]
    rxi[i,:] = traj[frame].particles.position[colloids,0]
    ryi[i,:] = traj[frame].particles.position[colloids,1]
    rzi[i,:] = traj[frame].particles.position[colloids,2]


  # do the same for all the colloid subpopulations
  rxi_c1 = np.zeros((nframes,ncolloid1))
  ryi_c1 = np.zeros((nframes,ncolloid1))
  rzi_c1 = np.zeros((nframes,ncolloid1)) 
  for i in range(nframes):
    frame = framechoice_ps[i]
    rxi_c1[i,:] = traj[frame].particles.position[colloid1,0]
    ryi_c1[i,:] = traj[frame].particles.position[colloid1,1]
    rzi_c1[i,:] = traj[frame].particles.position[colloid1,2]
  rxi_c2 = np.zeros((nframes,ncolloid2))
  ryi_c2 = np.zeros((nframes,ncolloid2))
  rzi_c2 = np.zeros((nframes,ncolloid2)) 
  for i in range(nframes):
    frame = framechoice_ps[i]
    rxi_c2[i,:] = traj[frame].particles.position[colloid2,0]
    ryi_c2[i,:] = traj[frame].particles.position[colloid2,1]
    rzi_c2[i,:] = traj[frame].particles.position[colloid2,2]


  # calculate the void_size for all the selected data
  print("...calculating void size for all colloids")
  pop_type = '-allc'
  void_size_calculation.void_size_calc(data_outpath,pop_type,ncolloids,nframes,framechoice_ps,nprobe,radii,dcell_init,rxi,ryi,rzi,box_length)

  void_size_calculation.deallocate_arrays()
  
  # do the same for all the colloid subpopulations
  print("...calculating void size for colloid type 1 only")
  pop_type = '-c1only'
  void_size_calculation.void_size_calc(data_outpath,pop_type,ncolloid1,nframes,framechoice_ps,nprobe,radii_c1,dcell_init,rxi_c1,ryi_c1,rzi_c1,box_length)

  void_size_calculation.deallocate_arrays()

  print("...calculating void size for colloid type 2 only")
  pop_type = '-c2only'
  void_size_calculation.void_size_calc(data_outpath,pop_type,ncolloid2,nframes,framechoice_ps,nprobe,radii_c2,dcell_init,rxi_c2,ryi_c2,rzi_c2,box_length)

  void_size_calculation.deallocate_arrays()

  print("...Void size distribution calculated for "+str(nprobe)+" probe points in each of "+str(nframes)+" frames")


  ## 2. CALCULATE THE AVERAGE VOID SIZE

  print("Calculating the average void size for "+str(len(framechoice_ps))+" frame(s)...")

  allc_filename = data_outpath+'/voidsize-allc.csv'
  c1_filename = data_outpath+'/voidsize-c1only.csv'
  c2_filename = data_outpath+'/voidsize-c2only.csv'

  filenames = [allc_filename, c1_filename, c2_filename]
  filetags = ['allc','c1only','c2only']

  frame_avg_list = []
  for f in range(len(filenames)):
    filepath = filenames[f]
    data_allframes_df = pd.read_csv(filepath)
    filetag = filetags[f]

    nvoids_list = []
    volumes_list = []
    T_avg_list = []
    G_avg_list = []
    frames_list = []
    for i in range(nframes):
      frame = framechoice_ps[i]

      data_df = data_allframes_df.loc[data_allframes_df['frame'] == frame]

      # calculate system volume and add to the dataframe
      L_X = int(np.ceil(max(data_df['probe_posx']))-np.floor(min(data_df['probe_posx'])))
      L_Y = int(np.ceil(max(data_df['probe_posy']))-np.floor(min(data_df['probe_posy'])))
      L_Z = int(np.ceil(max(data_df['probe_posz']))-np.floor(min(data_df['probe_posz'])))
      volumes_list.append(L_X*L_Y*L_Z)
      #print(f" - System volume: {[L_X, L_Y, L_Z]}")

      for voidsize_type in ['void_diameter_T','void_diameter_G']:
        #print(f' - {voidsize_type}')

        all_void_diameters = data_df[voidsize_type].to_numpy()

        #################################
        """ AVG PORE SIZE CALCULATION """
        #################################

        # get average
        void_diameters = all_void_diameters[np.isfinite(all_void_diameters)] # clean any NaNs and infs
        void_avg = np.mean(void_diameters)
        if voidsize_type == 'void_diameter_T':
          T_avg_list.append(void_avg)
        elif voidsize_type == 'void_diameter_G':
          G_avg_list.append(void_avg)
        else:
          print(f"ERROR: sizetype must be 'Torquato' or 'Gubbin', not: {sizetype}")
          exit(1)

        # calc nvoids for future normalization
        nvoids_dataset = len(void_diameters)

      # only add nvoids once, it's the same for Torquato and Gubbin
      nvoids_list.append(nvoids_dataset)

      #print("Void size average calculated for frame "+str(frame)+'...')
      frames_list.append(frame)

    avgs_df = pd.DataFrame(T_avg_list, columns=['T_avg'])
    avgs_df['G_avg'] = G_avg_list
    avgs_df['nvoids'] = nvoids_list
    avgs_df['volume'] = volumes_list
    avgs_df['population'] = filetag
    avgs_df.insert(loc=0,column='frame', value=frames_list)

    frame_avg_list.append(avgs_df)

  all_avgs_df = pd.concat(frame_avg_list, ignore_index=True)
  all_avgs_df.to_csv(data_outpath+'/voidsize-avg.csv', index=False)

  print("...Average void size calculated for "+str(len(filenames))+" populations in "+str(len(framechoice_ps))+" frame(s)")

#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':	
  # record params:
  print("Data Parameters")
  print(" - period:",period)
  print(" - Lbox_shortest:",Lbox_shortest)
  print(" - kT:",kT)
  print(" - R_C1:",R_C1)
  print(" - R_C2:",R_C2)
  print(" - eta0:",eta0)
  print(" - kappa:",kappa)
  print(" - colloid1_typeid:",colloid1_typeid)
  print(" - colloid2_typeid:",colloid2_typeid,'\n')

  # run analyses
  extract_properties_py(filepath)	
  coordination_number_py(filepath)
  msd_py(filepath)
  void_size_calc_py(filepath)

  print("\nAll analyses complete")
