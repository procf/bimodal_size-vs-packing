## Analyze the results of a DPD colloid simulation
## NOTE: requires matching Fortran module and solvopt module
## NOTE: this code assumes 1 solvent type (typeid=0) and 
##                           1 colloid types (typeid=1)
## NOTE: to select specific analyses, scroll to the bottom
##       of this file and comment out unwanted analyses in the
##       RUN CHECKS ON A SIMULATION section
"""
## This code performs the following analyses: 
##    - extracts temperature and pressure (the negative of stress) from all frames of a GSD file##    - calculates:
##        * colloid coordination number:
##           - coordination number distribution (Z counts) for each frame
##           - average coordination number (<Z>) for each frame
##        * mean squared displacement (MSD) of colloids and solvents for each frame
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
R_C = 1 #2
colloid_typeid = 1
kappa = 60


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
    KT=float(traj[i].log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0])
    # extract the transactios per secont (tps) for assessing speed/efficiency of the simulation
    tps=int(traj[i].log['Simulation/tps'][0])

    # write these values to the output file:
    # raw values
    f.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format(simframe, Pr, 
      Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, KT, tps))
		
    # rounded values
    #f.write('%f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %d %0.2f %d\n'%((i+1)*t1, Pr, 
    #  Vr[0], Vr[1], Vr[2], Vr[3], Vr[4], Pe, KT, tps))
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
  # use the last frame and typeid get the tags and total number colloids
  colloids = np.where(traj[-1].particles.typeid == [1])[0]
  ncolloids = len(colloids)
  # use last frame to get the radius of each colloid
  R_C = 0.5*traj[-1].particles.diameter[colloids]
  # use last frame to get the typeid of each colloid
  typeid = traj[-1].particles.typeid[colloids]

  ## 1. CALCULATE Z DISTRIBUTION

  # gather data for the whole simulation
  allpos_allframe = np.zeros((nframes,ncolloids,3))
  m_xys = np.zeros(nframes)
  for frame in range(nframes):
    # get all particle positions
    allpos_allframe[frame] = traj[frame].particles.position[colloids]
    # get the xy tilt factor (square=0.0, sheared-right=0.45)
    m_xys[frame] = traj[frame].configuration.box[3]

  print('Calculating Z-distribution data for '+str(nframes)+' times...')

  # run the fortran module
  Zs_array = module.coordination_number(nframes,Lbox,ncolloids,R_C,m_xys,allpos_allframe,cut_off)

  # convert array to data frame for easy saving
  allframes_Zs_df = []
  for frame in range(nframes):
    frame_df = pd.DataFrame(Zs_array[frame])
    frame_df.reset_index(inplace=True)
    frame_df = frame_df.rename(columns = {'index':'colloidID'})
    frame_df.insert(loc=1, column='typeid', value=typeid)
    frame_df.insert(loc=0, column='frame', value=frame)
    frame_df = frame_df.rename(columns = {0:'Z'})
    allframes_Zs_df.append(frame_df)

  total_Zsdf = pd.concat(allframes_Zs_df, ignore_index=True)

  total_Zsdf.to_csv(data_outpath+'/Z-counts.csv', index=False)

  print("...Coordination number calculated for all colloids in "+str(nframes)+" frames")


  ## 2. CALCULATE THE AVERAGE Z

  print("Calculating the average coordination number for "+str(nframes)+" frames...")

  # create an array to save the averages
  Zavgs_array = np.zeros(nframes)

  #caluclate the average coordination numbers
  for frame in range(nframes):
    all_sum = sum(Zs_array[frame,:])
    Zavgs_array[frame] = all_sum / ncolloids

  # convert array to data frame for easy saving
  Zavgs_df = pd.DataFrame(Zavgs_array)
  Zavgs_df.insert(loc=0,column='simframe', value=list(range(0,nframes)))
  for frame in range(nframes):
    Zavgs_df = Zavgs_df.rename(columns = {0:'Z_any'})

  Zavgs_df.to_csv(data_outpath+'/Zavg.csv', index=False)

  print("...Average coordination number calculated for "+str(nframes)+" frames")
#######

#######
## mean squared displacement (MSD)
"""
# calculate the mean squared displacement (MSD) for colloids and solvents in all frames,
# as well as the sample standard deviation of MSD. The MSD of solvents is only used to
# verify that the simulation is running with correct physics
#
# MSD measures how far a particle has moved from it's original position over a period of time 
# (if the data forms a diagonal line with a fixed slope, this indicates the particle is 
#   moving steadily (fluid); but if the data forms a flat/horizontal line, the particle 
#   has stopped moving (solid) -- MSD tells you about particle dynamics) 
"""

# calculate the time it takes a particle to diffuse half the shortest box length
d = 3 # dimension of the system (2D = 2, 3D = 3)
D = kT/(6*math.pi*eta0*R_C) # diffusion coefficient = r^2 / 2d*tau
tau_to_half = (Lbox_shortest/2)**2 / (2*d*D) # diffusion time to half-box (L/2)

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
    # get all the colloid particles in the last frame
    colloids = np.where(traj[-1].particles.typeid == [1])[0]
    # use this to count the total number of colloids
    ncolloids = len(colloids)

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


    # run the module to calculate MSD
    # (output file "msd.txt" is created in the Fortran module)      

    # only on colloids
    if ("_Colloids" in filepath) or ("BD" in filepath):
      module.msd_calculation_colloids(data_outpath,nframes,ncolloids,allcollpos)

    # colloids and solvent
    else:
      # get all the solvent particles in the last frame
      solvents = np.where(traj[-1].particles.typeid == [0])[0]
      # use this to count the total number of solvents
      nsolvents = len(solvents)

      # create an empty array for xyz positon of all solvent particles in all frames  
      allsolvpos = np.zeros((nframes,nsolvents,3))
      # get the initial positions from the first frame
      allsolvpos[0,:,:] = traj[0].particles.position[solvents]
      # correct the change in position for particles crossing a box boundary 
      for i in range(1,nframes):
        # calculate the change in position since the previous frame
        movement = traj[i].particles.position[solvents] - traj[i-1].particles.position[solvents]
        # if it is more than one box length (i.e. the particle crossed a 
        # boundary) correct the position to wrap inside the sim box
        # (matching our periodic boundaries)
        movement -= Lbox*np.rint(movement*inv_Lbox)
        # update and record the position of the particles in this frame
        allsolvpos[i,:,:] = allsolvpos[i-1,:,:] + movement

      module.msd_calculation_DPD(data_outpath,nframes,ncolloids,nsolvents,allcollpos,allsolvpos)

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
framechoice = [-1]

def void_size_calc_py(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  # get the number of frames
  nframes = len(traj)

  # convert negative framechoice values into specific frames 
  for i in range(len(framechoice)):
    if framechoice[i] < 0:
      framechoice[i] = nframes - abs(framechoice[i])
  # replace the total nframes with the desired nframes to analyze
  nframes = len(framechoice)

  # get the index of all type1 colloids
  colloids = np.where(traj[-1].particles.typeid == [1])[0]
  # calculate the number of type1 colloids
  ncolloids = len(colloids)
  # find the radii of every type1 colloid
  radii = 0.5*traj[-1].particles.diameter[colloids]

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
    frame = framechoice[i]
    box_length[i,:]=traj[frame].configuration.box[:3]
    rxi[i,:] = traj[frame].particles.position[colloids,0]
    ryi[i,:] = traj[frame].particles.position[colloids,1]
    rzi[i,:] = traj[frame].particles.position[colloids,2]

  # calculate the void_size for all the selected data
  void_size_calculation.void_size_calc(data_outpath,ncolloids,nframes,framechoice,nprobe,radii,dcell_init,rxi,ryi,rzi,box_length)

  print("Void size distribution calculated for "+str(nprobe)+" probe points in each of "+str(nframes)+" frames")
#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':	
  extract_properties_py(filepath)	
  coordination_number_py(filepath)
  msd_py(filepath)
  void_size_calc_py(filepath)
