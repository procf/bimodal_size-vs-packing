## Extract data to CSV and run Voronoi volume analysis on all colloids
## in a BD simulation
## NOTE: this code assumes 2 colloid types (typeid=0, typeid=1)
##
"""
## This code performs the following analyses:
##    - extracts colloid position data to from GSD to CSV files for all frames
##    - calculates (for all colloids and for each colloid sub-population):
##         * Voronoi volume around each colloid and the "free volume" 
##            (Voronoi volume - particle volume) in the simulation for each frame
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import gsd.hoomd
import math
import pyvoro
import pandas as pd
from scipy.stats import skew
from scipy.stats import kurtosis
import re
import os
import glob
import sys

##########################
""" INPUT PARAMETERS """
##########################

# manual simulation parameters
filepath = '../Gelation.gsd'
data_outpath = 'data'
colloid1_typeid = 1
colloid2_typeid = 2

# extra distance to make sure particles at the edge of the box are included
# range becomes (min+Lbuffer) to (max+Lbuffer) --> size+(2*Lbuffer) in each direction
Lbuffer = 0.001 

###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)

#######
## posCSV
"""
# create a CSV file with particle (i.e. node) position information for each frame
# (tag, x, y, z, typeID, radius)
"""

def posCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')

  # set path and filename  
  dir_path = data_outpath+'/frame-pos'
  pos_output = dir_path+'/positions_frame' # + <#>.csv in python loop

  # check for existing CSV data 
  if os.path.exists(dir_path) == False:
    os.mkdir(dir_path)
  if os.path.exists(dir_path) == True: 
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'positions_frame\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1

    if nframes == len(traj):  
      print('CSV files already seem to exist for all frames. Not creating new CSV files.')
      return 

  # if the CSV files didn't already exist, create them 
  nframes = len(traj)
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]

  for f in range(nframes):
    rpos = traj[f].particles.position[colloids]

    # make a data frame to export to CSV
    df_pos = pd.DataFrame()
    df_pos['tag'] = colloids
    df_pos['x'] = rpos[:,0]
    df_pos['y'] = rpos[:,1]
    df_pos['z'] = rpos[:,2] 
    df_pos['typeID'] = typeid
    df_pos['radius'] = radii

    df_pos.to_csv(pos_output+str(f)+'.csv', index=False)
  print("Position data saved to CSV for "+str(nframes)+" frames")
#######


#######
## Voronoi volumes
"""
# use pyvoro to calculate the volume surounding each colloid particle
# in each simulation frame:
# - for each particle pair, find the mid-point between them particle
# - use these midpoints to build polyhedra surrounding each particle
# - calculate the volume of each polyhedra and subtract particle volume
"""

def voronoi_volume_py(filename):
  traj = gsd.hoomd.open(filename, 'r')
  Lbox = traj[-1].configuration.box[:3]
  # set individual dimensions
  L_X = Lbox[0]
  L_Y = Lbox[1]
  L_Z = Lbox[2]

  # folder path
  dir_path = data_outpath+'/frame-pos'

  # get length of sim from the CSV data (match existing data not desired data!)
  # NOTE: expects data to be saved as "positions_frame#.csv" inside dir_path
  count = 0
  # Iterate directory
  for path in os.listdir(dir_path):
      # check if current path is a file
      if os.path.isfile(os.path.join(dir_path, path)):
          count += 1
  nframes = count

  # import all data into one dataframe
  all_frame_dfs = []
  for frame in range(nframes):
    #filepath = 'positions_frame99_plusframetag.csv'
    filepath = dir_path+'/positions_frame'+str(frame)+'.csv'
    # import CSV data
    df = pd.read_csv(filepath)
    # rename colums as needed
    df = df.rename(columns={"typeID":"typeid"})
    df.insert(loc=0, column='frame', value=frame)
    all_frame_dfs.append(df)
  allpos_df = pd.concat(all_frame_dfs, ignore_index=True)

  ###### GET SYSTEM PARAMS FROM LAST FRAME

  lastframe_df = allpos_df.loc[allpos_df['frame'] == nframes-1]

  # get the number of colloid types
  ntypes = len(np.unique(lastframe_df['typeid'].to_numpy()))

  # get the particle radii
  R_C = lastframe_df['radius']
  radii = np.unique(R_C)
  a1 = min(radii)
  a2 = max(radii)

  # number of particles
  ncolloids = len(np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0])
  ncolloid1 = len(lastframe_df.loc[lastframe_df['typeid'] == colloid1_typeid])
  ncolloid2 = len(lastframe_df.loc[lastframe_df['typeid'] == colloid2_typeid])
  if ncolloids != ncolloid1 + ncolloid2:
    print("ERROR: total number of colloids does not equal sum of subpopulations in voronoi volumes calc; Voronoi volumes NOT calculated")
    return

  # index numbers for all particles (extracted from first column of the csv file)
  colloid1 = lastframe_df.loc[lastframe_df['typeid'] == colloid1_typeid]["tag"]
  colloid2 = lastframe_df.loc[lastframe_df['typeid'] == colloid2_typeid]["tag"]
  colloids = lastframe_df["tag"] 

  # positions
  allpos1 = allpos_df.loc[allpos_df['typeid'] == colloid1_typeid][['frame','tag','typeid','x', 'y', 'z', 'radius']] # all colloid1 positions
  allpos2 = allpos_df.loc[allpos_df['typeid'] == colloid2_typeid][['frame','tag','typeid','x', 'y', 'z', 'radius']] # all colloid2 positions
  allpos = allpos_df[['frame','tag','typeid','x', 'y', 'z', 'radius']] # all colloid positions

  # create files to save the voronoi volumes
  f_all=open(data_outpath+'/voronoi-volumes.txt', 'w')
  f_all.write('frame particle typeid radius total-volume free-volume\n')

  f_c1=open(data_outpath+'/voronoi-volumes-c1only.txt', 'w')
  f_c1.write('frame particle typeid radius total-volume free-volume\n')

  f_c2=open(data_outpath+'/voronoi-volumes-c2only.txt', 'w')
  f_c2.write('frame particle typeid radius total-volume free-volume\n')

  #vorvols = np.zeros((nframes,ncolloids))

  for frame in range(0,nframes):

    frameradii = allpos.loc[allpos['frame'] == frame]['radius']
    frameradii_array = frameradii.to_numpy()
    frameradii_c1 = allpos1.loc[allpos1['frame'] == frame]['radius']
    frameradii_array_c1 = frameradii_c1.to_numpy()
    frameradii_c2 = allpos2.loc[allpos2['frame'] == frame]['radius']
    frameradii_array_c2 = frameradii_c2.to_numpy()

    frame_particle_volumes = (4/3)*np.pi*(frameradii_array**3)
    frame_particle_volumes_c1 = (4/3)*np.pi*(frameradii_array_c1**3)
    frame_particle_volumes_c2 = (4/3)*np.pi*(frameradii_array_c2**3)

    pos_c = allpos.loc[allpos['frame'] == frame][['x', 'y', 'z']]
    framepos_c = pos_c.to_numpy()
    pos_c1 = allpos1.loc[allpos1['frame'] == frame][['x', 'y', 'z']]
    framepos_c1 = pos_c1.to_numpy()
    pos_c2 = allpos2.loc[allpos2['frame'] == frame][['x', 'y', 'z']]
    framepos_c2 = pos_c2.to_numpy()  
 
    types_all = allpos.loc[allpos['frame'] == frame]['typeid'].to_numpy()
    types_1 = allpos1.loc[allpos['frame'] == frame]['typeid'].to_numpy()
    types_2 = allpos2.loc[allpos['frame'] == frame]['typeid'].to_numpy()

    cells_all = pyvoro.compute_voronoi(
      framepos_c, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array # matching list of particle radii -- required for polydisperse / radical tessellation 
    )

    cells_c1 = pyvoro.compute_voronoi(
      framepos_c1, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array_c1 # matching list of particle radii -- required for polydisperse / radical tessellation 
    )

    cells_c2 = pyvoro.compute_voronoi(
      framepos_c2, # point positions
      [[(-L_X/2)-Lbuffer, (L_X/2)+Lbuffer], [(-L_Y/2)-Lbuffer, (L_Y/2)+Lbuffer], [(-L_Z/2)-Lbuffer, (L_Z/2)+Lbuffer]], # limits 
      2.0, # block size
      radii=frameradii_array_c2 # matching list of particle radii -- required for polydisperse / radical tessellation 
    )
    
    volumes_all = np.zeros(len(cells_all)) 
    for j in range(len(volumes_all)):
      volumes_all[j] = cells_all[j].get("volume")

    free_volumes_all = np.zeros(len(volumes_all))
    for j in range(len(free_volumes_all)):
      free_volumes_all[j] = volumes_all[j] - frame_particle_volumes[j] 

    volumes_c1 = np.zeros(len(cells_c1))
    for j in range(len(volumes_c1)):
      volumes_c1[j] = cells_c1[j].get("volume")

    free_volumes_c1 = np.zeros(len(volumes_c1))
    for j in range(len(free_volumes_c1)):
      free_volumes_c1[j] = volumes_c1[j] - frame_particle_volumes_c1[j] 

    volumes_c2 = np.zeros(len(cells_c2))
    for j in range(len(volumes_c2)):
      volumes_c2[j] = cells_c2[j].get("volume")

    free_volumes_c2 = np.zeros(len(volumes_c2))
    for j in range(len(free_volumes_c2)):
      free_volumes_c2[j] = volumes_c2[j] - frame_particle_volumes_c2[j] 

 
    for particle in range(ncolloids):
      f_all.write('{0} {1} {2} {3} {4} {5}\n'.format(int(frame), particle, types_all[particle], frameradii_array[particle], volumes_all[particle], free_volumes_all[particle]))	

    for particle in range(ncolloid1):
      f_c1.write('{0} {1} {2} {3} {4} {5}\n'.format(int(frame), particle, types_1[particle], frameradii_array_c1[particle], volumes_c1[particle], free_volumes_c1[particle]))

    for particle in range(ncolloid2):
      f_c2.write('{0} {1} {2} {3} {4} {5}\n'.format(int(frame), particle, types_2[particle], frameradii_array_c2[particle], volumes_c2[particle], free_volumes_c2[particle]))

    #print(vorvols[1])

  print('Voronoi volume calculation complete for '+str(nframes)+' frames and '+str(ncolloids)+' colloid particles ('+str(ntypes)+'colloid types).')
#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':
  posCSV_calc(filepath)	
  voronoi_volume_py(filepath)

  print("\nAll Voronoi analyses complete.")
