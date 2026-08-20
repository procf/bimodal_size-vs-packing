## Extract data to CSV and run primary network analysis on the results 
## of a BD colloid simulation
## NOTE: requires matching Fortran module
## NOTE: this code assumes 2 colloid types (typeid=0, typeid=1)
##
"""
## This code performs the following analyses:
## - Extracts colloid position data to from GSD to CSV files for all frames
## - Extracts the list of all bonded colloid pairs (i.e network edges) from
##     GSD to CSV files for all frames
## - Calculates primary network analysis metrics for all frames
## 	- number of connected components 
##	- average degree (AKA average coordination number in gels)
##      - largest connected component (LCC)
##      - average clustering coefficient
##      - average square clustering coefficient 
"""
## (Rob Campbell)


########################
""" MODULE LIBRARY """
########################
import numpy as np
import gsd.hoomd
import math
import fortranmod as module 
import pandas as pd
import networkx as nx
import os
import re
import glob
import sys
from statistics import mean
from scipy.spatial.distance import pdist


##########################
""" INPUT PARAMETERS """
##########################

# manual simulation parameters
filepath = '../Gelation.gsd'
data_outpath = 'data'
colloid1_typeid = 1
colloid2_typeid = 2
kappa = 60

cut_off = round(3/kappa,2)

###########################
""" DEFINE SIM CHECKS """
###########################

# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


# define a function to recenter the dims on the center of mass, accounting for periodic boundaries in the cluster length
def pbc_centered_positions(positions, box):
    """
    Center the positions around their center of mass (COM),
    accounting for periodic boundaries using the minimum image convention.

    Returns:
    - centered positions (N, 3), shifted such that COM is at the origin (0,0,0)
    """
    # Convert positions to [0, L) for safety
    positions = positions % box

    # Compute COM using angles (Fourier trick for PBC-aware averaging)
    theta = 2 * np.pi * positions / box
    com_theta = np.arctan2(np.mean(np.sin(theta), axis=0),
                           np.mean(np.cos(theta), axis=0))
    com = (com_theta % (2 * np.pi)) * box / (2 * np.pi)

    # Shift all positions relative to COM using minimum image convention
    disp = positions - com
    disp -= box * np.round(disp / box)  # Apply minimum image convention
    return disp

#######
## posCSV
"""
# create a CSV of particle (i.e. node) position information for each frame
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
      print('position data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
  typeid = traj[-1].particles.typeid[colloids]
  radii = 0.5*traj[-1].particles.diameter[colloids]
  lbox = traj[-1].configuration.box[:3]

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
## edgelistCSV
"""
# create a CSV file of all the bonded particle pairs (i.e. edges of the network) for each frame
# i,j position for each bond/edge
"""

def edgelistCSV_calc(filename):
  # open the simulation GSD file as "traj" (trajectory)
  traj = gsd.hoomd.open(filename, 'r')
  
  # set path and filename
  edge_dir_path = data_outpath+'/frame-edges'
  edge_output = edge_dir_path+'/edgelist' # + <#>.csv in python loop
  edge1only_dir_path = data_outpath+'/frame-edges-c1only'
  edge1only_output = edge1only_dir_path+'/edgelist' # + <#>.csv in python loop
  edge2only_dir_path = data_outpath+'/frame-edges-c2only'
  edge2only_output = edge2only_dir_path+'/edgelist' # + <#>.csv in python loop
  
  # check for existing CSV data 
  if os.path.exists(edge_dir_path) == False:
    os.mkdir(edge_dir_path)
  if os.path.exists(edge1only_dir_path) == False:
    os.mkdir(edge1only_dir_path)
  if os.path.exists(edge2only_dir_path) == False:
    os.mkdir(edge2only_dir_path)

  if (os.path.exists(edge_dir_path) == True) and (os.path.exists(edge1only_dir_path) == True) and (os.path.exists(edge1only_dir_path) == True):
    # NOTE: this counts ALL CSV files in this directory
    nframes = 0
    # set the pattern for files ending in <number>.csv
    pattern = re.compile(r'edgelist\d+\.csv$')
    # Iterate directory
    for filename in os.listdir(edge_dir_path):
      # Check if the file has a CSV extension
      if pattern.match(filename):
        nframes += 1
    
    if nframes == len(traj):
      print('edgelist data CSV files already seem to exist for all frames. Not creating new CSV files.')
      return

  nframes = len(traj)
  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)
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

  rcut = cut_off 
  lbox = traj[-1].configuration.box[:3]
  inv_lbox = 1/lbox

  # create an array of xyz positon of all colloids in all frames    
  print("...calculating edgelist for all colloids")
  allpos = np.zeros((nframes,ncolloids,3))
  for i in range(0,nframes):
    allpos[i,:,:] = traj[i].particles.position[colloids] 
  module.edgelist_calc(nframes,ncolloids,radii,allpos,lbox,rcut,edge_output)

  # do the same for all the colloid subpopulations
  print("...calculating edgelist for colloid type 1 only")
  allpos1 = np.zeros((nframes,ncolloid1,3))
  for i in range(0,nframes):
    allpos1[i,:,:] = traj[i].particles.position[colloid1] 
  module.edgelist_calc(nframes,ncolloid1,radii_c1,allpos1,lbox,rcut,edge1only_output)
  print("...calculating edgelist for colloid type 2 only")
  allpos2 = np.zeros((nframes,ncolloid2,3))
  for i in range(0,nframes):
    allpos2[i,:,:] = traj[i].particles.position[colloid2] 
  module.edgelist_calc(nframes,ncolloid2,radii_c2,allpos2,lbox,rcut,edge2only_output)

  print("Edgelist calculation complete for "+str(nframes)+" frames")
#######


#######
## networkx analysis
"""
# use the networkx package to calculate for all frames:
#   - the number of connected components
#   - average degree (AKA average coordination number)
#   - largest connected component (LCC)
#   - average clustering coefficient
#   - average square clustering coefficient
"""

def primary_networkx_calc(filename):

  edge_output = data_outpath+'/frame-edges/edgelist' # + <#>.csv in f90
  edge1only_output = data_outpath+'/frame-edges-c1only/edgelist'
  edge2only_output = data_outpath+'/frame-edges-c2only/edgelist'
  edge_outs_list = [edge_output,edge1only_output,edge2only_output]

  tags = ["allc","c1only","c2only"]

  # get the number of particles
  traj = gsd.hoomd.open(filename, 'r')
  nframes = len(traj)

  colloids = np.where((traj[-1].particles.typeid == colloid1_typeid) | (traj[-1].particles.typeid == colloid2_typeid))[0]
  ncolloids = len(colloids)

  # do the same for all the colloid subpopulations
  colloid1 = np.where(traj[-1].particles.typeid == colloid1_typeid)[0]
  ncolloid1 = len(colloid1)
  colloid2 = np.where(traj[-1].particles.typeid == colloid2_typeid)[0]
  ncolloid2 = len(colloid2)
  if (ncolloid1 + ncolloid2) != ncolloids:
    print("ERROR: ncolloid1 + ncolloid2 != ncolloids in gofr calc; gofr was NOT calculated")
    return
  num_colloids = [ncolloids, ncolloid1, ncolloid2]

  index_arrays = [colloids, colloid1, colloid2]
  colloid_typeids = traj[-1].particles.typeid[colloids]

  # create empty dataframe for all network data  
  #network_frames_all_df = pd.DataFrame()
  #network_frames_c1only_df = pd.DataFrame()
  #network_frames_c2only_df = pd.DataFrame()
  #network_frames_dfs = [network_frames_all_df, network_frames_c1only_df, network_frames_c2only_df]
  network_frames_dfs = [[],[],[]]
  
  # import all data into one dataframe
  frame_dfs = [[], [], []]
  for frame in range(nframes):
    for i in range(len(edge_outs_list)): 
      data_path = edge_outs_list[i]
      # loop through all frames
      filepath = data_path+str(frame)+'.csv'
      # import CSV data
      df = pd.read_csv(filepath)
      # rename colums as needed
      df = df.rename(columns={"i": "source", "j": "target"})
      df.insert(loc=0, column='frame', value=frame)
      frame_dfs[i].append(df)

  frame_dfs_all = frame_dfs[0]
  frame_dfs_c1only = frame_dfs[1]
  frame_dfs_c2only = frame_dfs[2]

  alledge_df = pd.concat(frame_dfs_all, ignore_index=True) 
  alledge_c1only_df = pd.concat(frame_dfs_c1only, ignore_index=True) 
  alledge_c2only_df = pd.concat(frame_dfs_c2only, ignore_index=True) 
  all_edgelists = [alledge_df,alledge_c1only_df,alledge_c2only_df]

  # network analysis
  for frame in range(nframes):
    # get particle positions for percolation measurement
    pos_all = traj[frame].particles.position
    # all indices are in terms of the full list of colloids, so should always use the full list?
    pos_c = traj[frame].particles.position[colloids]
    pos_c1 = traj[frame].particles.position[colloid1]
    pos_c2 = traj[frame].particles.position[colloid2]
    pos_arrays = [pos_c, pos_c1, pos_c2]


    # box diagonal (for percolation)
    Lbox = traj[frame].configuration.box[:3]
    max_diameter = np.sqrt(Lbox[0]**2 + Lbox[1]**2 + Lbox[2]**2)

    for i in range(len(edge_outs_list)):

      edgelist = all_edgelists[i]

      ncolloids_curr = num_colloids[i]
      pos = pos_arrays[i]

      df = edgelist[edgelist['frame'] == frame][["source", "target"]]

      # create the network from edge list
      g = nx.from_pandas_edgelist(df)

      # if a node is not in the network, add it
      for particle in range(ncolloids_curr):
        if (  not(   g.has_node(particle)   )  ):
          g.add_node(particle)

      # calculate the total number of edges      
      nedges = nx.number_of_edges(g)

      # calculate the average degree of the network, i.e. avg contact number in gels
      avg_degree = 2 * nedges / ncolloids_curr

      # number of connected components
      n_cc = nx.number_connected_components(g)

      # return the indices of the nodes of the largest connected components
      lcc_nodes = max(nx.connected_components(g), key=len)

      # calculate the size of the largest connected component
      lcc_size = len(lcc_nodes)

      # physical diameter = maximum distance in the lcc
      lcc_pos = pos[list(lcc_nodes)]
      # use corrected positions centered on the center of mass
      centered_pos = pbc_centered_positions(lcc_pos, Lbox)
      physical_diameter = pdist(centered_pos).max() if len(centered_pos) > 1 else 0.0

      # --- per-axis extent of the LCC (crude anisotropy-aware span) ---
      # Reuses centered_pos so it shares a basis with physical_diameter.
      # extent/L -> 1 means the LCC fills that axis (1D-spanning proxy);
      # a compact blob has all three ~equal and < 1.
      if len(centered_pos) > 1:
          axis_extent = centered_pos.max(axis=0) - centered_pos.min(axis=0)  # (3,)
      else:
          axis_extent = np.zeros(3)
      axis_span_frac = axis_extent / Lbox            # per-axis, each in [0,1]
      max_axis_span_frac = float(axis_span_frac.max())
      n_axes_spanning = int((axis_span_frac > 0.85).sum())  # threshold tunable

      # compile outputs
      data={'network-type'           :[tags[i]],
            'frame'                  :[frame],
            'cut_off'                :[cut_off],
            'n_components'           :[n_cc], 
            'lcc_size'               :[lcc_size],
            'ncolloids'              :[ncolloids_curr],
            'f_lcc'                  :[lcc_size/ncolloids_curr],
            'lcc_span'               :[physical_diameter],
            'box_span'               :[max_diameter],
            'percolation'            :[physical_diameter/max_diameter],
            'lcc_extent_x'           :[float(axis_span_frac[0])],
            'lcc_extent_y'           :[float(axis_span_frac[1])],
            'lcc_extent_z'           :[float(axis_span_frac[2])],
            'max_axis_span'          :[max_axis_span_frac],
            'n_axes_spanning'        :[n_axes_spanning],
            'avg_degree'             :[avg_degree]
            }

      res_df = pd.DataFrame(data)
      network_frames_dfs[i].append(res_df) # = pd.concat([network_frames_dfs[i], res_df])
      #print(network_frames_dfs[i])

  network_frames_all_df = pd.concat(network_frames_dfs[0])
  network_frames_c1only_df = pd.concat(network_frames_dfs[1])
  network_frames_c2only_df = pd.concat(network_frames_dfs[2])
 
  # write all data to one CSV file
  network_frames_all_df.to_csv(data_outpath+'/networkx-allframes.csv',index = False)
  network_frames_c1only_df.to_csv(data_outpath+'/networkx-allframes-c1only.csv',index = False)
  network_frames_c2only_df.to_csv(data_outpath+'/networkx-allframes-c2only.csv',index = False)
  print("Primary network analysis complete")
#######


####################################
""" RUN CHECKS ON A SIMULATION """
####################################

if __name__ == '__main__':
  posCSV_calc(filepath)	
  edgelistCSV_calc(filepath)	
  primary_networkx_calc(filepath)

  print("\nAll network analyses complete.")	
