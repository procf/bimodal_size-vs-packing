########################
""" MODULE LIBRARY """
########################
import numpy as np
import pandas as pd
import igraph as ig
import networkx as nx
import os
import re
from collections import Counter
from collections import defaultdict


########################
""" DATA HANDLING """
########################

small_comps = [100, 75, 50, 25, 0]
compositions = ['100-0', '75-25', '50-50', '25-75' '0-100']
size_ratios = ['1-0', '1-2', '0-2']
volume_fractions = ['20']
attraction_strengths = ['12'] 
depletion_scaling = ['bysize'] 
attraction_range = ['60']
#clustering_style = ['standard','weighted']
clustering_style = ['standard']

eta0 = 0.3
kT = 0.1
L = 70
system_volume = L**3

data_outpath = "data"
# create "data" subfolder if it doesn't exit
if os.path.exists(data_outpath) == False:
  os.mkdir(data_outpath)


#################
# FUNCTIONS
#################

# =============================
# Find all tetrahedral cycles 
# =============================

def canonical_cycle(path):
    """
    Use canonical rotation to avoid duplicate cycle counts
    """
    L = len(path)
    rotations = [tuple(path[i:] + path[:i]) for i in range(L)]
    rev = path[::-1]
    rev_rotations = [tuple(rev[i:] + rev[:i]) for i in range(L)]
    return min(rotations + rev_rotations)

def find_cycles_upto_k(g, k=6):
    """
    Find simple cycles up to length k in an igraph graph.

    Returns list of cycles (as lists of vertex indices).
    """

    cycles = set()

    def dfs(start, current, visited, path):
        if len(path) > k:
            return

        for neighbor in g.neighbors(current):
            if neighbor == start and len(path) >= 3:
                # Found a cycle
                #cyc = tuple(sorted(path))
                cyc = canonical_cycle(path)
                cycles.add(cyc)
            elif neighbor not in visited:
                dfs(start, neighbor, visited | {neighbor}, path + [neighbor])

    for v in range(g.vcount()):
        dfs(v, v, {v}, [v])

    # convert back to list of lists
    return [list(c) for c in cycles]


# ========================================
# Get bond composition of all tetrahedra 
# ========================================

def loop_composition_with_bonds(g, max_len=8, normalize="nodes", edge_type_attr="bond_type"):
    """
    Compute loop composition and bond-type composition within loops.

    Parameters
    ----------
    g : igraph.Graph
    max_len : int
    normalize : str ("nodes", "edges", or None)
    edge_type_attr : str ("bond_type")
        Edge attribute storing bond type ("SS", "SL", "LL")

    Returns
    -------
    dict:
        {
          cycle_length: {
              "count": int,
              "density": float,
              "bond_counts": {"SS": int, "SL": int, "LL": int},
              "bond_fraction": {"SS": float, "SL": float, "LL": float}
          }
        }
    """

    #cycles = g.cycle_basis() # networkX command
    cycles = find_cycles_upto_k(g, k=6)

    # normalization
    if normalize == "nodes":
        norm = g.vcount()
    elif normalize == "edges":
        norm = g.ecount()
    else:
        norm = 1.0

    result = {}

    for cyc in cycles:
        L = len(cyc)
        if L > max_len:
            continue

        # initialize
        if L not in result:
            result[L] = {
                "count": 0,
                "bond_counts": {"SS": 0, "SL": 0, "LL": 0}
            }

        result[L]["count"] += 1

        # loop edges (wrap around!)
        for i in range(L):
            u = cyc[i]
            v = cyc[(i + 1) % L]

            eid = g.get_eid(u, v)
            btype = g.es[eid][edge_type_attr]

            result[L]["bond_counts"][btype] += 1

    # finalize densities + fractions
    for L in result:
        count = result[L]["count"]
        result[L]["density"] = count / norm

        total_bonds = sum(result[L]["bond_counts"].values())
        if total_bonds > 0:
            result[L]["bond_fraction"] = {
                k: v / total_bonds
                for k, v in result[L]["bond_counts"].items()
            }
        else:
            result[L]["bond_fraction"] = {"SS": 0, "SL": 0, "LL": 0}

    return result

def flatten_loop_results(loop_dict, prefix="loop"):
    """
    Flatten loop composition dictionary into CSV-friendly format.

    Parameters
    ----------
    loop_dict : dict
        Output of loop_composition_with_bonds
    prefix : str
        Prefix for column names

    Returns
    -------
    flat_dict : dict
    """

    flat = {}

    for L, data in loop_dict.items():
        base = f"{prefix}{L}"

        # basic metrics
        flat[f"{base}_count"] = data.get("count", 0)
        flat[f"{base}_density"] = data.get("density", 0.0)

        # bond fractions (preferred for comparison)
        bf = data.get("bond_fraction", {})
        flat[f"{base}_SS_frac"] = bf.get("SS", 0.0)
        flat[f"{base}_SL_frac"] = bf.get("SL", 0.0)
        flat[f"{base}_LL_frac"] = bf.get("LL", 0.0)

        # (optional) raw counts if you want them
        bc = data.get("bond_counts", {})
        flat[f"{base}_SS_count"] = bc.get("SS", 0)
        flat[f"{base}_SL_count"] = bc.get("SL", 0)
        flat[f"{base}_LL_count"] = bc.get("LL", 0)

    return flat


# =================================
# Classify tetrahedral aggregates 
# =================================

def build_tetrahedron_graph(tetrahedra):
    """
    Build adjacency graph of tetrahedra.

    Returns:
        tg: igraph.Graph where nodes = tetrahedra
    """
    tg = ig.Graph()
    tg.add_vertices(len(tetrahedra))

    # Convert to sets for fast overlap checks
    tet_sets = [set(t) for t in tetrahedra]

    edges = []

    for i in range(len(tetrahedra)):
        for j in range(i + 1, len(tetrahedra)):
            overlap = len(tet_sets[i].intersection(tet_sets[j]))

            # Face-sharing (3 nodes) or edge-sharing (2 nodes)
            if overlap >= 2:
                edges.append((i, j))

    tg.add_edges(edges)

    return tg

def extract_particle_subgraph(g, tetrahedra, component):
    nodes = set()
    for idx in component:
        nodes.update(tetrahedra[idx])

    subg = g.subgraph(list(nodes))
    return subg

def classify_clusters(g, tetrahedra, tg):
    comps = tg.components()
    cluster_info = []
    
    for comp in comps:
        size = len(comp)
        
        if size == 1:
            label = "ρ1T"
        elif size == 2:
            label = "ρ2T"  # always a single edge (face- or edge-share)
        else:
            # Get the subgraph of tg restricted to this component
            comp_subgraph = tg.subgraph(comp)
            n_edges = comp_subgraph.ecount()
            n_nodes = comp_subgraph.vcount()  # = size
            max_degree = max(comp_subgraph.degree())
            
            if size == 3:
                # Two possibilities: chain (2 edges) or triangle (3 edges)
                if n_edges == 2:
                    label = "ρ3T_chain"
                elif n_edges == 3:
                    label = "ρ3T_triangle"
                else:
                    label = "ρ3T_other"  # shouldn't happen
            
            elif size == 4:
                # Possibilities: linear chain (3 edges, max_deg=2),
                #                star/branch (3 edges, max_deg=3),
                #                cycle (4 edges, max_deg=2),
                #                triangle+pendant (4 edges, max_deg=3),
                #                K_4 (6 edges)
                if n_edges == 3:
                    if max_degree == 2:
                        label = "ρ4T_chain"
                    else:  # max_degree == 3 
                        label = "ρ4T_tree"
                elif n_edges == 4:
                    if max_degree == 2:
                        label = "ρ4T_cycle"
                    else:
                        label = "ρ4T_triangle_pendant"
                else:
                    label = "ρ4T_dense"  # >4 edges = denser than tree+1
            
            elif size == 5:
                # Pentagonal bipyramid is a specific particle-level pattern
                subg = extract_particle_subgraph(g, tetrahedra, comp)
                if subg.vcount() == 7 and subg.ecount() == 16:
                    label = "ρ5T_bipyramid"
                else:
                    # if not penta-bipyramid, classify by tg topology instead
                    if n_edges == 4 and max_degree == 2:
                        label = "ρ5T_chain"
                    elif n_edges == 4 and max_degree > 2:
                        label = "ρ5T_tree"
                    else:
                        label = "ρ5T_dense"
            
            else:  # size >= 6
                # For large components, summarize by topology rather than enumerate
                # cyclomatic = E - N + 1 (number of independent cycles)
                cyclomatic = n_edges - n_nodes + 1
                if cyclomatic == 0:
                    if max_degree == 2:
                        label = f"ρ{size}T_chain"
                    else:
                        label = f"ρ{size}T_tree"  # branched but acyclic
                else:
                    label = f"ρ{size}T_dense"  # has cycles
        
        cluster_info.append((comp, label))
    
    return cluster_info

def count_motifs(t_info):
    labels = [label for _, label in t_info]
    return Counter(labels)

def tetra_edges_from_nodes(g, tet):
    """
    Return set of edge IDs belonging to one tetrahedron (4 nodes).
    """
    edges = set()

    for i in range(len(tet)):
        for j in range(i+1, len(tet)):
            u, v = tet[i], tet[j]
            eid = g.get_eid(u, v, directed=False, error=False)
            if eid != -1:
                edges.add(eid)

    return edges

def motif_edge_set(g, tetrahedra, comp):
    """
    Exact union of tetrahedral edges for one motif component.
    """
    edges = set()

    for tidx in comp:
        tet = tetrahedra[tidx]
        edges |= tetra_edges_from_nodes(g, tet)

    return edges

def bond_type_counts_from_edges(g, edge_ids):
    counts = {"SS":0, "SL":0, "LL":0}

    for eid in edge_ids:
        bt = g.es[eid]["bond_type"]
        counts[bt] += 1

    return counts

def analyze_motifs(g, tetrahedra, tg, t_info):
    """Per-component analysis. Stores edge_ids for downstream overlap computations."""
    results = []
    for comp, label in t_info:
        edge_ids = motif_edge_set(g, tetrahedra, comp)
        edge_ids_set = set(edge_ids)

        nodes = set()
        for tidx in comp:
            nodes.update(tetrahedra[tidx])

        bond_counts = bond_type_counts_from_edges(g, edge_ids)
        total_bonds = sum(bond_counts.values())
        bond_fraction = {
            k: (v / total_bonds if total_bonds > 0 else 0)
            for k, v in bond_counts.items()
        }

        if len(comp) > 1:
            comp_subgraph = tg.subgraph(comp)
            n_adjacencies = comp_subgraph.ecount()
            cyclomatic = n_adjacencies - len(comp) + 1
            max_degree = max(comp_subgraph.degree())
        else:
            n_adjacencies = 0
            cyclomatic = 0
            max_degree = 0

        results.append({
            "label": label,
            "n_tetra": len(comp),
            "n_adjacencies": n_adjacencies,
            "cyclomatic": cyclomatic,
            "max_degree": max_degree,
            "n_nodes": len(nodes),
            "n_edges": len(edge_ids),
            "edge_ids": edge_ids_set,
            "bond_counts": bond_counts,
            "bond_fraction": bond_fraction,
        })
    return results


def parse_topology_class(label):
    """
    Map a fine-grained label to a coarse topology class.
    
    Returns one of: 'isolated', 'pair', 'chain', 'tree', 'cycle', 'dense', 
                    'bipyramid', 'triangle', 'triangle_pendant', 'other'
    """
    if label == "ρ1T":
        return "isolated"
    elif label == "ρ2T":
        return "pair"
    elif label == "ρ5T_bipyramid":
        return "bipyramid"
    elif "_chain" in label:
        return "chain"
    elif "_tree" in label or "_branch" in label:
        return "tree"
    elif "_cycle" in label:
        return "cycle"
    elif "_dense" in label:
        return "dense"
    elif "_triangle_pendant" in label:
        return "triangle_pendant"
    elif "_triangle" in label:
        return "triangle"
    else:
        return "other"


def build_tetrahedra_dataframe(results, config_label=None):
    """Build per-component dataframe. Auto-discovers overlap fields."""
    rows = []
    
    standard_keys = {
        'label', 'n_tetra', 'n_adjacencies', 'cyclomatic', 'max_degree',
        'n_nodes', 'n_edges', 'edge_ids', 'bond_counts', 'bond_fraction'
    }
    
    for comp_id, r in enumerate(results):
        row = {}
        if config_label is not None:
            row.update(config_label)
        
        row['component_id'] = comp_id
        row['label'] = r['label']
        row['topology_class'] = parse_topology_class(r['label'])
        row['n_tetra'] = r['n_tetra']
        row['n_adjacencies'] = r['n_adjacencies']
        row['cyclomatic'] = r['cyclomatic']
        row['max_degree'] = r['max_degree']
        row['n_nodes'] = r['n_nodes']
        row['n_edges'] = r['n_edges']
        
        bc = r['bond_counts']
        row['bond_count_SS'] = bc.get('SS', 0)
        row['bond_count_SL'] = bc.get('SL', 0)
        row['bond_count_LL'] = bc.get('LL', 0)
                
        rows.append(row)
    
    return pd.DataFrame(rows)

# ===============
# Main Analysis
# ===============

def analyze_topology(g):

    # TETRAHEDRA (exact K4)
    tetrahedra = [c for c in g.cliques(min=4, max=4)]  # already exact

    # EDGE SETS
    def has_edge(u, v):
        return g.are_adjacent(u, v)

    def get_edges_from_nodes(nodes):
        edges = set()
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                if has_edge(nodes[i], nodes[j]):
                    edges.add(g.get_eid(nodes[i], nodes[j]))
        return edges

    tetra_edges = set().union(*[get_edges_from_nodes(c) for c in tetrahedra]) if tetrahedra else set()

    # BOND COMPOSITION 
    def structure_composition(cliques):
        comp_counts = {'SS':0, 'SL':0, 'LL':0}
        total_edges = 0

        for c in cliques:
            edges = get_edges_from_nodes(c)
            for eid in edges:
                bond_type = g.es[eid]["bond_type"]
                comp_counts[bond_type] += 1
                total_edges += 1

        if total_edges == 0:
            return {k:0 for k in comp_counts}

        return {k: v/total_edges for k,v in comp_counts.items()}

    # comp in and out of tetrahedra
    def edge_set_composition(edge_set):
        comp_counts = {'SS':0, 'SL':0, 'LL':0}

        for eid in edge_set:
            bond_type = g.es[eid]["bond_type"]
            comp_counts[bond_type] += 1

        total = sum(comp_counts.values())

        if total == 0:
            return {k:0 for k in comp_counts}

        return {k: v/total for k,v in comp_counts.items()}

    # unique edges only
    tetra_comp = edge_set_composition(tetra_edges)

    ## TETRAHEDRAL AGGREGATES/MOTIFS
    tg = build_tetrahedron_graph(tetrahedra)
    t_info = classify_clusters(g, tetrahedra, tg)
    pt_motif_counts = count_motifs(t_info)
    pt_results = analyze_motifs(g, tetrahedra, tg, t_info)

    # edge -> ALL motif memberships
    edge_to_labels = {e.index: set() for e in g.es}

    for tidx, tet_nodes in enumerate(tetrahedra):

        label = tetra_idx_to_label.get(tidx, "other")

        tet_edges = get_edges_from_nodes(tet_nodes)

        for eid in tet_edges:
             edge_to_labels[eid].add(label)

    # ==========================================================
    # build dataframe
    # ==========================================================
    rows = []
    for e in g.es:
        eid = e.index
        i, j = e.tuple
        labels = edge_to_labels[eid]
        rows.append({
            "eid": eid,
            "i": i,
            "j": j,
            "type": e["bond_type"],
            # tetra motif flags
            "in_isolated": "ρ1T" in labels,
            "in_pair": "ρ2T" in labels,
            "in_chain": any("_chain" in lab for lab in labels),
            "in_tree": any("_tree" in lab or "_branch" in lab for lab in labels),
            "in_cycle": any("_cycle" in lab for lab in labels),
            "in_dense": any("_dense" in lab for lab in labels),
            "in_bipyramid": "ρ5T_bipyramid" in labels,
            "in_triangle": any("ρ3T_triangle" in lab for lab in labels),
            "in_triangle_pendant": any("_triangle_pendant" in lab for lab in labels),
            # Derived flags
            "multiT": len(labels) > 1,
            "noT": len(labels) == 0,
            # Readable membership list
            "motifs": sorted(labels)[0] if labels else "noT",
        })
    edge_data_df = pd.DataFrame(rows)

    return {
        "n_tetra": len(tetrahedra),
        'tetrahedra': tetrahedra,
        'tetra_comp': tetra_comp,

        't_info': t_info,
        'pt_motif_counts': pt_motif_counts,
        'pt_results': pt_results,

        'edge_data_df': edge_data_df,
    }
 


###########################
""" ANALYZE ALL SYSTEMS """
###########################
tetra_dfs = []
edge_data_dfs = []
system_dfs = []
for size in size_ratios:
  radii = size.split('-')
  R_C1 = float(radii[0])
  R_C2 = float(radii[1])
  for phi in volume_fractions:
      for comp_1 in small_comps:
          comp_2 = 100 - comp_1
          comp = f'{comp_1}-{comp_2}'
          # skip invalid fake mixtures for monodisperse case
          #print(size, comp, R_C1, R_C2, f'0-{int(R_C2)}', (size == f'0-{int(R_C2)}'))
          if size != '1-0' and comp == '100-0':
              #print('...skip bimodal size or mono_large size for mono_small comp')
              continue
          if size == '1-0' and comp != '100-0':
              #print('---skip mono_small size for bimodal comp or mono_large comp')
              continue
          if size != f'0-{int(R_C2)}' and comp == '0-100':
              #print('++++skip bimodal size or mono_small size for mono_large comp')
              continue
          if size == f'0-{int(R_C2)}' and comp != '0-100':
              #print('////skip mono_large size for bimodal comp or mono_small comp')
              continue
          if size == f'0-{int(R_C2)}' and comp == '0-100':
              size = f'{int(R_C2)}-0'
              #print('reformat size for mono_large size and mono_large comp')
              if R_C2 != float(int(R_C2)):
                  print('ERROR: large monomodal data processing script cannot accept float particle size: {R_C2}')
                  continue
          #print(size, comp)
          for D0 in attraction_strengths:
              for scaling in depletion_scaling:
                  for kappa in attraction_range:
                      if scaling == 'uniform':
                          if (size == '1-0') and (comp == '100-0'):
                             print(f'WARNING: monodisperse data only recognizes scaling "bysize", not {scaling}')
                             continue
                          elif (size == f'{int(R_C2)}-0') and (comp == '0-100'):
                             print(f'WARNING: monodisperse data only recognizes scaling "bysize", not {scaling}')
                             continue
                          elif comp != '100-0':
                             data_path = (f'/projects/props/Rob/colloids/bimodal/uniform-potential/r{size}'
                                f'/DPD/L70/poly0.0/seedNone/phi{phi}/{comp}/potential-morse/'
                                f'{D0}kT/kappa{kappa}/analysis-bi-colloids/data/') 
                      else:
                          if (size == '1-0') and (comp == '100-0'):
                             if scaling == 'bysize':
                                data_path = (f'/projects/props/Rob/colloids/bimodal/r{size}/DPD/L{L}/mono/seedNone/phi{phi}'
                                   f'/potential-morse-bysize/{D0}kT/kappa{kappa}/analysis-DPD/data/')
                                print(data_path)
                             else:
                                print(f'WARNING: monodisperse data only recognizes scaling "bysize", not {scaling}')
                                continue
                          elif (size == f'{int(R_C2)}-0') and (comp == '0-100'):
                             if scaling == 'bysize':
                                data_path = (f'/projects/props/Rob/colloids/bimodal/r{size}/DPD/L{L}/phi{phi}'
                                   f'/potential-morse-bysize/{D0}kT/kappa{kappa}/analysis-DPD/data/')
                             else:
                                print(f'WARNING: monodisperse data only recognizes scaling "bysize", not {scaling}')
                                continue
                          elif (comp != '100-0') and (comp != '0-100'):
                             data_path = (f'/projects/props/Rob/colloids/bimodal/r{size}'
                                f'/DPD/L70/poly0.0/seedNone/phi{phi}/{comp}/potential-morse-{scaling}/'
                                f'{D0}kT/kappa{kappa}/analysis-bi-colloids/data/') 
                          else:
                             continue

                      if os.path.exists(data_path) == False:
                          print(f"No data path: {data_path}")
                          continue

                      graph_path = f"{data_path}graph_lastframe.graphml"
                      if os.path.exists(graph_path) == False:
                          print(f"No GraphML file: {graph_path}")
                          print(f"...createing GraphML file")

                          gsd_path = f"{data_path}Gelation_Colloids.gsd"
                          if os.path.exists(gsd_path) == False:
                              print(f" - ERROR: Missing GSD file {gsd_path}")
                              continue

                          edgelist_path = f"{data_outpath}frame-edges"
                          if os.path.exists(edgelist_path) == False:
                              print(f" - ERROR: Missing edgelist data at")
                              print(f"    {edgelist_path}")
                              continue

                          else:
                              # make the graph
                              if "DPD" in data_path:
                                colloid1_typeid = 1
                                colloid2_typeid = 2
                              else:
                                colloid1_typeid = 0
                                colloid2_typeid = 1

                              cut_off = round(3/int(kappa),2)

                              # get the number of particles
                              traj = gsd.hoomd.open(gsd_path, 'r')
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
                                exit(1)

                              # import all data into one dataframe
                              frame_dfs_all = []
                              #pos_frame_dfs_all = []

                              edge_output = f'{data_outpath}/frame-edges/edgelist' # + <#>.csv in f90
                              #pos_output = f'{data_outpath}/frame-pos/positions_frame' # + <#>.csv in f90
                              for frame in range(nframes):
                                  # loop through all frames
                                  edge_file = edge_output+str(frame)+'.csv'
                                  # import CSV data
                                  edge_df = pd.read_csv(edge_file)
                                  # rename colums as needed
                                  edge_df = edge_df.rename(columns={"i": "source", "j": "target"})
                                  edge_df.insert(loc=0, column='frame', value=frame)
                                  frame_dfs_all.append(edge_df)

                              alledge_df = pd.concat(frame_dfs_all, ignore_index=True)

                              # network analysis on last frame only
                              frame = nframes-1

                              # get particle positions for percolation measurement
                              pos = traj[frame].particles.position
                              typeID = traj[frame].particles.typeid
                              radii = 0.5*traj[frame].particles.diameter

                              df = alledge_df[alledge_df['frame'] == frame][["source", "target"]]

                              # create the network from edge list
                              g = nx.from_pandas_edgelist(df)

                              # if a node is not in the network, add it
                              for particle in range(ncolloids):
                                if (  not(   g.has_node(particle)   )  ):
                                  g.add_node(particle)

                              # add node attributes
                              node_attrs = {
                                  i: {
                                      "type": int(typeID[i]),
                                      "radius": float(radii[i]),
                                      "x": float(pos[i][0]),
                                      "y": float(pos[i][1]),
                                      "z": float(pos[i][2]),
                                  }
                                  for i in range(ncolloids)
                              }

                              nx.set_node_attributes(g, node_attrs)

                              #graph_name = f"{data_outpath}graph_frame{frame}.graphml"
                              graph_name = f"{data_outpath}graph_lastframe.graphml"
                              nx.write_graphml(g, graph_name)
                              print(f"Saved GraphML → {graph_name}")

                      else:

                          # load the graph
                          g = ig.Graph.Read_GraphML(graph_path)

                          cut_off = round(3/int(kappa),2)

                          results = analyze_topology(g)
                          tetra_comp = results['tetra_comp']

                          #pt_motif_counts = results['pt_motif_counts']
                          #pt_motif_counts["ρ1T"] = pt_motif_counts.get("ρ1T", 0)
                          #pt_motif_counts["ρ2T"] = pt_motif_counts.get("ρ2T", 0)
                          #pt_motif_counts["ρ3T"] = pt_motif_counts.get("ρ3T", 0)
                          #pt_motif_counts["ρ5T"] = pt_motif_counts.get("ρ5T", 0)
                          #pt_motif_counts["other"] = pt_motif_counts.get("other", 0)

                          results_dict = {'size-ratio'          :[size],
                                          'phi'                 :[phi],
                                          'composition'         :[comp],
                                          'D0'                  :[D0],
                                          'scaled_D0'           :[scaling],
                                          'kappa'               :[kappa],

                                          'n_tetra': results['n_tetra'],
                                          'tetra_comp_SS': tetra_comp['SS'],
                                          'tetra_comp_SL': tetra_comp['SL'],
                                          'tetra_comp_LL': tetra_comp['LL'],

                                          #'p1T_counts': pt_motif_counts["ρ1T"],
                                          #'p2T_counts':pt_motif_counts["ρ2T"],
                                          #'p3T_counts':pt_motif_counts["ρ3T"],
                                          #'p5T_counts':pt_motif_counts["ρ5T"],
                                          #'pTother_counts':pt_motif_counts["other"],
                          }

                          # After computing pt_results = analyze_motifs(g, tetrahedra, t_info):
                          pt_results = results['pt_results']
                          tetrahedra_df = build_tetrahedra_dataframe(
                              pt_results,
                              config_label={
                                  'size-ratio': size,
                                  'phi': phi,
                                  'composition': comp,
                                  'D0': D0,
                                  'scaled_D0': scaling,
                                  'kappa': kappa,
                              }
                          )
                          tetra_dfs.append(tetrahedra_df)

                          res_df = pd.DataFrame(results_dict)
                          system_dfs.append(res_df)

                          edge_data_df = results['edge_data_df']

                          print(f' - topological analysis complete for: {R_C1}:{R_C2}, phi={phi}, {comp}, {D0}kT, kappa={kappa}')


all_results_df = pd.concat(system_dfs)
all_results_df.to_csv(f'{data_outpath}/topology.csv',index = False)
print('All results saved to "topology.csv"')

all_tetra_df = pd.concat(tetra_dfs, ignore_index=True)
all_tetra_df.to_csv(f'{data_outpath}/tetrahedra.csv',index = False)
print('Detailed tetrahedral structure results saved to "tetrahedra.csv"')

if len(edge_data_dfs) != 0:
  edge_labels_df = pd.concat(edge_data_dfs)
  edge_labels_df.to_csv(f'{data_outpath}/all_edge_data.csv',index = False)
  print('Labeled edge data saved to "all_edge_data.csv"')
