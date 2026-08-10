# Size matters more than packing in bimodal colloidal gel compositions
This is a simulation and analysis pipeline for evaluating the multi-scale structure of bimodal colloidal depletion gels.

Simulations can be run using the [sim-scripts](./sim-scripts) and the DPDMorse extension for HOOMD-blue v4.2.1, available as [hoomd4.2.1-mod](https://github.com/procf/hoomd4.2.1-mod).

## What to expect
For systems of ~10,000 colloidal particles, use [analysis-scripts](./analysis-scripts) to calculate:
- Coordination number distribution and average coordination number
- Void size distribution and average void size
- Network edgelist, the size of connected components, and the physical network diameter ("span" as a proportion of total box size)
- Angle distribution within the network
- Identify tetrahedral structures and classify their aggregates
- Gaussian Mixture Model (GMM) based mesoscale clustering
- Cauchy-Born estimate of the total elastic modulus from mesoscale (cluster) structure

Note: Standard analyses typically take less than 20min. GMM clustering can take ~9hrs. Cauchy-Born estimate takes <1min.

In bimodal depletion gels, colloid-colloid interactions scale with particle size. Therefore, the attraction strength for S-S, S-L, and L-L contacts will be different.
For size ratio 1:2 this scales roughly with the arithmetic mean, such that $D_0^{SL} = 1.5 D_0^{SS}$ and $D_0^{LL} = 2 D_0^{SS}$.

## Software/package requirements
In this project, the following packages are actively used:
1. GNU Fortran (GCC) 11.4.1
2. `python` v3.10.16 
3. `gsd` v3.2.1
4. `numpy` v1.24.40
5. `pandas` v2.2.3
6. `networkx` v3.4.2
7. `scipy` v1.15.3
8. `node2vec` v0.5.0
9. `umap-learn` v0.5.9
10. `scikit-learn` v1.7.2

## Hardware/OS tested
The program was tested on a single HPC-node running Rocky Linux 9.3 (kernel 5.14).

## Background

Colloidal gels have complex mechanical and rheological properties that emerge from the underlying structure of their particle network.
These structure-mechanics relationships are widely studied with models of uniform, monodisperse particles; however, most real colloidal 
applications use polydisperse and/or multi-modal formulations.
We use large-scale simulations to study how the addition of a larger particle population (size ratio 1:2) alters multi-scale gel structure. 
We vary the large-particle fraction at a fixed total volume fraction and fixed interaction condition. There is a significant increase in 
local packing as large particle fraction increases, but this does not produce enhanced rigidity percolation. Instead the mechanics of these 
systems appear to weaken due to changes in structural scale that are proportional to the mean particle size.


Here's the abstract: 

*Colloidal gels are frequently modeled as monodisperse particle networks, although practical formulations 
commonly contain particles with multiple characteristic sizes. Here, we use large-scale, hydrodynamically 
resolved simulations of colloidal depletion gels to isolate the effects of particle size and local packing 
in bimodal systems with a small-to-large size ratio of 1:2. Increasing the large-particle fraction introduces 
new heterotypic angular motifs and substantially increases the fraction of bonds participating in tetrahedral 
structures, with a maximum at intermediate composition. However, these additional rigid motifs do not 
reorganize into larger or more highly connected tetrahedral aggregates. The mean coordination and characteristic 
aggregate size remain nearly composition independent. By contrast, the void and cluster-size distributions coarsen
systematically as the large-particle fraction increases. These mesoscale distributions largely collapse when
normalized by a composition-dependent particle length scale, indicating that changes in composition primarily 
rescale gel architecture rather than producing distinct rigid-network topologies. An elastic modulus estimated 
using Cauchy–Born theory similarly follows this effective length scale more closely than the abundance of local 
tetrahedral motifs. These results show that, for moderate size disparity, particle size controls the structural
scale and predicted mechanical response of bimodal colloidal gels more strongly than enhanced local packing.*


## Contributors
This project was in collaboration with the [Soft Matter Engineering Laboratory](https://smel.eng.uci.edu/) at the University of California, Irvine. \
This work was done by [Rob Campbell](https://scholar.google.com/citations?user=i8S54zYAAAAJ&hl=en), Calvin (Ziye) Zhuang, 
[Ali Mohraz](https://scholar.google.com/citations?user=pW80NaAAAAAJ&hl=en), and [Safa Jamali](https://scholar.google.com/citations?user=D1asaYIAAAAJ&hl=en).\
Authors acknowledge support from the National Science Foundation (PMP-2025613) and NASA ROSES FINESST (80NSSC23K0015).
