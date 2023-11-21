# melanoma_dynamic_curvature
Analysis code for study: Multi-Scale Geometric Network Analysis Identifies Melanoma Immunotherapy Response Gene Modules. By K.A. Murgas et al. (under submission)

This repository contains the following files used for analysis:
1. prep_GSE91061.R was used to pre-process the publically available RNA-seq dataset for analysis. After downloading dataset GSE91061 from NCBI GEO (link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061), this script will prepare data including filtering and quantile normalization to reproduce the data analyzed in this study
2. dynamic_melanoma_parallel_nbrs.py was used to construct correlation network and determine dynamic Ollivier-Ricci curvature in the network, then save results for each edge
3. topo_stringdb11sparse.mat is the data file for network topology used in this study. The methodology was previously described in Banerji et al. 2013 SciRep (see supplemental information), and later adapted by Murgas et al. 2022 SciRep
4. melanoma_clustering_louvain.py was then used to determine gene modules and create several subfigures and statistical analysis
5. melanoma_cluster_analysis.R was used to create all other figures and statistical analyses
