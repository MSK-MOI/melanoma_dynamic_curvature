# script for examining melanoma immunotherapy dynamic curvature analysis results

# analysis code written by Kevin A Murgas
# Stony Brook University, Dept. of Biomedical Informatics
# advisor: Allen Tannenbaum PhD
# in collaboration with: Rena Elkin PhD, Nadeem Riaz PhD, Emil Saucan PhD, Joseph Deasy PhD

#%%
import os
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import networkx as nx

# load in orc values from dynamic curvature parallelized script
results_dir = "../Results/melanoma_qnorm_pearsonshiftinverse_PD1Bnbrs"
genes = pd.read_csv( os.path.join(results_dir, "genes.csv")).values.squeeze()
tau = pd.read_csv( os.path.join(results_dir, "tau.csv"))
orc_diff = pd.read_csv( os.path.join(results_dir, "orc_diff.csv"), index_col=[0,1])
corr_method = 'pearson'

#%%
# recompute graph distances
import scipy.io, time, scipy.sparse as sp
from scipy.sparse import csgraph as cg
expr_file = '../Data/GSE91061_rld_qnorm_HGNC.csv'
ppi_file = '../Data/topo_stringdb11sparse.mat'
genelist_file = '../gene_lists/WP_CANCER_IMMUNOTHERAPY_BY_PD1_BLOCKADE.v2023.1.Hs.grp'

# load expression data
print("\nLoad in gene expression data from melanoma immunotherapy study GSE91061...")
expr_data = pd.read_csv( expr_file, index_col=0)
expr_hgnc = expr_data.index.values
print(" Expression data loaded. %i genes x %i samples." % expr_data.shape )

# determine pre- and on- treatment groups from column labels
ids_data = expr_data.columns.values
pat_id, condition = [np.array(e) for e in zip(*[id.split("_")[0:2] for id in ids_data])]


# CONSTRUCT CANCER GENE NETWORK TOPOLOGY FROM DATA

# load STRINGdb protein-protein interaction network (sparsified)
print("Load in PPI topology (file: %s)" % ppi_file)
ppi_data = scipy.io.loadmat(ppi_file, simplify_cells = True)
adj = sp.csr_matrix(ppi_data['adj'].astype('bool'))
topoNames = ppi_data['topoNames'] # need to convert string array to char array in order for python to read
topoNames = [s.strip() for s in topoNames] # remove white space
print(" PPI loaded. %i vertices, %i edges" % (adj.shape[0] , adj.nnz))

# load gene list to subset
print("Load in gene list to subset expression data...")
genelist = pd.read_csv(genelist_file, skiprows=2, names=["Gene"])
genelist_hgnc = genelist['Gene'].values.astype(str)
print(" Gene list loaded. %i unique genes." % len(np.unique(genelist_hgnc)) )

# find matching genes in PPI topology genes + expression data + subset gene list
matching_genes = list(set(topoNames) & set(expr_hgnc) & set(genelist_hgnc))
matching_genes.sort()
print("\nMatched %i genes in expression, MGDB, and PPI." % len(matching_genes) )

# optional flag to include all 1-hop neighbors of subset genes
neighbors_flag = True
if (neighbors_flag == True):
    # find indices of matching genes in ppi topoNames
    ind_ppi = [topoNames.index(g) for g in matching_genes]

    # get indices of neighbors of matching genes
    nbrs_all = np.split(adj.indices, adj.indptr)[1:-1] # array of array of neighbors at each vertex/gene
    ind_nbrs = np.unique(np.concatenate(np.array(nbrs_all, dtype=object)[ind_ppi])) # get unique neighbors of matched genes
    
    # combine with the initial matching genes and get unique gene names
    ind_match_nbrs = np.unique(np.concatenate( (ind_ppi, ind_nbrs) ))
    match_nbrs_hgnc = np.array(topoNames)[ind_match_nbrs]
    matching_genes = list(set(expr_hgnc) & set(match_nbrs_hgnc)) # only those in expr_hgnc
    matching_genes.sort()
    print(" included neighbors for total of %i genes" % len(matching_genes))

# filter out genes with initial degree < 5
degree_orig = adj.A.sum(axis=1)
ind_ppi = np.array([topoNames.index(g) for g in matching_genes])
matching_genes = np.array(topoNames)[ind_ppi[degree_orig[ind_ppi]>=5]].tolist()
matching_genes.sort()
print(" Retain %i genes with degree >= 5." % len(matching_genes) )

# find indices of matching genes in expression gene names and ppi topoNames
ind_exp = [list(expr_hgnc).index(m) for m in matching_genes]
ind_ppi = [topoNames.index(m) for m in matching_genes]

# subset expression data and adjacency matrix
expr_match = expr_data.iloc[ind_exp,:]
adj_match = adj[ind_ppi,:][:,ind_ppi]

# take largest connected component and subset
cc = cg.connected_components(adj_match)[1]
ind_maxcc = np.where(cc == np.bincount(cc).argmax())[0]
expr_max = expr_match.iloc[ind_maxcc,:]
adj_max = adj_match[ind_maxcc,:][:,ind_maxcc]
geneNames_max = [expr_hgnc[i] for i in (ind_exp[j] for j in ind_maxcc)]

# determine graph topological properties:
nv = adj_max.shape[0] # number of vertices
ne = adj_max.nnz # number of edges
lv,rv = adj_max.nonzero() # get indices of left vertices (rows) and right vertices (columns)
print(" Largest connected component: %i vertices, %i edges" % (nv, ne))


#%%
# construct correlation network based on matched expression changes (on-treatment minus pre-treatment)

# find matched patient samples
# compute difference in expression in matched patients
# matched_pats[on] - match_pats[pre]
print('Compute difference of matched patients')
match_pats = np.intersect1d(pat_id[condition == 'Pre'],
                            pat_id[condition == 'On'])
pre_on_ind = [[np.where((pat_id == p) & (condition == 'Pre') )[0].item(),
            np.where((pat_id == p) & (condition == 'On') )[0].item()] for p in match_pats]
pre_ind, on_ind = [np.array(e) for e in zip(*pre_on_ind)]
data_diff = pd.DataFrame(expr_max.iloc[:, on_ind].values - expr_max.iloc[:, pre_ind].values,
                            index=expr_max.index, columns=match_pats)

# function to compute correlation matrix, graph laplacian
def data_to_graph_laplacian(data, corr_method = 'pearson'):
    
    # compute correlation matrix
    print(" Compute correlation matrix using "+corr_method)
    corr_mat = data.T.corr(method=corr_method).values

    # create weighted graph using the inverse square-root of absolute correlation as "distance" metric
    corr_g = adj_max._with_data( (1+corr_mat[lv,rv])/2 ) # shift correlation up 1 and divide by 2 as similarity metric
    wg = adj_max._with_data( 1 / corr_g.data ) # use inverse of shifted correlation as distance

    # shortest path distance via Djikstra's algorithm
    dist_mat = cg.shortest_path(wg)

    # manually compute graph laplacian: L = I - K^-1 * A
    rowsum = np.sum(corr_g.A,axis=0) # weighted degree
    L = np.eye(nv) - np.matmul(np.diag(1/rowsum), corr_g.A)
    return corr_mat, wg, dist_mat, L

# compute correlation-based distance matrix and graph Laplacian matrix
print("\nCompute correlation distance and graph Laplacian...")
t1 = time.time()
corr_diff, wg_diff, dist_diff, L_diff = data_to_graph_laplacian(data_diff, corr_method)
print(" Time for computation:", time.time()-t1)


#%%
# determine threshold tau as lowest when 99%ile of ORC exceeds 0.75
thresh = 0.75
pctile = 0.99
tau_crit_diff = orc_diff.quantile(pctile,axis=0).gt(thresh).idxmax()
print("Threshold tau:", tau[tau.tau_id == tau_crit_diff].value.values)

# extract curvature at critical threshold
ti = np.where(tau.tau_id == tau_crit_diff)[0].item()
orc_crit = orc_diff.iloc[:,ti].values.squeeze()

# integrate orc up to orc_crit (smoothed estimate)
dtau = np.diff(tau.value)
orc_int = orc_diff.iloc[:,:ti].multiply(dtau[:ti],axis=1).sum(axis=1).values


# %%
# Fig 1A:
# plot lines of each edge curvature over diffusion
plot_df = orc_diff.reset_index(drop=True).reset_index(names="edge_id") \
    .melt(id_vars="edge_id", var_name="tau_id", value_name="orc") \
        .merge(tau, on="tau_id").rename(columns={"value":"tau"})

# add orc_crit of each edge
orc_crit_df = pd.DataFrame({'orc_crit': orc_crit,
                            'orc_int': orc_int}).reset_index(names="edge_id")
plot_df = plot_df.merge(orc_crit_df, on="edge_id")

# plot all edges as lines, colored by threshold ORC
sns.set_style("ticks")
fig, ax = plt.subplots(figsize=(10,8))
sns.lineplot(ax=ax,
    data = plot_df,
    x="tau", y="orc", alpha=0.4,
    hue="orc_int", estimator = None,
    palette = "coolwarm_r", # RdBu coolwarm_r turbo_r
    legend=False)
ax.set_xscale('log') # log scale

# mark line at threshold tau
ax.axvline(tau.value[ti], color="black", linestyle=":")

# add colorbar
norm = plt.Normalize(orc_int.min(),
                     orc_int.max())
sm = plt.cm.ScalarMappable(cmap="coolwarm_r", norm=norm)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('$\\bar{\kappa}_{crit}$')

ax.set_xlabel('Scale Parameter $\\tau$')
ax.set_ylabel('Ollivier-Ricci curvature $\kappa$')
plt.show()

fig.savefig(os.path.join(results_dir,"kappa_tau_lineplot_int.png"), bbox_inches='tight')

#%%
# Table 1: top highest and lowest edges
print("Highest")
for e in np.argsort(-orc_int)[:5]:
    print(geneNames_max[lv[e]],
          geneNames_max[rv[e]],
          orc_crit[e],
          orc_int[e],
          corr_diff[lv[e],rv[e]])
print("Lowest")
for e in np.argsort(orc_int)[:5]:
    print(geneNames_max[lv[e]],
          geneNames_max[rv[e]],
          orc_crit[e],
          orc_int[e],
          corr_diff[lv[e],rv[e]])


#%% weighted louvain clustering method
# weighted modularity based on ORC_int
# using networkx louvain implementation

# create networkx graph object
g = nx.from_edgelist(orc_diff.index)

# scale orc_int to [0,1] and set as edge attribute
orc_int_scaled = (orc_int - min(orc_int)) / (max(orc_int) - min(orc_int))
attr_dict = dict(zip(*(orc_diff.index, orc_int_scaled.tolist())))
nx.set_edge_attributes(g, attr_dict, "orc_int_scaled")

# louvain clustering with weighted modularity
louvain_partition = nx.community.louvain_communities(g,
                                                     weight='orc_int_scaled',
                                                     resolution=1,
                                                     seed=0)
clusters_lv = np.zeros(nv).astype(int)
for ic in range(len(louvain_partition)):
    cluster_indices = np.array(list(louvain_partition[ic]))
    clusters_lv[cluster_indices] = ic

nclust = len(np.unique(clusters_lv))
print("Louvain clustering identified %i modules" % nclust)
print(np.bincount(clusters_lv))

## save cluster list
# include 'module' as 'cluster' + 1 for index to start at 1
pd.DataFrame({'Gene': geneNames_max,
              'cluster': clusters_lv,
              'module': clusters_lv+1}) \
            .sort_values('cluster') \
            .to_csv(os.path.join(results_dir,'louvain_weighted_clusters.csv'), index=False)

# load previous computed clusters
cluster_data = pd.read_csv(os.path.join(results_dir,'louvain_weighted_clusters.csv'), index_col=0)
clusters_lv = cluster_data.loc[geneNames_max].cluster

#%%
# Fig 1B:
# heatmap of correlation matrix
import PyComplexHeatmap as pch

plot_matrix = pd.DataFrame(
    corr_diff,
    index=geneNames_max,
    columns=geneNames_max)

# annotations for gene clusters
gene_df = pd.DataFrame({'Gene': geneNames_max,
                        'cluster_index': clusters_lv,
                        'cluster_name': np.char.add('module',(clusters_lv+1).astype(str))}, # add 1 to module label to start at 1 instead of 0
                        index=geneNames_max)
col_anno = pch.HeatmapAnnotation(cluster=pch.anno_simple(gene_df.cluster_name, colors={'module'+str(i+1):plt.cm.Set1(i) for i in range(nclust)},
                                                         height=5, add_text=True, text_kws={'rotation':0,'fontsize':10,'color':'black'},
                                                         legend=False))
row_anno = pch.HeatmapAnnotation(axis=0,
                                 cluster=pch.anno_simple(gene_df.cluster_name, colors={'module'+str(i+1):plt.cm.Set1(i) for i in range(nclust)},
                                                         height=5, add_text=True, text_kws={'rotation':90,'fontsize':10,'color':'black'},
                                                         legend=False))

# create heatmap with gene clusters, no gene names
fig = plt.figure(figsize=(10,9))
hm = pch.ClusterMapPlotter(data=plot_matrix, label = 'Pearson Corr.',
                           top_annotation=col_anno,
                           col_cluster=True, col_split = gene_df.cluster_index, col_split_order = list(range(nclust)),
                           left_annotation=row_anno, #row_dendrogram=False,
                           row_cluster=True, row_split = gene_df.cluster_index, row_split_order = list(range(nclust)),
                           show_rownames=False, show_colnames=False,
                           cmap='coolwarm', vmin=-1, vmax=1,
                           )
plt.show()
fig.savefig(os.path.join(results_dir,"correlation_"+corr_method+"_heatmap_lv.png"), bbox_inches='tight', dpi=600)


# %% partition layout functions
# partitional layout functions
# Source: https://stackoverflow.com/questions/65069624/networkx-cluster-nodes-in-a-circular-formation-based-on-node-color

NODE_LAYOUT = nx.spring_layout # circular_layout, spring_layout
COMMUNITY_LAYOUT = nx.circular_layout
def partition_layout(g, partition, ratio=0.3):
    """
    Compute the layout for a modular graph.

    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance
        network to plot

    partition -- dict mapping node -> community or None
        Network partition, i.e. a mapping from node ID to a group ID.

    ratio: 0 < float < 1.
        Controls how tightly the nodes are clustered around their partition centroid.
        If 0, all nodes of a partition are at the centroid position.
        if 1, nodes are positioned independently of their partition centroid.

    Returns:
    --------
    pos -- dict mapping int node -> (float x, float y)
        node positions

    """
    pos_communities = _position_communities(g, partition)

    pos_nodes = _position_nodes(g, partition)
    pos_nodes = {k : ratio * v for k, v in pos_nodes.items()}

    # combine positions
    pos = dict()
    for node in g.nodes():
        pos[node] = pos_communities[node] + pos_nodes[node]

    return pos


def _position_communities(g, partition, **kwargs):
    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(g, partition)

    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = COMMUNITY_LAYOUT(hypergraph, **kwargs)

    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]

    return pos


def _find_between_community_edges(g, partition):
    edges = dict()
    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]

        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]

    return edges


def _position_nodes(g, partition, **kwargs):
    # Positions nodes within communities.
    communities = dict()
    for node, community in partition.items():
        if community in communities:
            communities[community] += [node]
        else:
            communities[community] = [node]

    pos = dict()
    for community, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = NODE_LAYOUT(subgraph, **kwargs)
        pos.update(pos_subgraph)

    return pos

#%%
# Fig 1C: network colored and grouped by cluster

# create networkx graph object
g = nx.from_edgelist(orc_diff.index)

# set edge attribute for orc_int (integral-smoothed ORC at threshold tau)
attr_dict = dict(zip(*(orc_diff.index, orc_int.tolist())))
nx.set_edge_attributes(g, attr_dict, "orc_int")

# set node attribute for cluster label
nx.set_node_attributes(g, louvain_partition, "cluster_lv")

# plot network with clusters
fig = plt.figure(figsize=(10,8))

# compute partition layout (functions above)
layout = partition_layout(g, {n:clusters_lv[i] for i,n in enumerate(g.nodes())}, 0.2)

# draw nodes separately from edges to color properly
nx.draw_networkx_nodes(g, pos=layout,
                       node_color=plt.cm.Set1(clusters_lv),
                       node_size=10)
nx.draw_networkx_edges(g, pos=layout, alpha=0.3,
                       edge_color = [e[2] for e in g.edges(data="orc_int")],
                       edge_cmap = plt.cm.coolwarm_r,
                       edge_vmin=-1.2, edge_vmax=1.2 # for orc_int
                       )

# add ghost legend
for v in set(clusters_lv):
    plt.scatter([],[], c=[plt.cm.Set1(v)], label='module{}'.format(v+1)) # add 1 to module label
leg = plt.legend(loc = 'upper right', frameon=True)
frame = leg.get_frame()
frame.set_color('white')

# add colorbar
norm = plt.Normalize(-1.2, 1.2)
sm = plt.cm.ScalarMappable(cmap="coolwarm_r", norm=norm)
cbar = plt.colorbar(sm, ax=plt.gca())
cbar.set_label('$\\bar{\kappa}_{crit}$')
plt.axis('off')
#plt.tight_layout()
plt.show()

fig.savefig(os.path.join(results_dir,"clustergraph_louvain.png"), bbox_inches='tight', dpi=600)

#%%
# compute statistics:
# intra-cluster, inter-cluster average ORC_crit

# cluster coincidence matrix
nclust = len(np.unique(clusters_lv))
lv0,rv0 = zip(*list(orc_diff.index)) # get left vertices (rows) and right vertices (columns) from data MultiIndex

# for each cluster pair:
# get edge indices, counts, and average cross-cluster orc_crit
cross_cluster_inds = np.zeros((nclust,nclust), dtype=object)
cross_cluster_cts = np.zeros((nclust,nclust), dtype=int)
coincidence_matrix = np.zeros((nclust,nclust))
for ci in range(nclust):
    for cj in range(nclust):
        vi = np.where(clusters_lv==ci)[0]
        vj = np.where(clusters_lv==cj)[0]
        edges_f = np.where(np.isin(lv0, vi) & np.isin(rv0, vj))[0] # forward
        edges_r = np.where(np.isin(lv0, vj) & np.isin(rv0, vi))[0] # and backward connections
        edges = np.union1d(edges_f,edges_r)

        # store in matrix
        cross_cluster_inds[ci,cj] = edges
        cross_cluster_cts[ci,cj] = len(edges)
        coincidence_matrix[ci,cj] = np.mean(orc_int[edges])
#sns.heatmap(coincidence_matrix)

# Fig 1D:
# dotplot of cross-cluster edge curvature
from matplotlib.collections import PatchCollection
x, y = np.meshgrid(np.arange(nclust), np.arange(nclust))
s = np.sqrt(cross_cluster_cts) # radius proportional to sqrt(N) so area is proportional to N
c = coincidence_matrix

sns.set_style("white")
fig, ax = plt.subplots(figsize=(8,6.5))

R = s/s.max()/2
circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
col = PatchCollection(circles, array=c.flatten(), cmap="coolwarm_r")
col.set_clim([-0.6, 0.6])
ax.add_collection(col)

ax.set(xticks=np.arange(nclust), yticks=np.arange(nclust))
ax.set_xticks(np.arange(nclust+1)-0.5, minor=True)
ax.set_yticks(np.arange(nclust+1)-0.5, minor=True)
ax.set_xticklabels(np.arange(nclust)+1) # add 1 to module label
ax.set_yticklabels(np.arange(nclust)+1) # add 1 to module label
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, left=True)
ax.invert_yaxis()
ax.grid(which='minor')

cbar = fig.colorbar(col)
cbar.set_label('$\\bar{\kappa}_{crit}$')
plt.show()

fig.savefig(os.path.join(results_dir,"cross_cluster_dotplot_lv.pdf"), bbox_inches='tight')

#%%
# compare all within and between cluster edges
wi_edges = np.array([],dtype=int)
bw_edges = np.array([],dtype=int)
for ci in range(nclust):
    for cj in range(nclust):
        if ci==cj:
            wi_edges = np.union1d(wi_edges,cross_cluster_inds[ci,cj])
        else:
            bw_edges = np.union1d(bw_edges,cross_cluster_inds[ci,cj])

print( "Number edges: within=%i, between=%i" % (len(wi_edges), len(bw_edges)) )
print("Within cluster avg k_int: %0.4f\nBetween cluster avg k_crit: %0.4f" %
    (np.mean(orc_int[wi_edges]), np.mean(orc_int[bw_edges])) )

print("\nT-test of k_crit")
res = scipy.stats.ttest_ind(orc_crit[wi_edges], orc_crit[bw_edges])
print(res)
print("T-test of k_int")
res = scipy.stats.ttest_ind(orc_int[wi_edges], orc_int[bw_edges])
print(res)
