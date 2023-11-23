# analysis of acral melanoma pre-treatment vs on-treatment with nivolumab
# Dataset GEO accession code: GSE91061
# 109 RNASeq samples (58 On-treatment and 51 Pre-treatment) from 65 patients

# analysis code written by Kevin A Murgas
# Stony Brook University, Dept. of Biomedical Informatics
# advisor: Allen Tannenbaum PhD
# in collaboration with: Rena Elkin PhD, Nadeem Riaz PhD, Emil Saucan PhD, Joseph Deasy PhD

### HELPER FUNCTIONS ###
# function to compute data-correlation-based graph
# returns distance matrix and graph Laplacian
def data_to_graph_laplacian(data, corr_method='pearson'):
    # compute correlation matrix
    print(" Compute correlation matrix using "+corr_method)
    corr_mat = data.T.corr(method=corr_method).values

    # create weighted graph using the inverse square-root of absolute correlation as "distance" metric
    corr_g = adj_max._with_data( (1+corr_mat[lv,rv])/2 ) # shift correlation up 1 and divide by 2 as similarity metric
    wg = adj_max._with_data( 1/corr_g.data ) # use inverse of shifted correlation as distance

    # shortest path distance via Djikstra's algorithm
    dist_mat = cg.shortest_path(wg)

    # manually compute graph laplacian: L = I - K^-1 * A
    rowsum = np.sum(corr_g.A,axis=0) # weighted degree
    L = np.eye(nv) - np.matmul(np.diag(1/rowsum), corr_g.A)
    return dist_mat, L

# function given graph Laplacian matrix and a vector of tau to try
# returns a list of diffused distributions
def diffuse_for_tau(L, tau_range):
    diffused_list = []
    t0 = time.time()
    for ti in range( len(tau_range) ):
        tau_i = tau_range[ti]
        
        # scipy.linalg computes a matrix of diffused distributions starting from each vertex
        # resulting matrix: each row is a marginal for the diffusion of vertex i
        # i.e. p_i = diffused_i[i,:] being a diffusion from vertex i for time tau_i
        diffused_i = la.expm(-tau_i * L)
        diffused_list.append( np.ascontiguousarray(diffused_i) )
    print(" Done. Time for all diffusions:", time.time()-t0)
    return diffused_list

# function to compute Wasserstein distance between
# distributions of diffusion processes at nodes u and v
# precomputed, only need tau index
def orc_edge(ie, ti):
    u = lv[ie]; v = rv[ie]

    # pull from list of diffusions across all nodes
    m_u = diffusion_list[ti][u,:]
    m_v = diffusion_list[ti][v,:]

    # compute wasserstein distance between u and v marginals
    w = emd2(m_u, m_v, d0)
    return w

# function to compute ORC over all graph edges for each tau
# uses parallel processing pool to accelerate computation
# suggest running on high performance cluster
def compute_orc_all():
    edge_w_all = []
    edge_orc_all = []
    t_all = time.time()
    for ti in range( len(tau_range) ):
        print("Tau = ", tau_range[ti], " (", ti, "/", len(tau_range), ")")
        
        t0 = time.time()
        with mp.Pool(n_proc) as pool:
            res = pool.starmap(orc_edge, [(ie, ti) for ie in range(ne)])
        edge_w = list(res)
        edge_orc = 1 - ( edge_w / d0[lv,rv] )
        print(" Graph done. Time for all edges:", time.time()-t0)

        edge_w_all.append(edge_w)
        edge_orc_all.append(edge_orc)
    print("Complete. Time for all tau:", time.time()-t_all)
    return edge_w_all, edge_orc_all


### MAIN SCRIPT ###
if __name__ == "__main__":
    import sys, os, time, scipy.io
    import pandas as pd, numpy as np, scipy.sparse as sp
    from scipy.sparse import csgraph as cg, linalg as la
    from ot import emd2
    import multiprocessing as mp
    n_proc = mp.cpu_count()
    mp.set_start_method('fork')

    # get arguments from sys
    try:
        expr_file = sys.argv[1]
        ppi_file = sys.argv[2]
        genelist_file = sys.argv[3]
        save_dir = sys.argv[4]
        print("Recieved input arguments")
    except IndexError:
        # use defaults
        print("Too few input arguments, setting default")
        expr_file = 'melanoma_immunotherapy/GSE91061_rld_qnorm_HGNC.csv'
        ppi_file = 'Topology_python/topo_stringdb11sparse.mat'
        genelist_file = 'melanoma_immunotherapy/WP_CANCER_IMMUNOTHERAPY_BY_PD1_BLOCKADE.v2023.1.Hs.grp'
        save_dir = 'Results/melanoma_qnorm_pearsonshiftinverse_PD1Bnbrs'

    # load data
    print("Loading input files:")
    print(" Expression =", expr_file)
    print(" PPI =", ppi_file)
    print(" GeneList =", genelist_file)

    # load expression data
    print("\nLoad in gene expression data from melanoma immunotherapy study GSE91061...")
    expr_data = pd.read_csv( expr_file, index_col=0)
    expr_hgnc = expr_data.index.values
    print(" Expression data loaded. %i genes x %i samples." % expr_data.shape )

    # determine pre- and on- treatment groups from column labels
    ids_data = expr_data.columns.values
    pat_id, condition = [np.array(e) for e in zip(*[id.split("_")[0:2] for id in ids_data])]

    ### CONSTRUCT MELANOMA GENE NETWORK TOPOLOGY ###
    # load specified gene list to subset
    print("Load in gene list to subset expression data...")
    genelist = pd.read_csv(genelist_file, skiprows=2, names=["Gene"])
    genelist_hgnc = genelist['Gene'].values.astype(str)
    print(" Gene list loaded. %i unique genes." % len(np.unique(genelist_hgnc)) )

    # load STRINGdb protein-protein interaction network (sparsified)
    print("Load in PPI topology (file: %s)" % ppi_file)
    ppi_data = scipy.io.loadmat(ppi_file, simplify_cells = True)
    adj = sp.csr_matrix(ppi_data['adj'].astype('bool'))
    topoNames = ppi_data['topoNames'] # need to convert string array to char array in order for python to read
    topoNames = [s.strip() for s in topoNames] # remove white space
    print(" PPI loaded. %i vertices, %i edges" % (adj.shape[0] , adj.nnz))

    # find matching genes in PPI topology genes + expression data + subset gene list
    matching_genes = list(set(expr_hgnc) & set(genelist_hgnc) & set(topoNames))
    matching_genes.sort()
    print("\nMatched %i genes in expression, genelist, and PPI." % len(matching_genes) )

    # optional flag to include all 1-hop neighbors of subset genes
    include_neighbors = True
    if (include_neighbors == True):
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
    ind_exp = [list(expr_hgnc).index(m) for  m in matching_genes]
    ind_ppi = [topoNames.index(m) for  m in matching_genes]
    
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


    ### CORRELATION NETWORK ANALYSIS ###
    # construct correlation network based on matched expression changes (on-treatment minus pre-treatment)

    # find matched patient samples
    # compute difference in expression in matched patients
    # matched_pats[on] - match_pats[pre]
    match_pats = np.intersect1d(pat_id[condition == 'Pre'],
                                pat_id[condition == 'On'])
    pre_on_ind = [[np.where((pat_id == p) & (condition == 'Pre') )[0].item(),
                np.where((pat_id == p) & (condition == 'On') )[0].item()] for p in match_pats]
    pre_ind, on_ind = [np.array(e) for e in zip(*pre_on_ind)]
    data_diff = pd.DataFrame(expr_max.iloc[:, on_ind].values - expr_max.iloc[:, pre_ind].values,
                             index=expr_max.index, columns=match_pats)
    print("\nComputed difference of on-pre (n=%i pts)" % (data_diff.shape[1]))
    
    # compute correlation-based distance matrix and graph Laplacian matrix
    print("\nCompute correlation distance and graph Laplacian...")
    dist_diff, L_diff = data_to_graph_laplacian(data_diff)

    # compute diffusion of the graph Laplacian
    # at multiple points of scale parameter tau (time-step)
    print("\nCreate diffusion processes...")
    tau_range = np.logspace(-2, 2, 101, base=10) # here we select 21 steps on a log scale from 10^-2 to 10^2
    diffusion_diff = diffuse_for_tau(L_diff, tau_range)
    

    ### DYNAMIC CURVATURE ANALYSIS ###
    # for each tau:
    # compare Wasserstein distance of diffusion distributions
    # marginals == diffusion process around nodes i,j
    # optimal transport of distributon over ALL vertices (ALL distance pairs)
    # compute ORC over all edges in adj_max

    # check if adjacency matrix is symmetric
    # if so, take upper triangle edges (so as to not compute same edges twice)
    if (abs(adj_max-adj_max.T)>1e-10).nnz == 0:
        print("Symmetric adjacency matrix detected. Reducing edgelist...")
        lv,rv = sp.triu(adj_max).nonzero()
        ne = int(ne/2)
    else:
        print("Asymmetric adjacency matrix detected. Continuing...")

    print("\n-----")
    print("Analyze curvature in diffused network...")
    diffusion_list = diffusion_diff
    d0 = dist_diff
    edge_w_diff, edge_orc_diff = compute_orc_all()


    ### SAVE RESULTS ###

    print("\nSave with pandas...")

    # check if path exists, if not create a new directory
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # # Save gene names of vertices
    save_file = "genes.csv"
    save_data = pd.DataFrame({'Gene': geneNames_max})
    save_data.to_csv(os.path.join(save_dir, save_file), index=None)

    # save tau and names for columns of results
    tau_names = ['tau' + str(i) for i in range(len(tau_range))]
    save_file = "tau.csv"
    save_data = pd.DataFrame({'tau_id': tau_names, 'value': tau_range})
    save_data.to_csv(os.path.join(save_dir, save_file), index=None)

    # Save edge Ollivier-Ricci curvature values, with edge vertices as index
    save_file = "orc_diff.csv"
    save_data = pd.DataFrame(np.array(edge_orc_diff).T, columns=tau_names).astype('single')
    save_data.insert(0,'lv',lv); save_data.insert(1,'rv',rv)
    save_data.to_csv(os.path.join(save_dir, save_file), index=None)

    print("Script complete.")
