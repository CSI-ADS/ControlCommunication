import numpy as np
import copy
import networkx as nx
import scipy.sparse as sparse
from tqdm import tqdm

def is_zero(l, tol=1e-8):
	return np.abs(l) < tol

def all_zero(l, tol=1e-8):
	return int(np.sum(~is_zero(l, tol=tol))) == 0

def n_zeros(l, tol=1e-8):
	return int(np.sum(is_zero(l, tol=tol)))

def n_nonzeros(l, tol=1e-8):
	return int(np.sum(~is_zero(l, tol=tol)))

def get_nonzeros(l, tol=1e-8):
	return l[~is_zero(l, tol=tol)]

def get_idx_nonzeros(l, tol=1e-8):
	return np.where(~is_zero(l, tol=tol))[0]

def get_unweighted_adj(A, tol=1e-8):
    return A > tol

def get_root_nodes(A, tol=1e-8):
    A_conn = get_unweighted_adj(A, tol=tol)
    root_nodes = A_conn.any(axis=0).flatten()
    root_nodes = np.asarray(root_nodes).flatten() # weird behavior of numpy
    root_nodes = np.where(~root_nodes)[0]
    return root_nodes

def get_leaf_nodes(A, tol=1e-8):
    A_conn = get_unweighted_adj(A, tol=tol)
    leaf_nodes = A_conn.any(axis=1).flatten()
    leaf_nodes = np.asarray(leaf_nodes).flatten() # weird behavior of numpy
    leaf_nodes = np.where(~leaf_nodes)[0]
    return leaf_nodes

def get_parents(A, node, tol=1e-8):
    return get_idx_nonzeros(A[:, node], tol=tol)

def get_children(A, node, tol=1e-8):
    return get_idx_nonzeros(A[node,:], tol=tol)

def connected_component_subgraph_adj(G):
    assert isinstance(G, nx.Graph), "should be a graph"
    for nodes in nx.connected_components(G.to_undirected()):
        yield tuple(nodes), copy.deepcopy(nx.adjacency_matrix(G.subgraph(nodes)))

def get_all_disconnected(A, verbose=False):
    assert isinstance(A, sparse.csr_matrix), "adj should be sparse"
    G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
    Adict = {}
    for nodes, Asub in tqdm(connected_component_subgraph_adj(G), total=A.shape[0], disable=not verbose):
        yield nodes, Asub

def is_connected(A):
    assert isinstance(A, sparse.csr_matrix), "adj should be sparse"
    G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
    return nx.is_connected(G.to_undirected())
