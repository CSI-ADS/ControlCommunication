import numpy as np
from collections import deque
from tqdm import tqdm
import scipy.sparse as sparse
import networkx as nx

from .utils import get_children, n_nonzeros, get_all_disconnected, is_connected

def depth_first_pair_control(start_node, end_node, path, cum_control, Cnode, A):
    for child in get_children(A, end_node):
        if child in path:
            continue # no cycles
        # go deeper
        weight = A[end_node, child]
        control = cum_control[-1]*weight if len(cum_control)>0 else weight
        Cnode[child] += control
        cum_control.append(control)
        path.add(child)
        # print("FROM TO :", start_node+1, child+1, [p+1 for p in path], control, Cnode[child])
        depth_first_pair_control(start_node, child, path, cum_control, Cnode, A)
        cum_control.pop() # remove last element
        path.remove(child)
        # print("***")

# a wrapping function
def compute_pair_control(A, check_disonnected=True, verbose=False):
    assert isinstance(A, sparse.csr_matrix), "adj should be sparse"
    if check_disonnected:
        if not is_connected(A):
            if verbose:
                print("Found disconnected components, splitting into disconnected subgraphs")
            Cdict = {}
            for nodes, Asub in get_all_disconnected(A, verbose=verbose):
                Cdict[nodes] = compute_pair_control_(Asub, verbose=False)
            return Cdict
        else:
            if verbose:
                print("No disconnected components, proceeding as usual")
    if verbose:
        print("Computing pair control for full graph")
    Cdict = { # all together
        tuple(np.arange(A.shape[0])) : compute_pair_control_(A, verbose=verbose)
        }
    return Cdict

def compute_pair_control_(A, verbose=False):
    assert A.shape[0] == A.shape[1], "adj must be square adjacency matrix for now"
    assert isinstance(A, sparse.csr_matrix), "adj must be sparse"
    if A.shape[0] == 1: # shortcut!
        return sparse.csr_matrix(np.zeros((1,1)))
    if verbose:
        print("Converting to dense matrix (due to implementation)")
    A = np.asarray(A.todense()) # for now, converting anyways
    N = A.shape[0]
    C = np.zeros((N,N))
    if verbose:
        print("Looping over start nodes")
    for node in tqdm(range(N), disable=not verbose):
        cum_control = deque()
        path = {node}
        depth_first_pair_control(node, node, path, cum_control, C[node, :], A)
    assert n_nonzeros(np.diagonal(C)) == 0, "something went wrong, found elements on diagonal"
    return sparse.csr_matrix(C)
