#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import networkx as nx
from tqdm import tqdm
import scipy.sparse as sparse
import pickle

import warnings
from modules.utils import *
from modules.plotting import make_disc_plots

import multiprocessing
from joblib import Parallel, delayed
#
# import warnings
# from scipy.sparse import SparseEfficiencyWarning
# warnings.simplefilter('ignore',SparseEfficiencyWarning)

def unwrap_self(arg, **kwarg):
    return Network.compute_control(*arg, **kwarg)

class Network:

    def __init__(self, A, nodes=None, method='T'):
        assert isinstance(A, sparse.csr_matrix), "adjacency matrix must be scipy sparse matrix"
        self.N = A.shape[0]
        if self.N > 1: # avoid unnecessary calculations
            self.A = A
            self.set_control_matrix(A, method=method)
            self.G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
        self.nodes = nodes
        if self.nodes is None:
            self.nodes = np.arange(self.N)

    def set_control_matrix(self, A, method='W'):
        self.method = method
        self.C = A

    def get_descendants(self, i):
        return np.array(list(nx.descendants(self.G, i)))

    def get_Bsub(self, desc):
        return self.C[desc,:][:, desc]

    def get_d(self, i, desc):
        return self.C[i, desc].reshape(1, -1) # row vector

    def get_Isub(self, desc):
        return sparse.eye(len(desc), format='csr')

    def get_dtilde(self, i, desc):
        Ndesc = len(desc)

        if Ndesc == 0:
            return None

        Bsub = self.get_Bsub(desc)

        d = self.get_d(i, desc)

        Isub = self.get_Isub(desc)
        to_invert = Isub - Bsub

        to_invert = to_invert.tocsc()
        inverted = sparse.linalg.inv(to_invert)
        #print(type(to_invert), type(inverted))
        # to_invert = to_invert.todense()
        #inverted = np.linalg.inv(to_invert)
        #inverted = sparse.csr_matrix(inverted)

        if isinstance(inverted, np.ndarray):
            warnings.warn("inverted matrix is not sparse, converting projection to sparse", RuntimeWarning)
            d = d.todense()
            dtilde = d.dot(inverted)
            dtilde = sparse.csr_matrix(dtilde)
        else:
            inverted = inverted.tocsr()
            dtilde = d.dot(inverted)


        return dtilde

    def compute_control(self, i):
        if self.N == 1:
            return np.array([]), None
        try:
            desc = self.get_descendants(i)
            dtilde = self.get_dtilde(i, desc)
            return desc, dtilde
        except Exception as e:
            print(e)
            print("Possibly singular value found")
            return desc, None

    def compute_control_mat(self, n_jobs=-1):
        if self.N == 1:
            return None

        list_of_comps = Parallel(n_jobs=n_jobs, backend="threading")(delayed(unwrap_self)(i) for i in zip([self]*self.N, range(self.N)) )

        W = sparse.csr_matrix((self.N, self.N)).tolil()

        fail = False
        for i, (desc, dtilde) in enumerate(list_of_comps):
            #desc, dtilde = self.compute_control(i)
            if len(desc) == 0:
                continue
            if len(desc) > 0 and dtilde is None:
                return None # inversion error
            W[i,desc] = dtilde
        return W.tocsr()

    def connected_components(self):
        for nodes in nx.connected_components(self.G.to_undirected()):
            nodes = np.array(list(nodes))
            Asub = self.A[nodes,:][:,nodes]
            Asub = Asub.copy()
            yield Network(Asub, nodes=self.nodes[nodes], method=self.method)


def decompose_network(A, analyse=True, dump_disc=True):
    assert isinstance(A, sparse.csr_matrix), "adj should be sparse"

    print("Removing selfloops or empty elements in adj matrix if present")
    A = remove_selfloops(A)
    N = A.shape[0]
    Ndisc = A.shape[0] # just an estimate


    if analyse:
        print("Analyzing graph")
        Ndisc, hist = get_num_disconnected(A)
        print("Found {} disc. comp.".format(Ndisc))
        make_disc_plots(hist, Nmax=None, xlog=False)
        make_disc_plots(hist, Nmax=50, xlog=False)
        make_disc_plots(hist, Nmax=None, xlog=True)
        make_disc_plots(hist, Nmax=1000, xlog=True)

    g = Network(A)
    if Ndisc < 2:
        print("Single component found")
        conn_comps = [g]
    else:
        print("Getting disconnected components (upper bound estimation = number of nodes)")
        conn_comps = [gsub for gsub in tqdm(g.connected_components(), total=Ndisc)]
        print("Finished decomposing:", len(conn_comps))

    if dump_disc:
        print("Dumping disconnected components to pkl")
        with open("dump_disc.pkl", "wb") as f:
            pickle.dump(conn_comps, f)

    return conn_comps


def compute_control_vitali(A, n_jobs=-1, verbose=2, decompose=False, **kwargs):

    assert isinstance(A, sparse.csr_matrix), "adj should be sparse"

    def compute_sub(gsub):
        nodes = gsub.nodes
        Csub = None # if a single node: just None as C
        if len(nodes) > 1:
            Csub = gsub.compute_control_mat(n_jobs=n_jobs)
        return nodes, Csub

    if decompose:
        conn_comps = decompose_network(A, **kwargs)

        controls = Parallel(n_jobs=n_jobs, backend="threading", verbose=verbose)(delayed(compute_sub)(gsub) for gsub in  conn_comps)
        return controls
    else:
        g = Network(A)
        control = compute_sub(g)
        return control

