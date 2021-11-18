import numpy as np
import scipy as sp
import itertools
import copy
from sympy import Matrix
# import operator as op


from .utils import n_zeros, n_nonzeros, get_idx_nonzeros

#https://nl.mathworks.com/matlabcentral/fileexchange/49357-exact-controllability-of-complex-networks

def max_multipl(l, tol=1e-8):
	# return maximum multiplicity and number of counts
	# np.unique does not allow tolerance, so do it manually
	max_m = 1
	max_x = l[0]
	for x in l:
		# number of similar items, also itself incl.
		m = n_zeros(l-x, tol=tol)
		if m > max_m:
			max_m, max_x = m, x
	return max_x, max_m # value, multiplicity

def calc_muM(A, lambdaM):
	N = min(A.shape)
	return N-np.linalg.matrix_rank(lambdaM*np.eye(N)-A)

def calc_echelon_form(M, tol=1e-8):
	# bring to canonical form (row canonical function)
	can, _ = Matrix(M).rref()
	return np.array(can) # back to column form

def calc_independence(M, tol=1e-8):
	# Configurations of driver nodes would include all the dependent rows and
	# a node from each of the independent row if the multiplicity at row i is >1.
	E = calc_echelon_form(M.T, tol=tol).T
	li = [[] for _ in range(E.shape[0])]
	ld = []
	for i in range(E.shape[0]):
		if n_nonzeros(E[i,:]) == 1:
			idx = get_idx_nonzeros(E[i,:], tol=tol)[0]
			li[idx].append(i)
		else:
			ld.append(i)
	li = [l for l in li if len(l)>1] # larger than one to avoid 1 sets
	return li, ld

def calc_N_configs(li):
	if len(li) == 0:
		return 0
	N_configs = 1
	for idx in li:
		N_configs *= sp.special.comb(len(idx), len(idx)-1, exact=True)
	return N_configs
#
# def get_driver_node_permutations(li):
# 	N_configs = 1
# 	driver_configs = []
# 	for idx in li:
# 		if len(idx) <= 1: continue
# 		driver_configs.append(list(itertools.combinations(idx, len(idx)-1)))
# 		N_configs = N_configs*len(driver_configs[-1])
# 	if all(len(x)==0 for x in driver_configs):
# 		N_configs = 0
# 		driver_configs = []
# 	else:
# 		driver_configs = [x for x in driver_configs if len(x)>0] # reduce a bit, just to be sure
# 	return N_configs, driver_configs

# def generate_driver_node_set(driver_configs, ld):
# 	driver_nodes = copy.copy(ld)
# 	for x in driver_configs:
# 		if len(x) > 0:
# 			k = np.random.randint(len(x), size=1)
# 			driver_nodes += list(op.itemgetter(*k)(x))
# 	return driver_nodes

def yield_configs(l):
	if len(l) == 0:
		return l
	for x in itertools.combinations(l[0], len(l[0])-1):
		if len(l) == 1:
			yield list(x)
		else:
			for y in yield_configs(l[1:]):
				yield list(x) + list(y)

class DriverConfigGenerator:
	def __init__(self, A, tol=1e-8, verbose=False):
		self.A = A
		self.tol = tol
		self.verbose = verbose
		self.N = None
		self.lambdaM = None
		self.muM = None
		self.li = None
		self.ld = None
		self.N_configs = None
		self.driver_configs = None
		self.pre_compute(A)

	def pre_compute(self, A):
		self.N = min(A.shape)
		lambdas = np.linalg.eigvals(A)
		if self.verbose:
			print("lambdas: ", lambdas)
		self.lambdaM, self.muM = max_multipl(lambdas, tol=self.tol)
		# something to check still: was expecting this muM to be the same (unclear in the paper!!!)
		self.muM = calc_muM(A, self.lambdaM) # overwrite
		if self.verbose:
			print("(lambdaM, muM)=({},{})".format(self.lambdaM, self.muM))
		reg_mat = A - self.lambdaM*np.eye(self.N)
		canonical_M = calc_echelon_form(reg_mat)
		#if self.verbose:
		#	print("Canonical M.T")
		#	print(canonical_M.T)
		self.li, self.ld = calc_independence(canonical_M.T, tol=self.tol)
		if self.verbose:
			print("li = ", self.li)
			print("ld = ", self.ld)
		# self.N_configs, self.driver_configs = get_driver_node_permutations(self.li, self.N)
		self.N_configs = calc_N_configs(self.li)
		if self.verbose:
			print("N_configs:", self.N_configs)
			# print("driver_configs:", self.driver_configs)
	#
	# def get_driver_configs(self):
	# 	return self.driver_configs

	def get_min_driver_nodes(self):
		return self.muM

	def get_N_configs(self):
		return self.N_configs

	def generate_all(self):
		if len(self.li) > 0:
			for l in yield_configs(self.li):
				drivers = copy.copy(self.ld)
				drivers.extend(l)
				yield drivers
		else:
			yield self.ld

	def get_driver_node_set(self):
		drivers = next(iter(self.generate_all()))
		return drivers

	# def get_driver_node_set(self):
	# 	driver_nodes = generate_driver_node_set(self.driver_configs, self.ld, self.N)
	# 	if self.verbose:
	# 		print("driver sets:",driver_nodes)
	# 	return driver_nodes
