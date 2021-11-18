from graph_tool.all import *
import collections
import networkx as nx
import numpy as np
import copy
from tqdm.notebook import trange, tqdm

def network_to_structured(graph):
    if type(graph) is graph_tool.Graph or type(graph) is graph_tool.GraphView:
        structured = {}
        for v in graph.get_vertices():
            neighbors = graph.get_out_neighbors(v)
            structured[v] = list(neighbors)
        #print('Converted to structured matrix')
    elif type(graph) is nx.classes.graph.Graph:
        structured = nx.to_dict_of_lists(graph)
    elif type(graph) is nx.classes.digraph.DiGraph:
        structured = nx.to_dict_of_lists(graph)
    elif type(graph) is dict:
        structured = graph
        print('Reading in graph as given dict, make sure it has the form {[from]:[to]} for all nodes.')
    else:
        print('Unsupported graph type.')
        return 0
    return structured

#extension on implementation by David Eppstein, UC Irvine, 27 Apr 2002
def bipartiteMatchHK(network, matching = {}, iterate=False):
    '''Find maximum cardinality matching of a bipartite graph (U,V,E).
    The input format is a dictionary mapping members of U to a list
    of their neighbors in V. For a directed network this is done 
    by making a copy of all nodes, one representing the indegree (i-) of a node
    and one for the outdegree (i+). A selfloop would for instance thus go from (i+)->(i-). 
    The output is a tuple (M,A) where M is a dictionary mapping members of i- 
    to their matches in i+, A are the unmatched nodes in i+.
    The same object may thus occur in both U and V, and is treated as two distinct vertices if this happens.'''
    if iterate == True:
        graph = network
    else:
        graph = network_to_structured(network)
    all_nodes = set(list(graph.keys())).union(set([item for sublist in list(graph.values()) for item in sublist]))
    isolated_nodes = []
    U_nodes = list(graph)
    for i in range(0, len(U_nodes)):
        if graph[U_nodes[i]] == [None]:
            isolated_nodes.append(U_nodes[i])
            del graph[U_nodes[i]]
    # initialize greedy matching (redundant, but faster than full search)
    if matching == {}:
        for u in graph:
            for v in graph[u]:
                if v not in matching:
                    matching[v] = u
                    break
        #print('Performed greedy matching')

    while 1:
        # structure residual graph into layers
        # pred[u] gives the neighbor in the previous layer for u in U
        # preds[v] gives a list of neighbors in the previous layer for v in V
        # unmatched gives a list of unmatched vertices in final layer of V,
        # and is also used as a flag value for pred[u] when u is in the first layer
        preds = {}
        unmatched = []
        pred = dict([(u,unmatched) for u in graph])
        for v in matching:
            if matching[v] in list(pred): #ZEER TIJDSKOSTELIJK
                del pred[matching[v]]

        layer = list(pred)
        
        # repeatedly extend layering structure by another pair of layers
        while layer and not unmatched:
            newLayer = {}
            for u in layer:
                if u in graph: #ZEER TIJDSKOSTELIJK
                    for v in graph[u]:
                        if v not in preds:
                            newLayer.setdefault(v,[]).append(u)
            if newLayer == {}:
                for i in list(pred):
                    if pred[i] == []:
                        isolated_nodes.append(i)

            layer = []
            for v in newLayer:
                preds[v] = newLayer[v]
                if v in matching:
                    layer.append(matching[v])
                    pred[matching[v]] = v
                else:
                    unmatched.append(v)
        # did we finish layering without finding any alternating paths?
        if not unmatched:
            unlayered = {}
            for u in graph:
                for v in graph[u]:
                    if v not in preds:
                        unlayered[v] = None
            results= matching
            matching = {}
            return results, all_nodes.symmetric_difference(set(list(results.keys())))
#return (matching,list(pred),list(unlayered))

        # recursively search backward through layers to find alternating paths
        # recursion returns true if found path, false otherwise
        def recurse(v):
            if v in preds:
                L = preds[v]
                del preds[v]
                for u in L:
                    if u in pred:
                        pu = pred[u]
                        del pred[u]
                        if pu is unmatched or recurse(pu):
                            matching[v] = u
                            return 1
            return 0

        for v in unmatched: recurse(v)
            
def run_over_matching_configurations(network,matching):
    #finds set of matchings by initializing bipartiteMatchHK once with every edge in network.
    Graph = network_to_structured(network)
    matching_edges_over_configurations = {}
    all_matchings = [matching]
    all_drivers = [get_drivers_from_matching(Graph,matching)]
    for i in range(0, len(list(Graph))):
        selected_u = list(Graph)[i]
        for j in range(0,len(Graph[selected_u])):
            selected_v = Graph[selected_u][j]
            if selected_v != None:
                initializing_edge = {selected_v:selected_u}
                new_match, unmatched = bipartiteMatchHK(Graph, initializing_edge, iterate=True)
                if len(list(new_match)) == len(list(matching)): #is net_match of minimal size?
                    if new_match!=matching: #did we already find this matching?
                        all_matchings.append(new_match)
                        all_drivers.append(get_drivers_from_matching(Graph,new_match))
                        for index in range(0, len(list(new_match))):
                            if list(new_match)[index] in list(matching_edges_over_configurations):
                                if new_match[list(new_match)[index]] not in matching_edges_over_configurations[list(new_match)[index]]:
                                    matching_edges_over_configurations.setdefault(list(new_match)[index],[]).append(new_match[list(new_match)[index]])
                            else:
                                matching_edges_over_configurations.setdefault(list(new_match)[index],[]).append(new_match[list(new_match)[index]])
    return matching_edges_over_configurations, all_matchings, all_drivers
                                           
def get_drivers_from_matching(graph,matching):
    all_nodes = set(list(graph.keys())).union(set([item for sublist in list(graph.values()) for item in sublist]))
    return all_nodes.symmetric_difference(set(list(matching.keys())))