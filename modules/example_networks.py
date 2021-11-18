import networkx as nx
import pylab as plt

def bow_tie_vitali():
    # The Network of Global Corporate Control, Vitali et al.
    G = nx.DiGraph()

    G.add_node(0, pos=(0,0))
    G.add_node(1, pos=(1,0))
    G.add_node(2, pos=(2,-1))
    G.add_node(3, pos=(2, 1))
    G.add_node(4, pos=(3, 0))
    G.add_node(5, pos=(4, 0))

    G.add_edge(0, 1, weight=0.1)

    G.add_edge(1, 2, weight=0.5)
    G.add_edge(2, 1, weight=0.3)

    G.add_edge(1, 4, weight=0.2)
    G.add_edge(4, 1, weight=0.3)

    G.add_edge(1, 3, weight=0.5)
    G.add_edge(3, 1, weight=0.3)

    G.add_edge(2, 4, weight=0.2)
    G.add_edge(4, 2, weight=0.5)

    G.add_edge(3, 4, weight=0.6)
    G.add_edge(4, 3, weight=0.5)

    G.add_edge(4, 5, weight=1.0)
    return nx.adjacency_matrix(G), G

if __name__ == '__main__':
    A, G = bow_tie_vitali()
    #G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph())
    pos = nx.get_node_attributes(G, 'pos')
    plt.figure()
    nx.draw(G,pos,labels={node:node for node in G.nodes()}, arrows=True, connectionstyle='arc3, rad = 0.3')
    edge_labels={(u,v):"({},{}):{}".format(u, v, d['weight'])
             for u,v,d in G.edges(data=True)}
    print(edge_labels)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels,font_color='red', alpha=0.2)
    plt.axis('off')
    # plt.show()
