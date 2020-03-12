import numpy as np
import networkx as nx


cat_cortex = np.loadtxt("cat-cortex-adj.csv", unpack="True")
G = nx.DiGraph(cat_cortex)

ent = dict(G.in_degree(G.nodes()))

sai = dict(G.out_degree(G.nodes()))


