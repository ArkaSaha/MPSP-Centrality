
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys

"""
    python3 cc.py [.graph file] [mpsp output file]
"""

def parse_line(x):
    return int(x[0]), int(x[1]), int(x[2]), float(x[3])

def parse_line2(x):
    return int(x[0]), float(x[1]), float(x[2])

with open(sys.argv[1], 'r') as f:
    n, m = [int(x) for x in f.readline().strip().split()]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in tqdm(range(m)):
        u, v, l, p = parse_line(f.readline().strip().split()) 
        G.add_edge(u, v)

    pos = {}
    for i in tqdm(range(n)):
        u, lat, lon = parse_line2(f.readline().strip().split())
        pos[u] = (lat, lon)



graph_pos = {k : (pos[k][1], pos[k][0]) for k in pos} # latitude and longitude are flipped wrt x/y axis
cc = sorted(nx.connected_components(G), key=len, reverse=True)
print(len(cc))

plt.figure(figsize=(16, 9))
# big component in blue
nx.draw_networkx_nodes(G, pos=graph_pos, nodelist=cc[0], node_size=1, node_color="b")
colors = "grcmykw"

ci = 0
for comp in cc[1:]:
    nx.draw_networkx_nodes(G, pos=graph_pos, nodelist=comp, node_size=1, node_color=colors[ci])
    ci = (ci + 1) % len(colors)
    
if False:
    plt.show()
    sys.exit(0)

with open(sys.argv[2], 'r') as f:
    data = [x.strip().split() for x in f.readlines()]
    i = 0
    for _ in range(8):
        while len(data[i]) == 0 or data[i][0] != "With": 
            i = i + 1
        i = i + 1
        s, t = [int(x) for x in data[i-2]]
        print("{} -> {}".format(s, t))
        pruning_path = []
        while data[i][0] != "Without":
            u, v, l, p = parse_line(data[i])
            pruning_path.append((u,v))
            i = i + 1
        while data[i][0] != "SP:":
            i = i + 1
        i = i + 1
        shortest_path = []
        while data[i][0] != "Number":
            u, v, l, p = parse_line(data[i])
            shortest_path.append((u,v))
            i = i + 1

        plt.figure(figsize=(16, 9))
        # big component in blue
        nx.draw_networkx_nodes(G, pos=graph_pos, nodelist=cc[0], node_size=1, node_color="b")
        colors = "grcmykw"

        ci = 0
        for comp in cc[1:]:
            nx.draw_networkx_nodes(G, pos=graph_pos, nodelist=comp, node_size=1, node_color=colors[ci])
            ci = (ci + 1) % len(colors)


        # draw paths
        nx.draw_networkx_edges(G, pos=graph_pos, edgelist=pruning_path, width = 10, edge_color='r')
        nx.draw_networkx_edges(G, pos=graph_pos, edgelist=shortest_path, width=6, edge_color='g')

        plt.show()








