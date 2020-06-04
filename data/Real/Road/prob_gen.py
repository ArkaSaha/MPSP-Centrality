from geopy.distance import distance
from statistics import mean
import sys
from os import popen
import re

# sys.argv[1] is the .osm file and sys.argv[2] is the .log file

f = open(sys.argv[1],'r')
coord = {}

node_pattern = re.compile(".*<node id=\"(\d*)\".* lat=\"(-?\d*.\d*)\" lon=\"(-?\d*.\d*)\".*")
for line in f:
    m = node_pattern.match(line)
    if m:
        id_no = int(m.group(1))
        latitude = float(m.group(2))
        longitude = float(m.group(3))
        coord[id_no] = (latitude, longitude)

f.close()

f = open(sys.argv[2],'r')
vertices = {}
edges = {}
max_speed = 0
num = 0

for line in f:
    nodes = [int(x) for x in line.strip().split()]
    speeds = [float(x) for x in f.readline().strip().split()]
    for (u, v), s in zip(zip(nodes[:-1], nodes[1:]), speeds):
        if u not in vertices:
            vertices[u] = num
            num = num + 1
        if v not in vertices:
            vertices[v] = num
            num = num + 1
        if (u,v) not in edges:
            edges[(u,v)] = []
        edges[(u,v)].append(s)
        if max_speed < s:
            max_speed = s

lens = {}
probs = {}
print("{} {}".format(len(vertices), len(edges)))
for u, v in edges:
    #print("{} {}".format(u, v))
    lens[(u,v)] = round(distance(coord[u], coord[v]).km * 1000)
    probs[(u,v)] = 0.5 if sum(edges[(u,v)]) == 0 else sum(edges[(u,v)]) / (max_speed * len(edges[(u,v)]))

m_p = max(probs.values())
for u, v in edges:
    print("{} {} {:d} {:.6f}".format(vertices[u], vertices[v], lens[(u,v)], probs[(u,v)] / m_p))

f.close()

