from geopy.distance import distance
from statistics import mean
import sys
from os import popen
import re
import xml.etree.ElementTree as ET
from numpy.random import normal

# sys.argv[1] is the .osm file and sys.argv[2] is the .log file
assert sys.argv[3] in ['yes','no']

f = open(sys.argv[2],'r')
vertices = {}
edges = {}
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
f.close()
mean_speed = mean([max(l) for l in edges.values()])
max_speed = max([mean(l) for l in edges.values()])

coord = {}
if sys.argv[3] == 'no':
    f = open(sys.argv[1],'r')
    node_pattern = re.compile(".*<node id=\"(\d*)\".* lat=\"(-?\d*.\d*)\" lon=\"(-?\d*.\d*)\".*")
    for line in f:
        m = node_pattern.match(line)
        if m:
            id_no = int(m.group(1))
            latitude = float(m.group(2))
            longitude = float(m.group(3))
            coord[id_no] = (latitude, longitude)
    f.close()
else:
    tree = ET.parse(sys.argv[1])
    for child in tree.getroot():
        if child.tag == 'node':
            number = int(child.attrib['id'])
            coord[number] = (float(child.attrib['lat']),float(child.attrib['lon']))
        elif child.tag == 'way':
            nodes = []
            road = False
            oneway = False
            avg = mean_speed
            for ch in child:
                if ch.tag == 'nd':
                    nodes.append(int(ch.attrib['ref']))
                elif ch.tag == 'tag':
                    if ch.attrib['k'] == 'highway' and ch.attrib['v'] in ['motorway','motorway_link','trunk','trunk_link','primary','primary_link','secondary','secondary_link','tertiary','tertiary_link','unclassified','residential','service']:
                        road = True
                    if (ch.attrib['k'] == 'oneway' and ch.attrib['v'] == 'yes') or (ch.attrib['k'] == 'highway' and ch.attrib['v'] in ['motorway','motorway_link']):
                        oneway = True
                    if road and ch.attrib['k'][0:8] == 'maxspeed':
                        s = ch.attrib['v']
                        if s[-4:] == ' mph':
                            avg = float(s[0:-4]) * 1.609
                        elif s[-6:] == ' knots':
                            avg = float(s[0:-6]) * 1.852
                        elif not isinstance(s, float):
                            continue
                        else:
                            avg = float(s)
            speed = normal(avg, avg/4)
            while speed > avg or speed <= 0:
                speed = normal(avg, avg/4)
            if max_speed < speed:
                max_speed = speed
            for u,v in zip(nodes[:-1], nodes[1:]):
                if u in coord and v in coord:
                    if u not in vertices:
                        vertices[u] = num
                        num += 1
                    if v not in vertices:
                        vertices[v] = num
                        num += 1
                    if (u,v) not in edges:
                        edges[(u,v)] = [speed]
            if not oneway:
                nodes.reverse()
                for u,v in zip(nodes[:-1], nodes[1:]):
                    if u in coord and v in coord:
                        if u not in vertices:
                            vertices[u] = num
                            num += 1
                        if v not in vertices:
                            vertices[v] = num
                            num += 1
                        if (u,v) not in edges:
                            edges[(u,v)] = [speed]

lens = {}
probs = {}
print("{} {}".format(len(vertices), len(edges)))
for u, v in edges:
    lens[(u,v)] = round(distance(coord[u], coord[v]).km * 1000)
    probs[(u,v)] = 0.5 if sum(edges[(u,v)]) == 0 else sum(edges[(u,v)]) / (max_speed * len(edges[(u,v)]))

m_p = max(probs.values())
for u, v in edges:
    print("{} {} {:d} {:.6f}".format(vertices[u], vertices[v], lens[(u,v)], probs[(u,v)] / m_p))
