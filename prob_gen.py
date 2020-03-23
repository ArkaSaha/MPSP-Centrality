from geopy.distance import distance
from statistics import mean
import sys
from os import popen

f = open(sys.argv[1],'r')
coord = {}
for i in range(int(popen('wc -l '+sys.argv[1]).read().split()[0])):
	line = f.readline().lstrip('\t')
	if (line[0:6] == '<node '):
		l = line.split()
		id_no = l[1][4:-1]
		coord[id_no] = (float(l[2][5:-1]),float(l[3][5:-1]))
f.close()
f = open(sys.argv[2],'r')
vertices = {}
edges = {}
max_speed = 0
num = 0
for i in range(0,int(popen('wc -l '+sys.argv[2]).read().split()[0]),2):
	nodes = f.readline().rstrip('\t\n').split('\t')
	speeds = []
	l = f.readline().rstrip('\t\n').split('\t')
	for s in l:
		sp = float(s)
		speeds.append(sp)
		if sp > max_speed:
			max_speed = sp
	for j in range(len(speeds)-1):
		if nodes[j] not in vertices:
			vertices[nodes[j]] = num
			num += 1
		if nodes[j+1] not in vertices:
			vertices[nodes[j+1]] = num
			num += 1
		if (nodes[j],nodes[j+1]) not in edges:
			edges[(nodes[j],nodes[j+1])] = []
		edges[(nodes[j],nodes[j+1])].append(speeds[j])
f.close()
print(str(len(vertices))+' '+str(len(edges)))
for u, v in edges:
	print(str(vertices[u])+' '+str(vertices[v])+' '+str(round(distance(coord[u],coord[v]).km*10**3))+' '+str(mean([s/max_speed for s in edges[(u,v)]])))