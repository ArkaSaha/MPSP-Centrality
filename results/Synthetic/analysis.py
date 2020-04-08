
import sys
import re

def read_wise_file(filename):
    f = open(filename)
    data = {}
    summary = {}
    cur_hop =0 
    lines = f.read().splitlines()

    i = 0
    while i < len(lines):
        if "Number of hops" in lines[i]:
            cur_hop = int(lines[i][-1])
            data[cur_hop] = {}
            summary[cur_hop] = {}
        elif lines[i] != "":
            # process new block
            m = re.match("(\d*)\s(\d*)", lines[i])
            if m:
                u, v = int(m.group(1)), int(m.group(2))
                stats = {}
                if lines[i+1][0] == "*":
                    stats['timeout'] = True
                    i = i+1
                else:
                    stats['timeout'] = False
                stats['candidates'] = int(lines[i+1].split()[-1])
                i = i+3
                stats['edges'] = []
                stats['path'] = [int(lines[i].split()[0])]
                while lines[i][0] != 'L':
                    line_split = lines[i].split()
                    stats['edges'].append({"u" : int(line_split[0]),
                                           "v" : int(line_split[1]),
                                           "l" : int(line_split[2]),
                                           "p" : float(line_split[3])})
                    stats['path'].append(int(line_split[1]))
                    i = i + 1
                stats['path'] = tuple(stats['path'])
                assert(stats['path'][0] == u)
                assert(stats['path'][-1] == v)
                stats['length'] = int(lines[i].split()[-1])
                stats['prob'] = float(lines[i+1].split()[-1])
                stats['path_prob'] = float(lines[i+2].split()[-1])
                stats['candidate_time'] = float(lines[i+3].split()[-2])
                stats['prob_time'] = float(lines[i+4].split()[-2])

                data[cur_hop][(u,v)] = stats

                i = i + 4
            else:
                summary[cur_hop]['length'] =  float(lines[i].split()[-1])
                summary[cur_hop]['prob'] =  float(lines[i+1].split()[-1])
                summary[cur_hop]['candidate_time'] =  float(lines[i+2].split()[-2])
                summary[cur_hop]['prob_time'] =  float(lines[i+3].split()[-2])
                summary[cur_hop]['total_time'] =  float(lines[i+4].split()[-2])
                summary[cur_hop]['timeouts'] =  int(lines[i+5].split()[-1])
                i = i+5
        i = i+1
    f.close()
    return data, summary

def read_mpsp_file(filename):
    f = open(filename)
    data = {"pruning": {}, "no_pruning" : {}}
    summary = {"pruning" : {}, "no_pruning" : {}}
    cur_hop =0 
    lines = f.read().splitlines()

    i = 0
    while i < len(lines):
        if "Number of hops" in lines[i]:
            cur_hop = int(lines[i][-1])
            data["pruning"][cur_hop] = {}
            data["no_pruning"][cur_hop] = {}
            summary["pruning"][cur_hop] = {}
            summary["no_pruning"][cur_hop] = {}
        elif lines[i] != "":
            # process new block
            m = re.match("(\d*)\s(\d*)", lines[i])
            if m:
                u, v = int(m.group(1)), int(m.group(2))
                pruning, no_pruning = {}, {}
                if lines[i+1][0] == "N":
                    pruning["path_found"] = False
                    i = i+1
                else:
                    pruning["path_found"] = True
                    i = i+2
                    pruning['edges'] = []
                    pruning['path'] = [int(lines[i].split()[0])]
                    while lines[i][0] != 'W':
                        line_split = lines[i].split()
                        pruning['edges'].append({"u" : int(line_split[0]),
                                               "v" : int(line_split[1]),
                                               "l" : float(line_split[2]),
                                               "p" : float(line_split[3])})
                        pruning['path'].append(int(line_split[1]))
                        i = i + 1
                    pruning['path'] = tuple(pruning['path'])
                    assert(pruning['path'][0] == u)
                    assert(pruning['path'][-1] == v)

                    i = i+1
                    no_pruning['edges'] = []
                    no_pruning['path'] = [int(lines[i].split()[0])]
                    while lines[i][0] != 'N':
                        line_split = lines[i].split()
                        no_pruning['edges'].append({"u" : int(line_split[0]),
                                               "v" : int(line_split[1]),
                                               "l" : float(line_split[2]),
                                               "p" : float(line_split[3])})
                        no_pruning['path'].append(int(line_split[1]))
                        i = i + 1
                    no_pruning['path'] = tuple(no_pruning['path'])
                    assert(no_pruning['path'][0] == u)
                    assert(no_pruning['path'][-1] == v)

                pruning['dijkstra_runs'] = int(lines[i].split()[-1])
                no_pruning['dijkstra_runs'] = int(lines[i].split()[-1])
                pruning['pruned'] = int(lines[i+1].split()[-1])
                pruning['candidates'] = int(lines[i+2].split()[-1])
                no_pruning['candidates'] = int(lines[i+3].split()[-1])
                pruning['freq_modal_path'] = int(lines[i+4].split()[-1])
                pruning['mpsp_is_modal'] = int(lines[i+5].split()[-1])
                no_pruning['mpsp_is_modal'] = int(lines[i+6].split()[-1])
                pruning['length'] = int(lines[i+7].split()[-1])
                no_pruning['length'] = int(lines[i+8].split()[-1])
                pruning['prob_path'] = float(lines[i+9].split()[-1])
                no_pruning['prob_path'] = float(lines[i+10].split()[-1])
                pruning['prob'] = float(lines[i+11].split()[-1])
                no_pruning['prob'] = float(lines[i+12].split()[-1])
                pruning['candidate_time'] = float(lines[i+13].split()[-2])
                no_pruning['candidate_time'] = float(lines[i+14].split()[-2])
                pruning['prob_time'] = float(lines[i+15].split()[-2])
                no_pruning['prob_time'] = float(lines[i+16].split()[-2])

                data["pruning"][cur_hop][(u,v)] = pruning
                data["no_pruning"][cur_hop][(u,v)] = no_pruning

                i = i + 16
            else:
                summary["pruning"][cur_hop]['non_empty'] =  int(lines[i].split()[-1])
                summary["no_pruning"][cur_hop]['non_empty'] =  int(lines[i].split()[-1])
                summary["pruning"][cur_hop]['mathces'] =  int(lines[i+1].split()[-1])
                summary["no_pruning"][cur_hop]['matches'] =  int(lines[i+2].split()[-1])
                summary["pruning"][cur_hop]['length'] =  float(lines[i+3].split()[-1])
                summary["no_pruning"][cur_hop]['length'] =  float(lines[i+4].split()[-1])
                summary["pruning"][cur_hop]['prob'] =  float(lines[i+5].split()[-1])
                summary["no_pruning"][cur_hop]['prob'] =  float(lines[i+6].split()[-1])
                summary["pruning"][cur_hop]['dijkstra_runs'] =  float(lines[i+7].split()[-1])
                summary["no_pruning"][cur_hop]['dijkstra_runs'] =  float(lines[i+7].split()[-1])
                summary["pruning"][cur_hop]['samples_pruned'] =  float(lines[i+8].split()[-1])
                summary["no_pruning"][cur_hop]['candidates'] =  float(lines[i+9].split()[-1])
                summary["pruning"][cur_hop]['candidates'] =  float(lines[i+10].split()[-1])
                summary["no_pruning"][cur_hop]['candidate_time'] =  float(lines[i+11].split()[-2])
                summary["pruning"][cur_hop]['candidate_time'] =  float(lines[i+12].split()[-2])
                summary["no_pruning"][cur_hop]['prob_time'] =  float(lines[i+13].split()[-2])
                summary["pruning"][cur_hop]['prob_time'] =  float(lines[i+14].split()[-2])
                summary["no_pruning"][cur_hop]['total_time'] =  float(lines[i+15].split()[-2])
                summary["pruning"][cur_hop]['total_time'] =  float(lines[i+16].split()[-2])

                i = i+16
        i = i+1
    f.close()
    return data, summary

def convert_file_to_data(graph):
    graph_type = graph[:2]
    wise = {}

    filename = "{}/Converge/{}.output".format(graph_type, graph)
    mpsp = {}
    mpsp["data"], mpsp["summary"] = read_mpsp_file(filename)
    res = []

    for seconds in [1, 10, 60]:
        filename = "{}/Wise/{}_WISE_{}.output".format(graph_type, graph, seconds)
        wise[seconds] = {}
        wise[seconds]["data"], wise[seconds]["summary"] = read_wise_file(filename)

    return mpsp, wise

def compare(wise, mpsp, seconds):
    s = seconds 
    tol = 2 
    same, same_prob, mpsp_better, wise_better, no_mpsp = 0, 0, 0,0,0
    for h in [2,4,6,8,0]: # mpsp["data"]["pruning"]:
        for k in mpsp["data"]["pruning"][h]:
            if not mpsp["data"]["pruning"][h][k]["path_found"]:
                no_mpsp = no_mpsp + 1
            else:
                if mpsp["data"]["pruning"][h][k]["path"] == wise[s]["data"][h][k]["path"]:
                    same = same + 1
                elif mpsp["data"]["pruning"][h][k]["prob"] > tol * wise[s]["data"][h][k]["prob"]:
                    mpsp_better = mpsp_better + 1
                elif mpsp["data"]["pruning"][h][k]["prob"] * tol <  wise[s]["data"][h][k]["prob"]:
                    wise_better = wise_better + 1
                else:
                    same_prob = same_prob + 1
    return [same, same_prob, mpsp_better, wise_better, no_mpsp]

def compare_latex(mpsp, wise):
    res = []
    for seconds in [1, 10, 60]:
        res = res + compare(wise, mpsp, seconds)
    return " & ".join(str(x) for x in res)

def average_time(mpsp, wise):
    output = []
    #for h in mpsp["data"]["pruning"]:
    for h in [2,4,6,8,0]:
        mpsp_time = 0
        wise_time = {1: 0, 10: 0, 60: 0}
        nr = 0
        for k in mpsp["data"]["pruning"][h]:
            if mpsp["data"]["pruning"][h][k]["path_found"]:
                nr = nr + 1
                mpsp_time = mpsp_time + mpsp["data"]["pruning"][h][k]["prob_time"] +  \
                                        mpsp["data"]["pruning"][h][k]["candidate_time"]
                for s in [1, 10, 60]:
                    wise_time[s] = wise_time[s] + wise[s]["data"][h][k]["prob_time"] +  \
                                              wise[s]["data"][h][k]["candidate_time"]
        output = output + [mpsp_time/nr, wise_time[1]/nr, wise_time[10]/nr, wise_time[60]/nr]
    return output

def average_time_latex(mpsp, wise):
    return " & ".join("{:.5f}".format(x) for x in average_time(mpsp, wise))


def format_nodes(n):
    n = int(n)
    if n>= 1000000:
        return "{}M".format(n//1000000)
    else:
        return "{}K".format(n//1000)




BAgraphs = ["BA_10000_19997", "BA_20000_39997", "BA_50000_99997", "BA_100000_199997", 
            "BA_500000_999997", "BA_1000000_1999997", "BA_5000000_9999997", "BA_10000000_19999997"]

ERgraphs = ["ER_10000_20000", "ER_20000_40000", "ER_50000_100000", "ER_100000_200000", 
            "ER_500000_1000000", "ER_1000000_2000000", "ER_5000000_10000000", "ER_10000000_20000000"]

# compare in numbers
for graph_list in (BAgraphs, ERgraphs):
    print("\n\n")
    for g in graph_list:
        m = re.match(".._(\d*)_\d*", g)
        print("& {} & ".format(format_nodes(m.group(1))), end="")
        mpsp, wise = convert_file_to_data(g)
        print(compare_latex(mpsp, wise), end="")
        print(" \\\\")


# compare average times
for graph_list in (BAgraphs, ERgraphs):
    print("\n\n")
    for g in graph_list:
        m = re.match(".._(\d*)_\d*", g)
        print("& {} & ".format(format_nodes(m.group(1))), end="")
        mpsp, wise = convert_file_to_data(g)
        print(average_time_latex(mpsp, wise), end="")
        print(" \\\\")




