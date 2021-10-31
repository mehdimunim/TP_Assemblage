import networkx as nx
from kmer import *


def build_graph(kmer_dic):
    digraph = nx.DiGraph()
    for kmer in kmer_dic:
        sub_kmer1 = kmer[:-1]
        sub_kmer2 = kmer[1:]
        digraph.add_edge(sub_kmer1, sub_kmer2, weight=kmer_dic[kmer])
    return digraph


def get_starting_nodes(graph):
    res = []
    for node in graph:
        preds = [node for node in graph.predecessors(node)]
        if len(preds) == 0:
            res.append(node)
    return res


def get_sink_nodes(graph):
    res = []
    for node in graph:
        succ = [node for node in graph.successors(node)]
        if len(succ) == 0:
            res.append(node)
    return res


def get_contig_from_path(path):
    contig = ""
    for node_1, node_2 in zip(path[:-1], path[1:]):
        contig += node_1[0]
        contig += node_2[-1]
    return contig


def get_contigs(graph, starting_nodes, ending_nodes):
    res = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if nx.has_path(graph, start_node, end_node):
                paths = [get_contig_from_path(path) for path in nx.all_simple_paths(
                    graph, start_node, end_node)]
            res.append(paths)
    return res


def amax(data):
    """
    Return the index of the maximum in data. If several maxima, return the one with the smaller index.
    """
    max = 0
    amax = 0
    for i, item in enumerate(data):
        if max < item:
            max = item
            amax = i
    return amax


def main():
    mylist = [1, 2, 3, 4, 5]
    print(amax(mylist))
    mylist2 = [1, 6, 1, 1, 6]
    mylist.pop()
    print(amax(mylist2))


main()
