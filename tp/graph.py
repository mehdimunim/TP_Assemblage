import networkx as nx
from kmer import *


def build_graph(kmer_dic):
    digraph = nx.DiGraph()
    for kmer in kmer_dic:
        sub_kmer1 = kmer[:-1]
        sub_kmer2 = kmer[1:]
        digraph.add_edge(sub_kmer1, sub_kmer2, weight=kmer_dic[kmer])
    return digraph


def main():
    seq = "ATTCGGGGCCA"
    len = 2
    dic = build_kmer_dict(seq, len)
    graph = build_graph(dic)
    print(graph.edges)


main()
