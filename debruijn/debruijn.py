#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import statistics
from random import randint
import argparse
import os
import sys
from operator import itemgetter
import pickle
import random
import matplotlib.pyplot as plt
import networkx as nx

random.seed(9001)

__author__ = "Mehdi MUNIM"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Mehdi MUNIM"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mehdi MUNIM"
__email__ = "mehdi.munim@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """
    Read a multi-fastq file and yield the sequences.
    """
    with open(fastq_file, "r+") as file:
        for i, line in enumerate(file):
            if i % 4 == 1:
                # remove trailing whitespaces
                yield line.rstrip()


def cut_kmer(read, kmer_size):
    """
    Cut a sequence in kmers
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i: i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """
    Build a dictionary for kmers in sequence
    """
    dic = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in dic:
                dic[kmer] += 1
            else:
                dic.update({kmer: 1})
    return dic


def build_graph(kmer_dict):
    """
    Build a graph from kmer_dict
    """
    digraph = nx.DiGraph()
    for kmer in kmer_dict:
        sub_kmer1 = kmer[:-1]
        sub_kmer2 = kmer[1:]
        digraph.add_edge(sub_kmer1, sub_kmer2, weight=kmer_dict[kmer])
    return digraph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    Remove paths in graph

    Parameters:
    ----------
        graph: networks.Digraph
            Graph object
        path_list: list
            list of path to delete
        delete_entry_node: boolean
            delete entry nodes for each path
        delete_sink_node:
            delete sink nodes for each path
    Returns:
    --------
        graph: networks.Digraph
            graph with paths removed
    """
    start = None
    end = None
    if not delete_sink_node:
        end = -1
    elif not delete_entry_node:
        start = 1
    for path in path_list:
        graph.remove_nodes_from(path[start:end])
    return graph


def std(data):
    """
    Compute the standard deviation of data

    Parameters:
    ----------
        data: list
            numerical data

    Return:
    ------
        std: float
            standard deviation
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass


def path_average_weight(graph, path):
    """
    Get the average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    pass


def simplify_bubbles(graph):
    pass


def solve_entry_tips(graph, starting_nodes):
    pass


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    """
    Get the starting nodes of the graph
    """
    res = []
    for node in graph:
        preds = list(graph.predecessors(node))
        if len(preds) == 0:
            res.append(node)
    return res


def get_sink_nodes(graph):
    """
    Ge the sink nodes of the graph
    """
    res = []
    for node in graph:
        succ = list(graph.successors(node))
        if len(succ) == 0:
            res.append(node)
    return res


def get_contig_from_path(path):
    """
    Get the contig corresponding to the input path
    """
    contig = ""
    contig += path[0]
    for node in path[1:]:
        contig += node[-1]
    return contig


def get_contigs(graph, starting_nodes, ending_nodes):
    """
    Get all contigs of the graph
    """
    res = []
    # find all paths between starting and ending nodes
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            contigs = [get_contig_from_path(path) for path in nx.all_simple_paths(
                graph, start_node, end_node)]
            res += [(contig, len(contig)) for contig in contigs]
    return res


def save_contigs(contigs_list, output_file):
    """
    Save contigs to fasta file
    """
    with open(output_file, "w+", newline="\n") as out:
        for i, contig in enumerate(contigs_list):
            # header
            out.write(f">contig_{i} len={contig[1]}\n")
            # sequence
            out.write(fill(contig[0])+"\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i: i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    # print(elarge)
    esmall = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    # nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
