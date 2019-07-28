#!/usr/bin/env python3

import matplotlib.pyplot as plt
import networkx as nx
import parsers as ps
import argparse
import time
from graph import Graph
#import sys
from internal_coordinate import InternalCoordinate
#sys.setrecursionlimit(30)


def molecular_graph_test():
    paser = argparse.ArgumentParser()
    paser.add_argument("mol", type=str)

    try:
        arg = paser.parse_args()
        molFile = arg.mol

    except:
        print("Error >>> Need input")
        return 1

    mol1 = ps.read.mol2(molFile)
    mol2 = ps.read.mol2('M2.mol2')
    print('>>>> M1', mol1.natm, mol1.nbond, '>>>> M2', mol2.natm, mol2.nbond)

    mol1_graph = Graph(mol1)
    mol2_graph = Graph(mol2)

    """
    for node in mol1_graph.nodes:
        print(node.atom.nam, node.atom.num, node.rings)

    for node in mol2_graph.nodes:
        print(node.atom.nam, node.atom.num, node.rings)
    """
    mol2_graph.construct_tree(head=mol2_graph.nodes[0])

    start = time.time()
    print('(1-2) ')
    for node in mol1_graph.nodes:
        mol1_graph.construct_tree(head=node)
        print(mol1_graph.head.atom.num, mol1_graph.score(mol2_graph))
    finish =time.time()
    print(finish-start)

    for node in mol1_graph.nodes:
        print(node.atom.num, "--------------")
        for neighbor in node.neighbors:
            print(neighbor.atom.num, node._bond_order[neighbor])
    #ps.write.mol2(mol1_graph, 'M1_test')

if __name__ == '__main__':
    molecular_graph_test()
