#!/usr/bin/env python3

import matplotlib.pyplot as plt
import networkx as nx
import parsers as ps
import argparse
import time
from internal_coordinate import Graph
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

    mol = ps.read.mol2(molFile)
    mol2 = ps.read.mol2('M2.mol2')
    print('>>>> M1', mol.natm, mol.nbond, '>>>> M2', mol2.natm, mol2.nbond)
    mol.construct_graph()
    mol2.construct_graph()

    mol.tree(mol.selectbyAtomnum(1))
    mol2.tree(mol2.selectbyAtomnum(1))

    start = time.time()
    print(Graph.get_similarity(mol.graph.head, mol2.graph.head))
    finish =time.time()
    print(finish-start)
    exit()
    """
    conectivity = {}
    for iatom in mol.atoms:
        neighbor = []
        for bond in mol.selectbondbyAtom(iatom):
            for atom in bond.atoms:
                if atom != iatom:
                    neighbor.append(atom.num)
        conectivity[iatom.num] = neighbor

    for atom in conectivity:
        print(atom, '->', conectivity[atom])

    for atom in conectivity:
        for atom_neighbor in conectivity[atom]:
            if atom in conectivity[atom_neighbor]:
                conectivity[atom_neighbor].remove(atom)
    print('***********************************')
    for atom in conectivity:
        print(atom, '->', conectivity[atom])
    """
    mol.construct_graph()
    i = 6
    print('atom: ', i, '----------------------------------------')
    head = mol.selectbyAtomnum(i)
    mol.tree(head)
    #mol.graph.get_rings()
    #for i in range(1, 4):
    #    list = []
    #    for node in mol.graph.nodes:
    #        if i in node.ring:
    #            list.append(node.atom.num)
    #    print('ring ', i, 'contains', list )
    start = time.time()
    mol.graph.get_rings2()
    end = time.time()
    print('Ring Search took', end - start)
    print('----------------------------------------')

    exit()
    for ring_id, ring in mol.graph.rings.items():
        print(ring_id, [node.atom.num for node in ring])

    print('----------------------------------------')

    for node in mol.graph.nodes:
        print(node.atom.num, node.ring_id)

if __name__ == '__main__':
    molecular_graph_test()
