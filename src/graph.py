#!usr/bin/python3
"""
    module :: internal_coordinate
    :platform: Unix
"""

import numpy as np
import pickle
from atom import Atom
from molecule import Molecule
from atom_table import *
import itertools as it
from internal_coordinate import *


class Node(object):
    """
    This class is used represent the node in a graph representation of a molecule.

       - **Parm atom** (:class:`Atom`):      The input atom object.
       >>> print('Write the example')
    """

    def __init__(self, atom):
        self._atom = None
        self._parents = []
        self._children = []
        self._neighbors = []
        self.visited = False
        self.head_flag = False
        #self.color = 0
        #self.number_of_cycles = 0
        self._rings = []

        # Setting the atom by setter for tyoe checking
        self._atom = atom

    #def __hash__(self):
    #    return self.atom.num

    @property
    def atom(self):
        """
            Sets/Returns the atom member.
        """
        return self._atom

    @atom.setter

    def atom(self, atom):
        if type(atom) is Atom:
            self._atom = atom
        else:
            raise ValueError('Error >>> Pass an Atom object')

    @property
    def number_of_neighbor(self):
        """
        Returns the number of neighrs of the node.
        """
        return len(self._neighbors)

    @property
    def neighbors(self):
        return self._neighbors

    @neighbors.setter
    def neighbors(self, neighbor):
        if type(neighbor) is Node:
            self._neighbors.append(neighbor)
        else:
            raise ValueError('Error >>> Pass a Node object')

    @property
    def number_of_parents(self):
        return len(self._parents)

    @property
    def parents(self):
        return self._parents

    @parents.setter
    def parents(self, parents):
        if any([type(node) != Node for node in parents]):
            raise ValueError('Pass a list of Node objects')
        else:
            # TODO should be log
            # if self._parents:
            #    print('Warning >>> Overwriting the node\'s parents')
            self._parents = parents

    def add_parent(self, parent):
        if type(parent) is Node:
            self._parents.append(parent)
        else:
            raise ValueError('Error >>> Pass a Node object')

    @property
    def number_of_children(self):
        return len(self._children)

    @property
    def children(self):
        return self._children

    @children.setter
    def children(self, children):
        if any([type(node) != Node for node in children]):
            raise ValueError('Pass a list of Node objects')
        else:
            # TODO should be log
            # if self._children:
            #    print('Warning >>> Overwriting the node\'s children')
            self._children = children

    def add_children(self, child):
        if type(child) is Node:
            self._children.append(child)
        else:
            raise ValueError('Error >>> Pass a Node object')

    @property
    def rings(self):
        return self._rings

    @rings.setter
    def rings(self, ring):
        if type(ring) is int:
            self._rings.append(ring)
        elif ring == []:
            self._rings = []
        else:
            raise ValueError('Error >>> ring is should be an integer')

class Graph(object):
    """
    This class is used represent a molecule as a graph structure.

       - **Parm atom** (:class:`InternalCoordinate`):      The input atom object.
       >>> print('Write the example')
    """
    def __init__(self, internalCoord=None):
        self._nodes = None
        self._head = None
        #self.count = 0
        self.rings = {}
        #
        if internalCoord:
            self.construct_graph(internalCoord)
            self.set_rings()

    @property
    def head(self):
        return self._head

    @head.setter
    def head(self, head):
        if type(head) is Node:
            if head in self.nodes:
                # First rest the previous head node flag if already defined
                if self._head:
                    self._head.head_flag = False
            # Point to the new head
            self._head = head
            self._head.head_flag = True
        else:
            raise ValueError('head should be a Node object')

    @property
    def number_of_nodes(self):
        return len(self._nodes)

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, nodes):
        if any([type(node) != Node for node in nodes]):
            raise ValueError('Pass a list of Node objects')
        else:
            if self._nodes:
                print('Warning >>> Overwriting the graph Nodes')
            self._nodes = nodes

    def add_node(self, node):
        if type(node) is Node:
            self._nodes.append(node)
        else:
            raise ValueError('Error >>> Pass a Node object')

    def construct_graph(self, internalcoord):
        # Check the input. It is allowed for a graph to have one atom and no bond.
        if type(internalcoord) is not InternalCoordinate:
            print('Error >>> an InternalCoordinate is needed.')
            return None

        # overwrite the references
        self._nodes = []

        # Populate the nodes
        self._nodes = [Node(atom) for atom in internalcoord.atoms]

        # Find the neighbor list of each node based on the bonds
        for node_i in self._nodes:
            for node_j in self._nodes:
                if Bond(node_i.atom, node_j.atom) in internalcoord.bonds:
                    node_i.neighbors = node_j

    def reset_visited(self):
        # Reset the vested flag
        #self.count = 0
        for node in self.nodes:
            node.visited = False
            node.color = 0

    def reset_tree(self):
        """
        This method resets the all tree features (the tree head, the nodes head_flag, parents and childeren).
        the graph features (nodes neighbor list and rings list) are not affected.
        """
        # Reset the flag
        #self.count = 0
        self._head = None
        for node in self.nodes:
            node.visited = False
            node.head_flag = False
            node.parents = []
            node.children = []

    def construct_tree(self, head, method='BFS'):

        if head not in self.nodes:
            print('head atom not found in graph nodes')
            return None

        # Reset the tree for new tree structure
        self.reset_tree()

        # Set the new head
        self.head = head

        if method == 'DFS':
            self._tree_DFS(head)
        else:
            self._tree_BFS(head)

    def _tree_DFS(self, I):

        # Reset only for the first round
        if I.head_flag:
            self.reset_visited()
            I.visited = True

        for my_neighbor in I.neighbors:
            # print('call', self.count)
            # print('I', I.visited, I.atom.num, ('parent: ', [p.atom.num for p in I.parents], 'children: ', [c.atom.num for c in I.children]),
            #      ' neighbor: ', my_neighbor.visited, my_neighbor.atom.num, ('parent: ', [p.atom.num for p in my_neighbor.parents], 'children: ', [c.atom.num for c in my_neighbor.children]))
            #self.count += 1

            if (I not in my_neighbor.parents) and (I not in my_neighbor.children):
                my_neighbor.add_parent(I)

            if (my_neighbor not in I.parents) and (my_neighbor not in I.children):
                I.add_children(my_neighbor)

            if not my_neighbor.visited:
                my_neighbor.visited = True
                self._tree_DFS(my_neighbor)

    def _tree_BFS(self, head):

        # Search
        search_current = [head]
        search_next = []
        while search_current:
            # print(self.count)
            for I in search_current:
                if not I.visited:
                    I.visited = True
                    for my_neighbor in I.neighbors:
                        if (my_neighbor not in I.parents) and (my_neighbor not in I.children):
                            I.add_children(my_neighbor)
                        if (I not in my_neighbor.parents) and (I not in my_neighbor.children):
                            my_neighbor.add_parent(I)
                        search_next.append(my_neighbor)

                        # print('I', I.visited, I.atom.num, ('parent: ', [p.atom.num for p in I.parents], 'children: ', [c.atom.num for c in I.children]),
                        #    ' neighbor: ', my_neighbor.visited, my_neighbor.atom.num, ('parent: ', [p.atom.num for p in my_neighbor.parents], 'children: ', [c.atom.num for c in my_neighbor.children]))
            search_current = search_next
            search_next = []
            #self.count += 1

    def set_rings(self, depth=8):
        """
        This function set the ring membership of all nodes of the graph. This is a brute force DFS
        search that finds all rings for a atom up to a cut off which is defined by depth parameter.
        :param depth: The cut off for Max ring length. Default is XX
        """
        max_depth = depth
        # Find the shortest circuit for each node
        #rings = {}
        for node in self.nodes:
            if node.number_of_neighbor < 2:
                continue # Skip the leaves

            # Reset the node ring list
            node.rings = []

            # Initialize
            rings = set() # It is a set of sets, to save the unique rings.
            searched_path = []
            node_i = node
            node_reff = node
            parent_i = Node(Atom(num=0)) # Dummy

            # Find the rings
            self._find_cycle_DFS(node_reff, node_i, parent_i, searched_path, rings, max_depth)

            # Save only the rings length. Could be [6, 6], i.e. the atom belongs to two six-membered rings.
            for ring in rings:
                node.rings.append(len(ring))
            node.rings.sort()
        """
        rings_unique = set()
        node_superset = set()

        for node, node_rings in rings.items():
            # the previous cycle
            print('node_superset', [[node.atom.num for node in node_superset]])
            node_superset = set()
            l = min(map(len, node_rings))
            for iring in node_rings:
                if len(iring) <= l:
                    node_superset.update(iring)
                if iring.issubset(node_superset):
                    pass
                print(node.atom.num, '#r: ',l, [node.atom.num for node in iring])

            #rings_unique.add(frozenset(node for node in ring))

        exit()
        # Give the rings ID and save them in the graph.rings
        # Set the ring ID each node
        for i, ring in enumerate(rings_unique):
            self.rings[i] = ring
            for node in self.nodes:
                if node in ring:
                    node.ring_id.append(i)
            #print(len(ring), [node.atom.num for node in ring])
        """

    def _find_cycle_DFS(self, node_reff, node_i, parent_i, searched_path, rings, max_depth):

        # Keep track of search depth
        max_depth -= 1

        # Pass by search path by value to avoid changing the one of the caller function.
        current_searched = searched_path.copy()

        # Add current node to the search path of this function call.
        current_searched.append(node_i)

        for neighbor in node_i.neighbors:
            # Check if the reference node is found. If so add the ring length to the node ring list.
            if neighbor == node_reff and neighbor != parent_i:
                rings.add(frozenset(current_searched))
                return
            # Else go further.
            if neighbor not in current_searched and neighbor.number_of_neighbor > 1 and max_depth > 0:
                self._find_cycle_DFS(node_reff, neighbor, node_i, current_searched, rings, max_depth)

    def score(self, other):
        if not self.head or not other.head:
            print('Error >>> To get the similarity score, construct the tree structure')
            return None
        else:
            return Graph.get_score(self.head, other.head)

    @staticmethod
    def get_score(node1, node2):
        # TODO make sure it is a tree

        score_max = 0
        score_node_pair_matrix = dict()

        for node2_child in node2.children:
            for node1_child in node1.children:

                score_node_pair_matrix[(node2_child, node1_child)] = 0
                # Check nodes
                if node1_child.atom.typ == node2_child.atom.typ:
                    score_node_pair_matrix[(node2_child, node1_child)] +=1
                    #print([i for i in score_node_pair_matrix.items()])
                    #print(node2_child.atom.num, node1_child.atom.num)
                # Check the sub trees
                score_node_pair_matrix[(node2_child, node1_child)] += Graph.get_score(node1_child, node2_child)

        #print([(key[0].atom.nam, key[-1].atom.nam, value) for key, value in score_node_pair_matrix.items()])
        #print(len([(key[0].atom.num, key[-1].atom.num, value) for key, value in score_node_pair_matrix.items()]))
        # Get the combination of nodes that maximize the score
        """
        for permutation_node1 in it.permutations(node1.children):
            for permutation_node2 in it.permutations(node2.children):

                # Calculate the score with this mapping
                score_tm = 0
                for node_pair in zip(permutation_node2, permutation_node1):
                    score_tm += score_node_pair_matrix[node_pair]

                if score_tm > score_max:
                    score_max = score_tm

        """
        if node1.number_of_children >= node2.number_of_children:

            for permutation in it.permutations(node1.children):
                # Calculate the score with this mapping
                score_tm = 0
                for node_pair in zip(node2.children, permutation):
                    score_tm += score_node_pair_matrix[node_pair]

                if score_tm > score_max:
                    score_max = score_tm

        else:

            for permutation in it.permutations(node2.children):
                # Calculate the score with this mapping
                score_tm = 0
                for node_pair in zip(permutation, node1.children):
                    score_tm += score_node_pair_matrix[node_pair]

                if score_tm > score_max:
                    score_max = score_tm

        #print(score_max)
        return score_max

