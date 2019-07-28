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
        self._bond_order = dict()
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

    def bond_order(self, neighbor, value=None):
        """
        Returns the bond_order of the bond the node and its neighbor.
        :parameter neighbor
        """
        # Return bond order
        if value is None:

            if type(neighbor) is not Node:
                print("Error >>> Expecting a Node object, received a ", type(neighbor))
                return None

            if neighbor not in self.neighbors:
                print("Error >>> The neighbor was not found in the node's neighbor list ", type(neighbor))
                return None

            return self._bond_order[neighbor]

        # Set bond order
        else:
            try:
                value = float(value)
            except ValueError:
                print("Error >>> Bond order should be a float. Recieved ", value)
                return None
            self._bond_order[neighbor] = value

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
    def __init__(self, internalcoord=None):
        self._nodes = None
        self._head = None
        self._tree = False
        self._name = None
        #self.count = 0
        self.rings = {}
        #
        if internalcoord:
            self.construct_graph(internalcoord)
            self.set_rings()

    def __str__(self):
        return self.name

    @property
    def head(self):
        return self._head

    @property
    def tree(self):
        return self._tree

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
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

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
        #self._nodes = []
        self.name = internalcoord.nam
        # Populate the nodes
        self.nodes = [Node(atom) for atom in internalcoord.atoms]

        # Find the neighbor list of each node based on the bonds
        for node_i in self.nodes:
            for node_j in self.nodes:
                ibond = Bond(node_i.atom, node_j.atom)
                if ibond in internalcoord.bonds:
                    # Get the actual bond in a list format
                    ibond = internalcoord.selectbond(Bond(node_i.atom, node_j.atom))
                    if len(ibond) == 1:
                        node_i.neighbors = node_j
                        node_i._bond_order[node_j] = ibond[0].order
                        #print(ibond[0].order, "**")
                    else:
                        print("Error >>> atoms ", node_i.atom.num, node_j.atom.num, " seems to have repeated bond definitions")

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
        self._tree = False
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

        self._tree = True

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

    def select_by_atom_num(self, atom_num):

        # Check the input
        if type(atom_num) is not int:
            print("Error >>> an integer expected, a ", type(atom_num), " is given.")
            return None

        for node in self.nodes:
            if node.atom.num == atom_num:
                return node

        # If we get here the node was not found.
        print("Warn >>> the atom number ", atom_num, " was not found.")
        return None

    def prune(self, node):

        # Check the input
        if type(node) is not Node:
            print("Error >>> a node is expected. A ", type(node), " is given.")
            return None

        # Check the tree.
        if not self.tree:
            print("Error >>> pruning can only be performed on a tree")
            return None

        # Check the head
        if node == self.head:
            print("Error >>> the head atom can not be pruned.")
            return None

        # Check node is in the tree
        if node not in self.nodes:
            print("Error >>> the node was not found.")
            return None

        # initialize
        next_search = []
        current_search = node.children
        nodes_to_delete = [node]

        # Do a BFS to sweep the tree.
        while current_search:
            for cnode in current_search:
                if cnode not in nodes_to_delete:
                    nodes_to_delete.append(cnode)
                for child in cnode.children:
                    if child not in current_search and child not in next_search:
                        next_search.append(child)
            current_search = next_search
            next_search = []

        # delete the nodes
        for node in nodes_to_delete:
            self.delete(node)

    def delete(self, node):

        # Check the input
        if type(node) is not Node:
            print("Error >>> a node is expected. A ", type(node), " is given.")
            return None

        # Check node is in the graph
        if node not in self.nodes:
            print("Error >>> the node was not found.")
            return None

        # Remove the node from the neighbor list of the neighboring nodes.
        for node_neighbor in node.neighbors:
            node_neighbor.neighbors.remove(node)
            node_neighbor._bond_order.pop(node)

        # Remove the node from the children list of its parents.
        for node_parent in node.parents:
            node_parent.children.remove(node)

        # Remove the node from the parent list of its children.
        for node_child in node.children:
            node_child.parents.remove(node)

        # Remove the node itself
        self._nodes.remove(node)

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
        if not self.tree or not other.tree:
            print('Error >>> To get the similarity score, first construct the tree structure.')
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

    @property
    def bonds(self):
        bonds = []
        for node in self.nodes:
            for neighbor in node.neighbors:
                ibond = Bond(i=node.atom, j=neighbor.atom, order=node.bond_order(neighbor))
                if ibond not in bonds:
                    bonds.append(ibond)
        bonds.sort()
        return bonds


class Fragment(Graph):
    ForceFileds = ['UFF', 'OPLS', 'AMBER']

    def __init__(self, internalcoord, name):
        super().__init__(internalcoord)
        self._name = name
        self._atom_type = dict()

    def __str__(self):
        return self.name


    def atom_type(self, forcefield=None, value=None):

        # Return the atom type definition as a list of lists.
        if forcefield is None and value is None:
            return [[ff, typ] for ff, typ in self._atom_type.items()]

        # Return the atom type of the given force field
        elif forcefield is not None and value is None:
            return self._atom_type[forcefield]

        # Set the force field atom type
        elif forcefield is not None and value is not None:
            self._atom_type[forcefield] = value

        # Default
        else:
            print("Error >>> invalid arguments for the atom type function")
            return None

    @property
    def fragment(self):
        # Check whether the fragment meet the minimum requirement.
        finished = True

        # Check force field definition
        try:
            # Check the minimum FF definition (UUF) exists.
            if len(self._atom_type['UFF']) == 0:
                print("Warn >>> no UFF atom type is set.")
                finished = False
        except KeyError:
            print("Warn >>> No atom type is set.")
            finished = False

        # Check the structure
        if not self.tree:
            print("Warn >>> the head atom is not set.")
            finished = False

        return finished


class Fragment_Library(object):
    def __init__(self, name):
        self.name = str(name)
        self.element = dict()

    def __str__(self):
        return self.name

    def __iter__(self):
        for element, fragment in self.element.items():
            yield element, fragment
