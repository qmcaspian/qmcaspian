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


class Node(object):
    def __init__(self, atom):
        self._atom = None
        self._parents = []
        self._children = []
        self._neighbors = []
        self.visited = False
        self.head_flag = False
        self.color = 0
        self.number_of_cycles = 0
        self.ring_id = []

        self._atom = atom

    def __hash__(self):
        return self.atom.num

    @property
    def atom(self):
        return self._atom

    @atom.setter
    def atom(self, atom):
        if type(atom) is Atom:
            self._atom = atom
        else:
            raise ValueError('Error >>> Pass an Atom object')

    @property
    def number_of_neighbor(self):
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


class Graph(object):
    def __init__(self):
        self._nodes = None
        self._head = None
        self.count = 0
        self.rings = {}

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

    def reset_visited(self):
        # Reset the vested flag
        self.count = 0
        for node in self.nodes:
            node.visited = False
            node.color = 0

    def reset_tree(self):
        # Reset the vested flag
        self.count = 0
        for node in self.nodes:
            node.visited = False
            node.parents = []
            node.children = []

    def make_tree_DFS(self, I):

        # Reset only for the first round
        if I.head_flag:
            self.reset_visited()
            I.visited = True

        if I not in self.nodes:
            print('head atom not found in graph nodes')
            return None

        for my_neighbor in I.neighbors:
            # print('call', self.count)
            # print('I', I.visited, I.atom.num, ('parent: ', [p.atom.num for p in I.parents], 'children: ', [c.atom.num for c in I.children]),
            #      ' neighbor: ', my_neighbor.visited, my_neighbor.atom.num, ('parent: ', [p.atom.num for p in my_neighbor.parents], 'children: ', [c.atom.num for c in my_neighbor.children]))
            self.count += 1

            if (I not in my_neighbor.parents) and (I not in my_neighbor.children):
                my_neighbor.add_parent(I)

            if (my_neighbor not in I.parents) and (my_neighbor not in I.children):
                I.add_children(my_neighbor)

            if not my_neighbor.visited:
                my_neighbor.visited = True
                self.make_tree_DFS(my_neighbor)

    def make_tree_BFS(self, head):

        if head not in self.nodes:
            print('the head atom not found in graph nodes')
            return None

        # reset
        self.reset_tree()

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
            self.count += 1

    # TODO Not working
    def get_rings(self, node=None):

        if not self.nodes:
            print('Error >>> to find ring graph need at least three nodes, None fond')

        elif self.number_of_nodes < 3:
            print('Error >>> to find ring graph need at least three nodes,', self.number_of_nodes, ' found.')

        self.reset_visited()
        self.reset_tree()
        self.number_of_cycles = 0

        self.number_of_cycles
        node = self.head
        parent = Node(Atom(num=0)) # Dummy
        self.dfs_cycle(node, parent)

    # TODO Not working
    def dfs_cycle(self, node, parent):
        print('n', node.atom.num, 'with color', node.color ,'p', parent.atom.num)
        # already (completely) visited vertex.
        if node.color == 2:
            return

        # seen vertex, but was not completely visited -> cycle detected.
        # backtrack based on parents to find the complete cycle.
        if node.color == 1:

            self.number_of_cycles += 1
            cur = parent
            cur.ring.append(self.number_of_cycles)
            print('Start back track', cur.atom.num, cur.ring, 'parent', parent.atom.num)
            # backtrack the vertex which are in the current cycle that is found
            count = 1
            while cur != node:
                cur = cur.parents[-1]
                cur.ring.append(self.number_of_cycles)
                print('back track', cur.atom.num,'ring', cur.ring, 'p', [node.atom.num for node in cur.parents])
                count += 1
            return

        # Using overwrite function to avoid aappending
        node.add_parent(parent)

        #partially visited.
        node.color = 1


        for neighbor in node.neighbors:
            # if it has not been visited previously
            if neighbor not in node.parents:
                self.dfs_cycle(neighbor, node)

        # completely visited.
        node.color = 2

    def get_rings2(self, depth=8):
        max_depth = depth
        # Find the shortest circuit for each node
        rings = {}
        for node in self.nodes:
            if node.number_of_neighbor < 2:
                continue # Skip the leaves
            # initialize
            rings[node] = []
            searched_path = []
            node_i = node
            node_reff = node
            parent_i = Node(Atom(num=0)) # Dummy
            #print(node_i.atom.num, '-----------')
            #print('searchin ring for atom ', node_reff.atom.num )
            self.dfs_cycle2(node_reff, node_i, parent_i, searched_path, rings, max_depth)

        # Find the unique rings

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

    def dfs_cycle2(self, node_reff, node_i, parent_i, searched_path, rings, max_depth):
        max_depth -= 1
        # pass by value
        current_searched = searched_path.copy()
        current_searched.append(node_i)
        #print('n', node_i.atom.num, 'p', parent_i.atom.num, 'path', [node.atom.num for node in current_searched])
        for neighbor in node_i.neighbors:
            if neighbor == node_reff and neighbor != parent_i:
                found_ring = set(current_searched)
                if found_ring not in rings[node_reff]:
                    rings[node_reff].append(found_ring)
                    return

                """
                if len(current_searched) < len(rings[node_reff]):
                    #print('ring->', [node.atom.num for node in current_searched])
                    rings[node_reff] = current_searched
                    return
                """
            if neighbor not in current_searched and neighbor.number_of_neighbor > 1 and max_depth > 0:
                self.dfs_cycle2(node_reff, neighbor, node_i, current_searched, rings, max_depth)

    @staticmethod
    def get_similarity(mol1, mol2):
        # TODO make sure it is a tree

        score_max = 0
        score_node_pair_matrix = dict()

        for node_m2 in mol2.children:
            for node_m1 in mol1.children:

                score_node_pair_matrix[(node_m2, node_m1)] = 0
                # Check nodes
                if node_m1.atom.typ == node_m2.atom.typ:
                    score_node_pair_matrix[(node_m2, node_m1)] +=1
                    #print([i for i in score_node_pair_matrix.items()])
                    #print(node_m2.atom.num, node_m1.atom.num)
                # Check the sub trees
                score_node_pair_matrix[(node_m2, node_m1)] += Graph.get_similarity(node_m1, node_m2)

        #print([(key[0].atom.nam, key[-1].atom.nam, value) for key, value in score_node_pair_matrix.items()])
        #print(len([(key[0].atom.num, key[-1].atom.num, value) for key, value in score_node_pair_matrix.items()]))
        # Get the combination of nodes that maximize the score
        for permutation in it.permutations(mol1.children):
            # Calculate the score with this mapping
            score_tm = 0
            for node_pair in zip(mol2.children, permutation):
                score_tm += score_node_pair_matrix[node_pair]

            if score_tm > score_max:
                score_max = score_tm
        #print(score_max)
        return score_max


class atom_type_library(object):

    def __init__(self):
        # check the directory or path for updates
        pass

    def load(self):
        pass

    def write(self, library):
        pass

    def update(self):
        pass


class Bond(object):
    """
    A class for Bond terms
       - **Parm i** (:class:`Atom`):    atom i
       - **Parm j** (:class:`Atom`):    atom j
       - **Parm f** (:class:`float`):   force constant
       - **Parm r** (:class:`float`):  equilibrium distance

    >>> atom_i = Atom(nam='O1', typ='O', num=1)
    >>> atom_j = Atom(nam='H2', typ='H', num=2)
    >>> bondi = Bond(i=atom_i, j=atom_j, f=1200, r=1.4)
    >>> print('bond i: ', bondi.i.nam, bondi.j.nam, bondi.f, bondi.r)
    bond i:  O1 H2 1200.0 1.4
    """

    def __init__(self, i=Atom, j=Atom, f=0.0, r=0.0):
        self._i = None
        self._j = None
        # Call the setter to check the input type
        self.i = i
        self.j = j
        self._r = float(r)
        self._f = float(f)

    def __eq__(self, other):
        """
        Two bonds are equal if their atoms are equal. That is they have the same atom name and type.
        """
        if type(other) is Bond:
            return ((self._i == other.i) and (self._j == other.j)) or ((self._j == other.i) and (self._i == other.j))
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if type(other) is Bond:
            if self.i.num < other.i.num:
                return True
            elif (self.i.num == other.i.num) and (self.j.num < other.j.num):
                return True
        else:
            return NotImplemented

    def __gt__(self, other):
        return not self < other

    def __contains__(self, item):
        if type(item) is Atom:
            return (self.i == item) or (self.j == item)
        else:
            raise ValueError('Error >>> Only accept Atom object')

    @property
    def i(self):
        return self._i

    @i.setter
    def i(self, i):
        if type(i) is Atom:
            self._i = i
        else:
            raise ValueError('Pass an Atom object')

    @property
    def j(self):
        return self._j

    @j.setter
    def j(self, j):
        if type(j) is Atom:
            self._j = j
        else:
            raise ValueError('Pass an Atom object')

    @property
    def atoms(self):
        return [self._i, self._j]

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, k):
        self._f = float(k)

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        self._r = float(r)

    @property
    def show(self):
        return [self.i.num, self.j.num, self.f, self.r]


class Angle(object):
    """
    A class for Bond terms
       - **Parm i** (:class:`Atom`):    atom i
       - **Parm j** (:class:`Atom`):    atom j
       - **Parm k** (:class:`Atom`):   force constant
       - **Parm r** (:class:`float`):  equilibrium distance

    >>> atom_i = Atom(nam='O1', typ='O', num=1)
    >>> atom_j = Atom(nam='H2', typ='H', num=2)
    >>> atom_k = Atom(nam='H3', typ='H', num=3)
    >>> anglei = Angle(i=atom_i, j=atom_j, k=atom_k, f=1200, t=1.4)
    >>> print('Angle i: ', anglei.i.nam, anglei.j.nam, anglei.k.nam, anglei.f, anglei.t)
    Angle i:  O1 H2 H3 1200.0 1.4
    """

    def __init__(self, i=Atom, j=Atom, k=Atom, f=0.0, t=0.0):
        self._i = None
        self._j = None
        self._k = None
        # Call the setter to check the input type
        self.i = i
        self.j = j
        self._k = k
        self._t = float(t)
        self._f = float(f)

    def __eq__(self, other):
        """
        Two Angles are equal if their atoms are equal. That is, for each atom, they have the same atom number, name and type.
        """
        if type(other) is Angle:
            return ((self._i == other.i) and (self._j == other.j) and (self._k == other.k)) or \
                   ((self._k == other.i) and (self._j == other.j) and (self._i == other.k))
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if type(other) is Angle:
            if self.i.num < other.i.num:
                return True
            elif (self.i.num == other.i.num) and (self.j.num < other.j.num):
                return True
            elif (self.i.num == other.i.num) and (self.j.num == other.j.num) and (self.k.num < other.k.num):
                return True
        else:
            return NotImplemented

    def __gt__(self, other):
        return not self < other

    def __contains__(self, item):
        if type(item) is Atom:
            return (self._i == item) or (self._j == item) or (self._k == item)

        elif type(item) is Bond:
            return ((self._i == item.i) and (self._j == item.j) or (self._i == item.j) and (self._j == item.i) or
                    (self._j == item.i) and (self._k == item.j) or (self._j == item.j) and (self._k == item.i))
        else:
            raise ValueError('Error >>> Only accept Atom or Bond object')

    @property
    def i(self):
        return self._i

    @i.setter
    def i(self, i):
        if type(i) is Atom:
            self._i = i
        else:
            raise ValueError('Pass a Atom object')

    @property
    def j(self):
        return self._j

    @j.setter
    def j(self, j):
        if type(j) is Atom:
            self._j = j
        else:
            raise ValueError('Pass a Atom object')

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, k):
        if type(k) is Atom:
            self._k = k
        else:
            raise ValueError('Pass a Atom object')

    @property
    def atoms(self):
        return [self._i, self._j, self._k]

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, f):
        self._f = float(f)

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, t):
        self._t = float(t)

    @property
    def show(self):
        return [self.i.num, self.j.num, self.k.num, self.f, self.t]


class Torsion(object):
    def __init__(self, i=Atom, j=Atom, k=Atom, l=Atom, f1=0.0, f2=0.0, f3=0.0, t1=0.0, t2=0.0, t3=0.0, p1=0.0, p2=0.0,
                 p3=0.0):
        self._i = None
        self._j = None
        self._k = None
        self._l = None
        self._f = [0.0, 0.0, 0.0]
        self._t = [0.0, 0.0, 0.0]
        self._p = [0.0, 0.0, 0.0]

        # Call the setter to check the input type
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def __eq__(self, other):
        """
        Two Angles are equal if their atoms are equal. That is, for each atom, they have the same atom number, name and type.
        """
        if type(other) is Torsion:
            return ((self._i == other.i) and (self._j == other.j) and (self._k == other.k) and (self._l == other.l) or
                    (self._l == other.i) and (self._k == other.j) and (self._j == other.k) and (self._i == other.l))
        else:
            return NotImplemented

    def __ne__(self, other):
        return self == other

    def __lt__(self, other):
        if type(other) is Torsion:
            if self.i.num < other.i.num:
                return True
            elif (self.i.num == other.i.num) and (self.j.num < other.j.num):
                return True
            elif (self.i.num == other.i.num) and (self.j.num == other.j.num) and (self.k.num < other.k.num):
                return True
            elif (self.i.num == other.i.num) and (self.j.num == other.j.num) and (self.k.num == other.k.num) \
                    and (self.l.num < other.l.num):
                return True
        else:
            return NotImplemented

    def __gt__(self, other):
        return not self < other

    def __contains__(self, item):
        if type(item) is Atom:
            return (self._i == item) or (self._j == item) or (self._k == item) or (self._l == item)

        elif type(item) is Bond:
            return ((self._i == item.i) and (self._j == item.j) or (self._i == item.j) and (self._j == item.i) or
                    (self._j == item.i) and (self._k == item.j) or (self._j == item.j) and (self._k == item.i) or
                    (self._k == item.i) and (self._l == item.j) or (self._k == item.j) and (self._l == item.i))

        elif type(item) is Angle:
            return ((self._i == item.i) and (self._j == item.j) and (self._k == item.k) or
                    (self._i == item.k) and (self._j == item.j) and (self._k == item.i) or
                    (self._j == item.i) and (self._k == item.j) and (self._l == item.k) or
                    (self._j == item.k) and (self._k == item.j) and (self._l == item.i))
        else:
            raise ValueError('Error >>> Only accept Atom, Bond, or Angle object')

    @property
    def i(self):
        return self._i

    @i.setter
    def i(self, i):
        if type(i) is Atom:
            self._i = i
        else:
            raise ValueError('Pass a Atom object')

    @property
    def j(self):
        return self._j

    @j.setter
    def j(self, j):
        if type(j) is Atom:
            self._j = j
        else:
            raise ValueError('Pass a Atom object')

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, k):
        if type(k) is Atom:
            self._k = k
        else:
            raise ValueError('Pass a Atom object')

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, l):
        if type(l) is Atom:
            self._l = l
        else:
            raise ValueError('Pass a Atom object')

    @property
    def atoms(self):
        return [self._i, self._j, self._k, self._l]

    @property
    def n(self):
        if (len(self._f) == len(self._t)) and (len(self._f) == len(self._p)) and (len(self._p) == len(self._t)):
            return len(self._f)
        else:
            print("Warning >>> The number of torsion parameters (f, t, and phi) don't match")

    @property
    def f(self):
        return [f for f in self._f]

    @f.setter
    def f(self, f):
        try:
            for i in range(3):
                self._f[i] = float(f[i])
        except:
            raise ValueError('Pass a list of numbers')

    @property
    def f1(self):
        return self._f[0]

    @f1.setter
    def f1(self, f1):
        try:
            self._f[0] = float(f1)
        except:
            raise ValueError('Pass a number')

    @property
    def f2(self):
        return self._f[1]

    @f2.setter
    def f2(self, f2):
        try:
            self._f[1] = float(f2)
        except:
            raise ValueError('Pass a number')

    @property
    def f3(self):
        return self._f[2]

    @f3.setter
    def f3(self, f3):
        try:
            self._f[2] = float(f3)
        except:
            raise ValueError('Pass a number')

    @property
    def t(self):
        return [t for t in self._t]

    @t.setter
    def t(self, t):
        try:
            for i in range(3):
                self._t[i] = float(t[i])
        except:
            raise ValueError('Pass a list of numbers')

    @property
    def t1(self):
        return self._t[0]

    @t1.setter
    def t1(self, t1):
        try:
            self._t[0] = float(t1)
        except:
            raise ValueError('Pass a number')

    @property
    def t2(self):
        return self._t[1]

    @t2.setter
    def t2(self, t2):
        try:
            self._t[1] = float(t2)
        except:
            raise ValueError('Pass a number')

    @property
    def t3(self):
        return self._t[2]

    @t3.setter
    def t3(self, t3):
        try:
            self._t[2] = float(t3)
        except:
            raise ValueError('Pass a number')

    @property
    def p(self):
        return [p for p in self._p]

    @p.setter
    def p(self, p):
        try:
            for i in range(3):
                self._f[i] = float(p[i])
        except:
            raise ValueError('Pass a list of numbers')

    @property
    def p1(self):
        return self._p[0]

    @p1.setter
    def p1(self, p1):
        try:
            self._p[0] = float(p1)
        except:
            raise ValueError('Pass a number')

    @property
    def p2(self):
        return self._p[1]

    @p2.setter
    def p2(self, p2):
        try:
            self._p[1] = float(p2)
        except:
            raise ValueError('Pass a number')

    @property
    def p3(self):
        return self._p[2]

    @p3.setter
    def p3(self, p3):
        try:
            self._p[2] = float(p3)
        except:
            raise ValueError('Pass a number')

    @property
    def show(self):
        return [self.i.num, self.j.num, self.k.num, self.l.num, self.f, self.t, self.p]


class Improper(object):
    def __init__(self, i=Atom, j=Atom, k=Atom, l=Atom, f=0.0, t=0.0, phi=0.0):
        self._i = None
        self._j = None
        self._k = None
        self._l = None
        self._f = None
        self._t = None

        # Call the setter to check the input type
        self.i = i  # The pivot atom
        self.j = j
        self.k = k
        self.l = l
        self.f = f
        self.t = t

    def __eq__(self, other):
        """
        Two Improper are equal if they have the same pivot atoms (atom_i) that is connected to the same three atoms, regardless of their order.
        """
        if type(other) is Torsion:
            return ((self._i == other.i) and (self._j == other.j) and (self._k == other.k) and (self._l == other.l) or
                    (self._i == other.i) and (self._j == other.j) and (self._l == other.k) and (self._k == other.l) or
                    (self._i == other.i) and (self._k == other.j) and (self._j == other.k) and (self._l == other.l) or
                    (self._i == other.i) and (self._k == other.j) and (self._l == other.k) and (self._j == other.l) or
                    (self._i == other.i) and (self._l == other.j) and (self._k == other.k) and (self._j == other.l) or
                    (self._i == other.i) and (self._l == other.j) and (self._j == other.k) and (self._k == other.l))
        else:
            return NotImplemented

    def __ne__(self, other):
        return self == other

    def __lt__(self, other):
        if type(other) is Improper:
            if self.i.num < other.i.num:
                return True
            elif (self.i.num == other.i.num) and (self.j.num < other.j.num):
                return True
            elif (self.i.num == other.i.num) and (self.j.num == other.j.num) and (self.k.num < other.k.num):
                return True
            elif (self.i.num == other.i.num) and (self.j.num == other.j.num) and (self.k.num == other.k.num) \
                    and (self.l.num < other.l.num):
                return True
        else:
            return NotImplemented

    def __gt__(self, other):
        return not self < other

    def __contains__(self, item):
        if type(item) is Atom:
            return (self._i == item) or (self._j == item) or (self._k == item) or (self._l == item)

        elif type(item) is Bond:
            return ((self._i == item.i) and (self._j == item.j) or (self._i == item.j) and (self._j == item.i) or
                    (self._i == item.i) and (self._k == item.j) or (self._i == item.j) and (self._k == item.i) or
                    (self._i == item.i) and (self._l == item.j) or (self._i == item.j) and (self._l == item.i))

        elif type(item) is Angle:
            return ((self._j == item.i) and (self._i == item.j) and (self._k == item.k) or
                    (self._j == item.k) and (self._i == item.j) and (self._k == item.i) or
                    (self._j == item.i) and (self._i == item.j) and (self._l == item.k) or
                    (self._j == item.k) and (self._i == item.j) and (self._l == item.i) or
                    (self._k == item.i) and (self._i == item.j) and (self._l == item.k) or
                    (self._k == item.k) and (self._i == item.j) and (self._l == item.i))
        else:
            raise ValueError('Error >>> Only accept Atom, Bond, or Angle object')

    @property
    def i(self):
        return self._i

    @i.setter
    def i(self, i):
        if type(i) is Atom:
            self._i = i
        else:
            raise ValueError('Pass a Atom object')

    @property
    def j(self):
        return self._j

    @j.setter
    def j(self, j):
        if type(j) is Atom:
            self._j = j
        else:
            raise ValueError('Pass a Atom object')

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, k):
        if type(k) is Atom:
            self._k = k
        else:
            raise ValueError('Pass a Atom object')

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, l):
        if type(l) is Atom:
            self._l = l
        else:
            raise ValueError('Pass a Atom object')

    @property
    def atoms(self):
        return [self._i, self._j, self._k, self._l]

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, f):
        self._f = float(f)

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, t):
        self._t = float(t)

    @property
    def show(self):
        return [self.i.num, self.j.num, self.k.num, self.l.num, self.f, self.t]


class InternalCoordinate(Molecule):
    """
     A container class for internal coordinates
       - **Parm bonds**  (:class:`Bond`):
       - **Parm angles** (:class:`Angle`):
       - **Parm torsions** (:class:`Torsion`):
       - **Parm impropers** (:class:`Improper`):
    """
    Rad2Degree = 57.2957795
    ImproperAngleCutoff = 10.0  # degrees

    def __init__(self):
        Molecule.__init__(self, num=1)
        # TODO this should be in the Molecule class
        self._charge = None
        # TODO Fix the errors in FF_bonded
        # self._structure = Molecule(num=1)
        self._bonds = []
        self._angles = []
        self._torsions = []
        self._impropers = []
        self._graph = []

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = float(charge)

    # @property
    # def natm(self):
    #    return self._structure.natm

    # @property
    # def structure(self):
    #    return self._structure

    # @structure.setter
    # def structure(self, struc):
    #    if isinstance(struc, Molecule):
    #        self._structure = struc
    #    else:
    #        raise ValueError('Error >>> Pass a Molecule object')

    @property
    def nbond(self):
        return len(self._bonds)

    @property
    def bonds(self):
        return self._bonds

    @bonds.setter
    def bonds(self, bonds_list):
        if any([type(bond) != Bond for bond in bonds_list]):
            raise ValueError('Pass a list of bond objects')
        else:
            self._bonds = bonds_list

    def addbond(self, bond):
        if type(bond) is Bond:
            self._bonds.append(bond)
        else:
            raise ValueError('Error >>> Pass a bond object')

    def addbonds(self, bonds):
        for bond in bonds:
            if type(bond) is Bond:
                self._bonds.append(bond)
            else:
                raise ValueError('Error >>> Pass a list of bond objects')

    def deletebond(self, bond):
        if type(bond) is Bond:
            self._bonds.remove(bond)
        else:
            raise ValueError('Error >>> Pass a bond object')

    def deletebonds(self, bonds=None):
        if bonds == None:
            self._bonds = []
        elif type(bonds) is list:
            for bond in bonds:
                if type(bond) is Bond:
                    self._bonds.remove(bond)
                else:
                    raise ValueError('Error >>> Pass a list of bond objects')
        else:
            raise ValueError('Error >>> Pass a list of bond objects')

    def selectbond(self, bond):
        if type(bond) is Bond:
            return [ibond for ibond in self._bonds if ibond == bond]

    def selectbondbyAtom(self, atom):
        return [ibond for ibond in self._bonds if atom in ibond]

    def constructbondsbyCovalentRadius(self):
        if self.natm > 1:
            # Construct the bonds
            for i in range(1, self.natm + 1):
                # get atom i
                atom_i = self.structure.selectbyAtomnum(i)
                ui = np.array(atom_i.cord)

                for j in range(i + 1, self.natm + 1):

                    # get atom j
                    atom_j = self.structure.selectbyAtomnum(j)
                    uj = np.array(atom_j.cord)

                    # Calculate the bond distance
                    dij = np.linalg.norm(uj - ui)

                    # if the bond is not already in the list, add it.
                    if dij <= bond_lenght(atom_i.typ, atom_j.typ):
                        ibond = Bond(i=atom_i, j=atom_j, r=dij)
                        if ibond not in self.bonds:
                            self.addbond(ibond)

            # Sort the bonds
            self.bonds.sort()

        else:
            print("Error >>> at least 2 atoms are needed to construct bond terms. Number of atoms are: ", self.natm)

    def constructbondsbyHessian(self, hessian):
        if self.natm > 1:
            h = hessian
            # Construct the bonds
            for i in range(1, self.natm + 1):
                lbound_i = (i - 1) * 3
                hbound_i = (lbound_i + 3)

                for j in range(i + 1, self.natm + 1):
                    lbound_j = (j - 1) * 3
                    hbound_j = (lbound_j + 3)

                    # Get the partial hessian
                    hij = h[lbound_i:hbound_i, lbound_j:hbound_j] * (-1)

                    # Get the Eigen values
                    wij = np.linalg.eigvals(hij)

                    # assign bonds with positive Eigen values
                    if (np.isreal(wij).all() and (wij > 0).all()):

                        # Calculate the bond distance
                        uij = np.linalg.norm(np.array(self.structure.selectbyAtomnum(j).cord) - np.array(
                            self.structure.selectbyAtomnum(i).cord))
                        ibond = Bond(i=self.structure.selectbyAtomnum(i), j=self.structure.selectbyAtomnum(j), r=uij)
                        if ibond not in self.bonds:
                            self.addbond(ibond)

            # Sort the bonds
            self.bonds.sort()

        else:
            print("Error >>> at least 2 atoms are needed to construct bond terms. Number of atoms are: ", self.natm)

    @property
    def nangle(self):
        return len(self._angles)

    @property
    def angles(self):
        return self._angles

    @angles.setter
    def angles(self, angles_list):
        if any([type(angle) != Angle for angle in angles_list]):
            raise ValueError('Pass a list of bond objects')
        else:
            self._angles = angles_list

    def addangle(self, angle):
        if type(angle) is Angle:
            self._angles.append(angle)
        else:
            raise ValueError('Error >>> Pass a Angle object')

    def addangles(self, angles):
        for angle in angles:
            if type(angle) is Angle:
                self._angles.append(angle)
            else:
                raise ValueError('Error >>> Pass a list of Angle objects')

    def deleteangle(self, angle):
        if type(angle) is Angle:
            self._angles.remove(angle)
        else:
            raise ValueError('Error >>> Pass an angle object')

    def deleteangles(self, angles=None):
        if angles == None:
            self._angles = []
        elif type(angles) is list:
            for angle in angles:
                if type(angle) is Angle:
                    self._angles.remove(angle)
                else:
                    raise ValueError('Error >>> Pass a list of angle objects')
        else:
            raise ValueError('Error >>> Pass a list of angle objects')

    def selectangle(self, angle):
        if type(angle) is Angle:
            return [iangle for iangle in self._angles if iangle == angle]
        else:
            raise ValueError('Error >>> Pass a list of Angle object')

    def selectanglebyAtom(self, atom):
        return [iangle for iangle in self._angles if atom in iangle]

    def selectanglebyPivotAtom(self, atom):
        return [iangle for iangle in self._angles if atom == iangle._j]

    def selectanglebyBond(self, bond):
        return [iangle for iangle in self._angles if bond in iangle]

    def constructangles(self):
        if self.nbond > 1:
            self._bonds.sort()
            for atom_i in self.structure.atoms:
                for atom_j in self.structure.atoms:
                    for atom_k in self.structure.atoms:
                        if (Bond(i=atom_i, j=atom_j) in self.bonds) and (
                                    Bond(i=atom_j, j=atom_k) in self.bonds) and atom_i != atom_k:

                            angle = Angle(i=atom_i, j=atom_j, k=atom_k)
                            if angle not in self.angles:
                                # Get the unit vectors
                                uji = np.array(atom_i.cord) - np.array(atom_j.cord)
                                uji /= np.linalg.norm(uji)

                                ujk = np.array(atom_k.cord) - np.array(atom_j.cord)
                                ujk /= np.linalg.norm(ujk)

                                # Get the angle
                                t = np.arccos(np.dot(uji, ujk))
                                t *= self.Rad2Degree

                                angle.t = t
                                self.addangle(angle)
            self.angles.sort()
        else:
            print("Error >>> at least bonds are needed to construct angle terms. Number of bonds are: ", self.nbond)

    @property
    def ntorsion(self):
        return len(self._torsions)

    @property
    def torsions(self):
        return self._torsions

    @torsions.setter
    def torsions(self, torsion_list):
        if any([type(torsion) != Torsion for torsion in torsion_list]):
            raise ValueError('Pass a list of torsion objects')
        else:
            self._torsions = torsion_list

    def addtorsion(self, torsion):
        if type(torsion) is Torsion:
            self._torsions.append(torsion)
        else:
            raise ValueError('Error >>> Pass a Torsion object')

    def addtorsions(self, torsions):
        for torsion in torsions:
            if type(torsion) is Torsion:
                self._torsions.append(torsion)
            else:
                raise ValueError('Error >>> Pass a list of Angle objects')

    def deletetorsion(self, torsion):
        if type(torsion) is Torsion:
            self._torsions.remove(torsion)
        else:
            raise ValueError('Error >>> Pass a torsion object')

    def deletetorsions(self, torsions=None):
        if torsions == None:
            self._torsions = []
        elif type(torsions) is list:
            for torsion in torsions:
                if type(torsion) is Torsion:
                    self._torsions.remove(torsion)
                else:
                    raise ValueError('Error >>> Pass a list of torsion objects')
        else:
            raise ValueError('Error >>> Pass a list of torsion objects')

    def selecttorsion(self, torsion):
        if type(torsion) is Torsion:
            return [itorsion for itorsion in self._torsions if itorsion == torsion]
        else:
            raise ValueError('Error >>> Pass a list of Angle object')

    def selecttorsionbyAtom(self, atom):
        return [itorsion for itorsion in self._torsions if atom in itorsion]

    def selecttorsionbyPivotBond(self, bond):
        return [itorsion for itorsion in self._torsions if ((bond.i == itorsion._j and bond.j == itorsion._k) or
                                                            (bond.i == itorsion._k and bond.j == itorsion._l))]

    def selecttorsionbyBond(self, bond):
        return [itorsion for itorsion in self._torsions if bond in itorsion]

    def selecttorsionbyAngle(self, angle):
        return [itorsion for itorsion in self._torsions if angle in itorsion]

    def constructtorsions(self):

        if self.nangle > 1:
            self._angles.sort()
            for atom_i in self.structure.atoms:
                for atom_j in self.structure.atoms:
                    for atom_k in self.structure.atoms:
                        for atom_l in self.structure.atoms:
                            if (Angle(i=atom_i, j=atom_j, k=atom_k) in self.angles) and \
                                    (Angle(i=atom_j, j=atom_k, k=atom_l) in self.angles) and atom_i != atom_l:

                                torsion = Torsion(i=atom_i, j=atom_j, k=atom_k, l=atom_l)
                                if torsion not in self.torsions:
                                    # Get the unit vectors uij, uik, ujk, ujl
                                    uij = np.array(atom_j.cord) - np.array(atom_i.cord)
                                    uij /= np.linalg.norm(uij)

                                    uik = np.array(atom_k.cord) - np.array(atom_i.cord)
                                    uik /= np.linalg.norm(uik)

                                    ujk = np.array(atom_k.cord) - np.array(atom_j.cord)
                                    ujk /= np.linalg.norm(ujk)

                                    ujl = np.array(atom_l.cord) - np.array(atom_j.cord)
                                    ujl /= np.linalg.norm(ujl)

                                    # Get the unit vectors angles tkij, tljk
                                    tkij = np.arccos(np.dot(uik, uij))
                                    tljk = np.arccos(np.dot(ujl, ujk))

                                    # Get the normals
                                    uijk = np.cross(uij, uik) / np.sin(tkij)
                                    ujkl = np.cross(ujk, ujl) / np.sin(tljk)

                                    # Get the torsion angle
                                    tijkl = np.sign(np.dot(uij, ujkl)) * np.arccos(np.dot(uijk, ujkl))

                                    torsion.t1 = tijkl * self.Rad2Degree
                                    self.addtorsion(torsion)

            self.torsions.sort()
        else:
            print("Error >>> at least 2 angles are needed to construct torsion terms. Number of angles are: ",
                  self.nangle)

    @property
    def nimproper(self):
        return len(self._impropers)

    @property
    def impropers(self):
        return self._impropers

    @impropers.setter
    def impropers(self, improper_list):
        if any([type(improper) != Improper for improper in improper_list]):
            raise ValueError('Pass a list of improper objects')
        else:
            self._impropers = improper_list

    def addimproper(self, improper):
        if type(improper) is Improper:
            self._impropers.append(improper)
        else:
            raise ValueError('Error >>> Pass a Improper object')

    def addimpropers(self, impropers):
        for improper in impropers:
            if type(improper) is Improper:
                self._impropers.append(improper)
            else:
                raise ValueError('Error >>> Pass a list of Improper objects')

    def deleteimproper(self, improper):
        if type(improper) is Improper:
            self._impropers.remove(improper)
        else:
            raise ValueError('Error >>> Pass an improper object')

    def deleteimpropers(self, impropers=None):
        if impropers == None:
            self._impropers = []
        elif type(impropers) is list:
            for improper in impropers:
                if type(improper) is Improper:
                    self._impropers.remove(improper)
                else:
                    raise ValueError('Error >>> Pass a list of impropers objects')
        else:
            raise ValueError('Error >>> Pass a list of torsion objects')

    def selectimproper(self, improper):
        if type(improper) is Improper:
            return [iimproper for iimproper in self._impropers if iimproper == improper]
        else:
            raise ValueError('Error >>> Pass a list of Improper object')

    def selectimproperbyAtom(self, atom):
        return [iimproper for iimproper in self._impropers if atom in iimproper]

    def selectimproperbyPivotAtom(self, atom):
        return [iimproper for iimproper in self._impropers if (atom == iimproper.i)]

    def selectimproperbyBond(self, bond):
        return [iimproper for iimproper in self._impropers if bond in iimproper]

    def selectimproperbyAngle(self, angle):
        return [iimproper for iimproper in self._impropers if angle in iimproper]

    def constructimpropers(self):

        if self.nangle > 1:
            self._angles.sort()
            for atom_i in self.structure.atoms:

                # Get all angles with the same middle atom
                angle_list = self.selectanglebyPivotAtom(atom_i)

                # Get the potential impropers with the give pivot atom.
                if len(angle_list) == 3:
                    atom_list = []
                    for angle in angle_list:
                        for atom in angle.atoms:
                            if atom != atom_i and atom not in atom_list:
                                atom_list.append(atom)
                    atom_list.sort()

                    # Calculate the improper angles
                    # Get the unit vectors uij, uik, ujk, ujl
                    uij = np.array(atom_list[0].cord) - np.array(atom_i.cord)
                    uij /= np.linalg.norm(uij)

                    uik = np.array(atom_list[1].cord) - np.array(atom_i.cord)
                    uik /= np.linalg.norm(uik)

                    ujk = np.array(atom_list[1].cord) - np.array(atom_list[0].cord)
                    ujk /= np.linalg.norm(ujk)

                    ujl = np.array(atom_list[2].cord) - np.array(atom_list[1].cord)
                    ujl /= np.linalg.norm(ujl)

                    # Get the unit vectors angles tkij, tljk
                    tkij = np.arccos(np.dot(uik, uij))
                    tljk = np.arccos(np.dot(ujl, ujk))

                    # Get the normals
                    uijk = np.cross(uij, uik) / np.sin(tkij)
                    ujkl = np.cross(ujk, ujl) / np.sin(tljk)

                    # Get the torsion angle
                    tijkl = np.sign(np.dot(uij, ujkl)) * np.arccos(np.dot(uijk, ujkl))
                    t = tijkl * self.Rad2Degree

                    if np.abs(t) <= self.ImproperAngleCutoff:
                        self.addimproper(Improper(i=atom_i, j=atom_list[0], k=atom_list[1], l=atom_list[2], t=t))

        else:
            print("Error >>> at least 2 angles are needed to construct improper terms. Number of angles are: ",
                  self.nangle)

    @property
    def graph(self):
        return self._graph

    def construct_graph(self):
        # Initialize
        self._graph = None
        self._graph = Graph()

        # Populate the nodes
        self._graph.nodes = [Node(atom) for atom in self.atoms]

        # Find the neighbor list of each node
        for node_i in self._graph.nodes:
            for node_j in self._graph.nodes:
                if Bond(node_i.atom, node_j.atom) in self.bonds:
                    node_i.neighbors = node_j

    def tree(self, head):
        if not self.graph:
            print('Error, >>> To construct a tree, first a graph should be made')
            return None

        # Reset the head node
        for node in self.graph.nodes:
            if node.atom == head:
                self.graph.head = node
                break

        if not self.graph.head:
            print('Error, >>> head atom was not found')
            return None

            # Make the tree
            # for node in self.graph.nodes:
            # print(node.atom.num,'->', [n.atom.num for n in node.neighbors])

        self.graph.make_tree_BFS(self.graph.head)

        # Make the tree
        count_cycles = 0
        for node in self.graph.nodes:
            print([p.atom.num for p in node.parents], '->', node.atom.num, '->', [c.atom.num for c in node.children])
            if node.number_of_parents > 1 :
                count_cycles +=1
        print('------------------------------------','Number of rings', count_cycles, '------------------------------------')

        #self.graph.make_tree_DFS(self.graph.head)
        # Make the tree
        #for node in self.graph.nodes:
        #    print([p.atom.num for p in node.parents], '->', node.atom.num, '->', [c.atom.num for c in node.children])

if __name__ == '__main__':

    a = Atom(num=1, nam='H')
    b = Atom(num=2, nam='O')
    ab = Bond(i=a, j=b)
    abcd = Torsion(i=a, j=b, k=a, l=b, f1=0.0, f2=1.0, t1=0.0, t2=1.0, p1=0.0, p2=1.0)

    tsuv = Torsion(i=a, j=b, k=a, l=b, f1=5., t1=1., p1=180)

    if (ab.i == ab.j):
        print('YES')

    print(tsuv.f)
    print(tsuv.t)
    print(tsuv.p)
    print(tsuv.n)
    print(tsuv.i.show)
