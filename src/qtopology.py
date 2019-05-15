#!/usr/bin/env python3
"""
.. module:: qtopology
   :platform: Unix

"""
from atom import Atom
from molecule import Molecule
from macromolecule import Macromolecule


class Bond(object):
    """
    A class for Bond terms
       - **Parm i** (:class:`Atom`):    atom i
       - **Parm j** (:class:`Atom`):    atom j
       - **Parm f** (:class:`float`):   force constant
       - **Parm r0** (:class:`float`):  equilibrium distance

    >>> atom_i = Atom(nam='O1', typ='O', num=1)
    >>> atom_j = Atom(nam='H2', typ='H', num=2)
    >>> bondi = Bond(i=atom_i, j=atom_j, f=1200, r0=1.4)
    >>> print('bond i: ', bondi.i.nam, bondi.j.nam, bondi.f, bondi.r0)
    bond i:  O1 H2 1200.0 1.4
    """

    def __init__(self, i=None, j=None, f=None, r0=None):
        self._i = None
        self._j = None
        self._r0 = float(r0)
        self._f = float(f)
        if i is not None:
            self.i = i
        if j is not None:
            self.j = j

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

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self, k):
        self._f = float(k)

    @property
    def r0(self):
        return self._r0

    @r0.setter
    def r0(self, r0):
        return self._r0


class Angle(object):
    """
    A class for Bond terms
       - **Parm i** (:class:`Atom`):    atom i
       - **Parm j** (:class:`Atom`):    atom j
       - **Parm
       - **Parm k** (:class:`float`):   force constant
       - **Parm r0** (:class:`float`):  equilibrium distance

    >>> atom_i = Atom(nam='O1', typ='O', num=1)
    >>> atom_j = Atom(nam='H2', typ='H', num=2)
    >>> atom_k = Atom(nam='H3', typ='H', num=3)
    >>> anglei = Angle(i=atom_i, j=atom_j, k=atom_k, f=1200, t0=1.4)
    >>> print('Angle i: ', anglei.i.nam, anglei.j.nam, anglei.k.nam, anglei.f, anglei.t0)
    Angle i:  O1 H2 H3 1200.0 1.4
    """
    def __init__(self, i=None, j=None, k=None, f=None, t0=None):
        self._i = None
        self._j = None
        self._k = None
        self._f = float(f)
        self._t0 = float(t0)
        if i is not None:
            self.i = i
        if j is not None:
            self.j = j
        if k is not None:
            self.k = k


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
    def f(self):
        return self._f

    @f.setter
    def f(self, f):
        self._f = f

    @property
    def t0(self):
        return self._t0

    @t0.setter
    def t0(self, t0):
        self._t0 = t0