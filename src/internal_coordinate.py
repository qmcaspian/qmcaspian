#!usr/bin/python3
"""
    module :: internal_coordinate
    :platform: Unix
"""
from atom import Atom
from molecule import Molecule
import numpy as np
from atom_table import *


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


class InternalCoordinate(object):
    """
     A container class for internal coordinates
       - **Parm bonds**  (:class:`Bond`):
       - **Parm angles** (:class:`Angle`):
       - **Parm torsions** (:class:`Torsion`):
       - **Parm impropers** (:class:`Improper`):
    """

    def __init__(self):
        self._charge = None
        self._structure = Molecule(num=1, nam='Freq Calculation')
        self._bonds = []
        self._angles = []
        self._torsions = []
        self._impropers = []

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = float(charge)

    @property
    def natm(self):
        return self._structure.natm

    @property
    def structure(self):
        return self._structure

    @property
    def nbond(self):
        return len(self._bonds)

    @property
    def bonds(self):
        return self._bonds

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
                        uij = np.linalg.norm(np.array(self.structure.selectbyAtomnum(j).cord) - np.array(self.structure.selectbyAtomnum(i).cord))
                        ibond = Bond(i=self.structure.selectbyAtomnum(i), j=self.structure.selectbyAtomnum(j), r=uij)
                        if ibond not in self.bonds:
                            self.addbond(ibond)

            # Sort the bonds
            self.bonds.sort()

        else:
            print("Error >>> at least 2 atoms are needed to construct bond terms. Number of atoms are: ", self.natm)

    def constructbondsbyMol2(self, ifile):
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
                        uij = np.linalg.norm(np.array(self.structure.selectbyAtomnum(j).cord) - np.array(self.structure.selectbyAtomnum(i).cord))
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

    @property
    def ntorsion(self):
        return len(self._angles)

    @property
    def torsions(self):
        return self._torsions

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

    @property
    def nimproper(self):
        return len(self._impropers)

    @property
    def impropers(self):
        return self._impropers

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


if __name__ == '__main__':

    a = Atom(num=1, nam='H')
    b = Atom(num=2, nam='O')
    ab = Bond(i=a, j=b)
    abcd = Torsion(i=a, j=b, k=a, l=b, f=[0.0, 1.0], t=[0.0, 1.0], phi=[0.0, 1.0])

    tsuv = Torsion(i=a, j=b, k=a, l=b, f=5., t=1., phi=180)

    if (ab.i == ab.j):
        print('YES')

    print(tsuv.f)
    print(tsuv.t)
    print(tsuv.phi)
    print(tsuv.n)
    print(tsuv.i.show)
