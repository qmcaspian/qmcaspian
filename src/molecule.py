#!/usr/bin/python
"""
.. module:: atom
   :platform: Unix

"""

from atom import Atom
from copy import deepcopy


class Molecule(object):
    """
    A class to store a collection of :class:`Atom` objects. This class is iterable.

       - **Parm num** (:class:`int`):      Molecule number.
       - **Parm nam** (:class:`str`):      Molecule name.
       - **Parm atms** (:class:`Atom`):    List of atoms in the molecule.
       - **Parm natm** (:class:`int`):     Number of atoms in the molecule
       - **Parm head** (:class:`Atom`):    The head atom in the molecule. Used for connecting residues.
       - **Parm tail** (:class:`Atom`):    The tail atom in the molecule. Used for connecting residues.
       - **Parm pos** (:class:`str`):      The position of the molecule in a macromolecule (FIRST or LAST). Used for connecting residues.
       - **Parm chainnum** (:class:`int`): Molecule chain number.
       - **Parm chainnam** (:class:`str`): Molecule chain name

    .. note::
       Initializing a :class:`Molecule` object with arguments is optional. It is possible to initialize a molecule object
       by one or a list of atom object. Number of atoms (*natm*) is set automatically.

       >>> O  = Atom()
       >>> O4 = [Atom(num=1), Atom(num=2), Atom(num=3), Atom(num=4)]
       >>> molecule_O  = Molecule(num=1, atms=O)
       >>> molecule_O4 = Molecule(num=1, atms=O4)
       >>> for index, atom in enumerate(molecule_O4):
       ...  print(index, '-', atom.num)
       0 - 1
       1 - 2
       2 - 3
       3 - 4
    """

    def __init__(self, num=0, nam=None, atms=None, head=None, tail=None, pos=None, chainnum=None, chainnam=None):

        self._num = int(num)
        self._nam = str(nam)
        self._natm = 0
        self._pos = str(pos)
        self._chainnam = str(chainnam)
        self._chainnum = str(chainnum)
        self._head = None
        self._tail = None
        self._atms = []
        if head is not None: self.head = head
        if tail is not None: self.tail = tail
        if atms is not None: self.addatm(atms)

    def __iter__(self):
        for atom in self._atms:
            yield atom

    def __contains__(self, item):
        if item in self._atms: return True

    def __eq__(self, other):
        if isinstance(other, Molecule):
            return (self._num == other._num) and (self._nam == other._nam) and (self._natm == other._natm) and \
                   (self._chainnam == other._chainnam) and (self._chainnum == other._chainnum) and \
                   (self._atms == other._atms)
        else:
            return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def addatm(self, atm):
        """
        Appends atom(s) object to a molecule

        >>> H1, O = Atom(typ='H'), Atom(typ='O')
        >>> H2 = H1.copy
        >>> HOH = Molecule(num=1,nam='HOH')
        >>> HOH.addatm(H1)
        >>> HOH.addatm(H2)
        >>> HOH.show
        [1, 'HOH', 2]
        >>> HOH.addatm(O)
        >>> HOH.show
        [1, 'HOH', 3]

        .. note::
            *addatm* appends a hard copy of input variables.

        >>> H1 = Atom(num=1,typ='H')
        >>> H2 = Atom(num=1,typ='H')
        >>> HOH = Molecule(num=1,nam='HH')
        >>> HOH.addatm([H1, H2])
        >>> print([[atom.num, atom.typ] for atom in HOH])
        [[1, 'H'], [1, 'H']]
        >>> H1.num = 500
        >>> print([[atom.num, atom.typ] for atom in HOH])
        [[1, 'H'], [1, 'H']]
        """
        if type(atm) is Atom:
            self._atms.append(atm.copy)
            self._natm = self._natm + 1
        elif type(atm) is list:
            for x in atm:
                if type(x) is Atom:
                    self._atms.append(x.copy)
                    self._natm = self._natm + 1
                else:
                    raise ValueError("One element in the list is not Atom type")
        else:
            raise ValueError("Pass an Atom object")

    @property
    def num(self):
        """
        Sets/Returns Molecule number.

        >>> HOH = Molecule()
        >>> HOH.num
        0
        >>> HOH.num = 12
        >>> HOH.num
        12
        """
        return self._num

    @num.setter
    def num(self, num):
        try:
            self._num = int(num)
        except:
            raise ValueError("Molecule number could not be converted to integer")

    @property
    def nam(self):
        """
        Sets/Returns Molecule number.

        >>> HOH = Molecule()
        >>> HOH.nam

        >>> HOH.nam = 'HOH'
        >>> HOH.nam
        'HOH'
        """
        return self._nam

    @nam.setter
    def nam(self, nam):
        try:
            self._nam = str(nam)
        except:
            raise ValueError("Molecule name could not be converted to str")

    @property
    def natm(self):
        """
        Returns number of atoms in a molecule.

            >>> O4 = [Atom(), Atom(), Atom(), Atom()]
            >>> molecule_O4 = Molecule(num=1, atms=O4)
            >>> molecule_O4.natm
            4
        """
        return self._natm

    @property
    def show(self):
        """
        Returns num, nam, natm.

        >>> H, O = Atom(typ='H'), Atom(typ='O')
        >>> HOH = Molecule(num=1,nam='HOH', atms=[H, H, O])
        >>> HOH.show
        [1, 'HOH', 3]
        """
        return [self._num, self._nam, self._natm]

    def selectbyAtomnam(self, name):
        """
        Returns a :class:`Atom` or a list of :class:`Atom` with the given :class:`Atom.nam`. If only one :class:`Atom`,
        was found and :class:`Atom` object is returned. Otherwise a list of :class:`Atom` objects is returned.

        >>> ATOMS = [Atom(nam='HC'), Atom(nam='HC'), Atom(nam='O')]
        >>> HOH = Molecule(num=1, nam='HOH', atms=ATOMS)
        >>> select = HOH.selectbyAtomnam('HC')
        >>> [atom.show for atom in select]
        [[0, 'HC', 0.0, 0.0, 0.0, 'ND'], [0, 'HC', 0.0, 0.0, 0.0, 'ND']]
        """
        select = [atom for atom in self if atom.nam == name]
        if len(select) is 1:
            "If only one atom found returns the Atom"
            return select[0]
            " Return a list"
        else:
            return select

    def selectbyAtomnum(self, number):
        """
        Returns a :class:`Atom` or a list of :class:`Atom` with the given :class:`Atom.num`. If only one :class:`Atom`,
        was found and :class:`Atom` object is returned. Otherwise a list of :class:`Atom` objects is returned.

        >>> ATOMS = [Atom(num=1), Atom(num=2), Atom(num=3)]
        >>> HOH = Molecule(num=1, nam='HOH', atms=ATOMS)
        >>> select = HOH.selectbyAtomnum(1)
        >>> select.show
        [1, 'ND', 0.0, 0.0, 0.0, 'ND']
        """
        select = [atom for atom in self if atom.num == number]
        if len(select) is 1:
            return select[0]
        else:
            return select

    @property
    def head(self):
        """
        Sets/Returns Molecule head atom. The :class:`Molecule.head` accepts an :class:`Atom`, a *str*, or an *int*. If
        an :class:`Atom` object is given, it should exist in the molecule otherwise a ValueError is raised.
        If a *str* or an *int* is given, the atoms in the molecule are searched for a unique match and
        :class:`Molecule.head` points to that. If multiple values was found an ValueError exception is raised.

        >>> ATOMS = [Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='O')]
        >>> HOH = Molecule(num=1, nam='HOH', atms=ATOMS)
        >>> print(HOH.head)
        None
        >>> HOH.head = HOH.selectbyAtomnum(1)
        >>> HOH.head.show
        [1, 'H1', 0.0, 0.0, 0.0, 'ND']
        >>> print((HOH.head) is HOH.selectbyAtomnum(1))
        True

        Here the head function tries to find an :class:`Atom` that is equal to the input and assigns the head to that.
        In other word, the head :class:`Atom` always points to an atom that belongs to the molecule, even if, the input
        :class:`Atom` does not belong to the molecule.

        >>> HOH.head = ATOMS[0]
        >>> print((HOH.head) is HOH.selectbyAtomnum(1))
        True
        >>> print(HOH.head == ATOMS[0])
        True
        >>> print(HOH.head is ATOMS[0])
        False
        >>> HOH.head = 'H2'
        >>> HOH.head.show
        [2, 'H2', 0.0, 0.0, 0.0, 'ND']
        >>> HOH.head = 3
        >>> HOH.head.show
        [3, 'O', 0.0, 0.0, 0.0, 'ND']
        >>> try:
        ...     HOH.head = 23
        ... except Exception as error:
        ...     print(error)
        The input Atom number either was not found or has multiple occurrences
        """
        return self._head

    @head.setter
    def head(self, headAtom):

        if type(headAtom) is Atom:
            select = [atom for atom in self if atom == headAtom]
            if len(select) == 1:
                self._head = select[0]
            else:
                raise ValueError("The input Atom either was not found or has multiple occurrences")

        elif type(headAtom) is str:
            select = self.selectbyAtomnam(headAtom)
            if type(select) is Atom:
                self._head = select
            else:
                raise ValueError("The input Atom name either was not found or has multiple occurrences")

        elif type(headAtom) is int:
            select = self.selectbyAtomnum(headAtom)
            if type(select) is Atom:
                self._head = select
            else:
                raise ValueError("The input Atom number either was not found or has multiple occurrences")

        else:
            raise ValueError("Input should be either Atom, str, or int")

    @property
    def tail(self):
        """
        Sets/Returns Molecule tail atom. The :class:`Molecule.tail` accepts an :class:`Atom`, a *str*, or an *int*. If
        an :class:`Atom` object is given, it should exist in the molecule otherwise a ValueError is raised (see :class:`Molecule.head`).
        If a *str* or an *int* is given, the atoms in the molecule are searched for a unique match and
        :class:`Molecule.head` points to that. If multiple values was found an ValueError exception is raised.

        >>> ATOMS = [Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='O')]
        >>> HOH = Molecule(num=1, nam='HOH', atms=ATOMS)
        >>> print(HOH.tail)
        None
        >>> HOH.tail = HOH.selectbyAtomnum(1)
        >>> HOH.tail.show
        [1, 'H1', 0.0, 0.0, 0.0, 'ND']
        >>> HOH.tail = 'H2'
        >>> HOH.tail.show
        [2, 'H2', 0.0, 0.0, 0.0, 'ND']
        >>> HOH.tail = 3
        >>> HOH.tail.show
        [3, 'O', 0.0, 0.0, 0.0, 'ND']
        >>> try:
        ...     HOH.tail = 23
        ... except Exception as error:
        ...     print(error)
        The input Atom number either was not found or has multiple occurrences
        """
        return self._tail

    @tail.setter
    def tail(self, tailAtom):
        if type(tailAtom) is Atom:
            select = [atom for atom in self if atom == tailAtom]
            if len(select) == 1:
                self._tail = tailAtom
            else:
                raise ValueError("The input Atom should belong to the current Molecule")

        elif type(tailAtom) is str:
            select = self.selectbyAtomnam(tailAtom)
            if type(select) is Atom:
                self._tail = select
            else:
                raise ValueError("The input Atom name either was not found or has multiple occurrences")

        elif type(tailAtom) is int:

            select = self.selectbyAtomnum(tailAtom)
            if type(select) is Atom:
                self._tail = select
            else:
                raise ValueError("The input Atom number either was not found or has multiple occurrences")

        else:
            raise ValueError("Input should be either Atom, str, or int")

    @property
    def copy(self):
        """
        Returns a deepcopy of :class:`Molecule` object.

        >>> HOH1 = Molecule(num=1, nam='HOH')
        >>> HOH2 = HOH1.copy
        >>> HOH2.show
        [1, 'HOH', 0]
        >>> HOH2.num = 2
        >>> HOH2.show
        [2, 'HOH', 0]
        >>> HOH1.show
        [1, 'HOH', 0]
        """
        return deepcopy(self)

    @property
    def pos(self):
        """
        Sets/Returns the position of the molecule in a :class:`Macromolecule`. Could be set to "FIRST" or "LAST" for
        first or last residue in a chain.

        >>> HOH = Molecule(Atom(num=1), Atom(num=2), Atom(num=3))
        >>> HOH.pos
        >>>
        >>> HOH.pos = 'FIRST'
        >>> if HOH.pos is 'FIRST': print('First residue')
        First residue
        """
        return self._pos

    @pos.setter
    def pos(self, position):
        try:
            self._pos = str(position)
        except:
            raise ValueError("Molecule position could not be converted to string")

    @property
    def chainnam(self):
        """
        Returns/Sets chain name that molecule belongs to.

        >>> HOH = Molecule()
        >>> HOH.chainnam
        >>>
        >>> HOH.chainnam = 'A'
        >>> HOH.chainnam
        'A'
        """
        return self._chainnam

    @chainnam.setter
    def chainnam(self, chainName):
        try:
            self._chainnam = str(chainName)
        except:
            raise ValueError("Molecule chain name could not be converted to string")

    @property
    def chainnum(self):
        """
        Returns/Sets chain number that molecule belongs to.

        >>> HOH = Molecule()
        >>> HOH.chainnum
        >>>
        >>> HOH.chainnum = 10
        >>> HOH.chainnum
        10
        """
        return self._chainnum

    @chainnum.setter
    def chainnum(self, chainNumber):
        try:
            self._chainnum = int(chainNumber)
        except:
            raise ValueError("Molecule chain number could not be converted to integer")


"An example of class usage"
if __name__ == '__main__':
    mol1 = Molecule(num=1, nam='HOH', atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='O')])
    mol2 = Molecule(num=1, nam='PTP', atms=[Atom(num=1, nam='O1'), Atom(num=2, nam='C2'), Atom(num=3, nam='O')])
    mol3 = Molecule(num=1, nam='HOH', atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='O')])

    print(mol1 == mol2)
    print(mol1 is mol2)
    print(mol1 == mol3)
    print(mol1 is mol3)
    print(mol1.selectbyAtomnum(1))

    for atom in mol1:
        print(type(atom), id(atom), atom.show)
