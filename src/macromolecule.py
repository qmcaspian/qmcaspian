#!/usr/bin/python
"""
.. module:: atom
   :platform: Unix

"""

from atom import Atom
from molecule import Molecule
from copy import deepcopy


class Macromolecule(object):
    """
    A class to store a collection of :class:`Molecule` objects.
       - **Parm nam** (:class:`str`):      Molecule name.
       - **Parm mols** (:class:`Atom`):    List of atoms in the molecule.

    .. todo::
        Under construction
    .. note::
       Initializing a :class:`Macromolecule` object with arguments is optional. It is possible to initialize a :class:`Macromolecule` object
       by one or a list of :class:`molecule` object. Number of molecules (*nmol*) is set automatically.

    >>> reslist = [Molecule(num=1,nam='ARG'), Molecule(num=1,nam='ALA'), Molecule(num=1,nam='GLY'),]
    >>> peptide = Macromolecule(nam='peptide A1', mols=reslist)
    """

    def __init__(self, nam=None, mols=None):
        self._nmol = 0
        self._nam = str(nam)
        self._mols = []
        if mols != None:
            self.addmol(mols)

    def __iter__(self):
        for mol in self._mols:
            yield mol

    def __contains__(self, item):
        if type(item) is Molecule and item in self._mols:
            return True
        elif type(item) is Atom and item in self.atoms:
            return True

    def addmol(self, mol):
        """
        Appends a hard copy of :class:`Molecule` object(s) to a :class:`Macromolecule` object.

        >>> res1, res2, res3 = Molecule(num=1,nam='ARG'), Molecule(num=1,nam='ALA'), Molecule(num=1,nam='GLY')
        >>> peptide = Macromolecule(nam='peptide A1')
        >>> peptide.show
        ['peptide A1', 0]
        >>> peptide.addmol(res1)
        >>> peptide.addmol([res2, res3])
        >>> peptide.show
        ['peptide A1', 3]
        """

        if type(mol) is Molecule:
            self._mols.append(mol.copy)
            self._nmol = self._nmol + 1
        elif type(mol) is list:
            for x in mol:
                if type(x) is Molecule:
                    self._mols.append(x.copy)
                    self._nmol = self._nmol + 1
                else:
                    raise ValueError("One element in the list is not Atom type")
        else:
            raise ValueError("Pass an Atom object")

    @property
    def show(self):
        """
        Returns nam, nmol.

        >>> peptide = Macromolecule(nam='A2', mols=[Molecule(), Molecule(), Molecule()])
        >>> peptide.show
        ['A2', 3]
        """
        return [self._nam, self._nmol]

    @property
    def atoms(self):
        """
        Returns the list of atoms in a :class:`Macromolecule`

        >>> res1 = Molecule(num=1, atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='H3')])
        >>> res2 = Molecule(num=2, atms=[Atom(num=4, nam='O1'), Atom(num=5, nam='O2'), Atom(num=6, nam='O3')])
        >>> res3 = Molecule(num=3, atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
        >>> peptide = Macromolecule(nam='peptide A1', mols=[res1, res2, res3])
        >>>
        >>> [[atom.num, atom.nam] for atom in peptide.atoms]
        [[1, 'H1'], [2, 'H2'], [3, 'H3'], [4, 'O1'], [5, 'O2'], [6, 'O3'], [7, 'N1'], [8, 'N2'], [9, 'N3']]
        """

        'This is a nested for loop to generated a flat list of list'
        return [atm for mol in self for atm in mol]

    @property
    def nmol(self):
        """
        Returns number of :class:`Molecule` in a :class:`Macromolecule`

        >>> peptide = Macromolecule(nam='A1', mols=[Molecule(num=1), Molecule(num=2), Molecule(num=3)])
        >>> peptide.nmol
        3
        """
        return self._nmol

    @property
    def copy(self):
        """
        Returns a deepcopy of :class:`Macromolecule` object.
        """
        return deepcopy(self)


    def selectbyMolnum(self, number):
        """
        Returns a :class:`Molecule` or a list of :class:`Molecule` with the given :class:`Molecule.num`. If only one
        :class:`Molecule` was found, an :class:`Molecule` object is returned. Otherwise a list of :class:`Molecule`
        objects is returned.

        >>> res1 = Molecule(num=1, nam='res1', atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='H3')])
        >>> res2 = Molecule(num=2, nam='res1', atms=[Atom(num=4, nam='O1'), Atom(num=5, nam='O2'), Atom(num=6, nam='O3')])
        >>> res3 = Molecule(num=3, nam='res1', atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
        >>> peptide = Macromolecule(nam='peptide A1', mols=[res1, res2, res3])
        >>> peptide.selectbyMolnum(1).show
        [1, 'res1', 3]
        """
        select = [mol for mol in self if mol.num == number]
        if len(select) == 1:
            return select[0]
        else:
            return select


    def selectbyMolnam(self, name):
        """
        Returns a :class:`Molecule` or a list of :class:`Molecule` with the given :class:`Molecule.nam`. If only one
        :class:`Molecule` was found, an :class:`Molecule` object is returned. Otherwise a list of :class:`Molecule`
        objects is returned.

        >>> res1 = Molecule(num=1, nam='res1', atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='H3')])
        >>> res2 = Molecule(num=2, nam='res2', atms=[Atom(num=4, nam='O1'), Atom(num=5, nam='O2'), Atom(num=6, nam='O3')])
        >>> res3 = Molecule(num=3, nam='res3', atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
        >>> peptide = Macromolecule(nam='peptide A1', mols=[res1, res2, res3])
        >>> peptide.selectbyMolnam('res2').show
        [2, 'res2', 3]
        """
        select = [mol for mol in self if mol.nam == name]
        if len(select) == 1:
            return select[0]
        else:
            return select

    def selectbyAtomnum(self, number):
        """
        Returns a :class:`Atom` or a list of :class:`Atom` with the given :class:`Atom.num`. If only one
        :class:`Atom` was found, an :class:`Atom` object is returned. Otherwise a list of :class:`Atom`
        objects is returned.

        >>> res1 = Molecule(num=1, nam='res1', atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='H3')])
        >>> res2 = Molecule(num=2, nam='res2', atms=[Atom(num=4, nam='O1'), Atom(num=5, nam='O2'), Atom(num=6, nam='O3')])
        >>> res3 = Molecule(num=3, nam='res3', atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
        >>> peptide = Macromolecule(nam='peptide A1', mols=[res1, res2, res3])
        >>> peptide.selectbyAtomnum(9).show
        [9, 'N3', 0.0, 0.0, 0.0, 'ND']
        """
        select = [atm for mol in self for atm in mol if atm.num == number]
        if len(select) == 1:
            return select[0]
        else:
            return select

    def selectbyAtomnam(self, name):
        """
        Returns a :class:`Atom` or a list of :class:`Atom` with the given :class:`Atom.nam`. If only one
        :class:`Atom` was found, an :class:`Atom` object is returned. Otherwise a list of :class:`Atom`
        objects is returned.

        >>> res1 = Molecule(num=1, nam='res1', atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='H3')])
        >>> res2 = Molecule(num=2, nam='res2', atms=[Atom(num=4, nam='O1'), Atom(num=5, nam='O2'), Atom(num=6, nam='O3')])
        >>> res3 = Molecule(num=3, nam='res3', atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
        >>> peptide = Macromolecule(nam='peptide A1', mols=[res1, res2, res3])
        >>> peptide.selectbyAtomnam('O1').show
        [4, 'O1', 0.0, 0.0, 0.0, 'ND']
        """
        select = [atm for mol in self for atm in mol if atm.nam == name]
        if len(select) == 1:
            return select[0]
        else:
            return select


if __name__ == '__main__':
    res1 = Molecule(num=1, atms=[Atom(num=1, nam='H1'), Atom(num=2, nam='H2'), Atom(num=3, nam='H3')])
    res2 = Molecule(num=2, atms=[Atom(num=4, nam='O1'), Atom(num=5, nam='O2'), Atom(num=6, nam='O3')])
    res3 = Molecule(num=3, atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
    res4 = Molecule(num=4, atms=[Atom(num=7, nam='N1'), Atom(num=8, nam='N2'), Atom(num=9, nam='N3')])
    peptide = Macromolecule(nam='peptide A1', mols=[res1, res2, res3])

    print(res4 in peptide)
    print(res2 in peptide)
    print(res2 is peptide.selectbyMolnum(2))
    print(peptide.selectbyMolnum(2))
    print(res2)
