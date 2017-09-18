#!/usr/bin/python
"""
.. module:: atom
   :platform: Unix

"""
from copy import deepcopy

class Atom(object):
    """
    A class to store the coordinates (x,y,z) and other atom properties.

       - **Parm num** (:class:`int`):      Atom number.
       - **Parm nam** (:class:`int`):      Atom name.
       - **Parm x** (:class:`float`):        Atom x coordinate.
       - **Parm y** (:class:`float`):        Atom y coordinate.
       - **Parm z** (:class:`float`):        Atom name.
       - **Parm typ** (:class:`Atom`):     Atom type (element)

    .. note::
       Initializing an :class:`Atom` object with arguments is optional.

       >>> O = Atom(121, 'OH', 13.015, 20.145, 21.457, 'O')
    """
    "tolerance value for float comparison"
    tolerance = 1E-5

    def __init__(self, num=0, nam='ND', x=0.0, y=0.0, z=0.0, typ='ND'):
        """
        Args:
           num (int): atom number
           x (float): x coordinate
           y (float): y coordinate
           z (float): z coordinate
           typ (string): atom type
        """
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)
        self._num = int(num)
        self._nam = str(nam)
        self._typ = str(typ)
    def __eq__(self, other):
        if isinstance(other, Atom):
            return (self._num == other._num) and (self._nam == other._nam) and \
                   (abs(self._x - other._x) <= Atom.tolerance) and             \
                   (abs(self._y - other._y) <= Atom.tolerance) and             \
                   (abs(self._z - other._z) <= Atom.tolerance) and (self._typ == other._typ)
        else:
            return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    @property
    def cord(self):
        """
            Sets/Returns coordinates.

            >>> O = Atom()
            >>> O.cord = [3.015, 20.145, 21.457]
            >>> x, y, z = O.cord
            >>> print(x,y,z)
            3.015 20.145 21.457
            >>> O.cord = [13.000, 20.145, 21.457]
            >>> O.cord
            [13.0, 20.145, 21.457]
            >>> xyz = O.cord
            >>> xyz
            [13.0, 20.145, 21.457]
        """
        return [self._x, self._y, self._z]

    @cord.setter
    def cord(self, cord):
        try:
            self._x, self._y, self._z = float(cord[0]), float(cord[1]), float(cord[2])
        except ValueError:
            raise ValueError("Pass a list of [x, y, z] coordinates")

    @property
    def num(self):
        """
            Sets/Returns atom number.

            >>> O = Atom()
            >>> O.num
            0
            >>> O.num = '10'
            >>> O.num
            10
        """
        return self._num

    @num.setter
    def num(self, num):
        try:
            self._num = int(num)
        except ValueError:
            raise ValueError("Atom number could not be converted to integer")

    @property
    def nam(self):
        """
            Sets/Returns atom name.

            >>> O = Atom()
            >>> O.nam
            'ND'
            >>> O.nam = 'OG'
            >>> O.nam
            'OG'
        """
        return self._nam

    @nam.setter
    def nam(self, nam):
        try:
            self._nam = str(nam)
        except:
            raise ValueError("Atom name could not be converted to string")

    @property
    def typ(self):
        """
            Sets/Returns atom type.

            >>> O = Atom()
            >>> O.typ
            'ND'
            >>> O.typ = 'O'
            >>> O.typ
            'O'
            >>> print(O.typ)
            O
        """
        return self._typ

    @typ.setter
    def typ(self, typ):
        try:
            self._typ = str(typ)
        except:
            raise ValueError("Atom type could not be converted to string")

    @property
    def show(self):
        """
            Returns num, nam, x, y, z, typ.

            >>> O = Atom(num=121, x=13.015, y=20.145, z=21.457, typ='O')
            >>> O.show
            [121, 'ND', 13.015, 20.145, 21.457, 'O']
        """

        return [self._num, self._nam, self._x, self._y, self._z, self._typ]

    @property
    def copy(self):
        """
            Returns a deepcopy of Atom object.

            >>> O = Atom(num=121, x=13.015, y=20.145, z=21.457, typ='O')
            >>> H = O.copy
            >>> O.show
            [121, 'ND', 13.015, 20.145, 21.457, 'O']
            >>> H.show
            [121, 'ND', 13.015, 20.145, 21.457, 'O']
            >>> H = Atom(num=12, x=13.015, y=20.145, z=21.457, typ='H')
            >>> O.show
            [121, 'ND', 13.015, 20.145, 21.457, 'O']
            >>> H.show
            [12, 'ND', 13.015, 20.145, 21.457, 'H']

            .. Note::
                A variable assignment only creates another link to the same object. If you wish to create a new
                object, use *copy* function.

            >>> H = Atom(num=1,typ='H')
            >>> O = H
            >>> O.typ = 'O'
            >>> O.num = 2
            >>> H.show
            [2, 'ND', 0.0, 0.0, 0.0, 'O']
        """
        return deepcopy(self)


" An example of Atom class usage"
if __name__ == '__main__':
    O1 = Atom(num=121, nam='O1', x='13.015', y=20.145, z=21.457, typ='O')
    O2 = Atom(num=121, nam='O1', x=13.015, y=20.145, z=21.457, typ='O')
    H2 = Atom(num=122, nam='H', x=13.015, y=20.145, z=21.457, typ='H')
    ATOMS = [Atom(num=121, nam='O1', x=13.015, y=20.145, z=21.457, typ='O'),
             Atom(num=122, nam='H', x=13.015, y=20.145, z=21.457, typ='H')]
    ATOMS2 = [Atom(num=121, nam='O1', x=13.015, y=20.145, z=21.457, typ='O'),
              Atom(num=122, nam='H', x=13.015, y=20.145, z=21.457, typ='H')]
    print(O1 == O2)
    print(O1 is O2)
    print(O1 in ATOMS)
    print(O1 == H2)
    print('--')
    print(ATOMS == ATOMS2)
    print(ATOMS is ATOMS2)
    print(type(O1.cord[0]))
