#!/usr/bin/python


def atom_property(query, target, value):
    """
    returns the atomic property (query) based on the given atomic property (traget) with a given value (value)

    - **parm query**: needed atomic property
    - **parm target**: given atomic property
    - **parm value**: the value of given atomic property

    >>> print(atom_property(query='symbol', target='number', value=25))
    Mn
    >>> print(atom_property(query='mass', target='number', value=25))
    54.93805
    """
    innput = -1
    output = -1
    properties = ['name', 'symbol', 'number', 'mass']
    atom_table = [['Hydrogen', 'H', 1, 1.00794],
                    ['Helium', 'He', 2, 4.002602],
                    ['Lithium', 'Li', 3, 6.941],
                    ['Beryllium', 'Be', 4, 9.01218],
                    ['Boron', 'B', 5, 10.811],
                    ['Carbon', 'C', 6, 12.011],
                    ['Nitrogen', 'N', 7, 14.00674],
                    ['Oxygen', 'O', 8, 15.9994],
                    ['Fluorine', 'F', 9, 18.998403],
                    ['Neon', 'Ne', 10, 20.1797],
                    ['Sodium', 'Na', 11, 22.989768],
                    ['Magnesium', 'Mg', 12, 24.305],
                    ['Aluminum', 'Al', 13, 26.981539],
                    ['Silicon', 'Si', 14, 28.0855],
                    ['Phosphorus', 'P', 15, 30.973762],
                    ['Sulfur', 'S', 16, 32.066],
                    ['Chlorine', 'Cl', 17, 35.4527],
                    ['Argon', 'Ar', 18, 39.948],
                    ['Potassium', 'K', 19, 39.0983],
                    ['Calcium', 'Ca', 20, 40.078],
                    ['Scandium', 'Sc', 21, 44.95591],
                    ['Titanium', 'Ti', 22, 47.88],
                    ['Vanadium', 'V', 23, 50.9415],
                    ['Chromium', 'Cr', 24, 51.9961],
                    ['Manganese', 'Mn', 25, 54.93805],
                    ['Iron', 'Fe', 26, 55.847],
                    ['Cobalt', 'Co', 27, 58.9332],
                    ['Nickel', 'Ni', 28, 58.6934],
                    ['Copper', 'Cu', 29, 63.546],
                    ['Zinc', 'Zn', 30, 65.39],
                    ['Gallium', 'Ga', 31, 69.723],
                    ['Germanium', 'Ge', 32, 72.61],
                    ['Arsenic', 'As', 33, 74.92159],
                    ['Selenium', 'Se', 34, 78.96],
                    ['Bromine', 'Br', 35, 79.904],
                    ['Krypton', 'Kr', 36, 83.8]]

    for i,parm in enumerate(properties):
        if parm == target: innput = i
        if parm == query: output = i

    if innput == -1 or output == -1:
        raise ValueError("The selected properties were not found. query=" + query + "  target=" + target + "\n" + \
                         "The accepted keywords are: 'name', 'symbol', 'number', 'mass'")
    for row in atom_table:
        if row[innput] == value:
            return row[output]

if __name__ == '__main__':
    print(atom_property(query='name', target='symbol', value='S'))

