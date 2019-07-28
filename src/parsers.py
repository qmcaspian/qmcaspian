#!/usr/bin/env python3

from atom import Atom
from molecule import Molecule
from macromolecule import Macromolecule
from atom_table import atom_property
from graph import *
from resp import ESP
from internal_coordinate import *
import sys
import os
import re
import math
import numpy as np


# TODO create a separate opt object to be returned by the __call__ function of OPTG09

def parseQM(ifile, program):
    """
    This is a wrapper function for classes that parse the QM calculations output files


    """
    if program == 'G':
        results = ParseGaussian(ifile)
    else:
        results = None

    # run the __call__ function
    return results()


class ParseGaussian(object):
    """
    This class reads the header section of a Gaussian output and initiate a corresponding parser class.

    """

    def __init__(self, ifile):
        self._f = None
        self._calc = ''
        self._method = ''
        self._version = ''
        self._results = None
        self.getheader(ifile)

        if ifile.endswith('fchk'):
            if re.search(r'freq', self._method, re.IGNORECASE):
                print("frequency calculations - input file is a formated check point file")
                results = FreqFchkG09(ifile, self._method)
                self._results = results()

        else:
            # TODO create a separate opt object to be returned by the __call__ function of OPTG09
            if re.search(r'\sopt', self._method, re.IGNORECASE):
                print("Gaussian optimization file")
                self._results = OptG09(ifile)


            elif re.search(r'\sfreq', self._method, re.IGNORECASE):
                print("Error >>> Not implemented yet")
                raise AssertionError

    def __call__(self):
        return self._results

    def getheader(self, ifile):
        """
        Checks input file header for Gaussian version and command line options
        """
        try:
            self._f = open(ifile, 'r')
        except Exception as e:
            sys.exit("Input/Output error: %s" % (str(e)))

        version_found = False
        method_found = False

        # Read a formatted checkpoint file
        if ifile.endswith('.fchk'):
            for line in self._f:
                if re.match(r'Route', str(line)):
                    # the Route line  contain the number of characters (c variable) forming the method section.
                    # 1 C = 12 characters and each line in this section contains 60 characters.
                    # number of line to read after Route section is ceil(C * 12 / 60)
                    method_found = True
                    c = int(line.strip('\n').split()[3])
                    lines_to_read = math.ceil(c * 12 / 60)
                    break
            for i in range(lines_to_read):
                self._method = self._method + str(self._f.readline().strip('\n'))

        # try to figure out what is the input file
        else:
            # For cashing the previous line
            line_1 = ""
            for line in self._f:
                if re.match(r'^\s\*+\n$', str(line_1)) and 'Gaussian' in line:
                    version_found = True
                    self._version = str(line.split()[0]) + ' ' + str(line.split()[1])
                if version_found and re.match(r'^\s-+\n$', str(line_1)):
                    method_found = True
                if (method_found and re.match(r'^\s-+\n$', str(line))) or 'Leave Link    1' in line:
                    break
                if method_found:
                    self._method = self._method + ' ' + line.strip('\n')
                line_1 = line
            self._f.close()

        if not method_found:
            raise ValueError("Cant identify the file. No header section was found in " + str(ifile))


class OptG09(Macromolecule):
    """
    A class to read and save G09 optimization output. This class inherits from :class:`Macromolecule`. Here we are using
    the :class:`Macromolecule` to as a capsule to keep each structure in the opt. Thus the molecules in :class:`Macromolecule`
    corresponds to each optimized structure.
    This class assumes the given file exists and is G09 optimization output. Use the :class:`ParseGaussian` if error
    management is needed.
    """

    # TODO create a separate opt object to be returned by the __call__ function of OPTG09

    def __init__(self, ifile):
        Macromolecule.__init__(self, nam=ifile.split('.')[0])
        self.foundStructure = False
        self.converged = False
        self.charge = 0
        self.multiplicity = 0
        self.energy = [0.0000]
        self.maxForce = [0.0000]
        self.averageForce = [0.0000]
        self.maxDisplacement = [0.0000]
        self.averageDisplacement = [0.0000]
        self.parse_opt_g09(ifile)

    def parse_opt_g09(self, ifile):
        with open(ifile, mode='r') as f:
            # The variables to cash previous lines
            line_5 = ""
            line_4 = ""
            line_3 = ""
            line_2 = ""
            line_1 = ""
            # Keep track of structure numbers
            imolnum = 0
            # Read flag
            geometry_found = False

            for line in f:
                if re.match(r'^\sNumber\s+Number\s+Type\s+X\s+Y\s+Z\n$', str(line_2)) and \
                        re.match(r'\s-+\n$', str(line_1)) and re.match(r'\s*Standard\sorientation:\s*\n$', str(line_5)):
                    geometry_found = True
                    self.foundStructure = True
                    imolnum += 1
                    self.addmol(Molecule(num=imolnum))
                    imol = self.selectbyMolnum(imolnum)

                elif re.match(r'\s-+\n$', str(line)) and geometry_found is True:
                    geometry_found = False

                if geometry_found:
                    inum, inam, dummy, ix, iy, iz = line.split()
                    imol.addatm(
                        Atom(num=inum, nam=atom_property(query='symbol', target='number', value=int(inam)), x=ix, y=iy,
                             z=iz, ))

                # Read energies
                try:
                    if re.search(r'^\sSCF\s+Done', str(line)):
                        self.energy.append(float(line.split()[4]))

                    # Read Max Force
                    if re.search(r'^\sMaximum Force\s+', str(line)):
                        self.maxForce.append(float(line.split()[2]))

                    # Read average Force
                    if re.search(r'^\sRMS     Force\s+', str(line)):
                        self.averageForce.append(float(line.split()[2]))

                    # Read Max Displacement
                    if re.search(r'^\sMaximum Displacement\s+', str(line)):
                        self.maxDisplacement.append(float(line.split()[2]))

                    # Read average Displacement
                    if re.search(r'^\sRMS     Displacement\s+', str(line)):
                        self.averageDisplacement.append(float(line.split()[2]))
                except:
                    print("Warning >>> Structure" + str(imolnum) + "seems to have problem!")

                # Convergence flag
                if re.search(r'^\sOptimization\s+completed', str(line)):
                    self.converged = True

                # Read Charge and multiplicity
                if re.search(r'Charge =', str(line)) and re.search(r'Multiplicity =', str(line)):
                    self.charge = line.split()[2]
                    self.multiplicity = line.split()[5]

                # Update the cash lines -2 and -1
                line_5 = line_4
                line_4 = line_3
                line_3 = line_2
                line_2 = line_1
                line_1 = line

            # Check any structure was found in the file
            if not self.foundStructure:
                raise ValueError("No Structure was found in the " + ifile)
            elif self.nmol == 1:
                print("Warning, only one structure was found. Optimization doesnt seem to proceeded.")

            # Push the lists by one element to match the numbering of the structures
            # First structure is the initial structure
            self.energy.insert(0, self.energy[0])
            self.maxForce.insert(0, self.maxForce[0])
            self.averageForce.insert(0, self.averageForce[0])
            self.maxDisplacement.insert(0, self.maxDisplacement[0])
            self.averageDisplacement.insert(0, self.averageDisplacement[0])


class Freq(object):
    def __init__(self):
        self._method = None
        self.internal_coord = InternalCoordinate()
        self.hessian_cartesian = None


class FreqFchkG09():
    Hartree2Kcalmol = 627.509608
    Bohr2Angstrom = 0.529177249
    Rad2Degree = 57.2957795

    def __init__(self, ifile, method=None):
        # Sequence of these commands is important. The internal coordinates are created based on the "structure" object.
        self._natm = self.get_natm(ifile)  # To temporarily cash the number of atoms
        self.freq = Freq()
        self.freq.method = method
        self.get_chrge(ifile)
        self.get_atoms(ifile)
        self.get_bonded_coord(ifile)
        self.get_hessian_cartesian(ifile)
        """
        print(self.freq.internal_coord.natm)
        print(self.freq.internal_coord.charge)
        for atom in self.freq.internal_coord.structure:
            print(atom.show)

        for bond in self.freq.internal_coord.bonds:
            print(bond.show)

        for angle in self.freq.internal_coord.angles:
            print(angle.show)

        for torsion in self.freq.internal_coord.torsions:
            print(torsion.show)

        np.set_printoptions(precision=1, suppress=True, threshold=np.inf, linewidth=520)
        print(self.hessian_cartesian)
        """

    def __call__(self):
        return self.freq

    def get_natm(self, ifile):
        with open(ifile, mode='r') as f:
            for line in f:
                if re.search(r'^Atomic numbers', str(line)):
                    return int(line.split()[4])
                elif not line:
                    print('Warning >>> Could not find the number of atoms')
                    return None

    def get_chrge(self, ifile):
        # Read the charge
        with open(ifile, mode='r') as f:
            for line in f:
                if re.search(r'^Charge', str(line)):
                    self.freq.internal_coord.charge = int(line.split()[2])
                elif not line:
                    print('Warning >>> Could not find the total charge')

    def get_atoms(self, ifile):
        atoms_found = False
        # Find the atom data section
        with open(ifile, mode='r') as f:
            # read  the atomic numbers
            for line in f:
                if re.search(r'^Atomic numbers', str(line)):
                    atoms_found = True
                    break  # bring the pointer to the correct place

            # each line contain Max 6 atomic numbers, So the number of lines to be read is ceil(natm / 6)
            if atoms_found:
                self.freq.internal_coord.structure.num = 1
                self.freq.internal_coord.structure.nam = os.path.splitext(ifile)[0]
                count = 1
                for i in range(math.ceil(self._natm / 6)):
                    for atm_num in f.readline().split():
                        # Get the atom symbol by atomic number
                        atm_symbol = atom_property(query='symbol', target='number', value=int(atm_num))
                        self.freq.internal_coord.structure.addatm(Atom(num=count, nam=atm_symbol, typ=atm_symbol))
                        count += 1
            else:
                print('Warning >>> Atom numbers were not found')

            # Find the coordinates section
            atoms_found = False
            for line in f:
                if re.search(r'^Current cartesian coordinates', str(line)):
                    # Number of lines to read is the number of elements in the section / 5 (number of elements in each line)
                    lines_to_read = math.ceil(int(line.split()[5]) / 5)
                    atoms_found = True
                    break

            # Read the coordinate
            if atoms_found:
                coord_list = []
                for i in range(lines_to_read):
                    for coord in f.readline().split():
                        coord_list.append(self.Bohr2Angstrom * float(coord))

                # Deposit the coordinates to the molecule
                for atom in self.freq.internal_coord.structure:
                    lower_bound = ((atom.num - 1) * 3)
                    higher_bound = lower_bound + 3
                    atom.cord = coord_list[lower_bound:higher_bound]
            else:
                print('Warning >>> Atom coordinates were not found in ', ifile)

    def get_bonded_coord(self, ifile):

        with open(ifile) as f:

            # Get the numbers of the internal coord
            internal_coord_found = False
            for line in f:
                if re.search(r'Redundant internal dimensions', str(line)):
                    internal_coord_found = True
                    break

            if internal_coord_found:
                n_iternal_coords, n_bonds, n_angles, n_torsions = f.readline().split()
                n_iternal_coords, n_bonds, n_angles, n_torsions = int(n_iternal_coords), int(n_bonds), int(
                    n_angles), int(n_torsions)
            else:  # Get out
                print('Warning >>> The bonded terms were not found in ', ifile)
                return None

            # Find the internal coordinates indices
            internal_coord_found = False
            for line in f:
                if re.search(r'Redundant internal coordinate indices', str(line)):
                    internal_coord_found = True
                    n_indices = int(line.split()[6])
                    lines_to_read = math.ceil(n_indices / 6)  # There are 6 indices per line
                    break

            # Read internal coordinates indices into a list
            if internal_coord_found:
                internal_coord_indices = []
                for i in range(lines_to_read):
                    for indice in f.readline().split():
                        internal_coord_indices.append(int(indice))
            else:  # Get out
                print('Warning >>> The internal coordinates were not found')
                return None

            # Find the equilibrium values of the internal coordinates
            internal_coord_found = False
            for line in f:
                if re.search(r'Redundant internal coordinates', line):
                    internal_coord_found = True
                    lines_to_read = math.ceil(n_iternal_coords / 5)  # There are 5 values per line
                    break

            # Read the values into a list
            internal_coord_val = []
            if internal_coord_found:
                for i in range(lines_to_read):
                    for val in f.readline().split():
                        internal_coord_val.append(float(val))
            else:  # get out
                print('Warning >>> The internal coordinates were not found')
                return None

        # Generate the internal_coord object. Every 4 elements corresponds to one degree of freedom. That is,
        # the bonds have zeros in the last two elements, etc.
        if internal_coord_found and self.freq.internal_coord.natm > 1:
            count = 0
            while count < n_bonds:
                lower_bound = count * 4
                higher_bound = lower_bound + 4
                i_indice = internal_coord_indices[lower_bound:higher_bound]
                atom_i = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[0])
                atom_j = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[1])
                i_bond = Bond(i=atom_i, j=atom_j, r=(self.Bohr2Angstrom * internal_coord_val[count]))
                self.freq.internal_coord.addbond(i_bond)
                count += 1

            while count < n_angles + n_bonds:
                lower_bound = count * 4
                higher_bound = lower_bound + 4
                i_indice = internal_coord_indices[lower_bound:higher_bound]
                atom_i = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[0])
                atom_j = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[1])
                atom_k = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[2])
                i_angle = Angle(i=atom_i, j=atom_j, k=atom_k, t=(self.Rad2Degree * internal_coord_val[count]))
                self.freq.internal_coord.addangle(i_angle)
                count += 1

            while count < n_torsions + n_angles + n_bonds:
                lower_bound = count * 4
                higher_bound = lower_bound + 4
                i_indice = internal_coord_indices[lower_bound:higher_bound]
                atom_i = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[0])
                atom_j = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[1])
                atom_k = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[2])
                atom_l = self.freq.internal_coord.structure.selectbyAtomnum(i_indice[3])
                i_torsion = Torsion(i=atom_i, j=atom_j, k=atom_k, l=atom_l,
                                    t1=(self.Rad2Degree * internal_coord_val[count]))
                self.freq.internal_coord.addtorsion(i_torsion)
                count += 1

    def get_hessian_cartesian(self, ifile):

        # Read the hessian as a list.
        hessian_found = False
        with open(ifile) as f:
            for line in f:
                if re.search(r'Cartesian Force Constants', str(line)):
                    hessian_found = True
                    n_elements = int(line.split()[5])
                    # there are 5 elements per line
                    lines_to_read = math.ceil(n_elements / 5)
                    break

            # Read the lower triangle matrix into a list
            if hessian_found:
                hessian_list = []
                for i in range(lines_to_read):
                    for element in f.readline().split():
                        hessian_list.append(float(element))
            else:  # Get out
                print('Warning >>> The Cartesian Hessian matrix was not found')
                return None

        dim = self._natm * 3
        hessian = np.zeros(shape=(dim, dim), dtype=np.float64)
        i, j = np.tril_indices(dim)
        np.set_printoptions(precision=5, suppress=True, threshold=np.inf, linewidth=520)
        hessian[(i, j)] = hessian_list
        hessian[(j, i)] = hessian_list
        hessian *= self.Hartree2Kcalmol / (self.Bohr2Angstrom * self.Bohr2Angstrom)
        self.freq.hessian_cartesian = hessian


class parseMol2(object):
    def __init__(self, ifile):
        self.internal_coord = InternalCoordinate()
        self.get_atoms(ifile)
        self.get_bonds(ifile)
        self.getresult()

    def get_atoms(self, ifile):
        with open(ifile) as f:
            atom_found = False
            for line in f:
                if re.search(r'@<TRIPOS>ATOM', str(line)):
                    atom_found = True
                    continue
                if re.search(r'@<TRIPOS>BOND', str(line)):
                    break

                if atom_found:
                    inum, iname, ix, iy, iz, ityp = line.split()[0:6]
                    self.internal_coord.addatm(Atom(num=inum, nam=iname, x=ix, y=iy, z=iz, typ=ityp[0]))
        if not atom_found:
            print('Warning >>> Atom section was not found')
            return None

    def get_bonds(self, ifile):
        with open(ifile) as f:
            bond_found = False
            for line in f:
                if re.search(r'@<TRIPOS>BOND', str(line)):
                    bond_found = True
                    continue
                if re.search(r'@<TRIPOS>SUBSTRUCTURE', str(line)):
                    break

                if bond_found:
                    i = int(line.split()[1])
                    j = int(line.split()[2])
                    atom_i = self.internal_coord.selectbyAtomnum(i)
                    atom_j = self.internal_coord.selectbyAtomnum(j)
                    dij = np.linalg.norm(np.array(atom_j.cord) - np.array(atom_i.cord))
                    self.internal_coord.addbond(Bond(i=atom_i, j=atom_j, r=dij))
        if not bond_found:
            print('Warning >>> Bond section was not found')
            return None

    def getresult(self):
        return self.internal_coord


def read_esp(file, esp_ref=None):
    Hartree2Kcalmol = 627.509608
    Bohr2Angstrom = 0.529177249

    if esp_ref is None:
        esp = ESP()
    else:
        esp = esp_ref

    atom_sectiom_found = False
    with open(file) as f:
        for line in f:
            if re.search(r'CHARGE =', str(line)):
                esp.charge = int(line.split()[2])
            if re.search(r'ATOMIC COORDINATES AND ESP CHARGES. #ATOMS =', str(line)):
                atom_sectiom_found = True
                lines_to_read = int(line.split()[7])
                break

        if atom_sectiom_found:
            for i in range(lines_to_read):
                line = f.readline()
                atom, x, y, z, charge = str(line.split()[0]), \
                                        float(line.split()[1].replace('D', 'e')) * Bohr2Angstrom, \
                                        float(line.split()[2].replace('D', 'e')) * Bohr2Angstrom, \
                                        float(line.split()[3].replace('D', 'e')) * Bohr2Angstrom, \
                                        float(line.split()[4].replace('D', 'e'))
                esp.addatm(Atom(num=(i + 1), nam=atom, x=x, y=y, z=z, typ=atom, charge=charge))

    esp_matrix_found = False
    with open(file) as f:
        for line in f:
            if re.search(r'ESP VALUES AND GRID POINT COORDINATES', str(line)):
                esp_matrix_found = True
                numberofpoints = int(line.split()[8])
                break

        if esp_matrix_found:
            esp.points = np.zeros((numberofpoints, 4))
            for point in range(numberofpoints):
                e, x, y, z = f.readline().split()
                esp.points[point, :] = float(e.replace('D', 'e')) * Hartree2Kcalmol, \
                                       float(x.replace('D', 'e')) * Bohr2Angstrom, \
                                       float(y.replace('D', 'e')) * Bohr2Angstrom, \
                                       float(z.replace('D', 'e')) * Bohr2Angstrom

    if atom_sectiom_found and esp_matrix_found and esp_ref is None:
        return esp
    else:
        return None


class read(object):
    Hartree2Kcalmol = 627.509608
    Bohr2Angstrom = 0.529177249
    Rad2Degree = 57.2957795

    @staticmethod
    def mol2(ifile):
        nAtoms = 0
        nBonds = 0
        nSubStructures = 0
        internal_coord = InternalCoordinate() # Return object

        with open(ifile) as f:
            # Look for the molecule data
            molecule_info_found = False
            for line in f:
                if re.search(r'@<TRIPOS>MOLECULE', str(line)):
                    molecule_info_found = True
                    break

            if molecule_info_found:
                # Read molecule name
                internal_coord.nam = f.readline().strip()
                # Read number of atom, bond, and substructures
                line = f.readline().split()
                if len(line) == 1:
                    nAtoms = int(line[0])
                if len(line) == 2:
                    nAtoms, nBonds = int(line[0]), int(line[1])
                if len(line) >= 3:
                    nAtoms, nBonds, nSubStructures = int(line[0]), int(line[1]), int(line[2])

            # Make sure the file contains only one residue/molecule.
            if nSubStructures > 1:
                print("Error >>> The ", ifile, " file contain multiple residues/molecules")
                print("This is not supported, the mole2 are strictly used for structures with single substructure.")
                return None

            # Read the atom section
            f.seek(0)
            atom_found = False
            for line in f:
                if re.search(r'@<TRIPOS>ATOM', str(line)):
                    atom_found = True
                    break

            if atom_found:
                for atom in range(nAtoms):
                    line = f.readline().split()
                    # The first five records always exist. The atom type recored is split handel the SYBEL atom tye.
                    inum, iname, ix, iy, iz, ityp = int(line[0]), line[1], float(line[2]), float(line[3]), \
                                                                                float(line[4]), line[5].split(".")[0]
                    # The position of charge recored depends on the presence of residue record.
                    icharge = 0
                    if nSubStructures == 0:
                        # If no residue info is present. The charge (if exist) is at position 6 as float.
                        try:
                            icharge = float(line[6])
                        except:
                            pass # Then no charge record
                    else:
                        # The charge (if exist) is at position 8 as float.
                        try:
                            icharge = float(line[8])
                        except: # Then no charge record
                            pass

                    internal_coord.addatm(Atom(num=inum, nam=iname, x=ix, y=iy, z=iz, typ=ityp, charge=icharge))

            # Read the bond section
            bond_found = False
            f.seek(0)
            for line in f:
                if re.search(r'@<TRIPOS>BOND', str(line)):
                    bond_found = True
                    break

            if bond_found:
                for bond in range(nBonds):
                    line = f.readline()
                    i, j, bond_order = int(line.split()[1]), int(line.split()[2]), str(line.split()[3])
                    atom_i = internal_coord.selectbyAtomnum(i)
                    atom_j = internal_coord.selectbyAtomnum(j)
                    try:
                        bond_order = float(bond_order)
                    except ValueError:
                        # Keeping support for Mol2 Ar bond type
                        if bond_order == 'Ar':
                            bond_order = 1.5
                        else:
                            print('Warning >>> Bond order ', bond_order, ' between ', atom_i.num, '-', atom_j.num, 'is not supported, change it to numerical')
                    dij = np.linalg.norm(np.array(atom_j.cord) - np.array(atom_i.cord))
                    internal_coord.addbond(Bond(i=atom_i, j=atom_j, r=dij, order=bond_order))


        if atom_found and bond_found:
            return internal_coord
        elif atom_found and not bond_found:
            print('Warning >>> Bond section was not found')
            return internal_coord
        else:
            print('Warning >>> Atom section was not found')
            return None


    @staticmethod
    def esp(file, iESP=None):

        if iESP is None:
            esp = ESP()
        elif type(iESP) is ESP:
            esp = iESP
        else:
            raise ValueError('Error >>> iESP should be an instance of ESP class')

        atom_section_found = False
        lines_to_read = 0

        with open(file) as f:
            for line in f:
                if re.search(r'CHARGE =', str(line)):
                    esp.charge = int(line.split()[2])
                if re.search(r'ATOMIC COORDINATES AND ESP CHARGES. #ATOMS =', str(line)):
                    atom_section_found = True
                    lines_to_read = int(line.split()[7])
                    break

            if atom_section_found:
                for i in range(lines_to_read):
                    line = f.readline()
                    atom, x, y, z, charge = str(line.split()[0]), \
                                            float(line.split()[1].replace('D', 'e')) * read.Bohr2Angstrom, \
                                            float(line.split()[2].replace('D', 'e')) * read.Bohr2Angstrom, \
                                            float(line.split()[3].replace('D', 'e')) * read.Bohr2Angstrom, \
                                            float(line.split()[4].replace('D', 'e'))
                    esp.addatm(Atom(num=(i + 1), nam=atom, x=x, y=y, z=z, typ=atom, charge=charge))

        esp_matrix_found = False
        number_of_points = 0

        with open(file) as f:
            for line in f:
                if re.search(r'ESP VALUES AND GRID POINT COORDINATES', str(line)):
                    esp_matrix_found = True
                    number_of_points = int(line.split()[8])
                    break

            if esp_matrix_found:
                esp.points = np.zeros((number_of_points, 4))
                for point in range(number_of_points):
                    e, x, y, z = f.readline().split()
                    esp.points[point, :] = float(e.replace('D', 'e')) * read.Hartree2Kcalmol, \
                                           float(x.replace('D', 'e')) * read.Bohr2Angstrom, \
                                           float(y.replace('D', 'e')) * read.Bohr2Angstrom, \
                                           float(z.replace('D', 'e')) * read.Bohr2Angstrom

        if atom_section_found and esp_matrix_found and iESP is None:
            return esp
        elif atom_section_found and esp_matrix_found and iESP is not None:
            return None
        else:
            print("Error >>> The atom section or the ESP grid was not found.")
            return None

    @staticmethod
    def charge_restraints(file, iesp):

        try:
             # Check the inputs
            if type(iesp) is not ESP:
                raise ValueError('Error >>> read.charge_restraints() needs a ESP object as input')

        except ValueError as error:
            print(error.args)
            return None

        with open(file) as f:

            # find the charge restraints
            restraints_found = False
            lines_to_read = 0
            for line in f:
                if re.search(r'CONSTRAINTS', str(line), re.IGNORECASE):
                    restraints_found = True
                    for word in line.split():
                        if word.isdigit():
                            lines_to_read = int(word)
                            break
                if restraints_found:
                    break

            if restraints_found and lines_to_read > 0:
                restraints = []
                for i in range(lines_to_read):
                    line = f.readline()
                    try:
                        atom_num, charge = int(line.split()[0]), float(line.split()[1])
                        if atom_num not in [i[0] for i in restraints]:
                            restraints.append([atom_num, charge])
                    except ValueError as error:
                        print('Could not read the charge restraints at line: ' + str(i+1))
                iesp.restraints = restraints
            else:
                # Just give a notification and continue for the symmetry section
                print("Warning >>> No charge restraints were not found.")

            # find the charge symmetry.
            f.seek(0)
            symmetry_found = False
            lines_to_read = 0
            for line in f:
                if re.search(r'SYMMETRY', str(line), re.IGNORECASE):
                    symmetry_found = True
                    for word in line.split():
                        if word.isdigit():
                            lines_to_read = int(word)
                            break
                if symmetry_found:
                    break

            if symmetry_found and lines_to_read > 0:
                symmetry = []
                for i in range(lines_to_read):
                    line = f.readline()
                    row = []

                    for atom_num in line.split():
                        if not atom_num.isdigit():
                            break
                        num = int(atom_num)
                        if num in [atom for row in symmetry for atom in row] or num in row:
                            print("Error >>> Repeated atom in the symmetry section, line: " + str(i + 1))
                            break
                        else:
                            row.append(num)
                    if len(row) > 1:
                        symmetry.append(sorted(row))
                    else:
                        print("Error >>> Bad charge symmetries at line: " + str(i + 1))
                        break
                iesp.symmetry = symmetry
            else:
                # Just give a notification and continue for the symmetry section
                print('Warning >>> No symmetry restraints were not found.')

class write(object):

    @staticmethod
    def mol2(data, name):
        if type(data) is Graph or type(data) is Fragment:
            write._graph2mol2(data, name)
        else:
            print("Error >>> writing mol2 from ", type(data), " data type is not implemented." )

    @staticmethod
    def _graph2mol2(graph, name):

        # Number of atoms and bonds
        nAtom = graph.number_of_nodes
        bonds = graph.bonds
        nbond = len(bonds)
        nSubstructure = 0

        if nAtom == 0 :
            print("Error >>> number of atoms can not be 0.")
            return None


        with open(name + ".mol2", 'w') as f:
            # Write the MOLECULE section
            f.write('@<TRIPOS>MOLECULE \n')
            f.write("%s \n" % (graph.name))
            f.write( "%i   %i   %i \n" %(nAtom, nbond, nSubstructure))
            f.write('SMALL \n')
            f.write('USER_CHARGES \n')

            if graph.fragment:
                # Write out the head atom
                f.write('\n@<TRIPOS>HEAD \n')
                f.write('{0:3d}   {1:3s} \n'.format(graph.head.atom.num, graph.head.atom.nam))

                f.write('\n@<TRIPOS>TYPE \n')
                for ff, typ in graph.atom_type():
                    f.write('FF   {0:5s}   {1:5s} \n'.format(ff, typ))

            # Write The atom section with out residue/substructure definition.
            f.write('\n@<TRIPOS>ATOM \n')
            for node in graph.nodes:
                f.write('{0:3d} {1:5s} {2: .4f}  {3: .4f}  {4: .4f}  {5:4s} {6: .4f} \n'.format(node.atom.num,
                                node.atom.nam, node.atom.x, node.atom.y, node.atom.z, node.atom.typ, node.atom.charge))

            if nbond > 0:
                # Write the bond section
                f.write('\n@<TRIPOS>BOND \n')
                for i, bond in enumerate(bonds):
                    # Take care of aromatic bonds. The internal representation is 1.5f
                    order = bond.order
                    if order == 1.5:
                        order = 'Ar'
                    else:
                        order = str(int(order))

                    f.write('{0:3d} {1:3d} {2:3d}   {3:3s} \n'.format(i, bond.i.num, bond.j.num, order))






" An example of class usage"
if __name__ == '__main__':
    pass
