#!/usr/bin/python

from atom import Atom
from molecule import Molecule
from macromolecule import Macromolecule
from atom_table import atom_property
import sys
import re


def parseQM(ifile, program):
    """
    This is a wrapper function for classes that parse the QM calculations output files


    """
    if program == 'G':
        results = ParseGaussian(ifile)
    else:
        results = None

    return results()


class ParseGaussian(object):
    """
    This class reads the header section of a Gaussian output and initiate a corresponding parser class.


    """

    def __init__(self, ifile):

        self._calc = ''
        self._method = ''
        self._version = ''
        self._results = None
        self.getheader(ifile)

        if re.search(r'\sopt', self._method, re.IGNORECASE):
            self._results = OptG09(ifile)

        elif re.search(r'\sfrq', self._method, re.IGNORECASE):
            print("not implemented yet")

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
            raise ValueError("No header section was found in " + str(ifile))


class OptG09(Macromolecule):
    """
    A class to read and save G09 optimization output. This class inherits from :class:`Macromolecule`.
    This class assumes the given file exists and is G09 optimization output. Use the :class:`ParseGaussian` if error
    management is needed.


    """

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
            # The variables to cash two previous lines
            line_2 = ""
            line_1 = ""
            # Keep track of structure numbers
            imolnum = 0
            # Read flag
            geometry_found = False

            for line in f:
                if re.match(r'^\sNumber\s+Number\s+Type\s+X\s+Y\s+Z\n$', str(line_2)) and \
                        re.match(r'\s-+\n$', str(line_1)):
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

                # Convergence flag
                if re.search(r'^\sOptimization\s+completed', str(line)):
                    self.converged = True

                # Read Charge and multiplicity
                if re.search(r'Charge =', str(line)) and re.search(r'Multiplicity =', str(line)):
                    self.charge = line.split()[2]
                    self.multiplicity = line.split()[5]

                # Update the cash lines -2 and -1
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


" An example of class usage"
if __name__ == '__main__':
    a = parseQM('int1_EOAwWAT.log', 'G')
    print(a.show)
