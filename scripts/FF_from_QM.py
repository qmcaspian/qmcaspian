#!/usr/bin/python3.5

from parse_qm import *
from ff_bonded_from_qm import *
import argparse

def FF_bonded_from_QM():
    """
    Read a formatted G09 checkpoint file of a frequency calculation to extract the force constant
    of the bonded terms. These terms are according to the harmonic approximation and should be used
     with care. The QM calculation should be performed with "freq=intmodes" keyword in G09 for the
     hessian matrix to be printed in internal coordinates.

         * How to set up:

        export the library directories to PYTHONPATH. In my case it is

          print(mol4)  **>>> export PYTHONPATH=$PYTHONPATH:/home/masoud/qmcaspian/src/**

        export the scripts directory to PATH

        **>>> export PATH=$PATH:/home/masoud/qmcaspian/scripts/**

    * How to use:

        *>>> FF_bonded_from_QM.py <inputfile>.chk*
    """
    # Read the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input file: Formatted checkpoint file of "
                                                     "freq calculation G09 with \"freq=intmodes\" keyword")
    args = parser.parse_args()
    ifile = args.input_file

    # Read the input frequency calculation checkpoint file
    freq = parseQM(ifile, 'G')
    #print(type(freq))

    np.set_printoptions(precision=3, suppress=True, threshold=np.inf, linewidth=520)
    #print(freq.hessian_cartesian)

    # Calculate the force constances and assign them to the internal coordinates.
    GetBondedfromQM(freq)

    for atom in freq.internal_coord.structure:
        print(atom.show)

    for i, bond in enumerate(freq.internal_coord.bonds):
        print(i, bond.show)
    #for angle in freq.internal_coord.angles:
    #    print(angle.show)
    #for torsion in freq.internal_coord.torsions:
    #    print(torsion.show)
    print('------------------------------------------------------------')
    freq.internal_coord._bonds = []
    freq.internal_coord.constructbondsbyHessian(freq.hessian_cartesian)
    freq.internal_coord.constructbondsbyCovalentRadius()

    for i, bond in enumerate(freq.internal_coord.bonds):
        print(i, bond.show)

if __name__ == '__main__':
    FF_bonded_from_QM()
