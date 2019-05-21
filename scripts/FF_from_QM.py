#!/usr/bin/python3.5

from parsers import *
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
    parser.add_argument("freq", type=str, help="freq file: Formatted checkpoint file of "
                                                     "freq calculation G09 with \"freq=intmodes\" keyword")
    parser.add_argument("-s", type=str, help="structure file: an optional mol2 file from which the bonds will be read.")
    try:
        args = parser.parse_args()
        freq_file = args.freq
        mol2_file = args.s
    except:
        print('Use -h for flags and argument\n')
        return 1

    # Read the input frequency calculation checkpoint file
    freq = parseQM(freq_file, 'G')
    freq.internal_coord.bonds.sort()
    freq.internal_coord.angles.sort()
    freq.internal_coord.torsions.sort()

    print('------------------------------------------------------------')
    np.set_printoptions(precision=1, suppress=True, threshold=np.inf, linewidth=520)
    #print(freq.hessian_cartesian)
    #for i, atom in enumerate(freq.internal_coord.structure): print(atom.show)
    #for i, bond in enumerate(freq.internal_coord.bonds): print(i, bond.show)
    #for i, angle in enumerate(freq.internal_coord.angles): print(i, angle.show)
    #for i, torsion in enumerate(freq.internal_coord.torsions): print(i, torsion.show)
    for i, improper in enumerate(freq.internal_coord.impropers): print(i, improper.show)

    print('------------------------------------------------------------')

    freq.internal_coord.deletebonds()
    freq.internal_coord.deleteangles()
    freq.internal_coord.deletetorsions()
    freq.internal_coord.deleteimpropers()

    # If the structure file is given then read the bonds from mol2 file.
    if mol2_file:
        # Erease all bonded terms
        freq.internal_coord.deletebonds()
        freq.internal_coord.deleteangles()
        freq.internal_coord.deletetorsions()
        freq.internal_coord.deleteimpropers()

        # Read the mol2 file
        mol2 = parseMol2(mol2_file).getresult()
        # Check the atoms are the same
        if mol2.structure.atoms == freq.internal_coord.structure.atoms:
            # Overwrite the bonds
            freq.internal_coord.bonds = sorted(mol2.bonds)
        else:
            raise Warning('The atoms in mol2 file don\'t match the frequency file\n')

    # If the bonds are not available from file
    if freq.internal_coord.nbond == 0:
        # Attempt to generate the bonds
        print('Warning >>> No bond term was found, generating them internally. ')
        freq.internal_coord.constructbondsbyHessian(freq.hessian_cartesian)
        freq.internal_coord.constructbondsbyCovalentRadius()

    if freq.internal_coord.nangle == 0:
        # Attempt to generate the angles
        print('Warning >>> No angle term was found, generating them internally. ')
        freq.internal_coord.constructangles()

    if freq.internal_coord.ntorsion == 0:
        # Attempt to generate the torsions
        print('Warning >>> No torsion term was found, generating them internally. ')
        freq.internal_coord.constructtorsions()

    if freq.internal_coord.nimproper == 0:
        # Attempt to generate the impropers
        print('Warning >>> No improper term was found, generating them internally. ')
        freq.internal_coord.constructimpropers()

    GetBondedfromQM(freq)

    print('///////////////////////////////////////////////////////////////')
    #np.set_printoptions(precision=1, suppress=True, threshold=np.inf, linewidth=520)
    for i, atom in enumerate(freq.internal_coord.structure): print(atom.show)
    for i, bond in enumerate(freq.internal_coord.bonds): print(i, bond.show)
    for i, angle in enumerate(freq.internal_coord.angles): print(i, angle.show)
    for i, torsion in enumerate(freq.internal_coord.torsions): print(i, torsion.show)
    for i, improper in enumerate(freq.internal_coord.impropers): print(i, improper.show)
    print('///////////////////////////////////////////////////////////////')

    #freq.internal_coord._bonds = []
    #freq.internal_coord.constructbondsbyHessian(freq.hessian_cartesian)
    #freq.internal_coord.constructbondsbyCovalentRadius()

if __name__ == '__main__':
    FF_bonded_from_QM()
