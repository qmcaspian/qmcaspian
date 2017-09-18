#!/usr/bin/env python3

import argparse
from parse_qm import parseQM


def fh_correction():
    """
    This script reads a Gaussian optimization output and uses the last structure in the file to generate the input files
    for frq, big basis set and solvation corrections.

    * How to set up:

        export the library directories to PYTHONPATH. In my case it is

        **>>> export PYTHONPATH=$PYTHONPATH:/home/masoud/QMcaspian/src/**

        export the scripts directory to PATH

        **>>> export PATH=$PATH:/home/masoud/QMcaspian/scripts/**

    * How to use:

        *>>> fh_correction.py [-n (int) (structure number to be used)] <inputfile>*

    """
    # Structure number to be used for generation of the input files
    structure_num = None

    # Read the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input file")
    parser.add_argument("-n", type=int, metavar="", help="Structure number to be used for generating input")
    args = parser.parse_args()
    ifile = args.input_file
    structure_num = args.n

    # Read the input file
    opt = parseQM(ifile, 'G')

    # Check the convergence
    if not opt.converged:
      print("Warning, the optimization is not converged")

    # Print the structure number that is used
    if structure_num == None:
        print("Using the last structure for correction files: " + str(opt.nmol))
        structure_num = opt.nmol

    elif structure_num > opt.nmol or structure_num <= 0:
        print("Warning: The given structure number is out of range")
        print("         Using the last structure for correction files: " + str(opt.nmol))
        structure_num = opt.nmol

    else:
        print("Structure " + str(structure_num) + " is used for input generation.")


    # Write the frq
    out_filename = ifile.split('.')[0] + '_frq.inp'
    out_f = open(out_filename, 'w')
    out_f.write('%nprocshared=16\n' + '%mem=16GB\n' + '#p freq b3lyp/6-31g(d,p) empiricaldispersion=gd3\n')
    out_f.write(' \n' + 'Title Card Required\n' + ' \n' + str(opt.charge) + '  ' + str(opt.multiplicity) + '\n')
    for atom in opt.selectbyMolnum(opt.nmol):
        out_f.write('%2s     %10f   %10f   %10f \n' % (atom.nam, atom.cord[0], atom.cord[1], atom.cord[2]))
    out_f.write(' \n' + ' \n')
    out_f.close()

    # Write the big basis set
    out_filename = ifile.split('.')[0] + '_bbs.inp'
    out_f = open(out_filename, 'w')
    out_f.write('%nprocshared=16\n' + '%mem=16GB\n' + '#p b3lyp/6-311+g(2d,2p) empiricaldispersion=gd3\n')
    out_f.write(' \n' + 'Title Card Required\n' + ' \n' + str(opt.charge) + '  ' + str(opt.multiplicity) + '\n')
    for atom in opt.selectbyMolnum(opt.nmol):
        out_f.write('%2s     %10f   %10f   %10f \n' % (atom.nam, atom.cord[0], atom.cord[1], atom.cord[2]))
    out_f.write(' \n' + ' \n')
    out_f.close()

    # Write the solvation
    out_filename = ifile.split('.')[0] + '_solv.inp'
    out_f = open(out_filename, 'w')
    out_f.write('%nprocshared=16\n' + '%mem=16GB\n' + '#p b3lyp/6-31g(d,p) scrf=(solvent=water,read,smd) empiricaldispersion=gd3\n')
    out_f.write(' \n' + 'Title Card Required\n' + ' \n' + str(opt.charge) + '  ' + str(opt.multiplicity) + '\n')
    for atom in opt.selectbyMolnum(opt.nmol):
        out_f.write('%2s     %10f   %10f   %10f \n' % (atom.nam, atom.cord[0], atom.cord[1], atom.cord[2]))
    out_f.write(' \n')
    out_f.write('eps=4')
    out_f.write(' \n' + ' \n')
    out_f.close()
    return 0

if __name__ == '__main__':
    fh_correction()







































if __name__ == '__main__':
    pass