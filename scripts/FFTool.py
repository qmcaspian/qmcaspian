#!/usr/bin/env python3

"""
    module :: FFTool
    :platform: Unix
"""

import numpy as np
import pickle
import re
import os
from atom import Atom
from molecule import Molecule
from atom_table import *
import itertools as it
from internal_coordinate import *
from graph import *
import parsers as ps
from colors import *

# TODO add a filter for the page options to limit the access to current page
# TODO read the first argument for the page options "So comments can be added to the input file"


class FFTool(object):

    def __init__(self):
        self.active_library = None
        self.libraries = dict()
        self.active_fragment = None
        self.fragments = dict()
        self.active_structure = None
        self.structures = dict()

    def analyze(self, flag, parm=None):
        switcher = {
                    'f': self.fragment_library_msg,
                    'fc': self.create_fragment_library,
                    'fp': self.print_libraries,
                    'fs': self.save_library,
                    'fl': self.load_library,
                    'ff': self.create_fragment_msg,
                    'ffl': self.load_structure,
                    'ffh': self.set_head_atom,
                    'ffp': self.prune_tree,
                    'ffc': self.change_to_wild_card,
                    'ffd': self.display_fragment,
                    'ffa': self.set_atom_type,
                    'ffs': self.save_fragment,
                    'ffn': self.name_fragment
                    }
        func = switcher.get(flag, lambda: 'Invalid')

        if parm is None:
            return func()
        else:
            return func(parm)

    def main_msg(self):

        input_msg = None

        main_page_options = """
        -------------------Main Page-------------------
        To proceed, choose one of the bellow options
        -----------------------------------------------
        q: Quit FFTool.
        f: Manage Fragment Library.
        -----------------------------------------------
        """
        # Page filter avoid accessing the option flags of other page
        option_flags = ['f', 'q']

        while input_msg != 'q':
            input_msg = input(main_page_options)
            if input_msg not in option_flags:
                print("Error >>> Not a valid option, try again.")
                continue #skip it
            self.analyze(input_msg)

    def fragment_library_msg(self):

        input_msg = None

        fragment_library_options = """
---------------------------Fragment Library Options---------------------------
q:  Quit this page.
fc: Create empty fragment library.
fp: Print the loaded libraries.
ff: Create new fragment.
fs: Save library.
fl: Load library.
------------------------------------------------------------------------------"""

        # Page filter avoid accessing the option flags of other page
        option_flags = ['fc','fp','ff','fs','fl', 'q']

        while input_msg != 'q':

            print(fragment_library_options)
            input_msg = input("Enter your option key: ")

            # Check the flag
            if input_msg not in option_flags:
                print("Error >>> Not a valid option, try again.")
                continue

            self.analyze(input_msg)

    def load_library(self):
        if self.active_library:
            overwrite_library = input('Warning >>> a library is already loaded, do you want to overwrite it? (Y/N) ')
            if not re.match(r'Y', str(overwrite_library), re.IGNORECASE):
                print("Info >>> Canceling operation")
                return None
            else:
                print("Info >>> overwriting existing fragment library")

        library_name = input('Type the library name: ')
        try:
            self.active_library = pickle.load(open(library_name.strip(), 'rb'))
        except FileNotFoundError as error:
            print("Error >>>", error.args[1], "'", library_name.strip(), "'")
            return None

    def save_library(self):
        if self.active_library is None:
            print("Error >>> No library is loaded.")
            return None
        library_name = input('Type the library name: ')
        pickle.dump(self.active_library, open(library_name.strip() + '.dat', 'wb'))

    def print_libraries(self):
        if self.active_library:
            print("Active Library: ", self.active_library.name)
            for element, fragments in self.active_library:
                print(element, '->', fragments)
        else:
            print(self.active_library)

    def create_fragment_library(self):

        # Warn if a library is loaded in the active memory.
        if self.active_library:
            overwrite_library = input('Warning >>> A library is already loaded, do you want to overwrite it?\n')
            if not re.match(r'Y', str(overwrite_library), re.IGNORECASE):
                print("Info >>> Canceling operation")
                return None
            else:
                print("Info >>> overwriting existing fragment library")

        # Get the library name.
        library_name = input("Enter the fragment library name: ")
        self.active_library = Fragment_Library(name=library_name)

        # Populate the library with the elements symbol as the keys.
        # Each key has a list that accommodates the atom types,
        # e.g. {'C':[(graph, {UFF:'C_3', OPLSAA:'C3'}), (graph, {UFF:'C_2', OPLSAA:'CD'}), ...]
        for i in range(1, 37):
            try:
                atomic_symbol = atom_property(query='symbol', target='number', value=i)
                self.active_library.element[atomic_symbol] = []
            except:
                print('Warning >>> No symbol was found for atomic number ', i)
        print('Info >>> An empty library is created for atomic number up to ', i, " ,and is loaded to active library.")

    #TODO Constructing a new fragment interactively.
    def create_fragment_msg(self):

        input_msg = None

        # Check for existing fragment.
        if self.active_fragment is not None:
            input_msg = input(" A fragment already exists. Do you want to overwrite: (Y/N) ")
            if re.match(r'N',input_msg, re.IGNORECASE):
                return None

        help_msg = """To create a new fragment
    1) Load a stricture file (mol2 format).
    2) Define the head atom, which represents the atom type.
    3) Modify the structure if necessary.
        a) Prune the structure to remove unnecessary atoms.
        b) Change atoms and bonds to wild cards "*" and "-1"
           respectively. These wild cards match anything.
    4) Set the force field specific atom type.
    5) Test against the existing fragment library to make sure
       the new fragment is not identical to the other entries.
    6) Save to the library, file, or memory."""

        fragment_options = """
------------------------------Fragment Options--------------------------------
q:   Quit this page.
h:   Print help.
ffl: Load structure.
ffd: Display fragment.
ffn: Set name.
ffh: Set head atom.
ffp: Prune fragment.
ffc: Set wild cards.
ffa: Set atom types.
ffs: Save fragment.
------------------------------------------------------------------------------"""

        # Page filter to avoid accessing the option flags of other pages.
        option_flags = ['ffl', 'ffh', 'ffp', 'ffd', 'ffa', 'ffc', 'ffs', 'ffn', 'q', 'h']

        print(help_msg)

        while input_msg != 'q':

            # Get the inout
            print(fragment_options)
            input_msg = input("Enter your option key: ")

            if input_msg not in option_flags:
                print("Error >>> Not a valid option, try again.")

            # Catch the page help.
            if re.match(r'h', input_msg, re.IGNORECASE):
                print(help_msg)
                continue

            # Clean the active_fragment before exit.
            if re.match(r'q', input_msg, re.IGNORECASE):
                if not self.active_fragment.tree:
                    input_msg = input("Warning >>> a fragment is loaded. Exiting would delete it. " +
                                                                                        "Do you want to exit? (Y/N) ")
                    if re.match(r'Y', input_msg, re.IGNORECASE):
                        self.active_fragment = None
                        break
                    else:
                        input_msg = "X"
                        continue

            self.analyze(input_msg)

    def name_fragment(self):

        # Check active memory
        if self.active_fragment is None:
            print("Error >>> no fragment is loaded.")
            return None

        self.active_fragment.name = input("Enter the fragment name: ")


    def save_fragment(self):

        # Check active memory
        if self.active_fragment is None:
            print("Error >>> no fragment is loaded.")
            return None

        # Check the minimum requirements of Fragment.
        if not self.active_fragment.fragment:
            print("Error >>> loaded fragment seems to be incomplete. Can not save it as a valid fragment yet.")
            return None

        # Get the save mode (f: mol2 file, m: memory)
        save_mode = input("Enter save mode, (\"f\" for saving to file and \"m\" for saving to memory): ")

        if re.match(save_mode, 'm', re.IGNORECASE):

            # Check memory for duplicate names
            if self.active_fragment.name in self.fragments:
                print("Error >>> a fragment with the same name is already in the memory. Change the name and try again.")
                return None

            # Save to memory
            self.fragments[self.active_fragment.name] = self.active_fragment

        elif re.match(save_mode, 'f', re.IGNORECASE):

            name = str(self.active_fragment.name).replace(" ", "")

            # Check for existing files
            if os.path.isfile(name.strip() + '.mol2'):

                overwrite = input("Warn >>> a fragment file with the same name exist. Do you want to overwite if? (Y/N) ")

                if re.match(overwrite, 'Y', re.IGNORECASE):
                    ps.write._graph2mol2(self.active_fragment, name)

                else:
                    "Info >>> canceling saving to file."
                    return None
            else:
                ps.write.mol2(self.active_fragment, name)

        else:
            print("Error >>> Invalid input.")
            return None





    def set_atom_type(self):

        # Check the memory
        if type(self.active_fragment) is not Fragment:
            print("Error >>> no fragment is load a fragment first.")
            return None

        # Check the fragment
        if not self.active_fragment.head:
            print("Error >>> the head atom of the loaded fragment is not set. Set the head atom first.")
            return None

        # Print instructions
        print("The UFF definition is the minimum requirement.  ")
        print("To set atom type enter ", red("\"forcefield  atom-type\". "))
        print("To finish. Press Enter or type" + red(" \"q\" "))

        input_msg = ''
        while True:

            # Get input
            input_msg = input()
            if input_msg == 'q' or input_msg == 'Q' or input_msg == '':
                break

            input_msg = input_msg.split()

            # Check the args.
            if len(input_msg) == 2:

                ff, atom_type = input_msg
                if ff not in Fragment.ForceFileds:
                    print("Info >>> the ", ff, " is not an internal force field.")

                # Set the ff and atom type.
                self.active_fragment.atom_type(forcefield=ff, value=atom_type)

            # Default
            else:
                print("Error >>> invalid input.")
                continue

        self.fragment_summary(self.active_fragment)

    def display_fragment(self):
        self.fragment_summary(self.active_fragment)

    def change_to_wild_card(self):
        # Rename the atoms and bonds of a fragment to wild cards.
        # Check the memory
        if type(self.active_fragment) is not Fragment:
            print("Error >>> no fragment is load a fragment first.")
            return None

        # Check the fragment
        if not self.active_fragment.head:
            print("Error >>> the head atom of the loaded fragment is not set. Set the head atom first.")
            return None

        print("To modify an atom label enter" + red(" \"atom-number\" "))
        print("To modify a bond order enter" + red(" \"atom-number  atom-number\" "))
        print("To finish. Press Enter or type" + red(" \"q\" "))

        chang_input = ''
        while True:
            chang_input = input()
            if chang_input == 'q' or chang_input == 'Q' or chang_input == '':
                break

            chang_input = chang_input.split()

            # Change the atom label.
            if len(chang_input) == 1:

                # Check the input and
                try:
                    atom_num = int(chang_input[0])
                except ValueError:
                    print("Error >>> invalid input.")
                    continue

                # Get the correct node
                node = self.active_fragment.select_by_atom_num(atom_num)

                if node is None:
                    print("Error >>> atom ", atom_num, "  was not found.")
                    continue

                # Change the atom label and symbol.
                node.atom.typ = 'R'
                node.atom.nam = 'R'

            # Change the bond
            elif len(chang_input) == 2:

                # Check the input and
                try:
                    atom_i, atom_j = int(chang_input[0]), int(chang_input[1])
                except ValueError:
                    print("Error >>> invalid input.")
                    continue

                # Get the correct nodes
                node_i = self.active_fragment.select_by_atom_num(atom_i)
                node_j = self.active_fragment.select_by_atom_num(atom_j)

                if node_i is None:
                    print("Error >>> atom ", atom_i, " was not found.")
                    continue

                if node_j is None:
                    print("Error >>> atom ", atom_j, " was not found.")
                    continue

                if node_i not in node_j.neighbors:
                    print("Error >>>  bond ", node_i.atom.num, node_j.atom.num, " was not found.")
                    continue

                if node_j not in node_i.neighbors:
                    print("Error >>>  bond ", node_j.atom.num,  node_i.atom.num, " was not found.")
                    continue

                # Change the bond order
                node_i.bond_order(node_j, -1.0)
                node_j.bond_order(node_i, -1.0)

            # Default
            else:
                print("Error >>> invalid input.")
                continue

        # Show a summary
        self.fragment_summary(self.active_fragment)

    def prune_tree(self):

        # Check for a loaded graph
        if type(self.active_fragment) is not Fragment:
            print("Error >>> no structure is loaded. First, load a structure.")
            return None

        # Check if its tree
        if not self.active_fragment.tree:
            print("Error >>> The head atom of structure is not set. Set the head atom first.")
            return None

        # Get a list of atom-numbers.
        prune_list = input("Enter the atom-numbers to be deleted from fragment: ")
        prune_list = prune_list.split()

        # prune tree
        for atom_num in prune_list:
            # Get the node
            node = self.active_fragment.select_by_atom_num(int(atom_num))

            # Check the nude still exist. In the case of multiple atom pruning. The node might be already deleted.
            if node is None:
                print("Info >>> The atom ", atom_num, " seems to be already deleted.")
                continue # Skip

            # Check for head
            if node == self.active_fragment.head:
                print("Error >>> the head atom can not be pruned.")
                return None

            self.active_fragment.prune(node)

        # Show a summary
        self.fragment_summary(self.active_fragment)

    def set_head_atom(self):

        # Check for a loaded graph
        if type(self.active_fragment) is not Fragment:
            print("Error >>> no structure is loaded. First, load a structure file.")
            return None

        # Check for the previous tree
        if self.active_fragment.head:
            input_msg = input("Warning >>> a head atom is already defined. Do you want to overwrite it? Y/N ")
            if re.match(r'Y', input_msg):
                self.active_fragment.reset_tree()
            else:
                return None

        # Get the head atom number
        head_atom_num = input("Enter the atom number of the head atom: ")

        # Check input
        try:
            head_atom_num = int(head_atom_num)
        except ValueError:
            print("Error >>> invalid input. The input expected to be an integer.")
            return None

        # Find the the head atom and set it.
        found = 0
        for node in self.active_fragment.nodes:
            if node.atom.num == head_atom_num:
                self.active_fragment.construct_tree(node)
                found += 1

        # Check for repeated number
        if found > 1:
            print("Error >>> multiple atoms has the same atom-number. Revise the structure file.")
            self.active_fragment.reset_tree()
            return None
        elif found == 0:
            print("Error >>> atom-number out of range.")
            self.active_fragment.reset_tree()
            return None

        # Show a summary
        self.fragment_summary(self.active_fragment)

        print("Info >>> to proceed, prune the structure if necessary, or set the atom type definition.")

    def load_structure(self):

        # Get the load mode flag (f/g)
        load_mode = input("Enter the load mode (\"f\" to load from file and \"m\" to load from memory): ")

        # Check the input
        if load_mode.strip() not in ['f', 'F', 'm', 'M']:
            print("Error >>> invalid input, ", load_mode)
            return None

        if load_mode.strip() in ['f', 'F']:

            # Read from file
            file_name = input("Enter the file name (mol2 format): ")

            if not os.path.isfile(file_name):
                print("Error >>> could not find the file " + str(file_name))
                return None

            mol = ps.read.mol2(file_name)

        elif load_mode.strip() in ['m', 'M']:

            # Load from memory
            file_name = input("Enter the structure name: ")

            if file_name not in self.structures:
                print("Error >>> could not find the structure file.")
                return None

            mol = self.structures[file_name].deepcopy()

        else:
            print("Error >>> invalid input.")
            return None


        print("Info >>> the loaded structure name is: ", str(mol.nam))
        name = input("Enter an name if you want to change it of press Enter to keep the current name: ")
        mol.nam = str(name)

        # Convert it to Fragment and save it to active memory.
        self.active_fragment = Fragment(mol, mol.nam)

        # Display a summary of structure.
        self.fragment_summary(self.active_fragment)

        print("Info >>> to proceed choose the head atom. ")

    def fragment_summary(self, fragment):

        # Check the type
        if type(fragment) is not Fragment:
            print("Error >>> a graph object is expected, a ", type(fragment), " was found.")
            return None


        # print heading.
        if fragment.head:
            print("Info >>> Loaded Fragment:",)
            print("Info >>> Name: ", fragment.name)
            print(blue("Info >>> head atom number: " + str(fragment.head.atom.num) + ", name: " + str(fragment.head.atom.nam)))

            if len(fragment.atom_type()) > 0:
                print("Info >>> Atom type definition: ")
                print("         Force field        Atom type ")
                for ff, atom_type in self.active_fragment.atom_type():
                    print('             {0:5s}              {1:7s}'.format(ff, atom_type))
        else:
            print("Info >>> Loaded structure:")

        print('----------------------------------------------------------------------')
        print(' index    number    atom name   #neighbor    ring membership    atomic charge')
        for i, node in enumerate(fragment.nodes):
            if len(node.rings) > 0:
                ring = True
            else:
                ring = False
            if node == fragment.head:
                print(blue(' {0:3d}       {1:3d}         {2:4s}       {3:3d}             {4:5s}           {5: .4f}'.format(i,
                                            node.atom.num, node.atom.nam, node.number_of_neighbor, str(ring), node.atom.charge)))
            else:
                print(' {0:3d}       {1:3d}         {2:4s}       {3:3d}             {4:5s}           {5: .4f}'.format(i,
                                             node.atom.num, node.atom.nam, node.number_of_neighbor, str(ring), node.atom.charge))
        print('\n index  atom i  atom j  bond order')
        for i, bond in enumerate(self.active_fragment.bonds):
            print(' {0:3d}     {1:3d}  -->{2:3d}      {3: .1f}'.format(i, bond.i.num, bond.j.num, bond.order))
        print('----------------------------------------------------------------------')



if __name__ == '__main__':
    A = FFTool()
    A.main_msg()



