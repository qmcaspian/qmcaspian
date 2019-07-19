#!/usr/bin/env python3

"""
    module :: FFTool
    :platform: Unix
"""

import numpy as np
import pickle
from atom import Atom
from molecule import Molecule
from atom_table import *
import itertools as it
from internal_coordinate import *


class FFTool(object):
    def __init__(self):
        self.current_library = None
        self.temp_library = None

    def analyze(self, arg):
        switcher = {
                    'f': self.fragment_library_msg,
                    'fc': self.create_fragment_library,
                    'fp': self.print_libraries,
                    'fs': self.save_library,
                    'fl': self.load_library
                    }
        func = switcher.get(arg, lambda:'Invalid')
        func()

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
        while input_msg != 'q':
            input_msg = input(main_page_options)
            self.analyze(input_msg)

    def fragment_library_msg(self):
        input_msg = None
        current_library = None
        fragment_library_options = """
        -----------Fragment Library Options------------
        q:  Quit this page.
        fc: Create a new fragment library.
        fp: Print the loaded libraries.
        fa: Add new fragments to the library.
        fs: Save library.
        fl: Load library.
        -----------------------------------------------
        """
        while input_msg != 'q':
            input_msg = input(fragment_library_options)
            self.analyze(input_msg)

    def load_library(self):
        library_name = input('Type the library name:\n')
        self.current_library = pickle.load(open(library_name.strip(), 'rb'))

    def save_library(self):
        library_name = input('Type the library name:\n')
        pickle.dump(self.current_library, open(library_name.strip() + '.dat', 'wb'))

    def print_libraries(self):
        if type(self.current_library) is dict:
            for key, value in self.current_library.items():
                print(key, '->', value)
        else:
            print(self.current_library)

        if type(self.temp_library) is dict:
            for key, value in self.temp_library.items():
                print(key, '->', value)
        else:
            print(self.temp_library)

    def create_fragment_library(self):
        self.current_library = dict()
        # Populate the library with the elements symbol as the keys.
        # Each key has a list that accommodates the atom types,
        # e.g. {'C':[(graph, {UFF:'C_3', OPLSAA:'C3'}), (graph, {UFF:'C_2', OPLSAA:'CD'}), ...]
        for i in range(1, 37):
            try:
                atomic_symbol = atom_property(query='symbol', target='number', value=i)
                self.current_library[atomic_symbol] = []
            except:
                print('Warning >>> No symbol was found for atomic number ', i)
        print('Info >>> An empty library is created for atomic number up to ', i)
        print('Info >>> The new library is loaded to the main cash')

if __name__ == '__main__':
    A = FFTool()
    A.main_msg()



