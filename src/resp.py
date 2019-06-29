#!/usr/bin/env python3
from molecule import Molecule
import numpy as np
from scipy.spatial import distance
import atom_table as at


class ESP(Molecule):
    COULOMBCONSTANT = 332.0636
    CONVERGED = 1.0e-8

    def __init__(self):

        Molecule.__init__(self)
        self.max_cycle = 1000
        self.points = None
        self.charge = None
        self.restraints = None
        self.symmetry = None

    def calculate_charges2(self, norest=False):

        if not self._check_inputs():
            return None

        # Get the dimensions of the charges and ESP points matrices
        # q run over the atomic charges and p run over the ESP points
        dim_q = self.natm
        dim_p = self.points.shape[0]
        dim_c = 1
        a_w = 0.0005
        a_s = 0.01
        b = 0.1

        q0 = np.zeros((dim_q, 1))
        Rq = np.array([atom.cord for atom in self.atoms])
        Rp = self.points[:, 1:]
        Rqp = np.array(distance.cdist(Rp, Rq))
        Rqp = 1 / Rqp
        V = np.array(self.points[:, 0]).reshape((dim_p, 1)) / self.COULOMBCONSTANT

        A = np.dot(Rqp.transpose(), Rqp)
        B = np.dot(Rqp.transpose(), V)

        # A and B with total charge restraints for charges
        Awf = np.ones((dim_q + dim_c, dim_q + dim_c))
        Awf[:dim_q, :dim_q] = A
        Awf[-1, -1] = 0

        Bwf = np.zeros((dim_q + dim_c, 1))
        Bwf[:dim_q] = B
        Bwf[-1] = self.charge

        # Count restrains
        if self.restraints:
            dim_c += len(self.restraints)
        if self.symmetry:
            dim_c += len([atom for row in self.symmetry for atom in row[1:]])

        # A and B with all restraints for charges
        Asc = np.zeros((dim_q + dim_c, dim_q + dim_c))
        Asc[:dim_q + 1, :dim_q + 1] = Awf


        Bsc = np.zeros((dim_q + dim_c, 1))
        Bsc[:dim_q + 1] = Bwf

        # Get the free charges
        qwf = self._optimize(Awf, Bwf, q0, a_w, b, dim_q)

        # TODO using the free charges to add restrain

        # Total charge restraint is already added
        res_count = 1
        if self.restraints:
            for atom in self.restraints:
                # find the locations
                row_index = dim_q + res_count
                col_index = atom[0] - 1
                res_count += 1

                # Assign
                Asc[row_index, col_index] = 1
                Bsc[row_index] = atom[-1]

        if self.symmetry:
            for row in self.symmetry:
                atom_i = row[0]
                for atom_j in row[1:]:
                    # find the locations
                    row_index = dim_q + res_count
                    col_index1 = atom_i - 1
                    col_index2 = atom_j - 1
                    res_count += 1

                    # Assign
                    Asc[row_index, col_index1] =  1
                    Asc[row_index, col_index2] = -1
                    Bsc[row_index] = 0

        # Complete the matrix
        Asc[:, dim_q+1:] = Asc[dim_q+1:, :].transpose()

        # Get the free charges
        qsc = self._optimize(Asc, Bsc, q0, a_s, b, dim_q)

        np.set_printoptions(precision=7, suppress=True, threshold=np.inf, linewidth=520)
        for atom, charge_free, charge_rest in zip(self.atoms, qwf, qsc):
            print(atom.num, atom.nam, atom.charge, '->', charge_free, '::', charge_rest )

        coord = np.array([atom.cord for atom in self.atoms])
        mass = np.array([at.atom_property(query='mass', target='symbol', value=atom.nam) for atom in self.atoms]).reshape(dim_q, 1)
        
        # TODO dipole function can not be dependent on the mass. 
        # Change to accept a point (x,y,z) as the refference for dipole calculations
        dipole_wf = self._dipole(coord, qwf, mass)
        dipole_sc = self._dipole(coord, qsc, mass)

        qesp = np.array([atom.charge for atom in self.atoms])
        dipole_esp = self._dipole(coord, qesp, mass)


    def _optimize(self, A, B, q0, a, b, dim_q):

        q_current = np.zeros((dim_q, 1))
        q_previous = np.zeros((dim_q, 1))
        converged = False

        for i in range(self.max_cycle):

            # Update restraints
            b2 = b * b
            q2 = np.square(q_current)
            aq = a * q_current
            restraints = np.multiply(aq, 1 / np.sqrt(q2 + b2))

            # Update A with new restraints contributions
            A_rest = np.array(A)
            A_rest[range(dim_q), range(dim_q)] = A_rest[range(dim_q), range(dim_q)] + restraints.transpose()

            # Update B with new restraints contributions
            B_rest = np.array(B)
            B_rest[:dim_q, :] = B_rest[:dim_q, :] + np.multiply(q0, restraints)

            # Calculate charges
            q_previous = q_current

            try:
                q_current = np.dot(np.linalg.inv(A_rest), B_rest)[:dim_q]

            except np.linalg.LinAlgError:

                print('Warning >>> the A matrix is singular, revise the constraints')
                q_current = np.dot(np.linalg.pinv(A_rest), B_rest)[:dim_q]

            q_diff = np.linalg.norm(q_current - q_previous)

            if q_diff < self.CONVERGED:
                converged = True
                break

        if not converged:
            print("Warning >>> RESP calculations did not converged.")
        else:
            print("Info >>> RESP calculations converged in ", i, "steps.")

        return q_current

    def _check_inputs(self):
        input_flag = False

        # Check the minimum inputs requirements
        if (self.natm > 1) and (self.points is not None) and (self.charge is not None):
            input_flag = True
        else:
            print('Error >>> charges cant be calculated, set the atoms, ESP points, and the formal charge')
            input_flag = False

        # Check the charge restraints
        if self.restraints:
            for restraint in self.restraints:
                if restraint[0] not in [atom.num for atom in self.atoms]:
                    print('Error >>> atom ', restraint, 'specified in charged CONSTRAINTS section does not exits')
                    input_flag = False

        # Check the symmetry constraints
        if self.symmetry:
            for row in self.symmetry:
                for atom_num in row:
                    if atom_num not in [atom.num for atom in self.atoms]:
                        print('Error >>> atom ', atom_num, 'specified in charged SYMMETRY section does not exits')
                        input_flag = False

        # Check charge restraints on symmetric atoms
        if self.restraints and self.symmetry:
            restrained_atoms = [rest[0] for rest in self.restraints]

            for row in self.symmetry:
                charges_of_symmetric_atoms = []

                for atom in row:
                    if atom in restrained_atoms:
                        index = restrained_atoms.index(atom)
                        charges_of_symmetric_atoms.append(self.restraints[index][2])

                if not all([i == charges_of_symmetric_atoms[0] for i in charges_of_symmetric_atoms]):
                    print('Error >>> symmetric atoms ', row, 'can not have different target charge restraints')
        return input_flag

    def _dipole(self, coord, charge, mass=None):

        if mass is None:
            center = np.sum(coord, axis=0) / coord.shape[0]
            center = center.reshape(1, 3)
            #print(center)

        if mass is not None:
            center = np.dot(coord.transpose(), mass) / np.sum(mass)
            center = center.reshape(1, 3)
            #print(center)

        dipole = np.zeros((1,3))
        for i in range(coord.shape[0]):
            dipole += (coord[i, :] - center) * charge[i]
        dipole /= 0.20819434
        print(dipole)
        print(np.linalg.norm(dipole))

        return dipole
