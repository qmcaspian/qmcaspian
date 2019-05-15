#!/bin/bash/python3

from parse_qm import FreqFchkG09
from internal_coordinate import *
import numpy as np

class GetBondedfromQM(object):
    """
    This class accepts a freq object and calculates the force constances of the bonded terms
    """
    def __init__(self, freq_input):
        self.ifreq = freq_input

        # Calculate the force constances
        self.calculate_bond_forces()
        self.calculate_angle_forces()
        self.calculate_torsion_forces()

    def calculate_bond_forces(self):
        # The calculations for each bond will be done for both i-j and j-i)
        for bond in self.ifreq.internal_coord.bonds:

            # Get the indices of the partial hessian
            lbound_i = (bond.i.num - 1) * 3
            hbound_i = (lbound_i + 3)
            lbound_j = (bond.j.num - 1) * 3
            hbound_j = (lbound_j + 3)

            # Get the partial hessian.
            hij = self.ifreq.hessian_cartesian[lbound_i:hbound_i, lbound_j:hbound_j]
            hji = self.ifreq.hessian_cartesian[lbound_j:hbound_j, lbound_i:hbound_i]

            # Get the Eigen values and vectors
            wij, vij = np.linalg.eig(hij)
            wji, vji = np.linalg.eig(hji)

            #print('***', bond.i.num, bond.j.num, wij)
            # Get the bond unite vector
            uij = np.array(bond.j.cord) - np.array(bond.i.cord)
            uij /= np.linalg.norm(uij)

            uji = - uij

            # project the Eigen values onto the bond vector and sum over them.
            kij = 0.0
            for i in range(3):
                kij = kij + (wij[i] * np.abs(np.dot(uij, vij[:, i])))

            kji = 0.0
            for i in range(3):
                kji = kji + (wji[i] * np.abs(np.dot(uji, vji[:, i])))

            # Take the average of ij and ji calculations and convert to Oplsaa format (Not sure about the later!!!)
            k = -(kij + kji) / 2

            # assign the bond constant
            bond.f = k

            """
            np.set_printoptions(precision=5, suppress=True, threshold=np.inf, linewidth=520)
            print(bond.show)
            print(lbound_i, hbound_i)
            print(lbound_j, hbound_j)
            print('*************************************************************')
            print(hij)
            print('*************************************************************')
            print(hji)
            print('*************************************************************')
            print(wij)
            print(vij)
            print(wji)
            print(vji)
            print(uij)
            print(uji)
            print(kij)
            print(kji)
            print('k= ', k)
            print('--------------------------------------------------------------')
            """
            #print(kij)
            #print(kji)

    def calculate_angle_forces(self):
        # The calculations for each angle will be done for both i-j-k and k-j-i
        for angle in self.ifreq.internal_coord.angles:

            # Get the indices of the partial hessian
            lbound_i = (angle.i.num - 1) * 3
            hbound_i = (lbound_i + 3)
            lbound_j = (angle.j.num - 1) * 3
            hbound_j = (lbound_j + 3)
            lbound_k = (angle.k.num - 1) * 3
            hbound_k = (lbound_k + 3)

            # Get the partial hessian for i --> j <-- k, which gives uij, ukj
            hij = self.ifreq.hessian_cartesian[lbound_i:hbound_i, lbound_j:hbound_j]
            hkj = self.ifreq.hessian_cartesian[lbound_k:hbound_k, lbound_j:hbound_j]

            # Get the partial hessian  for i <-- j --> k ,which gives uji, ujk
            #hji = self.ifreq.hessian_cartesian[lbound_j:hbound_j, lbound_i:hbound_i]
            #hjk = self.ifreq.hessian_cartesian[lbound_j:hbound_j, lbound_k:hbound_k]


            # Get the Eigen values and vectors for i --> j <-- k, which gives uij, ukj
            wij, vij = np.linalg.eig(hij)
            wkj, vkj = np.linalg.eig(hkj)

            # Get the Eigen values and vectors for i <-- j --> k ,which gives uji, ujk
            #wji, vji = np.linalg.eig(hji)
            #wjk, vjk = np.linalg.eig(hjk)

            # Get the angle unit vectors for i --> j <-- k, which gives uij, ukj
            # Keep the bond lengths
            uij = np.array(angle.j.cord) - np.array(angle.i.cord)
            bondij = np.linalg.norm(uij)
            uij /= bondij

            ukj = np.array(angle.j.cord) - np.array(angle.k.cord)
            bondkj = np.linalg.norm(ukj)
            ukj /= bondkj

            # Get the angle unit vectors for i <-- j --> k ,which gives uji, ujk
            #uji = - uij
            #ujk = - ukj

            # Get the normal unit vector to the plane i --> j <-- k
            n_ijkj = np.cross(ukj, uij)
            n_ijkj /= np.linalg.norm(n_ijkj)

            # Get the normal unit vector to the plane i <-- j --> k
            #n_jijk = np.cross(ujk, uji)
            #n_jijk /= np.linalg.norm(n_jijk)

            # Get the projection unit vectors for i --> j <-- k
            pij = np.cross(n_ijkj, uij)
            pkj = np.cross(ukj, n_ijkj)

            # Get the projection unit vectors for i <-- j --> k
            #pji = np.cross(n_jijk, uji)
            #pjk = np.cross(ujk, n_jijk)

            # Get the force constances k projected on the p unit vectors
            kij = 0.0
            kkj = 0.0
            #kji = 0.0
            #kjk = 0.0
            for i in range(3):
                # for i --> j <-- k
                kij = kij + (wij[i]) * np.abs(np.dot(pij, vij[:, i]))
                kkj = kkj + (wkj[i]) * np.abs(np.dot(pkj, vkj[:, i]))

                # for i <-- j --> k
                #kji = kji + (wji[i]) * np.abs(np.dot(pji, vji[:, i]))
                #kjk = kjk + (wjk[i]) * np.abs(np.dot(pjk, vjk[:, i]))

            kij = kij * np.square(bondij)
            kkj = kkj * np.square(bondkj)
            #kji = kji * np.square(bondij)
            #kjk = kjk * np.square(bondkj)

            k_ijkj = (1 / kij) + (1 / kkj)
            k_ijkj = -1 / k_ijkj

            #k_jijk = (1 / kji) + (1 / kjk)
            #k_jijk = -1 / k_jijk

            # assign the angle force constant
            angle.f = k_ijkj

            """
            print(pij, pkj)
            print(pji, pjk)
            print('*************************************************************')
            print(kij)
            print(kkj)
            print(kji)
            print(kjk)
            print(bondij)
            print(bondkj)
            """
            #print(k_ijkj)
            #print(k_jijk)

    def calculate_torsion_forces(self):
        for torsion in self.ifreq.internal_coord.torsions:

            # Get the indices of the partial hessian
            lbound_i = (torsion.i.num - 1) * 3
            hbound_i = (lbound_i + 3)
            lbound_j = (torsion.j.num - 1) * 3
            hbound_j = (lbound_j + 3)
            lbound_k = (torsion.k.num - 1) * 3
            hbound_k = (lbound_k + 3)
            lbound_l = (torsion.l.num - 1) * 3
            hbound_l = (lbound_l + 3)


            # Get the partial hessian
            hij = self.ifreq.hessian_cartesian[lbound_i:hbound_i, lbound_j:hbound_j]
            hlk = self.ifreq.hessian_cartesian[lbound_l:hbound_l, lbound_k:hbound_k]

            # Get the Eigen values and vectors
            wij, vij = np.linalg.eig(hij)
            wlk, vlk = np.linalg.eig(hlk)


            # Get the angle unit vectors for i --> j <-- k, which gives uij, ukj
            #                                  and j -->,k <-- l, which gives ujk, ulk
            uij = np.array(torsion.j.cord) - np.array(torsion.i.cord)
            bondij = np.linalg.norm(uij)
            uij /= bondij

            ukj = np.array(torsion.j.cord) - np.array(torsion.k.cord)
            bondkj = np.linalg.norm(ukj)
            ukj /= bondkj

            ujk = np.array(torsion.k.cord) - np.array(torsion.j.cord)
            bondjk = np.linalg.norm(ujk)
            ujk /= bondjk

            ulk = np.array(torsion.k.cord) - np.array(torsion.l.cord)
            bondlk = np.linalg.norm(ulk)
            ulk /= bondlk

            ukl = np.array(torsion.l.cord) - np.array(torsion.k.cord)
            bondkl = np.linalg.norm(ukl)
            ukl /= bondkl

            # Get the normal unit vector to the plane for i --> j <-- k, which gives uij, ukj
            #                                         and j -->,k <-- l, which gives ujk, ulk
            n_ijk = np.cross(ukj, uij)
            n_ijk /= np.linalg.norm(n_ijk)

            n_jkl = np.cross(ulk, ujk)
            n_jkl /= np.linalg.norm(n_jkl)

            # Get the force constances k projected on the p unit vectors
            kijk = 0.0
            klkj = 0.0
            for i in range(3):
                # for i --> j <-- k
                kijk = kijk + (wij[i]) * np.abs(np.dot(n_ijk, vij[:, i]))

                # for j --> k <-- l
                klkj = klkj + (wlk[i]) * np.abs(np.dot(n_jkl, vlk[:, i]))

            # Calculate the final cross products
            c_ijjk = np.cross(uij, ujk)
            c_jkkl = np.cross(ujk, ukl)
            kijk = kijk * np.square(bondij) * np.dot(c_ijjk, c_ijjk)
            klkj = klkj * np.square(bondkl) * np.dot(c_jkkl, c_jkkl)


            k = (1 / kijk) + (1 / klkj)
            k = -1 / k
            #k_jijk = (1 / kji) + (1 / kjk)
            #k_jijk = -1 / k_jijk

            # assign the angle force constant
            torsion.f1 = float(k)
