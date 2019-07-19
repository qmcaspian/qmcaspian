import numpy as np

import time


class System(object):
    def __init__(self):
        self.natom = 10
        self.coord = np.random.rand(self.natom * 3) * 10

if __name__ ==  '__main__':

    A = System()
    start = time.time()

    for i in range(A.natom-1):
        for j in range(i+1, A.natom):
            r = A.coord[j:j+3] - A.coord[i:i+3]
            print(np.dot(r, r))

    finish = time.time()
    print('time: ', finish-start)
