import numpy as np
import sys


class RES(object):
    def __init__(self):
        self.coord = None
        self.name = None

class MOL(object):
    def __init__(self):
        self.coord = np.array(range(10))
        self.res = []
        self.frozen = []
        self.res.append(RES())
        self.res.append(RES())
        self.res.append(RES())

        self.res[0].name = "ALA"
        self.res[0].coord = self.coord[:2]
        self.res[1].name = "GLY"
        self.res[1].coord = self.coord[2:3]
        self.res[2].name = "PRO"
        self.res[2].coord = self.coord[3:4]
        self.frozen.append(self.res[0])
        self.frozen.append(self.res[2])

if __name__ == '__main__':
    A = MOL()
    print(A.coord)
    print(A.res[0].coord, A.res[1].coord)
    A.coord[0] = 10
    print(A.coord)
    print(A.res[0].coord, A.res[1].coord)
    print(np.may_share_memory(A.res[0].coord, A.coord))
    A.res[0].coord *= 5
    print(A.res[0].coord, A.res[1].coord)
    print(id(A.res[0].name), id(A.frozen[0].name))
    print(A.frozen[0].name, A.frozen[1].name)
    A.res[0].name = "ASP"
    print(A.frozen[0].name, A.frozen[1].name)
    A.frozen[1].name = "TRP"
    print(A.frozen[0].name, A.frozen[1].name)

    print(type(A.coord))

