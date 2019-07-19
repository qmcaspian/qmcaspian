import numpy as np
import sumc as s
import time

def sumP(M):
    tmp = 0
    for row in M:
        for j in row:
            tmp += j
    return tmp

def sumPMatrices(M1, M2):
    tmp_row = []
    tmp_matrix = []
    for i, j in zip(M1, M2):
      for r, s in zip(i, j):
        tmp_row.append((r+s))
      tmp_matrix.append(tmp_row)
    return tmp_matrix

def benchmark():
    # A numpy matrix
    A = np.array(np.random.random(size=(5000, 5000)), dtype=np.double)
    # Apy the corresponding python nested list
    Apy = A.tolist()
    time_NumpySum = 0
    time_NumpySumMatrix = 0
    time_sumL = 0
    time_sumE = 0
    time_sumEMatrices = 0
    time_PySUM = 0
    time_sumPMatrices = 0
    n = 5

    for _ in range(n):
        # Numpy. sum of elements using numpy
        sum = np.double(0.0)
        start = time.time()
        sum = np.sum(A)
        end = time.time()
        time_NumpySum += (end - start)

        # sumL. sum of elements using pointers arithmatic
        sum = np.double(0.0)
        start = time.time()
        sum = s.sumL(A)
        end = time.time()
        time_sumL += (end - start)
        
        # sumE. sun of elements using EIGEN sum
        sum = np.double(0.0)
        start = time.time()
        sum = s.sumE(A)
        end = time.time()
        time_sumE += (end - start)

        # Pyrhon. The matrix is row major.
        start = time.time()
        sum = sumP(Apy)
        end = time.time()
        time_PySUM += (end - start)
        
        # Numpy. sum of two matrix
        sum_matrix = np.zeros((5000, 5000), dtype=np.double)
        start = time.time()
        sum = A + A
        end = time.time()
        time_NumpySumMatrix += (end - start)

        # sumEMatrices. sum of two matrix using EIGEN native function
        sum_matrix = np.zeros((5000, 5000), dtype=np.double)
        start = time.time()
        sum = s.sumEMatrices(A, A)
        end = time.time()
        time_sumEMatrices += (end - start)
        
        # sumPMatrices. sum of two matrix using python native function
        start = time.time()
        sum = sumPMatrices(Apy, Apy)
        end = time.time()
        time_sumPMatrices += (end - start)

    print('---------------------------------------')
    print('Ave time for NumpySum : ', time_NumpySum / n)
    print('Ave time for sumL : ', time_sumL / n)
    print('Ave time for sumE : ', time_sumE / n)
    print('Ave time for PySUM : ', time_PySUM / n)
    
    print('Ave time for NumpySumMatrix : ', time_NumpySumMatrix / n)
    print('Ave time for sumEMatrices : ', time_sumEMatrices / n)
    print('Ave time for sumPMatrices : ', time_sumPMatrices / n)
    
    print('---------------------------------------')
    print('relative ave time NumpySum/NumpySum', int(time_NumpySum/time_NumpySum))
    print('relative ave time sumL/NumpySum', int(time_sumL/time_NumpySum))
    print('relative ave time sumE/NumpySum', int(time_sumE/time_NumpySum))
    print('relative ave time PySUM/NumpySum', int(time_PySUM/time_NumpySum))
    
    print('relative ave time NumpySumMatrix/NumpySumMatrix', int(time_NumpySumMatrix/time_NumpySumMatrix))
    print('relative ave time sumEMatrices/NumpySumMatrix', int(time_sumEMatrices/time_NumpySumMatrix))
    print('relative ave time sumPMatrices/NumpySumMatrix', int(time_sumPMatrices/time_NumpySumMatrix))
    print('---------------------------------------')


if __name__ == '__main__':

    # Run benchmark
    benchmark()


