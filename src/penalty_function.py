#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
x0 = 1
a = 0.1
b = 0.5
b2 = np.square(b)

x = np.arange(x0 - 2, x0 + 2, 0.1)


def um(x):
    x2 = np.square(x)
    dx = x0 - x
    dx2 = np.square(dx)
    return a * (np.sqrt(dx2 + b2) - b)

def dum(x):
    x2 = np.square(x)
    dx = x0 - x
    dx2 = np.square(dx)
    return -a * (dx / np.sqrt(dx2 + b2))


def u0(x):
    x2 = np.square(x)
    return a * (np.sqrt(x2 + b2) - b)

def du0(x):
    x2 = np.square(x)
    return a * (x / np.sqrt(x2 + b2))

def h(x):
    return a * np.square(x0 - x)

def dh(x):
    return -2 * a * (x0 - x)

print('um ', um(x0))
print('u0 ', u0(0))
print('h ', h(x0))


#for a in np.arange(0.0, 1, 0.005):
#    plt.plot(x, a * np.sqrt(np.square(x - x0) + b*b) - b)

#plt.plot(x, 0.0005 * np.sqrt(np.square(x - x0) + b*b) - b, 'b')
plt.plot(x, um(x), '.g')
plt.plot(x, dum(x), 'g')

plt.plot(x, u0(x), '.b')
plt.plot(x, du0(x), 'b')

plt.plot(x, h(x), '.r')
plt.plot(x, dh(x), 'r')

#plt.plot(x, 0.5 * np.square(x - x0))
plt.show()