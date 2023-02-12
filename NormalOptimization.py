import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
plt.close("all")

def N(x):
    return np.exp(-x ** 2 / 2)

def f(x, x0):
    a = 0
    c = 2 * x0
    B = 0
    
    Filter = (x < x0)
    Array = np.empty_like(x)
    Array[Filter] = 1 - a / 2 * x[Filter] ** 2
    Array[~Filter] = np.exp(-c * (x[~Filter] - x0) + B)
    return Array

def Target(x0):
    a = 0
    c = 2 * x0
    B = 0

    return x0 - a / 6 * x0 ** 3 + 1 / c * np.exp(B)

x0 = 1 / np.sqrt(2)

x = np.linspace(0, 3, 1000)

plt.figure()
plt.plot(x, N(x))
plt.plot(x, f(x, x0))

Count = 100000
xScan = np.linspace(0, 10, Count)
x0Scan = np.linspace(0, 2, 1000)
I = np.empty(1000)

for i, x0 in enumerate(x0Scan):
    I[i] = np.sum(f(xScan, x0)) * 10 / Count
    
plt.figure()
plt.plot(x0Scan, I)

print(x0Scan[np.argmin(I)])
print(np.min(I))

Value = opt.minimize(Target, 1.15)
print(opt.minimize(Target, 1.15))


def f(x, x0):
    a = 0
    c = x0
    B = -x0 ** 2 / 2
    
    Filter = (x < x0)
    Array = np.empty_like(x)
    Array[Filter] = 1 - a / 2 * x[Filter] ** 2
    Array[~Filter] = np.exp(-c * (x[~Filter] - x0) + B)
    return Array

x0 = 1.4

x = np.linspace(0, 3, 1000)

plt.figure()
plt.plot(x, N(x))
plt.plot(x, f(x, x0))

Count = 100000
xScan = np.linspace(0, 10, Count)
x0Scan = np.linspace(0, 2, 1000)
I = np.empty(1000)

for i, x0 in enumerate(x0Scan):
    I[i] = np.sum(f(xScan, x0)) * 10 / Count
    
plt.figure()
plt.plot(x0Scan, I)

print(x0Scan[np.argmin(I)])
print(np.min(I))