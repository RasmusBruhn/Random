import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as opt
plt.close("all")

def f(n, N, p):
    return np.exp(sp.loggamma(N + 1) - sp.loggamma(n + 1) - sp.loggamma(N - n + 1) + n * np.log(p) + (N - n) * np.log(1 - p))

def P(n, N, p, mu, a, b):
    n_m = np.ceil((N * p - np.exp(a) * (1 - p)) / (p + np.exp(a) * (1 - p)))
    n_p = np.ceil((N * p - np.exp(-b) * (1 - p)) / (p + np.exp(-b) * (1 - p)))
    A = a * np.exp(sp.loggamma(N + 1) - sp.loggamma(n_m + 1) - sp.loggamma(N - n_m + 1) - a * (n_m - mu) + np.log(p) * n_m + np.log(1 - p) * (N - n_m)) / (np.exp(a) - 1)
    B = b * np.exp(sp.loggamma(N + 1) - sp.loggamma(n_p + 1) - sp.loggamma(N - n_p + 1) + b * (n_p - mu) + np.log(p) * n_p + np.log(1 - p) * (N - n_p)) / (1 - np.exp(-b))

    Result = np.empty_like(n, dtype = float)
    
    if N * p < 1:
        return B / b * np.exp(-b * (n - mu)) * (1 - np.exp(-b))
    
    elif N * p > N - 1:
        return A / a * np.exp(a * (n - mu)) * (np.exp(a) - 1)
    
    F_m = n < np.floor(mu)
    F_p = n > np.floor(mu)
    F_l = (~F_m & ~F_p)
    Result[F_m] = A / a * np.exp(a * (n[F_m] - mu)) * (np.exp(a) - 1)
    Result[F_p] = B / b * np.exp(-b * (n[F_p] - mu)) * (1 - np.exp(-b))
    Result[F_l] = A / a * (1 - np.exp(a * (n[F_l] - mu))) + B / b * (1 - np.exp(-b * (n[F_l] + 1 - mu)))
    return Result

def Puse(n, N, p):
    if N * p < 1:
        return P(n, N, p, 1, 1 / (np.sqrt(1 - 1 / N) - np.log(1 / N)), 1 / (np.sqrt(1 - 1 / N) - np.log(1 - 1 / N)))
    
    if N * p > N - 1:
        return P(n, N, p, N - 1, 1 / (np.sqrt(1 - 1 / N) - np.log(1 - 1 / N)), 1 / (np.sqrt(1 - 1 / N) - np.log(1 / N)))
    
    else:
        return P(n, N, p, N * p, 1 / (np.sqrt(N * p * (1 - p)) - np.log(p)), 1 / (np.sqrt(N * p * (1 - p)) - np.log(1 - p)))

def I(Params, N, p, mu):
    a, b = Params
    n_m = np.ceil((N * p - np.exp(a) * (1 - p)) / (p + np.exp(a) * (1 - p)))
    n_p = np.ceil((N * p - np.exp(-b) * (1 - p)) / (p + np.exp(-b) * (1 - p)))
    A = a * np.exp(sp.loggamma(N + 1) - sp.loggamma(n_m + 1) - sp.loggamma(N - n_m + 1) - a * (n_m - mu) + np.log(p) * n_m + np.log(1 - p) * (N - n_m)) / (np.exp(a) - 1)
    B = b * np.exp(sp.loggamma(N + 1) - sp.loggamma(n_p + 1) - sp.loggamma(N - n_p + 1) + b * (n_p - mu) + np.log(p) * n_p + np.log(1 - p) * (N - n_p)) / (1 - np.exp(-b))
    
    if N * p < 1:
        return B / b * np.exp(b * mu) * (1 - np.exp(-b * (N + 1)))
    
    if N * p > N - 1:
        return A / a * (np.exp(a * (N + 1 - mu)) - np.exp(-a * mu))
    
    return A / a * (1 - np.exp(-mu * a)) + B / b * (1 - np.exp(-b * (N + 1 - mu)))

def Iuse(N, p):
    if N * p < 1:
        return I((1 / (np.sqrt(1 - 1 / N) - np.log(1 / N)), 1 / (np.sqrt(1 - 1 / N) - np.log(1 - 1 / N))), N, p, 1)
    
    if N * p > N - 1:
        return I((1 / (np.sqrt(1 - 1 / N) - np.log(1 - 1 / N)), 1 / (np.sqrt(1 - 1 / N) - np.log(1 / N))), N, p, N - 1)
    
    else:
        return I((1 / (np.sqrt(N * p * (1 - p)) - np.log(p)), 1 / (np.sqrt(N * p * (1 - p)) - np.log(1 - p))), N, p, N * p)

N = 10
p = 0.9999999999999999
mu = N * p
a = 1 / (np.sqrt(N * p * (1 - p)) - np.log(p))
b = 1 / (np.sqrt(N * p * (1 - p)) - np.log(1 - p))
x = np.arange(0, N, 1)

fValues = f(x, N, p)
PValues = Puse(x, N, p)

fig, ax = plt.subplots()
ax.bar(x, fValues, width = 1, fill = False, edgecolor = "red")
ax.bar(x, PValues, width = 1, fill = False, edgecolor = "green")
print(f"I = {Iuse(N, p)}")

Mask = fValues > PValues
if np.any(Mask.any() and ((fValues - PValues)[Mask] > 1e-8).any()):
    print(f"Error: {x[Mask]}, {(fValues - PValues)[Mask]}")
"""
a_list = np.linspace(0.5 / Std, 2 / Std, 10)
b_list = np.linspace(0.5 / Std, 2 / Std, 10)

a_values = np.empty_like(a_list)
b_values = np.empty_like(b_list)

for i, at in enumerate(a_list):
    a_values[i] = I((at, b), N, p, mu)

for i, bt in enumerate(b_list):
    b_values[i] = I((a, bt), N, p, mu)

fig, ax = plt.subplots()
ax.plot(a_list * Std, a_values)
print(f"a = {Std * a_list[np.argmin(a_values)]}")

fig, ax = plt.subplots()
ax.plot(b_list * Std, b_values)
print(f"b = {Std * b_list[np.argmin(b_values)]}")
"""
"""
p_list = np.linspace(0.001, 0.999, 100)
p_a_values = np.empty_like(p_list)
p_b_values = np.empty_like(p_list)
p_I_values = np.empty_like(p_list)
p_I0_values = np.empty_like(p_list)

for i, pt in enumerate(p_list):
    Std1 = np.sqrt(N * pt * (1 - pt)) - np.log(pt) 
    Std2 = np.sqrt(N * pt * (1 - pt)) - np.log(1 - pt)
    Result = opt.minimize(I, [1 / Std1, 1 / Std2], args = (N, pt, N * pt))#, bounds = ((0.1 / Std1, 10 / Std1), (0.1 / Std2, 10 / Std2)))
    p_a_values[i], p_b_values[i] = Result.x * np.array([Std1, Std2])
    p_I_values[i] = Result.fun
    p_I0_values[i] = I((1 / Std1, 1 / Std2), N, pt, N * pt)
    
fig, ax = plt.subplots()
ax.plot(p_list, p_I_values)
ax.plot(p_list, p_I0_values)

fig, ax = plt.subplots()
ax.plot(p_list, p_a_values)

fig, ax = plt.subplots()
ax.plot(p_list, p_b_values)
"""
