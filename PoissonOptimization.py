import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as opt
plt.close("all")

def f(n, lam):
    return st.poisson.pmf(n, lam)

def P(n, lam, mu, a, b):
    n_m = np.ceil(lam * np.exp(-a)) - 1
    n_p = np.ceil(lam * np.exp(b)) - 1
    A = a * lam ** n_m * np.exp(-(lam + a * (n_m - mu))) / (sp.gamma(n_m + 1) * (np.exp(a) - 1))
    B = b * lam ** n_p * np.exp(-(lam - b * (n_p - mu))) / (sp.gamma(n_p + 1) * (1 - np.exp(-b)))

    if lam <= 1:
        C = np.exp(1) / 2
        return C * np.exp(- n)
        
    Result = np.empty_like(n, dtype = float)
    F_m = n <= np.floor(mu) - 1
    F_p = n >= np.ceil(mu)
    F_l = (~F_m & ~F_p)
    Result[F_m] = A / a * np.exp(a * (n[F_m] - mu)) * (np.exp(a) - 1)
    Result[F_p] = B / b * np.exp(-b * (n[F_p] - mu)) * (1 - np.exp(-b))
    Result[F_l] = A / a * (1 - np.exp(a * (n[F_l] - mu))) + B / b * (1 - np.exp(-b * (n[F_l] + 1 - mu)))
    return Result

def I(Params, lam):
    mu, a, b = Params
    n_m = np.ceil(lam * np.exp(-a)) - 1
    n_p = np.ceil(lam * np.exp(b)) - 1
    A = a * lam ** n_m * np.exp(-(lam + a * (n_m - mu))) / (sp.gamma(n_m + 1) * (np.exp(a) - 1))
    B = b * lam ** n_p * np.exp(-(lam - b * (n_p - mu))) / (sp.gamma(n_p + 1) * (1 - np.exp(-b)))
    return A / a * (1 - np.exp(-mu * a)) + B / b

lam = 0.01
mu = lam
a = 1 / np.sqrt(lam)
b = 1 / np.sqrt(lam)
x = np.arange(np.ceil(3 * lam) + 2)

fig, ax = plt.subplots()
ax.bar(x, f(x, lam), width = 1, fill = False, edgecolor = "red")
ax.bar(x, P(x, lam, mu, a, b), width = 1, fill = False, edgecolor = "green")
"""
a_list = np.linspace(0.5 / np.sqrt(lam), 2 / np.sqrt(lam), 1000)
b_list = np.linspace(0.5 / np.sqrt(lam), 2 / np.sqrt(lam), 1000)
mu_list = np.linspace(0, 2 * lam, 1000)

a_values = np.empty_like(a_list)
b_values = np.empty_like(b_list)
mu_values = np.empty_like(mu_list)

for i, at in enumerate(a_list):
    a_values[i] = I((mu, at, b), lam)

for i, bt in enumerate(b_list):
    b_values[i] = I((mu, a, bt), lam)

for i, mut in enumerate(mu_list):
    mu_values[i] = I((mut, a, b), lam)

fig, ax = plt.subplots()
ax.plot(a_list, a_values)
print(f"a = {np.sqrt(lam) * a_list[np.argmin(a_values)]}")

fig, ax = plt.subplots()
ax.plot(b_list, b_values)
print(f"b = {np.sqrt(lam) * b_list[np.argmin(b_values)]}")

fig, ax = plt.subplots()
ax.plot(mu_list, mu_values)
print(f"mu = {mu_list[np.argmin(mu_values)] / lam}")
""""""
lam_list = np.linspace(1, 100, 100)
lam_mu_values = np.empty_like(lam_list)
lam_a_values = np.empty_like(lam_list)
lam_b_values = np.empty_like(lam_list)
lam_I_values = np.empty_like(lam_list)
lam_I0_values = np.empty_like(lam_list)

for i, lamt in enumerate(lam_list):
    Result = opt.minimize(I, [lamt, 1 / np.sqrt(lamt), 1 / np.sqrt(lamt)], args = (lamt,), bounds = ((0.5 * lamt, 2 * lamt), (0.5 / np.sqrt(lamt), 2 / np.sqrt(lamt)), (0.5 / np.sqrt(lamt), 2 / np.sqrt(lamt))))
    lam_mu_values[i], lam_a_values[i], lam_b_values[i] = Result.x * np.array([1 / lamt, np.sqrt(lamt), np.sqrt(lamt)])
    lam_I_values[i] = Result.fun
    lam_I0_values[i] = I((lamt, 1 / np.sqrt(lamt), 1 / np.sqrt(lamt)), lamt)
    
fig, ax = plt.subplots()
ax.plot(lam_list, lam_I_values)
ax.plot(lam_list, lam_I0_values)

fig, ax = plt.subplots()
ax.plot(lam_list, lam_a_values)

fig, ax = plt.subplots()
ax.plot(lam_list, lam_b_values)

fig, ax = plt.subplots()
ax.plot(lam_list, lam_mu_values)
"""