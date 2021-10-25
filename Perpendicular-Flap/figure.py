#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
import lotus

plt.style.use("Journal")

def xy_to_theta(x, y):
    theta = np.zeros_like(x)
    s = np.zeros_like(x)
    for i in range(1,len(x)):
        print((x[i]-x[i-1]))
        print((y[i]-y[i-1]))
        theta[i] = np.arctan((x[i]-x[i-1])/(y[i]-y[i-1]))
        s[i] = np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
    return theta, s

def theta_to_xy(theta, s):
    x = np.zeros_like(theta)
    y = np.zeros_like(theta)
    ds = s[1]-s[0]
    for i in range(1,len(theta)):
        x[i] = x[i-1]+np.sin(theta[i-1])*ds
        y[i] = y[i-1]+np.cos(theta[i-1])*ds
    return x, y

def effective_length(theta, s):
    ds = s[1]-s[0]
    return np.trapz(np.cos(theta)**3, x=s, dx=ds)

def effective_height(theta, s):
    ds = s[1]-s[0]
    return np.trapz(np.cos(theta), x=s, dx=ds)

def int_num(theta, theta_star, s_star, ds):
    func = np.cos(theta - theta_star) * np.cos(theta)**2
    return np.trapz(y=func, x=s_star, dx=ds)

def FluidForce(theta, s, ds):
    intgrl = np.empty_like(theta)
    for i in range(len(intgrl)):
        intgrl[i] = int_num(theta[i:], theta[i], s[i:], ds)
    return intgrl

def Linf(u): return max(abs(u))

def kappa(E, h, nu, L, rho_f=1, U=1):
    D = E*h**3/(12*(1-nu**2))
    return (rho_f*U**2*L**3)/D

def cm(x): return x/2.56

# data from simulation
rho = 1.; U =1.
L=64; b=1.; d=2.56
CD = 1.95

# numerical parameters
N = 2**5
s, ds = np.linspace(0, 1, N, retstep=True)
theta = np.zeros_like(s)

# build curvature matrix by hand (secnod order with BC)
ddx = diags([1, -2, 1], [-1, 0, 1], shape=(N,N)).toarray()/ds**2

# boundary condititons as explicit solution to the system
ddx[ 0,:]=0.; ddx[ 0, 0]=1.
ddx[-1,:]=0.; ddx[-1,-1]=1.

# stuff for plots
YoungsModulus = [2500.,5000.]

# plot results
fig = plt.figure(figsize=(cm(12.),cm(6.)))
gs = fig.add_gridspec(nrows=1, ncols=2)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
for i,E in enumerate(YoungsModulus):
    tip = []
    for i in range(1,250):
        tri,_,_,_  = lotus.read_vtu("res_"+str(int(E))+"/Fluid/dat/interF."+str(i)+".vtu")
        dat = tri[np.isclose(tri[:,2],0.5,atol=1e-2)]
        tip.append(np.sqrt(max(dat[:,0]/L)**2+max(dat[:,1]/L+0.5)**2))
    ax2.plot(tip)
    Cauchy = rho*U**2*L**3 / (E*d**3/12)
    ax1.plot(dat[:,0]/L, dat[:,1]/L+0.5, "-C0", alpha=0.5, lw=2,
             label=rf"$Ca=%.f$" % Cauchy)

tip = []
for i in range(1,250):
    tri,_,_,_  = lotus.read_vtu("Fluid/dat/interF."+str(i)+".vtu")
    dat = tri[np.isclose(tri[:,2],0.5,atol=1e-2)]
    tip.append(np.sqrt(max(dat[:,0]/L)**2+max(dat[:,1]/L+0.5)**2))
ax2.plot(tip)
Cauchy = rho*U**2*L**3 / (E*d**3/12)
ax1.plot(dat[:,0]/L, dat[:,1]/L+0.5, "-C0", alpha=0.5, lw=2,
         label=rf"$Ca=%.f$" % 1250)

for i,Ca in enumerate([8,16]):
   
    print("Cauchy Number (Gosselin) : %.4f" % Ca)
    E = 0.5*CD*b*U**2*L**3/(Ca*b*d**3/12)
    print("Young\'s modulus: (E) : %.2f" % E)

    resid = 1.; k=0
    while resid>1e-5:

        # compute fluid loading
        integral = Ca*FluidForce(theta, s, ds)

        # apply boundary conditions, as explicit solution
        integral[0] = 0.            # Homogeneous Dirichlet at root
        integral[-1] = -theta[-2]   # Homogeneous Neumann at tip

        # solve linear system
        res = np.linalg.solve(ddx, integral)

        # prepare next iteration
        resid = Linf(theta+res)
        theta = -res
        k+=1

    # plot once we converged
    if i==0:
        ax1.plot(*theta_to_xy(theta, s), "xk", lw=0.5, ms=3, label="Analytical")
    else:
        ax1.plot(*theta_to_xy(theta, s), "xk", lw=0.5, ms=3)

    print("Number of iteration for convergence : ", k)
    print("Effective blade length              : %.2f" % effective_length(theta, s))
    print("Effective blade height              : %.2f" % effective_height(theta, s))

ax1.legend(title="Numerical", loc=4)
ax1.set_xlim(-0.1, 1.); ax1.set_ylim(0.,1.1)
ax1.set_aspect("equal", adjustable="box")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax2.set_xlabel("time")
ax2.set_ylabel("deflection")
plt.tight_layout()
plt.savefig("deflection.png")
plt.close()
