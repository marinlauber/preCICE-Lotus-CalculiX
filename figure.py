#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import vtk
import matplotlib.pyplot as plt
from scipy.sparse import diags
import lotus

plt.style.use("Journal")

def vorticity(u,v,x,y):
    omega = np.zeros_like(u)
    omega[1:-1,1:-1] = 0.25*((v[2:,1:-1]+v[2:,2:]-v[:-2,1:-1]-v[:-2,2:])-
                             (u[1:-1,2:]+u[2:,2:]-u[1:-1,:-2]-u[2:,:-2]))
    return omega

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

def cm(x): return x/2.56

# data from simulation
rho = 1.; U =1.
L=64; b=6.4; d=3.2
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

# plot results
plt.figure(figsize=(cm(6.),cm(6.)))

for i,E in enumerate([40000,20000]):

    tri,_,_,_  = lotus.read_vtu("interF."+str(E)+".vtu")
    dat = tri[np.isclose(tri[:,2],0.0,atol=1e-2)]
    plt.plot(dat[:,0]/L, dat[:,1]/L, '-C0', alpha=0.5+(0.5*i), lw=2,
            label=rf"$Ca=%.1f$" % (1.5*(i+1)))

for i,E in enumerate([40000,20000]):
    # solid
    I = b*d**3/12
    # Ca = rho*CD*b*U**2*L**3 / (2*E*I)
    Ca = rho*U**2*L**4 / (16*E*I)
    print('Cauchy Number (Gosselin) : %.4f' % Ca)

    resid = 1.; k=0
    while resid>1e-5:

        # cpmpute fluid loading
        integral = Ca*FluidForce(theta, s, ds)

        # apply boundary conditions, as explicit solution
        integral[0] = 0.          # Homogeneous Dirichlet at root
        integral[-1] = -theta[-2] # Homogeneous Neumann at tip

        # solve linear system
        res = np.linalg.solve(ddx, integral)

        # prepare next iteration
        resid = Linf(theta+res)
        theta = -res
        k+=1
    if i==0:
        plt.plot(*theta_to_xy(theta, s), 'xk', lw=0.5, ms=3, label="Analytical")
    else:
        plt.plot(*theta_to_xy(theta, s), 'xk', lw=0.5, ms=3)


print("Number of iteration for convergence : ", k)
print("Effective blade length              : %.2f" % effective_length(theta, s))
print("Effective blade height              : %.2f" % effective_height(theta, s))
plt.legend(title="Numerical", loc=4)
plt.axis([-0.1, 0.8, 0, 1.1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("deflection.png")
plt.close()
