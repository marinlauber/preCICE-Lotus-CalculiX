#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import vtk
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.sparse import diags

plt.style.use("Journal")

def read_file(fname):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()
    data = reader.GetOutput()
    vtkData = data.GetPoints().GetData()
    Npoints = data.GetPointData().GetNumberOfTuples()

    p = np.array(data.GetPointData().GetArray("Pressure")).reshape((Npoints))
    v = np.array(data.GetPointData().GetArray("Velocity")).reshape((Npoints,3))
    coords = np.array([vtkData.GetTuple3(x) for x in range(Npoints)])
    connectiv = []

    for i in range(data.GetNumberOfCells()):
        connectiv.append([data.GetCell(i).GetPointIds().GetId(j)
            for j in range(data.GetCell(i).GetPointIds().GetNumberOfIds())])

    connectiv = np.array(connectiv, dtype=np.float64)

    return coords, connectiv, p, np.transpose(v)

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
L=1.; b=L/64.; d=L/20
CD = 1.

# solid
E = 20000; I = b*d**3/12
Ca =rho*CD*b*U**2*L**3 / (2*E*I)
print('Cauchy Number (Gosselin) : %.4f' % Ca)

# numerical parameters
N = 2**5
s, ds = np.linspace(0, 1, N, retstep=True)
theta = np.zeros_like(s)

# build curvature matrix by hand (secnod order with BC)
ddx = diags([1, -2, 1], [-1, 0, 1], shape=(N,N)).toarray()/ds**2

# boundary condititons as explicit solution to the system
ddx[ 0,:]=0.; ddx[ 0, 0]=1.
ddx[-1,:]=0.; ddx[-1,-1]=1.

NN = 100
tri, connectivity, p, v  = read_file("Fluid/dat/interF."+str(NN)+".vtu")

# select on edge only at z=6.4
dat = tri[np.isclose(tri[:,2],1.0,atol=1e-2)]

plt.figure(figsize=(cm(8),cm(8)))
plt.plot(dat[:,0]/64.,dat[:,1]/64., lw=0.75, label="Lotus-CalculiX")

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

plt.plot(*theta_to_xy(theta, s), 'xk', lw=0.1, ms=2.5, label="Analytical")
print("Number of iteration for convergence : ", k)
print("Effective blade length              : %.2f" % effective_length(theta, s))
print("Effective blade height              : %.2f" % effective_height(theta, s))
plt.legend()
plt.axis([-0.2, 1, 0, 1.1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("deflection.png")
