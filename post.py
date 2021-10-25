import numpy as np
import matplotlib.pyplot as plt

cm = lambda x : x/2.56

conv = np.genfromtxt("Solid/precice-Calculix-convergence.log", skip_header=1)
iter = np.genfromtxt("Solid/precice-Calculix-iterations.log", skip_header=1)

x = np.array([])
resDispl = np.array([])
resForces = np.array([])

for it in range(1,iter.shape[0]+1):
    array = conv[conv[:,0]==it,2:]
    x = np.hstack((x,it+np.linspace(0,1,array.shape[0])))
    resDispl = np.hstack((resDispl,array[:,0]))
    resForces = np.hstack((resForces,array[:,1]))

fig1 = plt.figure("FSI iterations",figsize=(cm(12),cm(10)), constrained_layout=True)
gs = fig1.add_gridspec(nrows=1, ncols=1)
ax = fig1.add_subplot(gs[0,0])
ax.semilogy(x, resDispl, label="Displacements")
ax.semilogy(x, resForces, label="Forces")
ax.legend()
ax.set_ylabel("Relative Residuals")
ax.set_xlabel("Time-steps")
ax.set_xlim(0,iter.shape[0])
ax2 = ax.twinx()
ax2.hist(conv[:,0]-1,bins=iter.shape[0]+1,histtype='stepfilled',
         align="mid",color="gray",alpha=0.3)
ax2.set_ylabel("Iterations")
fig1.savefig("iterations.png", dpi=600)

mg = np.genfromtxt("Fluid/fort.8", skip_header=1)

fig2 = plt.figure("Fluid Solver",figsize=(cm(12),cm(10)), constrained_layout=True)
gs = fig2.add_gridspec(nrows=1, ncols=1)
ax = fig2.add_subplot(gs[0,0])
x = np.linspace(0,iter.shape[0],mg.shape[0])
ax.semilogy(x, mg[:,1], label="res0", lw=0.5, alpha=0.8)
ax.semilogy(x, mg[:,2], label="res", lw=0.5, alpha=0.8)
ax.semilogy(x, mg[:,3], label="inf", lw=0.5, alpha=0.8)
ax.legend()
ax.set_xlim(0,iter.shape[0])
ax.set_ylabel("Pressure Residuals")
ax.set_xlabel("Time-steps")
ax2 = ax.twinx()
ax2.plot(x, mg[:,0],"-k",label="iter",lw=0.5,alpha=0.3)
ax2.set_ylabel("MG Iterations")
fig2.savefig("fluid.png", dpi=600)

plt.show()