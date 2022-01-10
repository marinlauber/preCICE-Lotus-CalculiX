import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def get_rho_s(Mratio, L, h):
    return Mratio*L/h

def get_E(kappa, L, h, nu=0.):
    return kappa*L**3*12*(1-nu**2)/(h**3)

def to_float(ls):
    "converts the node id and its coordinates from string to floats"
    return int(ls[0]),float(ls[1]),float(ls[2]),float(ls[3])
def to_int(ls):
    return int(ls[0]),int(ls[1]),int(ls[2]),int(ls[3])

# get the legnth scale if provided
parser = argparse.ArgumentParser(description='Generates the Finite-lement mesh..')
parser.add_argument('-L','--lengthscale', help='Length scale of the body')
parser.add_argument('-Kb','--bending', help='bending cauchy number')
parser.add_argument('-M','--mass', help='mass ratio')
args = parser.parse_args()
if (not args.lengthscale):
    L = 64.
else:
    L = float(args.lengthscale)
if (not args.bending):
    K_B = 0.15
else:
    K_B = float(args.bending)
if (not args.mass):
    m_r = 0.5
else:
    m_r = float(args.mass)

# generate mesh file
os.system("gmsh geom.geo -order 1 -2 -o init.inp")

# open a file and remove lines
input = open("init.inp", "r")

# read eachline and check
lines = input.readlines()[3:]

# check where the good stuff is
finish = len(lines)
for k in range(len(lines)):
    # read until we are at the start of the element line
    if(len(lines[k])>19 and lines[k].strip()[15:19]==" E N"):  i = k
    # then read until we arrive at the shell element
    if(len(lines[k])>19 and lines[k].strip()[15:18]=="CPS"): j = k
    # read until next definition
    if(lines[k][:7]=="*ELSET,"): finish=k; break

# write to new file
output = open("geom.inp", "w")
output.write("*NODE, NSET=Ninterface\n")
for l in range(i):
    # get floats
    ids,x,y,z = to_float(lines[l].strip().replace(",","").split())
    line = "%d, %.8f, %.8f, %.8f \n" % (ids+1,x*L,y*L,z)
    output.write(line)
output.write("*ELEMENT, type=S3, ELSET=PLATE\n")
for l in range(j+1,finish):
    el,n1,n2,n3 = to_int(lines[l].strip().replace(",","").split())
    line =  "%d, %d, %d, %d \n" % (el,n1+1,n2+1,n3+1)
    output.write(line)

# prepare boundary node file
fixed = open("fixed.nam", "w")
fixed.write("*NSET,NSET=Nfixed\n")
# read each line and extract coordinates
for l in range(i):
    ids,x,y,z = to_float(lines[l].strip().replace(",","").split())
    r = abs(x-0.008)
    if r<1e-2:
        fixed.write(str(ids+1)+",\n")

# get the parameters
thickness=2.56; nu=0.0
# effective stiffness is less
E = get_E(K_B, L/2, thickness)
density = get_rho_s(m_r, L, thickness)

# generate calcilux.inp file
calculix = open("calculix.inp", "w")
calculix.write("*INCLUDE, INPUT=geom.inp\n*INCLUDE, INPUT=fixed.nam\n")
calculix.write("*BOUNDARY\n Nfixed, 1, 3, 0.0\n Ninterface, 3, 4, 0.0\n")
calculix.write("*MATERIAL,NAME=MEMBRANE\n")
calculix.write("*ELASTIC\n %.1f, %.2f\n" % (E, nu))
calculix.write("*DENSITY\n %.2f\n" % density)

calculix.write("*DAMPING, ALPHA=0.0, BETA=0.0\n")
calculix.write("*SHELL SECTION, MATERIAL=MEMBRANE, ELSET=PLATE\n %.4f\n" % thickness)

calculix.write("*STEP, NLGEOM, INC=10000000000\n*DYNAMIC, ALPHA=0.0\n")
calculix.write(" 0.25, 100000\n")

calculix.write("*CONTROLS,PARAMETERS=TIME INCREMENTATION\n ,,,,,,,10,,\n")

calculix.write("*CLOAD\n Ninterface, 1, 0.0\n Ninterface, 2, 0.0\n Ninterface, 3, 0.0\n")

calculix.write("*NODE FILE,OUTPUT=3D\n")
calculix.write(" U\n")

calculix.write("*END STEP\n")

# clean old file
os.system("rm -f init.inp")
# Done
print("Done")
