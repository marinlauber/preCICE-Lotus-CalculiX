import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def get_rho_s(Mratio, L, h):
    return Mratio*L/h

def get_E(kappa, L, h, nu=0.):
    return kappa*L**3*12*(1-nu**2)/(h**3)

def to_float(ls):
    return [float(ls[i]) for i in range(len(ls))]
def to_int(ls):
    return [int(ls[i]) for i in range(len(ls))]

# get the legnth scale if provided
parser = argparse.ArgumentParser(description='Generates the Finite-element mesh..')
parser.add_argument('-L','--lengthscale', help='Length scale of the body')
parser.add_argument('-Kb','--bending', help='bending cauchy number')
parser.add_argument('-M','--mass', help='mass ratio')
args = parser.parse_args()
if (not args.lengthscale):
    L = 64.
else:
    L = float(args.lengthscale)
if (not args.bending):
    K_B = 0.35
else:
    K_B = float(args.bending)
if (not args.mass):
    m_r = 0.5
else:
    m_r = float(args.mass)

for file in ["geom_1","geom_2"]:
    # generate mesh file
    os.system("gmsh "+file+".geo -order 1 -2 -o init.inp")

    # open a file and remove lines
    input = open("init.inp", "r")

    # read eachline and check
    lines = input.readlines()[3:]
    input.close()

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
    output = open(file+".inp", "w")
    output.write("*NODE, NSET=Ninterface\n")
    for l in range(i):
        # get floats
        ids,x,y,z = to_float(lines[l].strip().replace(",","").split())
        line = "%d, %.8f, %.8f, %.8f\n" % (ids+1,x*L,y*L,z)
        output.write(line)
    output.write("*ELEMENT, type=S3, ELSET=PLATE\n")
    for l in range(j+1,finish):
        el,n1,n2,n3 = to_int(lines[l].strip().replace(",","").split())
        line =  "%d, %d, %d, %d\n" % (el,n1+1,n2+1,n3+1)
        output.write(line)
    output.close()
    # clean old file
    os.system("rm -f init.inp")

# merge the two files
geom_1 = open("geom_1.inp","r")
line_1 = geom_1.readlines()[1:]
geom_2 = open("geom_2.inp","r")
line_2 = geom_2.readlines()[1:]

# check where the good stuff is
finish_1 = len(line_1)
for k in range(len(line_1)):
    # read until we are at the start of the element line
    if(line_1[k][:4]=="*ELE"): i_1 = k
# check where the good stuff is
finish_2 = len(line_2)
for k in range(len(line_2)):
    # read until we are at the start of the element line
    if(line_2[k][:4]=="*ELE"):  i_2 = k

# open new file
output = open("geom.inp","w")
output.write("*NODE, NSET=Ninterface\n")
for l in range(i_1):
    # get floats
    ids,x,y,z = to_float(line_1[l].strip().replace(",","").split())
    line = "%d, %.8f, %.8f, %.8f\n" % (ids,x,y,z)
    output.write(line)
    step = int(ids)
print("End of the geom_1 lines %d" % step)
for l in range(i_2):
    # get floats
    ids,x,y,z = to_float(line_2[l].strip().replace(",","").split())
    line = "%d, %.8f, %.8f, %.8f\n" % (ids+step-1,x,y,z)
    output.write(line)
el0 = ids+step
output.write("*ELEMENT, type=S3, ELSET=PLATE\n")
for l in range(i_1+1,finish_1):
    el,n1,n2,n3 = to_int(line_1[l].strip().replace(",","").split())
    line =  "%d, %d, %d, %d\n" % (el0,n1,n2,n3)
    output.write(line)
    el0 += 1
for l in range(i_2+1,finish_2):
    el,n1,n2,n3 = to_int(line_2[l].strip().replace(",","").split())
    line =  "%d, %d, %d, %d\n" % (el0,n1+step-1,n2+step-1,n3+step-1)
    output.write(line)
    el0 += 1
output.close()

# prepare boundary node file
fixed = open("fixed.nam", "w")
fixed.write("*NSET,NSET=Nfixed\n")
# read each line and extract coordinates
for l in range(i_1):
    ids,x,y,z = to_float(line_1[l].strip().replace(",","").split())
    r = abs(x-0.0*L)
    if r<1e-2: fixed.write("%d,\n" % ids)
for l in range(i_2):
    ids,x,y,z = to_float(line_2[l].strip().replace(",","").split())
    r = abs(x-2.*L)
    if r<1e-2: fixed.write("%d,\n" % (ids+step-1))
fixed.close()