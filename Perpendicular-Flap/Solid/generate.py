import numpy as np
import matplotlib.pyplot as plt
import os

# generate mesh file
os.system("gmsh geom.geo -order 1 -2 -o init.inp")

# open a file and remove lines
input = open("init.inp", "r")

# read eachline and check
lines = input.readlines()[3:]

# check where the good stuff is
for k in range(len(lines)):
    if(len(lines[k])>19 and lines[k].strip()[15:19]==" E N"):
        i = k
    if(len(lines[k])>19 and lines[k].strip()[15:18]=="CPS"):
        j = k

# write to new file
output = open("geom.inp", "w")
output.write("*NODE, NSET=Ninterface\n")
for l in range(i):
    output.write(lines[l])
output.write("*ELEMENT, type=S3, ELSET=PLATE\n")
for l in range(j+1,len(lines)):
    output.write(lines[l])

# clean old file
os.system("rm -f init.inp")
# Done
print("Done")
