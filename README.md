## Coupled FSI Lotus-CalculiX simulation using the preCICE library

### Pre-processing

First, the structural model needs to be build. This is done by changing the `geom.geo` and the `calculix.inp` files. Once you have set the correct geometric dimensions, you can generate the mesh using 
```bash
cd Solid
python generate.py
```
which produces the mesh file (`geom.inp`) required in the `calculix.inp` file.

#### Boundary conditions for the structural sub-problem

You then have to specify the boundary nodes in the `fixed.nam` file. The corresponding boundary conditions are set in the `calculix.inp` file. Here we have a clamped boundary condition on one edge of the plate.

This can be done through the `cgx` pre-processor.
```bash
cgx -c geom.inp
```
which opens the `cgx` editor. You can then plot all the nodes with their number with `plot na all`. To generate a `fixed.nam` file you have two options, either add by hand the node number that are part of the boundary, or select them by drawing a box around it. This is done by first creating a node list
```bash
qadd fixed nam
```
then pressing `a` to enter the add mode. To draw a rectangle from which the nodes will be press `r` on each of the opposite corners. Move the rectangle such that all boundary nodes are inside and the press `n+a`. This adds all the nodes to the fixed set.
You can then exit the selection with `q`. To write the `fixed.nam` file type in
```bash
send fixed abqs nam
```
This should write the corresponding file.

### Running the Simulation

To run the simulations in parallel, default is np=1 (serial). 

```bash
./Allrun -parallel np
```

To clean all the files generated, simply clean the repo using

```bash
./Allclean
```

### Post-processing

The fluid is gradually accelerated from zero to one using an hyperbolic profile. The resulting fluid field once the structural motion has settled is shown below.

![Result 1](fluid_render.png)

An analytical expression for the blade deflection has been derived in [Luhar and Nepf](https://doi.org/10.4319/lo.2011.56.6.2003)


<a href="https://www.codecogs.com/eqnedit.php?latex=-\frac{d^2\theta}{d\hat{s}^2}\biggr\rvert_{\hat{s}^*}&space;=&space;Ca&space;\int_{\hat{s}^*}^{1}\cos(\theta-\theta^*)\cos^2\theta&space;d\hat{s}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?-\frac{d^2\theta}{d\hat{s}^2}\biggr\rvert_{\hat{s}^*}&space;=&space;Ca&space;\int_{\hat{s}^*}^{1}\cos(\theta-\theta^*)\cos^2\theta&space;d\hat{s}" title="-\frac{d^2\theta}{d\hat{s}^2}\biggr\rvert_{\hat{s}^*} = Ca \int_{\hat{s}^*}^{1}\cos(\theta-\theta^*)\cos^2\theta d\hat{s}" /></a>

with the Cauchy number $Ca$ defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=Ca&space;=&space;\frac{1}{2}\frac{\rho&space;C_dbU^2l^3}{EI}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Ca&space;=&space;\frac{1}{2}\frac{\rho&space;C_dbU^2l^3}{EI}" title="Ca = \frac{1}{2}\frac{\rho C_dbU^2l^3}{EI}" /></a>


The results and the analytical solution are shown below

<img src="deflection.png" alt="drawing" width="600"/>