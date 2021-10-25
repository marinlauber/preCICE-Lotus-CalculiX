#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
import os
import shutil
import vtk
import numpy as np
import matplotlib.pyplot as plt

def run(n_proc=0,run_folder='test',read_folder=None):
    "setup and run Lotus using lotus.f90 and the files in postproc"

    print('Number of proccessors :{}'.format(n_proc))
    print('Run folder            :{}'.format(run_folder))

    if read_folder:
        print('Read folder           :{}'.format(read_folder))
    else:
        print('No read folder')

    if os.path.exists(run_folder):
        print('Folder '+run_folder+' exists!')
        if read_folder=='./':
            print('Resuming in place')
        else:
            print('Moving contents to trash')
            subprocess.call('trash '+run_folder+'/*', shell=True)
    else:
        print('Creating '+run_folder)
        os.makedirs(run_folder)

    print('Setting up in '+run_folder)
    if os.path.isdir('postproc'):
        for file_name in os.listdir('postproc'):
            full_file_name = os.path.join('postproc', file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, run_folder)
    if os.path.isfile('lotus.f90'):
        shutil.copy('lotus.f90', run_folder)

    print('Making executable ')
    os.chdir(run_folder)
    subprocess.call('make -C $MGLHOME/geom_lib/ libgeom.a', shell=True)
    subprocess.call('make -C $MGLHOME/oop/ libfluid.a', shell=True)
    subprocess.call('make -f $MGLHOME/oop/Makefile lotus', shell=True)
    print('Finished executable ')

    print('Running executable ')
    if read_folder:
        if n_proc==0:
            subprocess.call('time ./lotus '+read_folder, shell=True)
        else:
            subprocess.call('mpirun -n '+str(n_proc)+' ./lotus '+read_folder, shell=True)
    else:
        if n_proc==0:
            subprocess.call('time ./lotus', shell=True)
        else:
            subprocess.call('mpirun -n '+str(n_proc)+' ./lotus', shell=True)

    print('Run all python files for postprocessing')
    for file_name in os.listdir('.'):
        if file_name.endswith('.py'):
            print(file_name)
            subprocess.call('python3 '+file_name, shell=True)

    print('Popping back up')
    os.chdir('../.')


def replace(template,dic):
    """
    Write lotus.f90 file by replacing the dic on the template.f90
    and return a potential folder name
    """
    f1 = open(template,'r')
    f2 = open('lotus.f90','w')
    for text in f1:
        for i, j in dic.items():
            text = text.replace(i, j)
        f2.write(text)
    f1.close()
    f2.close()
    return '_'.join([i+'_'+j for i,j in dic.items()])


def show_grid(fname, save, L=1):
    (u,v,w), (p), (x,y,z) = read_vtr(fname)
    x=(x-0.5)/L; y=(y-0.5)/L; z=(z-0.5)/L
    if len(z)==1:
        plot_2D_Grid(x, y, p[:,:,0], L, every=2)
    else:
        plot_3D_Gird(x, y, z, p, L, every=2)
    plt.tight_layout(pad=1.08, h_pad=1.2)
    if save: plt.savefig("grid.png", dpi=600)
    plt.show()


def plot_XY_grid(ax, x, y, dist, L, every, xlab=r"x/L", ylab=r"y/L"):
    ax.vlines(x[::every], ymin=min(y), ymax=max(y), color='gray', lw=0.2)
    ax.hlines(y[::every], xmin=min(x), xmax=max(x), color='gray', lw=0.2)
    ax.contour(x, y, dist.T, colors='k', levels=[0.])
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax.set_xlim(min(x-50/L),max(x+50/L))
    ax.set_ylim(min(y-50/L),max(y+50/L))
    ax.set_aspect('equal', adjustable='box')
    return ax


def plot_table(ax, *args):
    tot = []; cell_text = []
    for i,arg in enumerate(args):
        dx = arg[1:]-arg[:-1]
        mn, mx = min(abs(dx)), max(abs(dx))
        cell_text.append([str(len(dx)),"%.4f" % mn,"%.4f" % mx, "%.1f" % (mx/mn)] )
        tot.append(len(dx))
    cell_text.append( [str(np.product(np.array(tot))),"n-a","n-a","n-a"] )
    r=["X","Y","Z"]; rowlab = [r[j]+"-axis" for j in range(i+1)]; rowlab.append(r"$\sum$")
    the_table = ax.table(cellText=cell_text,
                         rowLabels=rowlab,
                         colLabels=[r"$N_i$",r"$\Delta x_i$ min",r"$\Delta x_i$ max",r"$r$ max"],
                         loc='center')
    the_table.scale(1, 1.5)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.axis("off")
    return ax


def plot_2D_Grid(x, y, dist, L, every):
    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot2grid((3,3), (0,0), colspan=3,rowspan=2, fig=fig)
    t  = plt.subplot2grid((3,3), (2,1), fig=fig)
    ax = plot_XY_grid(ax, x, y, dist, L, every)
    t  = plot_table(t, x, y)


def plot_3D_Gird(x, y, z, dist, L, every):
    """
    1st angle projection
    """
    fig = plt.figure(figsize=(10,6))
    xz = plt.subplot2grid((2,3), (0,0), colspan=2)
    xy = plt.subplot2grid((2,3), (1,0), colspan=2)
    yz = plt.subplot2grid((2,3), (0,2))
    t  = plt.subplot2grid((2,3), (1,2)) 
    xy = plot_XY_grid(xy, x, y, dist[:,:,np.argmin(abs(z))], L, every)
    yz = plot_XY_grid(yz, -y[::-1], z, dist[np.argmin(abs(x)),::-1,::], L, every, xlab=r"y/L", ylab=r"z/L")
    xz = plot_XY_grid(xz, x, z, dist[:,np.argmin(abs(y)),:], L, every, ylab=r"z/L")
    t  = plot_table(t, x, y, z)


def read_vti(fname):
	reader = vtk.vtkXMLPImageDataReader()
	reader.SetFileName(fname)
	reader.Update()
	data = reader.GetOutput()
	pointData = data.GetPointData()

	sh = data.GetDimensions()[::-1]
	ndims = len(sh)
    
	# get vector field
	v = np.array(pointData.GetVectors("Velocity")).reshape(sh + (ndims,))
	vec = []
	for d in range(ndims):
		a = v[..., d]
		vec.append(a)
	# get scalar field
	sca = np.array(pointData.GetScalars('Pressure')).reshape(sh + (1,))

	# Generate grid
	# nPoints = data.GetNumberOfPoints()
	(xmin, xmax, ymin, ymax, zmin, zmax) = data.GetBounds()
	grid3D = np.mgrid[xmin:xmax + 1, ymin:ymax + 1, zmin:zmax + 1]

	return np.transpose(np.array(vec), (0,3,2,1)), np.transpose(sca, (3,2,1,0))[0,:,:,:], grid3D


def read_vtr(fname):
	reader = vtk.vtkXMLPRectilinearGridReader()
	reader.SetFileName(fname)
	reader.Update()
	data = reader.GetOutput()
	pointData = data.GetPointData()

	sh = data.GetDimensions()[::-1]
	ndims = len(sh)

	# get vector field
	v = np.array(pointData.GetVectors("Velocity")).reshape(sh + (ndims,))
	vec = []
	for d in range(ndims):
		a = v[..., d]
		vec.append(a)
	vec = np.array(vec)

	# get scalar field
	sca = np.array(pointData.GetScalars('Pressure')).reshape(sh + (1,))

	# get grid
	x = np.array(data.GetXCoordinates())
	y = np.array(data.GetYCoordinates())
	z = np.array(data.GetZCoordinates())

	return np.transpose(vec, (0,3,2,1)), np.transpose(sca, (3,2,1,0))[0,:,:,:], np.array([x, y, z], dtype=object)


def read_vtu(fname):
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