from calculix import Calculix

# generate calculix.inp file
calculix = Calculix()
calculix.set_bc("fixed.nam", 2, 2, 0.0)
calculix.set_material(10000, 0.3, 10.0)
calculix.set_thickness(2.56)
calculix.set_amplitude("UFUNC_10","UFUNC_11","UFUNC_12")
calculix.set_bc("Nfixed", 1, 2, value=1.25, func="UFUNC_10")
calculix.set_bc("Nfixed", 3, 3, value=1.25, func="UFUNC_11")
calculix.set_bc("Nfixed", 2, 2, value=1.0, func="UFUNC_12")
calculix.set_bc("Ninterface", 3, 5, value=0.0)
calculix.set_solver(0.25)
calculix.close()


calculix = Calculix()
calculix.set_bc("fixed.nam", 1, 6, 0.0)
calculix.set_bc("Ninterface", 3, 4, 0.0)
calculix.set_material(10000, 0.00, 16.25)
calculix.add_node([[1001,0.,0.,0.]], nset="Nref")
calculix.set_bc("Nref", 1, 5, 0.0)
calculix.add_surface("S1", [130,131,134,257])
calculix.add_orientation([0.,0.,0.,],[0.,0.,1.])
calculix.add_coupling(1001, "S1", "OR1", "CN1", kinematic=[1,2])
calculix.set_thickness(2.56)
calculix.set_solver(0.25)
calculix.close()
