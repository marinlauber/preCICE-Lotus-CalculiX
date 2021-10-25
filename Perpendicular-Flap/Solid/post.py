import ccx2paraview
c = ccx2paraview.Converter('calculix.frd', ['vtu'])
c.run()