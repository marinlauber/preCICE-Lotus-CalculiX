#import logging
import ccx2paraview
#logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
c = ccx2paraview.Converter('pressure.frd', ['vtu'])
c.run()