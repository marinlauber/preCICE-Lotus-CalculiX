#
# builder for calculix input file: calculix.inp
#

EOL = "\n"

class Calculix(object):

    def __init__(self, mesh="geom.inp"):
        self.file = open("calculix.inp", "w")
        self.write("*INCLUDE, INPUT="+mesh+"\n")

    def write(self, args):
        self.file.write(args)

    def include(self, fname):
        for name in fname:
            self.write("*INCLUDE, INPUT="+name+"\n")
        self.write(EOL)

    def set_bc(self, fname, bc_1, bc_2, value, func=""):
        txt = "*BOUNDARY"
        if(func != ""): txt += ",AMPLITUDE="+func
        self.write(txt+"\n")
        self.set_constraints(fname, bc_1, bc_2, value)

    def set_constraints(self, fname, bc_1, bc_2, value):
        self.write(" %s, %s, %s, %5.6f" % (fname, bc_1, bc_2, value))
        self.write(EOL)

    def set_amplitude(self, *uamplitude):
        for amplitude in uamplitude:
            self.write("*AMPLITUDE,NAME="+amplitude+",USER\n")

    def set_material(self, E, nu, density, alpha=0.0, beta=0.0):
        self.write("*MATERIAL,NAME=MEMBRANE\n")
        self.write("*ELASTIC\n %.1f, %.2f\n" % (E, nu))
        self.write("*DENSITY\n %.4f\n" % density)
        self.write("*DAMPING, ALPHA=%.2f, BETA=%.2f\n" % (alpha, beta))
        self.write(EOL)

    def set_thickness(self, thickness, nodal=False):
        tmp = "*SHELL SECTION, MATERIAL=MEMBRANE, ELSET=PLATE"
        if(nodal):
            tmp+=",NODAL THICKNESS\n %.4f\n*INCLUDE, INPUT=thickness.nam\n" % thickness
        else:
            tmp += "\n %.4f\n" % thickness
        self.write(tmp)
        self.write(EOL)

    def add_node(self, nodes, nset="Nodes"):
        self.write("*NODE,NSET="+nset+"\n")
        for node in nodes:
            tmp = "%d" % node[0]
            for i in range(1,len(node)): tmp += ", %.4f" % node[i] 
            self.write(tmp)
        self.write(EOL)

    def add_surface(self, name, nodes):
        self.write("*SURFACE,NAME="+name+",TYPE=NODE\n")
        for node in nodes:
            self.write(" %d,\n" % node)
        self.write(EOL)

    def add_orientation(self, A, B, system="CYLINDRICAL", name="OR1"):
        self.write("*ORIENTATION,NAME="+name+",SYSTEM="+system+"\n")
        tmp = ""
        for i in range(3): tmp += "%.4f," % A[i]
        for i in range(3): tmp += "%.4f," % B[i]
        self.write(tmp[:-1])
        self.write(EOL+EOL)
    
    def add_coupling(self, ref, surface, orientation, name, kinematic=[]):
        self.write("*COUPLING,REF NODE="+str(ref)+",SURFACE="+surface)
        self.write(",ORIENTATION="+orientation+",CONSTRAINT NAME="+name+"\n ")
        if len(kinematic)>0:
            for i in range(len(kinematic)): self.write("%d," % kinematic[i])
            self.write(EOL)

    def set_solver(self, dt, T=1e5, Ninc=1e10, alpha=-0.05):
        self.write("*STEP, NLGEOM, INC=%d\n" % Ninc)
        self.write("*DYNAMIC, ALPHA=%-3f\n" % alpha)
        self.write(" %.4f, %.1f\n" % (dt, T))
        self.write("*CONTROLS,PARAMETERS=TIME INCREMENTATION\n ,,,,,,,10,,\n\n")

    def close(self, *args, frequency=1):
        self.write("*CLOAD\n")
        for dof in [1,2,3]:
            self.write(" Ninterface, %d, 0.0\n" % dof)
        self.write("*NODE FILE,OUTPUT=3D,FREQUENCY=%d\n" % frequency)
        out = " U"
        for arg in args:
            out += ", "+arg
        self.write(out+"\n")
        self.write("*END STEP\n")
        self.file.close()
