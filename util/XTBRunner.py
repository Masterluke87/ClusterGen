from ase.units import Hartree
from ase.io import read,write
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
import os,stat
import numpy as np

class XTB(FileIOCalculator):

    implemented_properties = ['energy']
    ignored_changes= {"positions"}
    command = "./runxtb.sh"

    default_parameters = dict(
            charge = 0, mult = 1, cycles=1000
            )


    def __init__(self, restart=None, label ="tmp/XTB", atoms=None, **kwargs):
        FileIOCalculator.__init__(self, restart=restart,label=label, atoms=atoms, **kwargs)
        FileIOCalculator.set_label(self, label)  

    def write_input(self,atoms,properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms)
        write(self.directory+"/xtbinput.xyz",atoms)
        f = open(self.directory+"/runxtb.sh","w")
        f.write("#!/bin/bash \n")
        f.write("~/Downloads/xtb-6.4.1/bin/xtb xtbinput.xyz --gfn 1 --uhf {} --chrg {} --opt normal --cycles {} > XTBOUT".format(self.parameters["mult"]-1,self.parameters["charge"], self.parameters["cycles"]))
        f.close()

        st = os.stat(self.directory+"/runxtb.sh")
        os.chmod(self.directory+"/runxtb.sh", st.st_mode | stat.S_IEXEC )

    def read_results(self):
        newatoms = read(self.directory+"/xtbopt.xyz")
        lines = open(self.directory+"/XTBOUT").readlines() 

        converged = [x for x in lines if "GEOMETRY OPTIMIZATION CONVERGED" in x]
        if len(converged)==0:
            raise(Exception("Did not converge"))

      
        energies = [float(x.split()[3]) for x in lines if "TOTAL ENERGY" in x]
        self.results["energy"] =  (energies[-1] * Hartree)

    def calculate(self,atoms=None, properties=['energy'], system_changes=['']):
        FileIOCalculator.calculate(self,atoms,properties,system_changes)
        newatoms = read(self.directory+"/xtbopt.xyz")
        atoms.positions = newatoms.positions.copy()
