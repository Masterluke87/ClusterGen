import argparse
from ase.io import read
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
import numpy as np
from ase.visualize import view
from ase.atoms import Atoms
from ase.ga.data import DataConnection,PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.io import read,write
from ase.atoms import Atoms
import os,sys
from util.XTBRunner import XTB
import pdb

logG = lambda x: print("\033[92m"+x+"\033[00m")
logR = lambda x: print("\033[91m"+x+"\033[00m")


def runTheParser():
    parser = argparse.ArgumentParser(description="""
    ClusterGen. A tool collection indended to generate low energy cluster. 
    The protocol combines the concepts of genetic algorithm (implemented in ase) with the XTB method (see, 
    https://github.com/grimme-lab/xtb).
    """)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-c','--charge', metavar='Charge',default=0, type=int,
                help='Defines the charge of your cluster')       
    optional.add_argument('-m','--mult', metavar='Mult',default=1, type=int,
                help='Set the multiplicty of the target cluster')
    optional.add_argument('-p','--population', metavar='nPOP',default=20, type=int,
                help='Size of the Population')
    optional.add_argument('-db','--database', metavar='nPOP',default="cluster.db", type=str,
                help='Sqlite database to handle the structures')          
    optional.add_argument('-l','--length', metavar=' [Angstrom] ',default="12", type=str,
                help='Size of the cubic simulation cell in Angstrom')          
                   
    optional.add_argument('-F','--formula', metavar='Formula',default="", type=str, required=True,
                help='Set empirical formula for the target cluster, eg. Ag12 or Pd16FeCl2 or anything reasonable')
    return parser.parse_args()

def generateNewDatabase(args):
    print("\n-->generateNewDatabase\n"+"-"*len("-->generateNewDatabase"))

    try:
        cluster = Atoms(args.formula)
    except Exception as e:
        print("ASE cannot understand your formular!")
        print(f"Error with {e} in '{args.formula}'") 
        sys.exit()
        
    print(f"Atomic list       : {cluster.get_atomic_numbers()}")


    Cell = Atoms("",pbc=False)
    Cell.set_cell([args.length,args.length,args.length])

    unique_atom_types = get_all_atom_types(Cell, cluster.get_atomic_numbers())
    blmin = closest_distances_generator(atom_numbers=unique_atom_types,ratio_of_covalent_radii=0.7)

    print(f"Unique Atom Types : {unique_atom_types}")
    sg = StartGenerator(Cell, cluster.get_atomic_numbers(), blmin)

    population_size = args.population

    d = PrepareDB(db_file_name=args.database,
                    population_size=population_size,
                    simulation_cell=Cell,
                    mult = args.mult,
                    charge = args.charge,
                    stoichiometry=cluster.get_atomic_numbers())
    
    
    da = DataConnection(args.database)


    """
    Let's start
    with some more candidates
    """
    for a in [sg.get_new_candidate() for i in range(int(population_size*1.5))]:
        d.add_unrelaxed_candidate(a)

    print(f"--| {int(population_size*1.5)} random start candidates added")

    counter=0
    for a in [x for x in da.get_all_unrelaxed_candidates() if x.info["key_value_pairs"]["extinct"]==0]:
        a.set_calculator(XTB(charge=args.charge,mult=args.mult))
        try:
            a.info["key_value_pairs"]["raw_score"] = - a.get_potential_energy()
            da.add_relaxed_step(a)
            counter+=1
            logG(f"--| Candidate #{counter} relaxed, and added to the database.")
        except Exception as e:
            logR(f"--| Relaxation stopped due to {e}; {repr(e)}")
            da.kill_candidate(a.info["confid"])
        
        if counter>=population_size:
            """
            Deactivate the candidates that we additionally added.
            """
            for i in da.get_all_unrelaxed_candidates():
                da.kill_candidate(i.info["confid"])
            return True

    return False
    
if __name__ == "__main__":
    args = runTheParser()

    print("\nClusterGen\n"+"="*len("ClusterGen")+"\n")
    print(f"Empirical Formula : {args.formula} ")
    print(f"Charge            : {args.charge} ")
    print(f"Multiplicity      : {args.mult} ")
    print(f"Population Size   : {args.population} ")
    print(f"Database          : {args.database} ")
    print(f"Cell Size         : {args.length} [Ang]  ")
    

    if os.path.isfile(args.database):
        print("Cluster file exists, we will just add candidates")
        if mutateAndAdd(args):

    else: 
        if generateNewDatabase(args):
            logG("--| New Database sucessfully created!")
        else:
            logR("--| Something went wrong!")

    

    



    


