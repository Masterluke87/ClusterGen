import argparse
from ase.io import read
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError
import numpy as np
from ase.visualize import view
from ase.atoms import Atoms
from ase.ga.data import DataConnection
from ase.optimize import BFGS
from ase.io import read,write

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
                   
    required.add_argument('-F','--formula', metavar='Formula',default="", type=str, required=True,
                help='Set empirical formula for the target cluster, eg. Ag12 or Pd16FeCl2 or anything reasonable')
    return parser.parse_args()


if __name__ == "__main__":
    args = runTheParser()
    print("ClusterGen\n"+"="*len("ClusterGen")+"\n")
    print(f"Empirical Formula : {args.formula} ")
    print(f"Charge            : {args.charge} ")
    print(f"Multiplicity      : {args.mult} ")
    print(f"Population Size   : {args.population} ")
    print(f"Database          : {args.database} ")
    print(f"Cell Size         : {args.length} [Ang]  ")
    
    



    


