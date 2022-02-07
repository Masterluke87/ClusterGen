import argparse
from random import random
from ase.io import read
import numpy as np
from ase.atoms import Atoms
from ase.ga.population import Population
from ase.ga.data import DataConnection,PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation,RattleMutation,PermutationMutation
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.io import read,write
from ase.atoms import Atoms
import os,sys
from util.XTBRunner import XTB

logG = lambda x: print("\033[92m"+x+"\033[00m")
logR = lambda x: print("\033[91m"+x+"\033[00m")

logH2 = lambda x: print("\n" +x+"\n"+"-"*len(x)+"\n")


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
    optional.add_argument('-db','--database', metavar='DB',default="cluster.db", type=str,
                help='Sqlite database to handle the structures') 
    optional.add_argument('-g','--generations', metavar='nGens',default=0, type=int,
                help='Number of Generations that should be generated.')                      
    optional.add_argument('-l','--length', metavar=' [Angstrom] ',default="12", type=str,
                help='Size of the cubic simulation cell in Angstrom')          
    optional.add_argument('-F','--formula', metavar='Formula',default="", type=str, required=True,
                help='Set empirical formula for the target cluster, eg. Ag12 or Pd16FeCl2 or anything reasonable')
    return parser.parse_args()

def generateNewDatabase(args):
    logH2("-->generateNewDatabase")

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
    

def mutateAndAdd(args):
    logH2("-->mutateAndAdd")

    da = DataConnection(args.database)

    mult   = da.get_param("mult")
    charge = da.get_param("charge")
    population_size = da.get_param("population_size")

    if (mult != args.mult):
        print(f"--| Note, Input mult: {args.mult} != database mult: {mult}. Using database value")
    if (charge != args.charge):
        print(f"--| Note, Input charge: {args.charge} != database charge: {charge}. Using database value")
    if (population_size != args.population):
        print(f"--| Note, Input pop: {args.population} != database charge: {population_size}. Using database value")

    mutation_probability = 0.3

    atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
    n_to_optimize = len(atom_numbers_to_optimize)
    slab = da.get_slab()
    all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
    blmin = closest_distances_generator(all_atom_types,ratio_of_covalent_radii=0.7)

    comp = InteratomicDistanceComparator(n_top=n_to_optimize,
                                        pair_cor_cum_diff=0.015,
                                        pair_cor_max=0.7,
                                        dE=0.02,
                                        mic=False)
    pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
    mutations = OperationSelector([1., 1., 1.],
                                    [MirrorMutation(blmin, n_to_optimize),
                                    RattleMutation(blmin, n_to_optimize),
                                    PermutationMutation(n_to_optimize)])


    for a in [x for x in da.get_all_unrelaxed_candidates() if x.info["key_value_pairs"]["extinct"]==0]: 
        a.set_calculator(XTB(charge=charge,mult=mult))
        try:
            a.info["key_value_pairs"]["raw_score"] = - a.get_potential_energy()
            logG("Candidate relaxed!")
            da.add_relaxed_step(a)
        except:
            logR("Relaxation stopped")
            da.kill_candidate(a.info["confid"])
    
    
    population = Population(data_connection=da,
           population_size=population_size,
                           comparator=comp)


    n_to_test = population_size * args.generations

    print(f"--| Will try to optimize {n_to_test} new candidates, correponding to {args.generations} new generation(s)!")

    counter = 1
    while(counter<= n_to_test):
        print(f'--| Now starting candidate # {counter}')
        a1, a2 = population.get_two_candidates()
        a3, desc = pairing.get_new_individual([a1, a2])
        if a3 is None:
            continue
        da.add_unrelaxed_candidate(a3, description=desc)
        if random() < mutation_probability:
            a3_mut, desc = mutations.get_new_individual([a3])
            if a3_mut is not None:
                da.add_unrelaxed_step(a3_mut, desc)
                a3 = a3_mut

        try:
            a3.set_calculator(XTB(charge=charge,mult=mult))
            a3.info["key_value_pairs"]["raw_score"] = - a3.get_potential_energy()
            da.add_relaxed_step(a3)
            population.update()
            logG(f"--| Optimization of candidate #{counter} was successful")
            counter+=1
            
        except:
            logR(f"--| Optimization failed..")
            da.kill_candidate(a3.info["confid"])

    logG(f"Done adding {counter-1} new candidates!")
    return True

def getStatistics(args):
    logH2("--> getStatistics")


    da = DataConnection(args.database)
    popSize = da.get_param("population_size")
    
    comp = InteratomicDistanceComparator(n_top=len(da.get_atom_numbers_to_optimize()),
                                        pair_cor_cum_diff=0.015,
                                        pair_cor_max=0.7,
                                        dE=0.02,
                                        mic=False)
    population = Population(data_connection=da,
                            population_size=popSize,
                            comparator=comp)


    print(f"# Generations         : {da.get_generation_number()}")
    print(f"# relaxed candidates  : {len(da.get_all_relaxed_candidates())}")
    print(f"Populations size:     : {popSize} ")

    if popSize<5:
        printSize=popSize
    else:
        printSize=5

    for i in range(1,da.get_generation_number()+1):
        print(str("{:3d} "+" {:3d}"*printSize+" {:12.6f}"*printSize).format(i,
                                *[x.info["confid"] for x in population.get_population_after_generation(i)][:printSize],
                                *[-x.info["key_value_pairs"]["raw_score"] for x in population.get_population_after_generation(i)][:printSize]))    


if __name__ == "__main__":
    args = runTheParser()

    
    print("\nClusterGen\n"+"="*len("ClusterGen")+"\n")
    logH2("User Input:")
    print(f"Empirical Formula : {args.formula} ")
    print(f"Charge            : {args.charge} ")
    print(f"Multiplicity      : {args.mult} ")
    print(f"Population Size   : {args.population} ")
    print(f"Database          : {args.database} ")
    print(f"Cell Size         : {args.length} [Ang]  ")
    print(f"# of Generations  : {args.generations}")
    print("\n")

    if not os.getenv("XTB_PATH"):
        logR("Please set the XTB_PATH to point to the xtb directory")
        sys.exit()

    if not os.path.isfile(args.database):
        if generateNewDatabase(args):
            logG("--| New Database sucessfully created!")
            if mutateAndAdd(args):
                getStatistics(args)
        else:
            logR("--| Something went wrong!")
            sys.exit()
    else:
        print("--| Cluster file exists, we will just add candidates")
        if mutateAndAdd(args):
            getStatistics(args)
            



    



    


