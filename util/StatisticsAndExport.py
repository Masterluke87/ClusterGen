from ase.ga.data import DataConnection
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.population import Population
from ase.visualize import view

logG = lambda x: print("\033[92m"+x+"\033[00m")
logR = lambda x: print("\033[91m"+x+"\033[00m")
logH2 = lambda x: print("\n" +x+"\n"+"-"*len(x)+"\n")


def write_xyz(fileobj, images, comment='', fmt='%22.15f'):
    """
    Shameless copy from ASE project to add the energy to the comment line
    """
    comment = comment.rstrip()
    if '\n' in comment:
        raise ValueError('Comment line should not have line breaks.')
    for atoms in images:
        natoms = len(atoms)
        fileobj.write('%d\n%s\n' % (natoms, -atoms.info["key_value_pairs"]["raw_score"]))
        for s, (x, y, z) in zip(atoms.symbols, atoms.positions):
            fileobj.write('%-2s %s %s %s\n' % (s, fmt % x, fmt % y, fmt % z))


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
    if args.number is not None:
        printSize=int(args.number)



    for i in range(1,da.get_generation_number()+1):
        print(str("{:4d} "+" {:4d}"*printSize+" {:12.6f}"*printSize).format(i,
                                *[x.info["confid"] for x in population.get_population_after_generation(i)][:printSize],
                                *[-x.info["key_value_pairs"]["raw_score"] for x in population.get_population_after_generation(i)][:printSize]))    


def exportOrViewCandidates(args):
    logH2("--> Export")
    if (args.export):   
        print(f"#Exporting {int(args.number)} candidates.")
    if (args.export):   
        print(f"#Viewing {int(args.number)} candidates.")
    


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

    
    candidates = population.get_current_population()[:int(args.number)]
    if (args.export):  
        with open(args.database.split(".")[0]+".xyz","w") as f:
            write_xyz(f, candidates)
        logG("Success.")
        
    if (args.view):
        view(candidates)
    return candidates