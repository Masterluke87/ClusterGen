from ase.ga.data import DataConnection
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.population import Population
from ase.visualize import view
from templates.QM import QMTemplates
import os

logG = lambda x: print("\033[92m"+x+"\033[00m")
logR = lambda x: print("\033[91m"+x+"\033[00m")
logH2 = lambda x: print("\n" +x+"\n"+"-"*len(x)+"\n")


def write_xyz(fileobj, atoms, charge=0, mult=1, fmt='%22.15f'):
    natoms = len(atoms)
    fileobj.write('{:d}\nE={:16.8f};c={:d};m={:d}\n'.format(natoms, -atoms.info["key_value_pairs"]["raw_score"],charge,mult))
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


    print(f"--| # Generations         : {da.get_generation_number()}")
    print(f"--| # relaxed candidates  : {len(da.get_all_relaxed_candidates())}")
    print(f"--| Populations size:     : {popSize} ")
 
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
    if (args.view):   
        print(f"#Viewing {int(args.number)} candidates.")
    


    da = DataConnection(args.database)
    popSize = da.get_param("population_size")
    mult   = da.get_param("mult")
    charge = da.get_param("charge")
    
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
        """
        First we create a folder 
        """ 
        PathToDatabase = os.path.dirname(args.database)
        Prefix = os.path.basename(args.database.split(".")[0])
        
        if PathToDatabase.strip() == "":
            PathToDatabase = "./"
        
        calcTemplate=args.template
        calcProg = args.prog
        
        print(f"--| # Generations         : {da.get_generation_number()}")
        print(f"--| # relaxed candidates  : {len(da.get_all_relaxed_candidates())}")
        print(f"--| Populations size:     : {popSize} ")
        print(f"--| PathToDatabase        : {PathToDatabase}")
        print(f"--| Prefix                : {Prefix}")
        print(f"--| QM Program            : {calcProg} ")
        print(f"--| QM Template           : {calcTemplate} ")
        

        if not os.path.exists(PathToDatabase+"/"+Prefix):
            os.mkdir(PathToDatabase+"/"+Prefix)

        for c,i in enumerate(candidates):
            if not os.path.exists(PathToDatabase+"/"+Prefix+f"/{c}/"):
                os.mkdir(PathToDatabase+"/"+Prefix+f"/{c}/")
            with open(PathToDatabase+"/"+Prefix+f"/{c}/{c}.xyz","w") as f:
                write_xyz(f, i, charge=charge,mult=mult)
            QMTemplates(calcProg,
                        calcTemplate,
                        PathToDatabase+"/"+Prefix+f"/{c}/",
                        coord=f"{c}.xyz",
                        charge=charge,
                        mult=mult)

            
        logG("Success.")
        
    if (args.view):
        view(candidates)
    return candidates