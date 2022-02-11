from util.StatisticsAndExport import getStatistics,exportOrViewCandidates
import argparse

logG = lambda x: print("\033[92m"+x+"\033[00m")
logR = lambda x: print("\033[91m"+x+"\033[00m")
logH2 = lambda x: print("\n" +x+"\n"+"-"*len(x)+"\n")

def runTheParser():
    parser = argparse.ArgumentParser(description="""
    ClusterGen -- Export. Export low energy, unique candidates from the population.
    """)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-S","--statistics",action='store_true', help="Show statistics")
    optional.add_argument("-e","--export",action='store_true', help="Export the lowest n candidates of a population")
    optional.add_argument("-v","--view",action='store_true', help="View lowest n candidates in ASE gui window")
    required.add_argument('-db','--database', metavar='DB',default="cluster.db", type=str,
                help='Sqlite database to handle the structures') 
    required.add_argument('-n','--number', metavar='num',default=5, type=str, required=True,
                help='How many candidates should be exported?')
    return parser.parse_args()



if __name__ == "__main__":
    args = runTheParser()

    
    print("\nClusterGen -- Export\n"+"="*len("ClusterGen -- Export")+"\n")
    args = runTheParser()

    if args.statistics:
        getStatistics(args)


    if (args.export) or (args.view):
        X = exportOrViewCandidates(args)