# ClusterGen
Script collection to generate low energy clusters via Genetic Algorithm and XTB. 

## Dependencies

## Configuration

## Examples

### Transition metal doped gold cluster, $Au12Fe$

```bash
python GenAndMutate.py -db cluster.db -F Au12Fe -p 20 -c 0 -m 5 -l 15.0 -f 2.0 
```
This command will generate and initialize a new $Au_{12}Fe$ population (with charge 0 and a multiplicity of 5), which is saved in the *cluster.db* database file. The size of the population is 20 and the size of the initial simulation box (where the atoms are placed) is 15.0 $\AA$. The *-f x* parameter is used to generate x times more starting candidates than the population size (mainly because not all starting candidates converge to a local minumum).   

