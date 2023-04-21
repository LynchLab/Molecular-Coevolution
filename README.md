# Molecular-Coevolution
The computer program TwoSites.cpp is structured to determine the long-term steady-state distributions for a two-locus biallelic system in a discrete-generation-time setting incorporating consecutive generations of reversible mutation, recombination, selection, and random genetic drift. After a sufficient burn-in period, the program compiles the long-term steady-state probability distribution of the four haplotypes as well as the rates of transition among them.

Underlying details of the population-genetic model can be found at the top of the code.

On lines 49-71 of the code, the user enters the mutation rates, selection coefficients, recombination rates, and relative effective population sizes of the two genetic loci.

Lines 102, 108, and 216 can be edited to name the output files.

The code is written so that runs can be made in parallel for 19 population sizes, which are entered on lines 224-242. If it is desired to alter the number of runs, edit line 107; and if >19 are to be run, the population size array and lines 249-267 need to be expanded.

Lines 249-267 can be edited to define the numbers of generations over which to simulated the runs for each population size. 

On lines 296-327, an effort is made to reduce the population sizes and increase the mutation rates and selection coefficients in coordination to yield the desired ratios of mutation and selection to random genetic drift. This speeds up the runtimes considerably, generally with the same results, but it is desired not to adhere to such rescaling, enter “kfac = 1;” on line 313.

There will be two sets of output files: “slurm” files containing the cumulative statistics allow the user to determine whether the runs have proceeded for sufficient time to equilibrate; the “dataout” files give the final results for each of the 19 population sizes. 



To run the program, enter on the unix command line:

module load intel/2019.4
icc -o TwoSites TwoSites.cpp -lm -lgsl
sbatch –array=1-19 TwoSites.sh

(The first line may need to be modified, depending on the system involved. The 19 in the final line needs to be edited if a different number of runs is being made).



A batch shell file must be provided in the local file space to set up the series of runs:

#! /bin/bash

#SBATCH -A mlynch11
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH -n 1
#SBATCH -t 11-4:00

echo "Running the script in parallel"
./TwoSites $SLURM_ARRAY_TASK_ID


Here, the #SBATCH lines will need to be modified to the users specifications.

If all is operating properly, upon submission of the job, the slurm and dataout files should immediately appear, and after the user-set burn-in time, the slurm files will begin to be periodically updated with the cumulative statistics, with the dataout files becoming populated after each run completes. The completion times will vary depending on the user-defined run lengths.  


Summary file of results:

Upon completion of all runs, an sbatch file called concatenate.sh (for example; lines below) can be submitted, which will create an output called summary.txt that contains the comma-delimited stacked set of final results, which can be imported to a spreadsheet. “dataoutxxx” needs to be edited to give the prefix of the output file names, which will come out in parallel as dataoutxxx_1.txt to dataoutxxx_19.txt.

To run this file, type on the command line: sbatch concatenate.sh


#! /bin/bash

#SBATCH -A mlynch11

#SBATCH -n 1
#SBATCH -t 0-4:00

myFiles1=`ls dataoutxxx_?.txt`
myFiles2=`ls dataoutxxx_??.txt`
cat $myFiles1 $myFiles2 > summary.txt

 
