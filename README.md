# DenovoSim
Simulate de novo mutations using model based on trinucleotide mutation rate model introduced in Samocha et. al, 2014. Users can condition on number of patients in a study, and number of mutations observed. Updates to the model planned and tracked in issues section.

# Running the simulations
DenovoSim takes the following command line arguments:
# running unsupervised simulation from ~/software/DenovoSim/scripts
cd /software/DenovoSim/scripts

bsub -J "coding_sim[1-1000:2]" -R'select[mem>3000] rusage[mem=3000]' -M3000

/software/R-3.2.2/bin/Rscript /nfs/users/nfs_p/ps14/software/DenovoSim/scripts/DenovoSim.R 

--verbose

--n_probands=4294
Number of probands in your cohort.

--iterations=10
Number of iterations per file.

--iteration_start=\$LSB_JOBINDEX
Which iteration this script should start at - it will also be the digit that is added onto 'base_name'.

--n_chunks=1 
This describes how many chunks this particular run will split into. For instance, if iterations = 100 and n_chunks = 10, you will get 10 files output with 10 simulations each. I recommend just keeping it at one and using bjobs array to generate multiple files.

--base_name=~/experiments/simulated_data/coding_sim  
This will save as coding_sim.n.txt where n is the iteration number

--regions=~/reference_data/gencode_exons.txt
Set of regions in which simulations should be mutated. Required columns are chr, start, stop.

Full run would look like this:
bsub -J "coding_sim[1-1000:2]" -R'select[mem>3000] rusage[mem=3000]' -M3000 /software/R-3.2.2/bin/Rscript /nfs/users/nfs_p/ps14/software/DenovoSim/scripts/DenovoSim.R--verbose --n_probands=4294 --iterations=10 --iteration_start=\$LSB_JOBINDEX --n_chunks=1  --base_name=~/experiments/simulated_data/coding_sim --regions=~/reference_data/gencode_exons.txt

Note that 'by' digit in the job array should match --iterations. E.g. if you run coding_sim[1-1000:5] then --iterations=5. This would generate files coding_sim.1.txt, coding_sim.6.txt, coding_sim.11.txt, ...

Where coding_sim.1.txt has iterations 1,2,3,4,5 and coding_sim.6.txt has iterations 6,7,8,9,10, etc.

You can write a quick script in R or bash to load these in one by one and rbind the data frames to get one massive dataframe with your simulations (and simulation number will be in column $iterations)



