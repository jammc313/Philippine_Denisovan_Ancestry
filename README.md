# Philippine Denisovan Ancestry
A repository containing scripts used to simulate data for assessing models of Denisovan ancestry in Philippine populations.
- Each of the simulation python scripts (e.g. msp_sims_NULL1.py) generate data under the specified model. 
- The SLURM submission script "submit_msp_sims_arrayjob.sh" submits each of those simulations as a job array of 500 replicates per model. Be careful: even though the tree sequence output for each simulation is compressed, this will take A LOT of memory when using the current length of simulated sequence (L=100M)! 
- The script "allModels_consistency_checks.py" takes the n=500 compressed tree sequences for each model, and calculates various summary statistics (f4 and D estimates of Denisovan & Neanderthal ancestry, introgressing tract lengths), outputting the results as numpy arrays.
- The jupyter notebook "allModels_consistency_checks_plotting.ipynb", loads in those results arrays for a specified model and plots the distributions of results.
