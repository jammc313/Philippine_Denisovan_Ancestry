import msprime as msp
import sys
import tszip

# # msprime simulations
# 
# Recombination = 1.25 e-8 (equal to mu)
# We have 1 run with replicated simulations. 
# L=100M
# 10 samples per population. The full haplotype matrix is returned.
# Simulates the coalescent with recombination under the specified model parameters and returns the resulting :class:`tskit.TreeSequence`. Note that Ne is the effective diploid population size (so the effective number of genomes in the population is 2*Ne), but ``sample_size`` is the number of (monoploid) genomes sampled.
# 
# Source pops related by ...
# nS vector of sample sizes (monoploid genomes) in each S, i.e. (nS1, nS2, nS3)
# tS vector ages of samples in each S, i.e. (tS1, tS2, tS3)
# t12 time of S1,S2 split (first split)
# t123 time of ((S1, S2), S3) split
# ta vector of admixture times from S2 and S3 to S1
# f vector of admixture fractions from from S2 and S3 to S1
# N vector of effective sizes in order of population numbers
# L is length of simulated region in bases (set to 100M here)
# 
# Here we are looking at admixed segment length distributions, f4 tests, and D tests results.
# The model we are testing is ALT2B, with an admixture event from Neanderthal into the shared ancestor of all Non-Africans, a private Denisovan admixture event into Papuans, and a Denisovan admixture event private to Ayta. The difference from ALT2 is that here includes an admixture event from East Asian into Ayta.

# ## Enter msprime functions

model="ALT2B"
simrep=sys.argv[1]

# enter functions

DEN0, DEN1, DEN2, DEN3, AFR, CEU, EAS, PAP, AYT, NEA, CHM = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10


def set_up_pops(nS, tS):
    samples = [msp.Sample(population=DEN3, time=tS[1])]*(1*nS[0]) # Denisovan 3 (Altai)
    samples.extend([msp.Sample(population=AFR, time=tS[0])]*(1*nS[0])) # Africa
    samples.extend([msp.Sample(population=CEU, time=tS[0])]*(1*nS[0])) # European
    samples.extend([msp.Sample(population=EAS, time=tS[0])]*(1*nS[0])) # East Asian
    samples.extend([msp.Sample(population=PAP, time=tS[0])]*(1*nS[0])) # Papuan
    samples.extend([msp.Sample(population=AYT, time=tS[0])]*(1*nS[0])) # Negrito (Ayta)
    samples.extend([msp.Sample(population=NEA, time=tS[2])]*(1*nS[0])) # Neanderthal
    samples.extend([msp.Sample(population=CHM, time=tS[0])]*(1*nS[0])) # Chimp
    return samples

def set_up_demography(t78, t68, t85, t54, t10, t20, t03, t93, t34, t410, ta1, ta2, ta3, ta4, f):
    #divergence of source populations (topology of tree)
    source_divergence = [msp.MassMigration(time=t78, source=PAP, destination=AYT, proportion=1),
                        msp.MassMigration(time=t68, source=EAS, destination=AYT, proportion=1),
                        msp.MassMigration(time=t85, source=AYT, destination=CEU, proportion=1), 
                        msp.MassMigration(time=t54, source=CEU, destination=AFR, proportion=1), 
                        msp.MassMigration(time=t10, source=DEN1, destination=DEN0, proportion=1), 
                        msp.MassMigration(time=t20, source=DEN2, destination=DEN0, proportion=1),
                        msp.MassMigration(time=t03, source=DEN0, destination=DEN3, proportion=1),
                        msp.MassMigration(time=t93, source=NEA, destination=DEN3, proportion=1),
                        msp.MassMigration(time=t34, source=DEN3, destination=AFR, proportion=1),
                        msp.MassMigration(time=t410, source=AFR, destination=CHM, proportion=1)] 
    #admixture times and proportions
    admix = [msp.MassMigration(time=ta1, source=AYT, destination=DEN2, proportion=f[1]), #fraction from Denisovan 2 to Ayta
            msp.MassMigration(time=ta2, source=PAP, destination=DEN1, proportion=f[2]), #fraction from Denisovan 1 to Papuan
            msp.MassMigration(time=ta3, source=CEU, destination=NEA, proportion=f[3]), #fraction from Neanderthal to OoA pops
            msp.MassMigration(time=ta4, source=AYT, destination=EAS, proportion=f[0])] #fraction from East Asian into Ayta
    #population parameter changes
    N_change = [msp.PopulationParametersChange(time=2400, initial_size=2000, growth_rate=0, population_id=CEU),
               msp.PopulationParametersChange(time=2800, initial_size=5000, growth_rate=0, population_id=CEU)]
    #combine and sort the demography
    demography = source_divergence + admix + N_change
    return sorted(demography, key = lambda x: x.time)


# Enter data
L=10000000
mu=1.25e-8
r=1.25e-8

generation_time=25
T78=46000/generation_time  # PAP joins AYT 
T68=50000/generation_time  # EAS joins AYT  
T85=55000/generation_time  # AYT joins CEU 
T54=70000/generation_time  # CEU joins AFR
T10=200000/generation_time # Denisovan 1 joins Denisovan 0
T20=200000/generation_time # Denisovan 2 joins Denisovan 0 
T03=300000/generation_time # Denisovan 0 joins Denisovan 3 (Altai) 
T93=400000/generation_time # Neanderthal joins Denisovan 3 (Altai) 
T34=600000/generation_time # Denisovan 3 (Altai) joins AFR
T410=4000000/generation_time # AFR joins Chimp
 
TA1=35000/generation_time   
TA2=45000/generation_time   
TA3=68000/generation_time   
TA4=2500/generation_time   

TS_NEA=60000/generation_time
TS_DEN3=40000/generation_time

NumSamples=80
nS=[10]
tS=[0,TS_DEN3,TS_NEA]
f=[0.10, 0.06, 0.04, 0.02]
N=[1500,1500,1500,1500,15000,5000,3500,3500,3500,2000,30000]
seed=None

samples = set_up_pops(nS,tS)
demography = set_up_demography(T78, T68, T85, T54, T10, T20, T03, T93, T34, T410, TA1, TA2, TA3, TA4, f)
pops = [msp.PopulationConfiguration(initial_size = n) for n in N]

ts = msp.simulate(samples=samples, Ne=N[0], population_configurations=pops, demographic_events=demography, mutation_rate=mu, length=L, recombination_rate=r, record_migrations=True, random_seed=seed)

# output resulting tree sequences to compressed .tsz file
tszip.compress(ts, "tree_seq_files/{}_model_{}.tsz".format(model, simrep))
