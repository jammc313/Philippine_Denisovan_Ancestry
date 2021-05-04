#!/usr/bin/env python
# coding: utf-8

# ## In this notebook we:
# ### a) show consistency in our f4 estimates using simulations. We visualize the distribution of results in statistics across simulation replicates. These distributions could be compared against observed values.
# ### b) we are not performing formal model selection, because we are not fitting our models to data. Rather, this is a simple qualitative way of showing some expectations under simulated models and comparing them to observed data. 
# ### c) Results are saved from simulations and used to construct PCs onto which we can project observed values to compare model consistency with observed data in a PCA

# In[1]:


import msprime as msp
import tskit
import tszip
import random
import numpy as np
import pandas as pd
import sys
import os
import math
import allel


# In[2]:


# Global toggle to save / not-save figures
isSaveFigures = False

# In[3]:


# some common variables from all models sims
DEN0, DEN1, DEN2, DEN3, AFR, CEU, EAS, PAP, AYT, NEA, CHM = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10

L=100000000
mu=1.25e-8
r=1.25e-8
generation_time=25


# In[4]:


# some statistic functions

def make_sample_list(ts):
    sample_list = []
    for pop in range(DEN3, ts.num_populations):
        sample_list.append(ts.samples(pop).tolist())
    return sample_list

# input sims genotype matrix, returns allele count array of pops sampled (10 haplotypes per population)
def pop_sample_ac(geno_mat):
    haplo_arr = allel.HaplotypeArray(geno_mat)
    ac_one = haplo_arr[:,0:10].count_alleles()
    ac_two = haplo_arr[:,10:20].count_alleles()
    ac_three = haplo_arr[:,20:30].count_alleles()
    ac_four = haplo_arr[:,30:40].count_alleles()
    ac_five = haplo_arr[:,40:50].count_alleles()
    ac_six = haplo_arr[:,50:60].count_alleles()
    ac_seven = haplo_arr[:,60:70].count_alleles()
    ac_eight = haplo_arr[:,70:80].count_alleles()
    # stack arrays with frames = population allele counts for all SNPs
    arrays = [ac_one,ac_two,ac_three,ac_four,ac_five,ac_six,ac_seven,ac_eight]
    ac_All = np.stack(arrays, axis=0)
    return ac_All

# function to estimate introgression using f4 on both Papuan and Ayta 
def f4_DenAncEst(sims, option):
    f4_res = []
    for testpop in (4,5): # do for both Papuan and Ayta
        if option=="a":  # f4(AFR, NEA, CEU, X) / f4(AFR, NEA, CEU, DEN)
            num_idx = [1,6,2,testpop] # supply indexing list
            den_idx = [1,6,2,0]   
        elif option=="b":
            num_idx = [1,6,3,testpop] # supply indexing list
            den_idx = [1,6,3,0]  
        else:
            raise ValueError("f4 test options must be a or b")
        test_sample_list = make_sample_list(sims)
        sample_set_num = [test_sample_list[i] for i in num_idx] # get sample list for numerator
        sample_set_den = [test_sample_list[i] for i in den_idx] # get sample list for denominator
        num = ts.f4(sample_sets=sample_set_num) # estimate f4 for numerator
        den = ts.f4(sample_sets=sample_set_den) # estimate f4 for denominator
        f4_res.append(num/den)
    return f4_res

# function to estimate introgression using f4
def f4_AncEst(sims, option):
    f4_res = []
    if option=="c": # f4(AFR, PAP, EAS, AYT) Australasian ancestry in Ayta
        idx = [1,4,3,5] # supply indexing list
    elif option=="d": # f4(AFR, EAS, PAP, AYT) East Asian ancestry in Ayta
        idx = [1,3,4,5] # supply indexing list
    elif option=="e": # f4(AFR, DEN3, PAP, AYT) Direct comparison of Denisovan ancestry in Papuan & Ayta
        idx = [1,0,4,5] # supply indexing list
    else:
        raise ValueError("f4 test options must be c, d or e")
    test_sample_list = make_sample_list(sims)
    sample_set = [test_sample_list[i] for i in idx] # get indexed sample list
    f4_res.append(float(ts.f4(sample_sets=sample_set)))
    return f4_res

# D test function, given archaic sample, will estimate D for PAP & AYT
def D_AncEst(sims, archaic):
    if archaic=="DEN":
        D_idx = np.array([[1,0,2,4],[1,0,2,5]])   # D{AFR, DEN3, CEU, X} with X={PAP,AYT}
    elif archaic=="NEA":
        D_idx = np.array([[7,6,1,4],[7,6,1,5]])   # D{CHM, NEA, AFR, X} with X={PAP,AYT}
    else:
        raise ValueError("Archaic population must be DEN or NEA")
    ac_AllSamples = pop_sample_ac(sims.genotype_matrix())  # convert tree seqs to allele count arrays
    # D test returning d,se,z for PAP & AYT
    D_est_arr = np.array([(allel.average_patterson_d(ac_AllSamples[x[0]], ac_AllSamples[x[1]], ac_AllSamples[x[2]], ac_AllSamples[x[3]], 1000)[0]) for x in D_idx])
    return D_est_arr

# get migrating tracts
def get_mig_segs(sims):
    # Get all tracts that migrated into the archaic populations
    mig_segs_DEN_ayt, mig_segs_DEN_pap, mig_segs_NEA_ayt, mig_segs_NEA_pap = [], [], [], []
    den_nd, nea_nd = [], []
    for migration in ts.migrations():
            if migration.dest == DEN1 or migration.dest == DEN2:
                den_nd.append((migration.left, migration.right, migration.node))
            elif migration.dest == NEA:
                nea_nd.append((migration.left, migration.right, migration.node))
    for seg in den_nd: 
        seg_node = seg[2]
        if sims.tables.nodes[seg_node].population == 7: # if birth population is PAP
            mig_segs_DEN_pap.append(seg[1]-seg[0])
        elif sims.tables.nodes[seg_node].population == 8: # if birth population is AYT
            mig_segs_DEN_ayt.append(seg[1]-seg[0])
    for seg in nea_nd: 
        seg_node = seg[2]
        if sims.tables.nodes[seg_node].population == 7:
            mig_segs_NEA_pap.append(seg[1]-seg[0])
        elif sims.tables.nodes[seg_node].population == 8:
            mig_segs_NEA_ayt.append(seg[1]-seg[0])
    return np.array(mig_segs_DEN_ayt), np.array(mig_segs_DEN_pap), np.array(mig_segs_NEA_ayt), np.array(mig_segs_NEA_pap)


# In[5]:


tajimas_pi = []
D_DenAnc = [] # D(AFR, DEN3, CEU, X) for X={PAP,AYT} returning d,se,z
D_NeaAnc = [] # D(CHM, NEA, AFR, X) for X={PAP,AYT} returning d,se,z
f4_DenAnc_a = [] # f4(AFR, NEA, CEU, X) / f4(AFR, NEA, CEU, DEN3) for X={PAP,AYT}
f4_DenAnc_b = [] # f4(AFR, NEA, EAS, X) / f4(AFR, NEA, EAS, DEN3) for X={PAP,AYT}
f4_Est_c = [] # f4(AFR, PAP, EAS, AYT) Australasian ancestry in Ayta
f4_Est_d = [] # f4(AFR, EAS, PAP, AYT) East Asian ancestry in Ayta
f4_Est_e = [] # f4(AFR, DEN3, PAP, AYT) Direct comparison of Denisovan ancestry in Papuan & Ayta
mean_tractL_DEN = []
mean_tractL_NEA = []
prop_tractL_DEN = []
prop_tractL_NEA = []


# ### Assign model and get list of .ts files for it

# In[12]:


arg_list = sys.argv
model = arg_list[1] + '_'
list_ts_files = [x for x in os.listdir('tree_seq_files/') if x.startswith(model)]
n = len(list_ts_files)


# ## 1) Compute statistics from list of .ts files

# In[13]:


# iterate through that list, getting stats for every file object
for file in list_ts_files:
    ts = tszip.decompress('tree_seq_files/'+file)
    sample_list = make_sample_list(ts)
    tajimas_pi.append(ts.diversity(sample_sets=sample_list))
    D_DenAnc.append(D_AncEst(ts, "DEN"))
    D_NeaAnc.append(D_AncEst(ts, "NEA"))
    f4_DenAnc_a.append(f4_DenAncEst(ts, "a"))
    f4_DenAnc_b.append(f4_DenAncEst(ts, "b"))
    f4_Est_c.append(f4_AncEst(ts, "c"))
    f4_Est_d.append(f4_AncEst(ts, "d"))
    f4_Est_e.append(f4_AncEst(ts, "e"))
    mig_segs_DEN_AYT, mig_segs_DEN_PAP, mig_segs_NEA_AYT, mig_segs_NEA_PAP  = get_mig_segs(ts)
    mean_tractL_DEN.append(np.array([np.mean(mig_segs_DEN_PAP), np.mean(mig_segs_DEN_AYT)]))
    mean_tractL_NEA.append(np.array([np.mean(mig_segs_NEA_PAP), np.mean(mig_segs_NEA_AYT)])) 
    prop_tractL_DEN.append(np.array([np.sum(mig_segs_DEN_PAP)/L, np.sum(mig_segs_DEN_AYT)/L]))
    prop_tractL_NEA.append(np.array([np.sum(mig_segs_DEN_PAP)/L, np.sum(mig_segs_DEN_AYT)/L]))


# ### Save / load numpy arrays

# In[20]:


# save down appropriately named numpy files
tajimas_pi_arr = np.array(tajimas_pi)
np.save(model+'tajimas_pi.npy', tajimas_pi_arr)
D_DenAnc_arr = np.array(D_DenAnc)
np.save(model+'D_DenAnc.npy', D_DenAnc_arr)
D_NeaAnc_arr = np.array(D_NeaAnc)
np.save(model+'D_NeaAnc.npy', D_NeaAnc_arr)
f4_DenAnc_a_arr = np.array(f4_DenAnc_a)
np.save(model+'f4_DenAnc_a.npy', f4_DenAnc_a_arr)
f4_DenAnc_b_arr = np.array(f4_DenAnc_b)
np.save(model+'f4_DenAnc_b.npy', f4_DenAnc_b_arr)
f4_Est_c_arr = np.array(f4_Est_c)
np.save(model+'f4_Est_c.npy', f4_Est_c_arr)
f4_Est_d_arr = np.array(f4_Est_d)
np.save(model+'f4_Est_d.npy', f4_Est_d_arr)
f4_Est_e_arr = np.array(f4_Est_e)
np.save(model+'f4_Est_e.npy', f4_Est_e_arr)
mean_tractL_DEN_arr = np.array(mean_tractL_DEN)
np.save(model+'mean_tractL_DEN.npy', mean_tractL_DEN_arr)
mean_tractL_NEA_arr = np.array(mean_tractL_NEA)
np.save(model+'mean_tractL_NEA.npy', mean_tractL_NEA_arr)
prop_tractL_DEN_arr = np.array(prop_tractL_DEN)
np.save(model+'prop_tractL_DEN.npy', prop_tractL_DEN_arr)
prop_tractL_NEA_arr = np.array(prop_tractL_NEA)
np.save(model+'prop_tractL_NEA.npy', prop_tractL_NEA_arr)




