# Python script which simulate data from the 3 different scenarios 

import msprime
import math
import numpy as np
from multiprocessing import Pool
import time
import concurrent.futures
import tskit 


###################################################################################################
# Scenario with a unique divergence for the 3 populations : simu pool of concomitent_div scenario #
###################################################################################################
# Definition of the scenario and parameters
def concomitent_div(Ne_ancpp=100000,Ne_fu=20000,Ne_wp=1000,Ne_sp=500,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_found=200,Ne_bot_sp=50,Ne_bot_wp=50,Ne_bot_fu=50,
                 tsplit_PP=10000,tsplit_L=500,m_spwp_anc=0.00005,m_spwp_rec=0.00001,m_spfu_anc=0.000005,m_wpfu_anc=0.0005,m_spfu_rec=0.000005,m_wpfu_rec=0.0005,
                 n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100):
# Definition demo
    # Setting exponential parameter for the evolution of the Ne of the population
    rexpsp_avbot=math.log(Ne_sp/Ne_sp_found)/(tsplit_PP-tsplit_L) #Setting the evolution of the SP population size before the bottleneck
    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/tsplit_L #Setting the evolution of the SP population size after the bottleneck
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_PP-tsplit_L) #Setting the evolution of the WP population size before the bottleneck
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/tsplit_L #Setting the evolution of the WP population size after the bottleneck
    rexpfu_avbot=math.log(Ne_fu/Ne_fu_found)/(tsplit_PP-tsplit_L) #Setting the evolution of the FU population size before the bottleneck
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/tsplit_L #Setting the evolution of the FU population size after the bottleneck

    # Setting the populations name and their initial size 
    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    # Setting the migration rate
    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], rate=m_spwp_rec)

    # Setting the time of modifications in population size and migration flow
    demoIND.add_population_split(time=tsplit_PP, derived=["WP","FU","SP"], ancestral="ANCPP")
    demoIND.add_population_parameters_change(time=tsplit_L,population="SP", initial_size=Ne_sp,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=tsplit_L,population="WP", initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=tsplit_L,population="FU", initial_size=Ne_fu,growth_rate=rexpfu_avbot)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L,populations=["SP","WP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L,populations=["FU","WP"],rate=m_wpfu_anc)
    demoIND.sort_events()
    return demoIND

# msprime simulate and output a genotype matrix
def worker(seeds,nind,demoIND,chr_length,nchr):
    ancestry_seed, mutation_seed = seeds
    ts = msprime.sim_ancestry(
        samples=nind,
        demography=demoIND,
        sequence_length=chr_length,
        recombination_rate=4e-8,
        discrete_genome=True,
        random_seed=int(ancestry_seed))
    ts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=int(mutation_seed))
    geno_mat=ts.genotype_matrix()
    return(geno_mat)

# Setting parameters and do msprime simulation to output a genotype matrix
def simu_concomitent_div(Ne_ancpp=100000,Ne_fu=20000,Ne_wp=1000,Ne_sp=500,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_found=200,Ne_bot_sp=50,Ne_bot_wp=50,Ne_bot_fu=50,
                 tsplit_PP=10000,tsplit_L=500,m_spwp_anc=0.00005,m_spwp_rec=0.00001,m_spfu_anc=0.000005,m_wpfu_anc=0.0005,m_spfu_rec=0.000005,m_wpfu_rec=0.0005,
                 n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100,n_cpu=6):
    demoIND=concomitent_div(Ne_ancpp,Ne_fu,Ne_wp,Ne_sp,Ne_sp_found,Ne_wp_found,Ne_fu_found,Ne_bot_sp,Ne_bot_wp,Ne_bot_fu,
                 tsplit_PP,tsplit_L,m_spwp_anc,m_spwp_rec,m_spfu_anc,m_wpfu_rec,m_spfu_rec,m_wpfu_anc,
                 n_sp,n_wp,n_fu,nchr,chr_length,seed)
    nind={"SP":n_sp, "WP":n_wp, "FU":n_fu}
    indv_names=[]
    deb=time.time()
    for k, v in nind.items():
        for i in range(v):
            indv_names.append(k+"_"+str(i+1))
    # params
    num_replicates = nchr
    rng = np.random.RandomState(int(seed))
    seeds = rng.randint(1, 2**31, size=(num_replicates, 2))
    geno_matrices = []
    print("start")
    with concurrent.futures.ProcessPoolExecutor(max_workers=int(n_cpu)) as executor:
        for geno_mat in executor.map(worker, seeds, [nind] * num_replicates, [demoIND] * num_replicates, [chr_length] * num_replicates, [nchr] * num_replicates):
            geno_matrices.append(geno_mat)
    end=time.time()
    print("Time = ",end-deb)
    geno_final = geno_matrices[0]
    for i in range (1,len(geno_matrices)):
        geno_final=np.concatenate((geno_final,geno_matrices[i]),axis=0)
    return geno_final

##################################################################################################
# Scenario with a unique divergence for the 3 populations : simu ind of concomitent_div scenario #
##################################################################################################

# Definition of scenario and parameters for simulate individual data
def simu_concomitent_div_ind(Ne_ancpp=100000,Ne_fu=20000,Ne_wp=1000,Ne_sp=500,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_found=200,Ne_bot_sp=50,Ne_bot_wp=50,Ne_bot_fu=50,
                 tsplit_PP=10000,tsplit_L=500,m_spwp_anc=0.00005,m_spwp_rec=0.00001,m_spfu_anc=0.000005,m_wpfu_anc=0.0005,m_spfu_rec=0.000005,m_wpfu_rec=0.0005,
                 n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100,seed_dr=1):
   
    rexpsp_avbot=math.log(Ne_sp/Ne_sp_found)/(tsplit_PP-tsplit_L) 
    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/tsplit_L 
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_PP-tsplit_L)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/tsplit_L
    rexpfu_avbot=math.log(Ne_fu/Ne_fu_found)/(tsplit_PP-tsplit_L)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/tsplit_L

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], rate=m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","FU","SP"], ancestral="ANCPP")
    demoIND.add_population_parameters_change(time=tsplit_L,population="SP", initial_size=Ne_sp,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=tsplit_L,population="WP", initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=tsplit_L,population="FU", initial_size=Ne_fu,growth_rate=rexpfu_avbot)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L,populations=["SP","WP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L,populations=["FU","WP"],rate=m_wpfu_anc)
    demoIND.sort_events()

    # generating data
    nind={"SP":n_sp, "WP":n_wp, "FU":n_fu}
    indv_names=[]
    for k, v in nind.items():
        for i in range(v):
            indv_names.append(k+"_"+str(i+1))

    all_mutations = []
    for i,ts in enumerate(msprime.sim_ancestry(samples=nind, demography=demoIND, random_seed=seed,
                          recombination_rate=4e-8, sequence_length=chr_length,
                          discrete_genome=True, num_replicates=nchr)):
        output_file = f"concomitent_div/dr_{seed_dr}/all_chromosomes_seed_{seed}_{i}.vcf"
        mts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=seed)
        for mutation in mts.mutations():
            all_mutations.append(mutation)
        with open(output_file, "w") as vcf_file:
            mts.write_vcf(vcf_file)

# Definition for simulate specific data for some summary statistics 
def simu_concomitent_div_Uhl(Ne_ancpp=100000, Ne_fu=20000, Ne_wp=1000, Ne_sp=500, Ne_sp_found=50, Ne_wp_found=50, Ne_fu_found=200,
                   Ne_bot_sp=50, Ne_bot_wp=50, Ne_bot_fu=50, tsplit_PP=10000, tsplit_L=500, m_spwp_anc=0.00005,
                   m_spwp_rec=0.00001, m_spfu_anc=0.000005, m_wpfu_anc=0.0005, m_spfu_rec=0.000005, m_wpfu_rec=0.0005,
                   n_sp=50, n_wp=50, n_fu=40, nchr=50, chr_length=1e7, seed=100):

    rexpsp_avbot = math.log(Ne_sp / Ne_sp_found) / (tsplit_PP - tsplit_L)
    rexpsp_apbot = math.log(Ne_sp / Ne_bot_sp) / tsplit_L
    rexpwp_avbot = math.log(Ne_wp / Ne_wp_found) / (tsplit_PP - tsplit_L)
    rexpwp_apbot = math.log(Ne_wp / Ne_bot_wp) / tsplit_L
    rexpfu_avbot = math.log(Ne_fu / Ne_fu_found) / (tsplit_PP - tsplit_L)
    rexpfu_apbot = math.log(Ne_fu / Ne_bot_fu) / tsplit_L

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp, growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp, growth_rate=rexpwp_apbot)
    demoIND.add_population(name="FU", initial_size=Ne_fu, growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU", "SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU", "WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP", "WP"], rate=m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP", "FU", "SP"], ancestral="ANCPP")
    demoIND.add_population_parameters_change(time=tsplit_L, population="SP", initial_size=Ne_sp, growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=tsplit_L, population="WP", initial_size=Ne_wp, growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=tsplit_L, population="FU", initial_size=Ne_fu, growth_rate=rexpfu_avbot)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L, populations=["SP", "WP"], rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L, populations=["SP", "FU"], rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_L, populations=["FU", "WP"], rate=m_wpfu_anc)
    demoIND.sort_events()

    # Simulate haplotype data
    nind = {"SP": n_sp, "WP": n_wp, "FU": n_fu}
    all_haplotypes = {"SP": [], "WP": [], "FU": []}  # Initialise les clés pour chaque population

    # Simulate haplotypes for each chromosome
    for i, ts in enumerate(msprime.sim_ancestry(samples=nind, demography=demoIND, random_seed=seed,
                                                recombination_rate=4e-8, sequence_length=chr_length,
                                                discrete_genome=True, num_replicates=nchr)):
        mts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=seed)

        if not mts.variants():
            print(f"No variants found for replicate {i}.")
            continue

        # Haplotypes for each population
        for variant in mts.variants():
            haplotypes = variant.genotypes

            idx_sp = list(range(2 * n_sp))
            idx_wp = list(range(2 * n_sp, 2 * (n_sp + n_wp)))
            idx_fu = list(range(2 * (n_sp + n_wp), 2 * (n_sp + n_wp + n_fu)))

            sp_haplotypes = haplotypes[idx_sp]
            wp_haplotypes = haplotypes[idx_wp]
            fu_haplotypes = haplotypes[idx_fu]

            all_haplotypes["SP"].append(sp_haplotypes)
            all_haplotypes["WP"].append(wp_haplotypes)
            all_haplotypes["FU"].append(fu_haplotypes)

    # If no variants are found, return empty matrices
    if not all_haplotypes["SP"]:
        all_haplotypes["SP"] = np.array([])
    if not all_haplotypes["WP"]:
        all_haplotypes["WP"] = np.array([])
    if not all_haplotypes["FU"]:
        all_haplotypes["FU"] = np.array([])

    # Convert lists into numpy table
    all_haplotypes["SP"] = np.array(all_haplotypes["SP"])
    all_haplotypes["WP"] = np.array(all_haplotypes["WP"])
    all_haplotypes["FU"] = np.array(all_haplotypes["FU"])

    return all_haplotypes

###############################################################################################################################################
# Initial divergence occurs between the WP and FU populations, followed by a divergence between SP and WP. : simu pool of recent_div scenario #
###############################################################################################################################################

def scen_recent_div(Ne_ancpp=100000,Ne_fu_ancfound=200,Ne_wp_ancfound=10000,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100):

    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/t_bot 
    rexpsp_avbot=math.log(Ne_sp/Ne_sp_found)/(tsplit_WP-t_bot)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/t_bot
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_WP-t_bot)
    rexpwp_avsplit=math.log(Ne_wp_anc/Ne_wp_ancfound)/(tsplit_PP-tsplit_WP)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/t_bot
    rexpfu_avbot=math.log(Ne_fu_anc/Ne_fu_ancfound)/(tsplit_PP-t_bot)

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot,initially_active=True)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","FU"], ancestral="ANCPP")
    demoIND.add_population_split(time=tsplit_WP, derived=["SP"], ancestral="WP")
    demoIND.add_population_parameters_change(time=t_bot,population="SP",initial_size=Ne_sp,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="WP",initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="FU",initial_size=Ne_fu_anc,growth_rate=rexpfu_avbot)
    demoIND.add_population_parameters_change(time=tsplit_WP,population="WP",initial_size=Ne_wp_anc,growth_rate=rexpwp_avsplit)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","FU"],rate=m_wpfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","SP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_WP,populations=["WP","FU"],rate=m_ancpp)
    demoIND.sort_events()
    return demoIND

def simu_recent_div(Ne_ancpp=100000,Ne_fu_ancfound=200,Ne_wp_ancfound=10000,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100,n_cpu=6):
    demoIND=recent_div(Ne_ancpp,Ne_fu_ancfound,Ne_wp_ancfound,Ne_sp_found,Ne_wp_found,Ne_fu_anc,Ne_wp_anc,
               Ne_bot_fu,Ne_bot_sp,Ne_bot_wp,Ne_sp,Ne_wp,Ne_fu,tsplit_PP,tsplit_WP,t_bot,m_ancpp,
                m_spwp_anc,m_spwp_rec,m_spfu_anc,m_spfu_rec,m_wpfu_anc,m_wpfu_rec,
                n_sp,n_wp,n_fu,nchr,chr_length,seed)
    nind={"SP":n_sp, "WP":n_wp, "FU":n_fu}
    indv_names=[]
    deb=time.time()
    for k, v in nind.items():
        for i in range(v):
            indv_names.append(k+"_"+str(i+1))

    num_replicates = nchr
    rng = np.random.RandomState(int(seed))
    seeds = rng.randint(1, 2**31, size=(num_replicates, 2))
    geno_matrices = []
    print("start!")
    with concurrent.futures.ProcessPoolExecutor(max_workers=int(n_cpu)) as executor:
        for geno_mat in executor.map(worker, seeds, [nind] * num_replicates, [demoIND] * num_replicates, [chr_length] * num_replicates, [nchr] * num_replicates):
            geno_matrices.append(geno_mat)
    end=time.time()
    print("Time = ",end-deb)
    geno_final = geno_matrices[0]
    for i in range (1,len(geno_matrices)):
        geno_final=np.concatenate((geno_final,geno_matrices[i]),axis=0)
    return geno_final

###############################################################################################################################################
# Initial divergence occurs between the WP and FU populations, followed by a divergence between SP and WP. : simu ind  of recent_div scenario #
###############################################################################################################################################

def simu_recent_div_ind(Ne_ancpp=100000,Ne_fu_ancfound=200,Ne_wp_ancfound=10000,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100,seed_dr=1):

    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/t_bot 
    rexpsp_avbot=math.log(Ne_sp/Ne_sp_found)/(tsplit_WP-t_bot)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/t_bot
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_WP-t_bot)
    rexpwp_avsplit=math.log(Ne_wp_anc/Ne_wp_ancfound)/(tsplit_PP-tsplit_WP)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/t_bot
    rexpfu_avbot=math.log(Ne_fu_anc/Ne_fu_ancfound)/(tsplit_PP-t_bot)

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot,initially_active=True)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","FU"], ancestral="ANCPP")
    demoIND.add_population_split(time=tsplit_WP, derived=["SP"], ancestral="WP")
    demoIND.add_population_parameters_change(time=t_bot,population="SP",initial_size=Ne_sp,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="WP",initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="FU",initial_size=Ne_fu_anc,growth_rate=rexpfu_avbot)
    demoIND.add_population_parameters_change(time=tsplit_WP,population="WP",initial_size=Ne_wp_anc,growth_rate=rexpwp_avsplit)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","FU"],rate=m_wpfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","SP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_WP,populations=["WP","FU"],rate=m_ancpp)
    demoIND.sort_events()

    nind={"SP":n_sp, "WP":n_wp, "FU":n_fu}
    indv_names=[]
    for k, v in nind.items():
        for i in range(v):
            indv_names.append(k+"_"+str(i+1))

    all_mutations = []
    for i,ts in enumerate(msprime.sim_ancestry(samples=nind, demography=demoIND, random_seed=seed,
                          recombination_rate=4e-8, sequence_length=chr_length,
                          discrete_genome=True, num_replicates=nchr)):
        output_file = f"recent_div/dr_{seed_dr}/all_chromosomes_seed_{seed}_{i}.vcf"
        mts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=seed)
        for mutation in mts.mutations():
            all_mutations.append(mutation)
        with open(output_file, "w") as vcf_file:    
            mts.write_vcf(vcf_file)

def simu_recent_div_Uhl(Ne_ancpp=100000,Ne_fu_ancfound=200,Ne_wp_ancfound=10000,Ne_sp_found=50,Ne_wp_found=50,Ne_fu_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100):

    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/t_bot 
    rexpsp_avbot=math.log(Ne_sp/Ne_sp_found)/(tsplit_WP-t_bot)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/t_bot
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_WP-t_bot)
    rexpwp_avsplit=math.log(Ne_wp_anc/Ne_wp_ancfound)/(tsplit_PP-tsplit_WP)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/t_bot
    rexpfu_avbot=math.log(Ne_fu_anc/Ne_fu_ancfound)/(tsplit_PP-t_bot)

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot,initially_active=True)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","FU"], ancestral="ANCPP")
    demoIND.add_population_split(time=tsplit_WP, derived=["SP"], ancestral="WP")
    demoIND.add_population_parameters_change(time=t_bot,population="SP",initial_size=Ne_sp,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="WP",initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="FU",initial_size=Ne_fu_anc,growth_rate=rexpfu_avbot)
    demoIND.add_population_parameters_change(time=tsplit_WP,population="WP",initial_size=Ne_wp_anc,growth_rate=rexpwp_avsplit)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","FU"],rate=m_wpfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","SP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_WP,populations=["WP","FU"],rate=m_ancpp)
    demoIND.sort_events()
    
    nind = {"SP": n_sp, "WP": n_wp, "FU": n_fu}
    all_haplotypes = {"SP": [], "WP": [], "FU": []}  # Initialise les clés pour chaque population

    for i, ts in enumerate(msprime.sim_ancestry(samples=nind, demography=demoIND, random_seed=seed,
                                                recombination_rate=4e-8, sequence_length=chr_length,
                                                discrete_genome=True, num_replicates=nchr)):
        mts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=seed)

        if not mts.variants():
            print(f"No variants found for replicate {i}.")
            continue

        for variant in mts.variants():
            haplotypes = variant.genotypes

            idx_sp = list(range(2 * n_sp))
            idx_wp = list(range(2 * n_sp, 2 * (n_sp + n_wp)))
            idx_fu = list(range(2 * (n_sp + n_wp), 2 * (n_sp + n_wp + n_fu)))

            sp_haplotypes = haplotypes[idx_sp]
            wp_haplotypes = haplotypes[idx_wp]
            fu_haplotypes = haplotypes[idx_fu]

            all_haplotypes["SP"].append(sp_haplotypes)
            all_haplotypes["WP"].append(wp_haplotypes)
            all_haplotypes["FU"].append(fu_haplotypes)

    if not all_haplotypes["SP"]:
        all_haplotypes["SP"] = np.array([])
    if not all_haplotypes["WP"]:
        all_haplotypes["WP"] = np.array([])
    if not all_haplotypes["FU"]:
        all_haplotypes["FU"] = np.array([])

    all_haplotypes["SP"] = np.array(all_haplotypes["SP"])
    all_haplotypes["WP"] = np.array(all_haplotypes["WP"])
    all_haplotypes["FU"] = np.array(all_haplotypes["FU"])

    return all_haplotypes
    
############################################################################################################################################
# Initial divergence occurs between the SP and WP populations, followed by a divergence between WP and FU : simu pool of old_div scenario #
############################################################################################################################################

def old_div(Ne_ancpp=100000,Ne_sp_ancfound=200,Ne_wp_ancfound=10000,Ne_fu_found=50,Ne_wp_found=50,Ne_sp_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100):

    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/t_bot 
    rexpsp_avbot=math.log(Ne_sp_anc/Ne_sp_ancfound)/(tsplit_PP-t_bot)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/t_bot
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_WP-t_bot)
    rexpwp_avsplit=math.log(Ne_wp_anc/Ne_wp_ancfound)/(tsplit_PP-tsplit_WP)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/t_bot
    rexpfu_avbot=math.log(Ne_fu/Ne_fu_found)/(tsplit_WP-t_bot)

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot,initially_active=True)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","SP"], ancestral="ANCPP")
    demoIND.add_population_split(time=tsplit_WP, derived=["FU"], ancestral="WP")
    demoIND.add_population_parameters_change(time=t_bot,population="SP",initial_size=Ne_sp_anc,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="WP",initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="FU",initial_size=Ne_fu,growth_rate=rexpfu_avbot)
    demoIND.add_population_parameters_change(time=tsplit_WP,population="WP",initial_size=Ne_wp_anc,growth_rate=rexpwp_avsplit)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","FU"],rate=m_wpfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","SP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_WP,populations=["WP","SP"],rate=m_ancpp)
    demoIND.sort_events()
    return demoIND

def simu_old_div(Ne_ancpp=100000,Ne_sp_ancfound=200,Ne_wp_ancfound=10000,Ne_fu_found=50,Ne_wp_found=50,Ne_sp_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100,n_cpu=6):
    demoIND=old_div(Ne_ancpp,Ne_fu_found,Ne_wp_ancfound,Ne_sp_ancfound,Ne_wp_found,Ne_sp_anc,Ne_wp_anc,
               Ne_bot_fu,Ne_bot_sp,Ne_bot_wp,Ne_sp,Ne_wp,Ne_fu,tsplit_PP,tsplit_WP,t_bot,m_ancpp,
                m_spwp_anc,m_spwp_rec,m_spfu_anc,m_spfu_rec,m_wpfu_anc,m_wpfu_rec,
                n_sp,n_wp,n_fu,nchr,chr_length,seed)
    nind={"SP":n_sp, "WP":n_wp, "FU":n_fu}
    indv_names=[]
    deb=time.time()
    for k, v in nind.items():
        for i in range(v):
            indv_names.append(k+"_"+str(i+1))

    num_replicates = nchr
    rng = np.random.RandomState(int(seed))
    seeds = rng.randint(1, 2**31, size=(num_replicates, 2))
    geno_matrices = []
    print("let's go")
    with concurrent.futures.ProcessPoolExecutor(max_workers=int(n_cpu)) as executor:
        for geno_mat in executor.map(worker, seeds, [nind] * num_replicates, [demoIND] * num_replicates, [chr_length] * num_replicates, [nchr] * num_replicates):
            geno_matrices.append(geno_mat)
    end=time.time()
    print("Les params = ",end-deb)
    geno_final = geno_matrices[0]
    for i in range (1,len(geno_matrices)):
        geno_final=np.concatenate((geno_final,geno_matrices[i]),axis=0)
    return geno_final

#############################################################################################################################################
# Initial divergence occurs between the SP and WP populations, followed by a divergence between WP and FU : simu ind of recent_div scenario #
#############################################################################################################################################

def simu_old_div_ind(Ne_ancpp=100000,Ne_sp_ancfound=200,Ne_wp_ancfound=10000,Ne_fu_found=50,Ne_wp_found=50,Ne_sp_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100,seed_dr=1):

    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/t_bot 
    rexpsp_avbot=math.log(Ne_sp_anc/Ne_sp_ancfound)/(tsplit_PP-t_bot)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/t_bot
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_WP-t_bot)
    rexpwp_avsplit=math.log(Ne_wp_anc/Ne_wp_ancfound)/(tsplit_PP-tsplit_WP)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/t_bot
    rexpfu_avbot=math.log(Ne_fu/Ne_fu_found)/(tsplit_WP-t_bot)

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot,initially_active=True)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","SP"], ancestral="ANCPP")
    demoIND.add_population_split(time=tsplit_WP, derived=["FU"], ancestral="WP")
    demoIND.add_population_parameters_change(time=t_bot,population="SP",initial_size=Ne_sp_anc,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="WP",initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="FU",initial_size=Ne_fu,growth_rate=rexpfu_avbot)
    demoIND.add_population_parameters_change(time=tsplit_WP,population="WP",initial_size=Ne_wp_anc,growth_rate=rexpwp_avsplit)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","FU"],rate=m_wpfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","SP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_WP,populations=["WP","SP"],rate=m_ancpp)
    demoIND.sort_events()

    nind={"SP":n_sp, "WP":n_wp, "FU":n_fu}
    indv_names=[]
    for k, v in nind.items():
        for i in range(v):
            indv_names.append(k+"_"+str(i+1))

    all_mutations = []
    for i,ts in enumerate(msprime.sim_ancestry(samples=nind, demography=demoIND, random_seed=seed,
                          recombination_rate=4e-8, sequence_length=chr_length,
                          discrete_genome=True, num_replicates=nchr)):
        output_file = f"old_div/dr_{seed_dr}/all_chromosomes_seed_{seed}_{i}.vcf"
        mts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=seed)
        for mutation in mts.mutations():
            all_mutations.append(mutation)
        with open(output_file, "w") as vcf_file:    
            mts.write_vcf(vcf_file)
            
def simu_old_div_Uhl(Ne_ancpp=100000,Ne_sp_ancfound=200,Ne_wp_ancfound=10000,Ne_fu_found=50,Ne_wp_found=50,Ne_sp_anc=20000,Ne_wp_anc=10000,
               Ne_bot_fu=200,Ne_bot_sp=50,Ne_bot_wp=500,Ne_sp=500,Ne_wp=500,Ne_fu=500,tsplit_PP=10000,tsplit_WP=800,t_bot=500,m_ancpp=0.005,
                m_spwp_anc=0.00005,m_spwp_rec=0.000005,m_spfu_anc=0.0005,m_spfu_rec=0.0005,m_wpfu_anc=0.0005,m_wpfu_rec=0.0005,
                n_sp=50,n_wp=50,n_fu=40,nchr=50,chr_length=1e7,seed=100):

    rexpsp_apbot=math.log(Ne_sp/Ne_bot_sp)/t_bot 
    rexpsp_avbot=math.log(Ne_sp_anc/Ne_sp_ancfound)/(tsplit_PP-t_bot)
    rexpwp_apbot=math.log(Ne_wp/Ne_bot_wp)/t_bot
    rexpwp_avbot=math.log(Ne_wp/Ne_wp_found)/(tsplit_WP-t_bot)
    rexpwp_avsplit=math.log(Ne_wp_anc/Ne_wp_ancfound)/(tsplit_PP-tsplit_WP)
    rexpfu_apbot=math.log(Ne_fu/Ne_bot_fu)/t_bot
    rexpfu_avbot=math.log(Ne_fu/Ne_fu_found)/(tsplit_WP-t_bot)

    demoIND = msprime.Demography()
    demoIND.add_population(name="SP", initial_size=Ne_sp,growth_rate=rexpsp_apbot)
    demoIND.add_population(name="WP", initial_size=Ne_wp,growth_rate=rexpwp_apbot,initially_active=True)
    demoIND.add_population(name="FU", initial_size=Ne_fu,growth_rate=rexpfu_apbot)
    demoIND.add_population(name="ANCPP", initial_size=Ne_ancpp)

    demoIND.set_symmetric_migration_rate(["FU","SP"], m_spfu_rec)
    demoIND.set_symmetric_migration_rate(["FU","WP"], m_wpfu_rec)
    demoIND.set_symmetric_migration_rate(["SP","WP"], m_spwp_rec)

    demoIND.add_population_split(time=tsplit_PP, derived=["WP","SP"], ancestral="ANCPP")
    demoIND.add_population_split(time=tsplit_WP, derived=["FU"], ancestral="WP")
    demoIND.add_population_parameters_change(time=t_bot,population="SP",initial_size=Ne_sp_anc,growth_rate=rexpsp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="WP",initial_size=Ne_wp,growth_rate=rexpwp_avbot)
    demoIND.add_population_parameters_change(time=t_bot,population="FU",initial_size=Ne_fu,growth_rate=rexpfu_avbot)
    demoIND.add_population_parameters_change(time=tsplit_WP,population="WP",initial_size=Ne_wp_anc,growth_rate=rexpwp_avsplit)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","FU"],rate=m_wpfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["WP","SP"],rate=m_spwp_anc)
    demoIND.add_symmetric_migration_rate_change(time=t_bot,populations=["SP","FU"],rate=m_spfu_anc)
    demoIND.add_symmetric_migration_rate_change(time=tsplit_WP,populations=["WP","SP"],rate=m_ancpp)
    demoIND.sort_events()

    nind = {"SP": n_sp, "WP": n_wp, "FU": n_fu}
    all_haplotypes = {"SP": [], "WP": [], "FU": []}  # Initialise les clés pour chaque population

    for i, ts in enumerate(msprime.sim_ancestry(samples=nind, demography=demoIND, random_seed=seed,
                                                recombination_rate=4e-8, sequence_length=chr_length,
                                                discrete_genome=True, num_replicates=nchr)):
        mts = msprime.sim_mutations(ts, rate=2.9e-9, random_seed=seed)

        if not mts.variants():
            print(f"No variants found for replicate {i}.")
            continue

        for variant in mts.variants():
            haplotypes = variant.genotypes

            idx_sp = list(range(2 * n_sp))
            idx_wp = list(range(2 * n_sp, 2 * (n_sp + n_wp)))
            idx_fu = list(range(2 * (n_sp + n_wp), 2 * (n_sp + n_wp + n_fu)))

            sp_haplotypes = haplotypes[idx_sp]
            wp_haplotypes = haplotypes[idx_wp]
            fu_haplotypes = haplotypes[idx_fu]

            # Ajouter les haplotypes pour chaque population
            all_haplotypes["SP"].append(sp_haplotypes)
            all_haplotypes["WP"].append(wp_haplotypes)
            all_haplotypes["FU"].append(fu_haplotypes)

    if not all_haplotypes["SP"]:
        all_haplotypes["SP"] = np.array([])
    if not all_haplotypes["WP"]:
        all_haplotypes["WP"] = np.array([])
    if not all_haplotypes["FU"]:
        all_haplotypes["FU"] = np.array([])

    all_haplotypes["SP"] = np.array(all_haplotypes["SP"])
    all_haplotypes["WP"] = np.array(all_haplotypes["WP"])
    all_haplotypes["FU"] = np.array(all_haplotypes["FU"])

    return all_haplotypes
