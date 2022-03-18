import pandas as pd
import numpy as np
from collections import Counter

colnames = ["gene"]
all_genes =  pd.read_csv("/home/sgok603/sreemol_personal_eRF/project_three/nafld_R_project/data/gene_reference/gencode_gene_list.txt", header = None, names = colnames, sep = "\t")

liver_genes = pd.read_csv("/home/sgok603/sreemol_personal_eRF/project_three/nafld_R_project/data/gene_reference/codes3d_liver_genes.txt", header = 0, sep = "\t")

def iteration(all_genes, liver_genes):

    rand = all_genes.sample(n=31)
    common = pd.merge(rand, liver_genes, how='inner', on=['gene'])
    
    return(common)

def repeat(times, f, all_genes, liver_genes):
    
    emp = []
    for i in range(times):
        common = iteration(all_genes, liver_genes)
        emp.append(len(common))
    
    return(emp)

num = 10000 
random_success = repeat(num, iteration, all_genes, liver_genes)

hist_data = pd.DataFrame({'overlaps':random_success})

hist_data.to_csv("/home/sgok603/sreemol_personal_eRF/project_three/nafld_R_project/results/risk_gene_analysis_results/stat_test_results/random_gene_occurance_values_codes3d.txt", header = True, index = False)

counts = dict((i, random_success.count(i)) for i in random_success)

counts_df = pd.DataFrame(counts.items())

counts_df.columns = ['rand_intersect', 'count']

counts_df.to_csv("/home/sgok603/sreemol_personal_eRF/project_three/nafld_R_project/results/risk_gene_analysis_results/stat_test_results/test_for_random_gene_occurance_in_codes3d_liver_res.txt", sep = "\t", header = True, index= False)

prob_3 = []

[prob_3.append(x) for x in random_success if x > 13]

print(prob_3)

pval_prob_3 = float(len(prob_3)/float(num))

print(pval_prob_3)

