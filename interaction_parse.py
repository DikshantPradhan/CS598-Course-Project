#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Load required modules
import numpy as np
#from sklearn import tree
#from sklearn.metrics import accuracy_score
import pandas as pd


# In[5]:


def import_data(filename):
    data = pd.read_csv(filename, sep="\t", header=0)
    return(data)

directory = '/home/dikshantpradhan/Documents/'
data = import_data(directory + 'Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExE.txt')
single_mutant_data = import_data(directory + 'Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/strain_ids_and_single_mutant_fitness.csv')
#single_mutant_data
GO_data = import_data('/home/dikshantpradhan/GitHub/CS598-Course-Project/GOTermGeneListForClassificationYeast.txt')

ExE_data = import_data(directory + 'Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExE.txt')
NxN_data = import_data(directory + 'Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_NxN.txt')
NxE_data = import_data(directory + 'Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_ExN_NxE.txt')


# In[6]:


#data


# In[7]:


#print(ExE_data.shape)
#print(NxN_data.shape)
#print(NxE_data.shape)
full_data = pd.concat([ExE_data, NxN_data, NxE_data])
full_data.index = range(full_data.shape[0])
#print(full_data.shape[0])
#full_data


# In[8]:


#full_data['Double mutant fitness'][list([2,3,4])]


# In[17]:


def parse_fitness_helper(i, j, ij, pred_tol = 0.1, class_tol = 0.1):
    pred_ij = i*j
    pred_match = np.isclose(pred_ij, ij, pred_tol)
    if (pred_match):
        return('independent')
    i_j_match = np.isclose(i,j,class_tol)
    i_ij_match = np.isclose(i,ij,class_tol)
    j_ij_match = np.isclose(j,ij,class_tol)

    if (i_j_match):
        if (j_ij_match and i_ij_match):
            return('coequal')
        elif (j < ij):
            return('synergistic')
    if (i < j):
        if (j <= ij):
            return('masking')
        elif ((i <= ij) and (ij <= j)):
            return('suppressive')
    if ((ij < i) and (ij < j)):
        return('antagonistic')

    #print(i); print(j); print(ij)
    return('unknown')

def parse_fitness(i, j, ij, pred_tol = 0.1, class_tol = 0.1):
    classification = parse_fitness_helper(i, j, ij, pred_tol, class_tol)
    if (classification == 'unknown'):
        classification = parse_fitness_helper(j, i, ij, pred_tol, class_tol)
    return(classification)


# In[10]:


# build matrix
from itertools import chain

def search_list(mylist, search):
    idx = [i for i, s in enumerate(mylist) if search in s]
    return(idx)

def search_gene(mylist, gene):
    cats = ['_', '-']
    idxs = []
    for cat in cats:
        search = gene + cat
        indices = search_list(mylist, gene) #[i for i, s in enumerate(mylist) if search in s]
        idxs.append(indices)

    idxs = list(chain.from_iterable(idxs))
    return(idxs)

def append_list(mylist, add):
    mylist.append(add)
    mylist = list(chain.from_iterable(mylist))
    return(mylist)

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def extract_single_mutant_fitness(df, gene, gene_column = 'Allele/Gene name',
                                  fitness_column = 'Single mutant fitness (26°)'):
    gene = gene.lower()
    idxs = search_gene(df[gene_column], gene)
    mean = np.mean(df[fitness_column][idxs])
    return(mean)

#print(extract_single_mutant_fitness(single_mutant_data, 'rpc82'))

def extract_double_mutant_fitness(df, query_gene, array_gene,
                                 query_column = 'Query allele name', array_column = 'Array allele name',
                                 fitness_column = 'Double mutant fitness'):
    q_idxs = search_gene(df[query_column], query_gene)
    a_idxs = search_gene(df[array_column], array_gene)

    dbl_mutant_idxs = intersection(q_idxs, a_idxs)
    dbl_mutant_idxs = list(np.unique(dbl_mutant_idxs))
    #print(dbl_mutant_idxs)
    #print(df[fitness_column][dbl_mutant_idxs])
    mean = np.mean(df[fitness_column][dbl_mutant_idxs])
    return(mean)

#extract_double_mutant_fitness(full_data,'tfc3','stu1')


# In[11]:


reference_genes = ['TFC3', 'ARH1']
query_genes = ['MCM2', 'LSM2', 'STU1']

#reference_genes = reference_genes.lower()
#query_genes = query_genes.lower()


def build_interaction_mtx(reference_genes, query_genes, double_mutant_df, single_mutant_df):
    m = len(reference_genes)
    n = len(query_genes)
    mtx = np.zeros([n,m])

    for i in range(0,m):
        gene_i = reference_genes[i].lower()
        for j in range(0,n):
            gene_j = query_genes[j].lower()
            fitness = 0
            fitness = extract_double_mutant_fitness(double_mutant_df, gene_i, gene_j)
            mtx[j,i] = fitness
    return(mtx)

#dbl_mutant_fit = build_interaction_mtx(reference_genes, query_genes, full_data, single_mutant_data)
#print(dbl_mutant_fit)
#print(dbl_mutant_fit[1][0])


# In[12]:


def build_fitness_mtx(reference_genes, query_genes, double_mutant_df, single_mutant_df):
    m = len(reference_genes)
    n = len(query_genes)
    interaction = pd.DataFrame(index=query_genes, columns=reference_genes)

    double_mutant_fitness = build_interaction_mtx(reference_genes, query_genes, double_mutant_df, single_mutant_df)

    for i in range(0,m):
        gene_i = reference_genes[i]
        reference_fitness = extract_single_mutant_fitness(single_mutant_df, gene_i)
        for j in range(0,n):
            gene_j = query_genes[j]
            query_fitness = extract_single_mutant_fitness(single_mutant_df, gene_j)
            dbl_fitness = double_mutant_fitness[j][i]
            #print(query_fitness); print(dbl_fitness)
            interact_label = parse_fitness(reference_fitness, query_fitness, dbl_fitness)
            interaction[gene_i][gene_j] = interact_label


    return(interaction)
    # fitness = extract_single_mutant_fitness(single_mutant_df, gene_i)

#mtx = build_fitness_mtx(reference_genes, query_genes, full_data, single_mutant_data)
#mtx
#mtx['MCM2']['TFC3']


# In[13]:


def diff(first, second):
        second = set(second)
        return [item for item in first if item not in second]

def sample_categories(GO_data, nsamples = 5):
    categories = np.unique(GO_data.GOTermDescription)
    reference_genes = []
    query_genes = []
    for category in categories:
        idxs = np.where(GO_data.GOTermDescription == category)
        genes = GO_data.GeneID[idxs[0]]
        ref = np.random.choice(genes, size=nsamples, replace=False)
        for gene in ref:
            reference_genes.append(gene)
        #print(ref)
        query = diff(genes, ref)
        for gene in query:
            query_genes.append(gene)
        #reference_genes = reference_genes.append(ref)
        #query_genes = append_list(query_genes, query)
    reference_genes = np.unique(reference_genes)
    query_genes = np.unique(query_genes)
    return([reference_genes, query_genes])

genes = sample_categories(GO_data)
#genes[0]


# In[16]:


def interaction_dataframe(reference_genes, query_genes, double_mutant_data, single_mutant_data):
    df = pd.DataFrame(columns=reference_genes, index=query_genes)
    mtx = build_fitness_mtx(reference_genes, query_genes, double_mutant_data, single_mutant_data)

    for query in query_genes:
        for ref in reference_genes:
            df.loc[query][ref] = mtx.loc[query][ref]

    return(df)

#df = interaction_dataframe(genes[0][range(5)], genes[1][range(20)], full_data, single_mutant_data)
#print(df)


# In[20]:


#mtx.to_csv('test.csv')
all_genes = GO_data.GeneID
df = interaction_dataframe(all_genes, all_genes, full_data, single_mutant_data)
df.to_csv('full_interaction_data.csv')


# In[ ]:
