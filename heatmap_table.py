#05/11/2018
# The meaning of this script is to make a heatmap table of jaccard index of shared protein families between  different organisms
# It takes two_files: OrthoMCL output file and a taxonomy file in which taxonomy name is separated with tabulator from protein identifier (id'\t'name) 
# Run command heatmap_table.py taxonomy_file OrthoMCLfile

import sys
import matplotlib.pyplot as plt
import matplotlib.cm, matplotlib.colors
import numpy as np
import math
from collections import defaultdict


#reads taxonomy file and returns dictionary. Protein sequence Identifier is the key and value is taxonomy name

def get_taxonomy():
    taxonomy = sys.argv[1]
    taxonomy_dictionary = {}
    f = open(taxonomy, 'r')
    for a in f:
        a_split = a.rstrip().split('\t')
        taxonomy_dictionary[a_split[0]] = a_split[1]
    f.close()
    return taxonomy_dictionary



#Identifies all the taxonomic pairs in the data set (all-against-all) and creating dictionary, takes a taxonomy_list from 
#get_taxonomy function and returns a dictionary containing all the taxonomical pairs

def make_taxonomy_pairs(taxonomy_dictionary):

    taxonomy_pairs_dictionary = {}
    taxonomy_list = list(set(taxonomy_dictionary.values()))

    for tax in range(len(taxonomy_list)):
        taxonomy_pairs_dictionary[taxonomy_list[tax]] = 0
        for tax2 in range(tax, len(taxonomy_list)):
            if taxonomy_list[tax]>taxonomy_list[tax2]:   #sorting of the taxonomies for key labels
                key_to_taxdict = taxonomy_list[tax], taxonomy_list[tax2]
                taxonomy_pairs_dictionary[key_to_taxdict] = 0
            elif taxonomy_list[tax]<taxonomy_list[tax2]:
                key_to_taxdict = taxonomy_list[tax2], taxonomy_list[tax]
                taxonomy_pairs_dictionary[key_to_taxdict] = 0

    return taxonomy_pairs_dictionary

# Reads clustering file and counts shared protein families between taxonomies

def read_clustering(taxonomy_dictionary):

    clustering_file = sys.argv[2]
    f = open(clustering_file, 'r')

    cluster_tax_set = defaultdict(set) # Collects all the protein ids into sets of protein families
    protein_family = 1
    for prot_fam_line in f:

        prot_fam = prot_fam_line.rstrip().split('\t')
        # The identifiers in the clustering file are genome_id|protein_id, only protein_id is needed. Each line in clustering file forms a protein family 
        for tax_gi in prot_fam:
            gi_id = tax_gi.split('|')[1]
            if gi_id in taxonomy_dictionary:
                taxonomy_for_gi = taxonomy_dictionary.get(gi_id)
                cluster_tax_set[protein_family].add(taxonomy_for_gi)
        protein_family +=1
    return cluster_tax_set


# Calculating all the protein families in one genome

def count_protein_families_per_organism(cluster_tax_set):

    total_number_of_protein_families_in_gen = {}

    for aa in cluster_tax_set:
        for a in cluster_tax_set[aa]:
            if a in total_number_of_protein_families_in_gen:
                total_number_of_protein_families_in_gen[a] += 1
            elif a not in total_number_of_protein_families_in_gen:
                total_number_of_protein_families_in_gen[a] = 1
    return total_number_of_protein_families_in_gen

#Calculates protein families shared between two genomes. Takes all_against_all pair dictionary, taxonomy dictionary as well as clustering information

def count_shared_protein_families_between_two_organism(taxonomy_dictionary, taxonomy_pairs_dict, cluster_tax_set):

    #tax_clust_list_dict.setdefault(taxonomy_for_gi,[]).append(number)

    for aa in cluster_tax_set:
        intermediate_list = sorted(list(cluster_tax_set[aa]))

        for a in range(len(intermediate_list)-1):
            for b in range(a, len(intermediate_list)):
                if intermediate_list[a] > intermediate_list[b]:
                    key = intermediate_list[a], intermediate_list[b]
                    taxonomy_pairs_dict[key] += 1
                elif intermediate_list[b] > intermediate_list[a]:
                    key = intermediate_list[b], intermediate_list[a]
                    taxonomy_pairs_dict[key] += 1

    return taxonomy_pairs_dict


# Calculating jaccard_index

def Jaccard_index(taxonomy_pairs_dict, taxonomy_dictionary, total_number_of_protein_families_in_gen):

    jaccard_index_dictionary = {}
    name_list = set(taxonomy_dictionary.values())
    for d in name_list:
        #print (d)
        for e in name_list:
            if d > e:
                key = d, e
            elif e>d:
                key = e, d
            else:
                key = e, d

            total_number_of_protein_families_in_1_organism = float(total_number_of_protein_families_in_gen.get(d))
            total_number_of_protein_families_in_2_organism = float(total_number_of_protein_families_in_gen.get(e))
            if d !=e:
                shared_protein_families = float(taxonomy_pairs_dict.get(key))
            else:
                shared_protein_families = float(total_number_of_protein_families_in_gen.get(d))
            union_of_protein_families = (total_number_of_protein_families_in_1_organism + total_number_of_protein_families_in_2_organism) - shared_protein_families
            jaccard_index = shared_protein_families / union_of_protein_families
            jaccard_index_dictionary[key] = jaccard_index
    return jaccard_index_dictionary
#Making heatmap, takes jaccard index dictionary and taxonomy dictionary. Function opens a heatmap on a screen

def make_heatmap(jaccard_index_dictionary, taxonomy_dictionary):

    labels = sorted(list(set(taxonomy_dictionary.values())))  #labels for the heatmap
    table = []

    # Creates a data matrix
    for d in labels:
        line_list = []
        for e in labels:
            if d<e:
                key = e, d
                jaccard_number = jaccard_index_dictionary.get(key)
            elif e<d:
                key = d, e
                jaccard_number = jaccard_index_dictionary.get(key)

            else:
                key = e, d
                jaccard_number = jaccard_index_dictionary.get(key)

            line_list.append(jaccard_number)

        table.append(line_list)
    data = np.array(table)  # Creating data array

    #Parameters for heatmap
    fig, ax = plt.subplots()
    draw = ax.matshow(data, norm=matplotlib.colors.LogNorm(), aspect = 'auto')
    ax.set(xticks=np.arange(len(labels)), xticklabels=labels, yticks=np.arange(len(labels)), yticklabels=labels)
    ax.yaxis.set_ticks_position('right')
    ax.xaxis.set_ticks_position('top')
    ax.set_xticklabels(labels, rotation = 90,  minor=False)
    ax.set_yticklabels(labels, minor=False)
    ax.yaxis.set_label_position('right')
    fig.colorbar(draw, pad = 0.3)
    plt.show()



# Main


tax = get_taxonomy()
pairs =  make_taxonomy_pairs(tax)
clustering = read_clustering(tax)
total_count  =   count_protein_families_per_organism(clustering)
shared_family_counts = count_shared_protein_families_between_two_organism(tax, pairs, clustering)
jaccard_dict =  Jaccard_index(shared_family_counts, tax, total_count)
make_heatmap(jaccard_dict, tax)

