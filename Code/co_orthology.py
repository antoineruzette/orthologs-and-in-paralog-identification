# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 18:30:42 2022

@author: Antoine Ruzette
Comparative Genomics, MSc in Bioinformatics, KU Leuven, 2021-2022
"""
import pandas as pd
import re

##DATA FRAME INITIALIZATION
early_match = []
toremove3 = []
toremove4 = []

def initializationDataFrame():
        query_gene = []
        alignment_gene = []
        score = [] 

        return query_gene, alignment_gene, score
    
##EXTRACT QUERY AND ALIGNES GENE NAMES AND SCORES FROM THE BLAST FILES
def extractGenes(filename, dict_align): 
    
    with open(filename, 'r') as align1: 
        for ln in align1: 
            #get the query gene name
            if ln.startswith("Query="):
                ID_query = re.search('[XN]P_\d{5,}\.\d{1}', ln).group()
                (dict_align[list(dict_align.keys())[0]]).append(ID_query)#7:21 corresponds to the query gene ID
            #get the best target gene name & score
                test = align1.readline()
            
                while (not re.search('[XN]P_\d{5,}\.\d{1}', test)):
                    test = align1.readline()#continue to read
                ID_target = re.search('[XN]P_\d{5,}\.\d{1}', test).group()
                (dict_align[list(dict_align.keys())[1]]).append(ID_target)
                (dict_align[list(dict_align.keys())[2]]).append(test[70:74])
    return dict_align; 


##DEFINE ORTHOLOGY
def orthology(dict_align1, dict_align2): 
    query_gene2 = (dict_align2[list(dict_align2.keys())[0]])
    query_gene1 = (dict_align1[list(dict_align1.keys())[0]])
    alignment_gene2 = (dict_align2[list(dict_align2.keys())[1]])
    alignment_gene1 = (dict_align1[list(dict_align1.keys())[1]])
    score2 = (dict_align2[list(dict_align2.keys())[2]])
    score1 = (dict_align1[list(dict_align1.keys())[2]])
            
    for j in range(len(query_gene1)) :
        for i in range(len(query_gene2)) :
            
            if ( query_gene2[i] == alignment_gene1[j] and alignment_gene2[i] == query_gene1[j] ) :
                early_match.append([query_gene2[i], score2[i], query_gene1[j], score1[j]])#Order: Danio, Latimeria
                #query_gene2 -> Danio
                #alignment_gene2 -> Latimeria
                #query_gene1 -> Latimeria
                #alignment_gene1 -> Danio 
                
def in_paralogy(filename, dict_align, toremove): 

    with open(filename,"r") as align_self:
        for ln in align_self:
            #get the query gene name
            if ln.startswith("Query="):
                ID_query = re.search('[XN]P_\d{5,}\.\d{1}', ln).group()
                (dict_align[list(dict_align.keys())[0]]).append(ID_query)
            #get the best alignment gene name & score 
                test = align_self.readline()
                
                while (not re.search('[XN]P_\d{5,}\.\d{1}', test)) :
                    test = align_self.readline()#continue to read
                    
                tmp = align_self.readline()#goes to second best alignment
                if (re.search('[XN]P_\d{5,}\.\d{1}', tmp)) :
                    ID_target = re.search('[XN]P_\d{5,}\.\d{1}', tmp).group()
                    (dict_align[list(dict_align.keys())[1]]).append(ID_target)
                    (dict_align[list(dict_align.keys())[2]]).append(tmp[70:74])
                    for z in range(len(early_match)) :
                        if (dict_align == dict_align3): #for latimeria_self
                            if (early_match[z][2] == ID_query and tmp[70:74] > early_match[z][3]) :#for query_gene1 = latimeria
                                toremove.append(ID_query)
                        else:                           #for danio_self
                            if (early_match[z][0] == ID_query and tmp[70:74] > early_match[z][1]) :#for query_gene2 = danio
                                toremove.append(ID_query)
    return toremove; 


def co_orthology(): 
    #TO BE DEFINED, as a pair of in-paralogous sequences that are both orthologous to 
    return

#data frame initialization
for i in 1,2,3,4: 
    globals() ['dict_align' + str(i)] = {('query_gene' + str(i)) : initializationDataFrame()[0], 
                                         ('alignment_gene' + str(i)) : initializationDataFrame()[1], 
                                         ('score' + str(i)) : initializationDataFrame()[2]}
    i += 1

#query and aligned genes extraction

#extractGenes("alignment_test.txt", dict(dict_align1))
#extractGenes("alignment_test.txt", dict(dict_align2))

extractGenes("danio_latimeria-2-alignment.txt", dict(dict_align1))#DB: danio
extractGenes("latimera_rerio-2-alignment.txt", dict(dict_align2))#DB: latimera


#compute orthologous genes according to BBH
orthology(dict_align1, dict_align2)

#save results of orthology to file
column_names = ["Gene1 (Latimeria)", "Score1", "Gene2 (Danio)", "Score2"]
full_ortho = pd.DataFrame(columns = column_names)
for i in range(len(early_match)):
    full_ortho = full_ortho.append({"Gene1 (Latimeria)": early_match[i][2], "Score1": early_match[i][3], "Gene2 (Danio)":early_match[i][0], "Score2":early_match[i][1]}, ignore_index=True)
#send it to a text file
full_ortho.to_csv(r'full_ortho.txt', header = True, sep = '\t')


in_paralogy("latimeria_self-alignment.txt", dict(dict_align3), toremove3)#DB and query: latimeria
in_paralogy("danio_self-alignment.txt", dict(dict_align4), toremove4)#DB and query: danio
#in_paralogy("alignment_test.txt", dict(dict_align3), toremove3)#DB and query: latimeria


#save results of in-paralogy to file
column_paralogs = ["Paralogs"]
paralogs = pd.DataFrame(columns = column_paralogs)
for i in range(len(toremove3)):
    paralogs = paralogs.append({"Paralogs": toremove3[i]}, ignore_index=True)
for i in range(len(toremove4)):
    paralogs = paralogs.append({"Paralogs": toremove4[i]}, ignore_index=True)
paralogs.to_csv(r'paralogs.txt', header = True, sep = '\t')
