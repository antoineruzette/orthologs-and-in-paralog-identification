# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 21:52:58 2022

@author: Antoine Ruzette
Comparative Genomics, MSc in Bioinformatics, KU Leuven, 2021-2022
"""

import re
import pandas as pd
import numpy as np
import math
from collections import Counter
import matplotlib.pyplot as plt


def split(word): #split a string into its characters
    return [char for char in word]  


def shannon(msa_list): #compute the shannon entropy of a list
    score=0
    count = Counter(msa_list)
    freq = []
    
    for key in count:
        p = count[key]
        p = p/len(msa_list)
        freq.append(p)
        
    for m in range(len(freq)):
        new_score = freq[m]*math.log2(freq[m])
        shannon = (score + new_score)*-1
    
    return shannon


msa = []
msa_tmp = []

#retrieve IDs and sequences from the .clw file
file_in = open("Original_MSA_Muscle.clw")
for line in file_in:        
    if line.startswith('CLUSTAL'): 
        continue #skip the header of the .clw file
    elif re.match(r'\s', line):
           continue #skip empty lines
    else:
        msa.append(line[0:31])#IDs
        msa_tmp.append(line[37:98])#protein sequences
            
msa_list = list(dict.fromkeys(msa))#get the IDs, remove duplicates

#merge the sequences from different lines
master = []
for i in range(len(msa_list)): 
    master.append(''.join([msa_tmp[i], msa_tmp[i+len(msa_list)]]))
    

output = []
for grp in master: #split each sequence into its residues and store it in a nested list
    output.append(split(''.join(grp)))
msa_table = pd.DataFrame(data=output, index= msa_list)


shannon_list = []
for i in range(len(output[0])):
    shannon_tmp = shannon(msa_table[i])
    shannon_list.append(shannon_tmp)
    
    
#plot 1D list of shannon entropy of each site
x = np.array(range(0, len(shannon_list)))
x_tick = np.arange(0, len(shannon_list), 10)
y = shannon_list

plt.scatter(x, y, color = "darkblue", label = "residue's shannon entropy", s = 8, alpha = 1)
plt.xlabel("Residue position in the MSA")
plt.ylabel("Shannon entropy")
plt.title("Shannon entropy of each position in the MSA")


