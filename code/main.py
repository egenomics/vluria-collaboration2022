#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:26:41 2022

@author: jlvillanueva
"""
#module load Python/2.7.11
import os, subprocess
import pandas as pd

main_f="/home/jlvillanueva/Documents/analysis/vluria-collaboration2022"
#Function that returns the loop length at the cterminal, nterminal and the lengths of the middle regions of Loops as a list
def calc_loop_lengths(secondary_structure_string, letter):
     #Remove starting and ending L's
     nterminal=0
     cterminal=0
     for position in secondary_structure_string:
         if position=='L':
             cterminal+=1
         else:
             break
     for position in secondary_structure_string[::-1]:
         if position=='L':
             nterminal+=1
         else:
             break
     middle=secondary_structure_string.strip('L')
     distances=[]
     i=0
     for position in middle:
         if position == letter:
             i+=1
         elif i!=0:
             distances.append(i)
             i=0
     return cterminal, nterminal, distances

def calc_loop_lengths_no_strip(secondary_structure_string, letter):
    i=0
    distances=[]
    if secondary_structure_string.count(letter)==0:
        return [0]
    else:
        if len(secondary_structure_string)==1:
            distances.append(secondary_structure_string.count(letter))
        else:
            for position in secondary_structure_string:
                if position == letter:
                    i+=1
                elif i!=0:
                    distances.append(i)
                    i=0
            if i!=0:
                distances.append(secondary_structure_string.count(letter))
    return distances


#Calculates the inter-pfam-domain distances. It returns a list of lists. Each value tells the number of the largest stretch of L.
def inter_domains_loop_length(struct_list):
    start_domain=False
    n=0
    inter_lengths=[]
    inter_seqs=[]
    inter_seq=''
    for position in struct_list:
        if position=='X':
            start_domain=True
            if n!=0:
                inter_seqs.append(inter_seq)
                inter_lengths.append(max(calc_loop_lengths_no_strip(inter_seq, 'L')))
                inter_seq=''
                n=0
        else:
            if start_domain==True:
                n+=1
                inter_seq+=position
    return inter_lengths


#Reads the hmmsearch tbl format file and returns a pandas dataframe
def load_domtable(domtbl):
     all_lines = []
     names = ['target', 't_accession','tlen', 'queryname','q_accession',
          'qlen', 'full_e-val','full_score','full_bias','dom_number', 
          'dom_of','dom_c-Evalue','dom_i-Evalue','dom_score','dom_bias',
          'hmmcoord_from','hmmcoord_to','alncoord_from','alncoord_to','envcoord_from', 
          'envcoord_to','acc','target_decsription']
     with open(domtbl, 'r') as f:
         lines = (line for idx,line in enumerate(f.xreadlines()) if idx > 2 )
         for line in lines:
             cols = [x for x in line.strip().split(' ') if x]
             # if len(cols) > 23, we have extra description columns
             # combine them all into one string in the 19th column
             if len(cols) > 23:
                 cols[22] = ' '.join(cols[22:])
                 del(cols[23:])
                 #assert len(cols) == 23
             elif len(cols) < 23:
                 cols.append('')
                 #assert len(cols) == 23
             all_lines.append(cols)
         return pd.DataFrame(all_lines, columns=names)

def create_word(structure, pfam): 
    '''
    This function divides a structure if it has no pfam domains and has stretches longer than 30 loops.
    Afterwards for each substructure it returns words in proteins. Only takes into account L's if the stretch is longer thant 5.
    '''
    words=[]
    structures=[]
    start_ends_loops=[]
    start=0
    end=0
    loop_length=0
    if pfam==False:
        for index, position in enumerate(list(structure)):
            if position=='L':
                if loop_length==0:
                    start=index
                loop_length+=1
            else:
                end=index
                if loop_length>=30:
                    start_ends_loops.append([start, end, loop_length])
                loop_length=0
        start=0
        for coord in start_ends_loops:
            struc=structure[start:coord[0]]
            start=coord[1]
            if struc!='':
                structures.append(struc)
        structures.append(structure[start:])
    else:
        structures.append(structure)
    for struct in structures:
        word=''
        L=0
        letter=''
        for position in struct:
            if letter=='':
                letter=position
            if position=='L':
                L+=1
            if position!=letter:
                if letter=='L' and L>5:
                    word+=letter
                elif letter=='H' or letter=='E':
                    word+=letter
                letter=''
                L=0
        if word =='':
            if letter=='L' and L>5:
                word+=letter
            elif letter=='H' or letter=='E':
                word+=letter
            letter=''
            L=0
        elif word[-1]!=position:
            if position=='L' and L>5:
                word+=position
            elif position=='H' or position=='E':
                word+=position
            letter=''
            L=0
        words.append(word)
    return [words, structures]

#run xlsparser
from subprocess import Popen, PIPE, STDOUT
a=open(main_f+"/logs/errors_parser.txt", 'w') # log error for parser
b=open(main_f+"/logs/errors_prots_failed.txt", 'w') # log error for prots

mypath=main_f+"/Hs_PP_XMLs"
folders = [f for f in os.listdir(mypath)]
for folder in folders:
    proteins = [f for f in os.listdir(mypath+'/'+folder+'/XMLs_'+folder)]
    for prot in proteins:
        if 'processed' not in prot: #avoid the ones already processed
            protein_file = mypath+'/'+folder+'/XMLs_'+folder+'/'+prot
            cmd= 'xsltproc /home/jlvillanueva/Documents/analysis/vluria-collaboration2022/code/xml_pp_parser.xsl ' + protein_file +' > ' + protein_file + '_processed'
            p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            output = p.stdout.read()
            if output!='':
                print>>a, output
                print>>b, protein_file
a.close()
b.close()

#run xlsparser
word_parts=['E','H','L']
structure_dictionary={}

mypath=main_f+"/Hs_PP_XMLs"
folders = [f for f in os.listdir(mypath)]
for folder in folders:
    proteins = [f for f in os.listdir(mypath+'/'+folder+'/XMLs_'+folder)]
    for prot in proteins:
        if 'processed' not in prot: #avoid the ones already processed
            protein_file = mypath+'/'+folder+'/XMLs_'+folder+'/'+prot
            a=open(protein_file+'_processed', 'r')
            word=''
            structure=''
            for line in a:
                info=line.rstrip().split('\t')
                if len(info)>2:
                    if info[1]=='helix' or info[1]=='strand':
                        #print info
                        word+=info[1][0].upper() #We generate the word summarizing the secondary structures.
                if line[0] in word_parts:
                    structure+=line.rstrip()
            structure_dictionary[prot.split('.xml')[0]]=[structure, folder]
            print word
       # elif 'processed' in prot:
            #protein_file = mypath+'/'+folder+'/XMLs_'+folder+'/'+prot
            #os.system('rm '+protein_file)

#Another dictionary for basic statistics about domains across age groups-
summary_domains={}
f0=open(main_f+"/results/summary_domains.txt", 'w')
print>>f0,'protein\tage\tnum_domains\ttotal_domain_length\tprotein_length'
for prot in structure_dictionary:
    summary_domains[prot]={}
    summary_domains[prot]['protein_length']=len(structure_dictionary[prot][0])
    summary_domains[prot]['age']=structure_dictionary[prot][1]
    summary_domains[prot]['num_domains']=0
    summary_domains[prot]['total_domain_length']=0
    summary_domains[prot]['struct']=structure_dictionary[prot][0]
    summary_domains[prot]['struct_list']=list(structure_dictionary[prot][0]) #We will change positions with domains for an X and count them.
    summary_domains[prot]['struct_list_original']=list(structure_dictionary[prot][0])
    summary_domains[prot]['words']=[]



from os import listdir
from os.path import isfile, join
mypath=main_f+"/Hs_PfamA_raw"
onlyfiles = [f for f in listdir(mypath)]
dfs=[]

for file_ in onlyfiles:
    try:
        df=load_domtable(mypath+'/'+file_+'/LP_'+file_+'/PfamA-'+file_+'.raw')
        dfs.append(df)
    except:
        print file_ #Hs_rand raw file appears to be empty

a=pd.concat(dfs) #Put together the different pandas dataframes into a single one
i=0
f1=open(main_f+"/results/pfam_structure_dataset.txt", 'w')
f3=open(main_f+"/results/internal_loop_lengths.txt", 'w')

#print headers
print>>f1,'protein\tpfam_domain\tpfam_length\tage_group\tstart\tend\tstructure\tnum_L\tnum_middle_L\tnum_H\tnum_E\tdomain_description\tword'
print>>f3,'protein\tgroup\tloop_length'

for index, row in a.iterrows(): #Iterate the single pandas dictionary and generate some variables for the different datasets
    if row['target']!='#': #Ignore commented lines / description of the dataset
        if row['target'] in structure_dictionary:
            if structure_dictionary[row['target']][0] != '':
                struc = structure_dictionary[row['target']][0][int(row['alncoord_from'])-1:int(row['alncoord_to'])-1]
                group = structure_dictionary[row['target']][1]
                for n in range(int(row['alncoord_from'])-1, int(row['alncoord_to'])-1):
                    summary_domains[row['target']]['struct_list'][n]='X'
                summary_domains[row['target']]['num_domains']+=1
                summary_domains[row['target']]['total_domain_length']=summary_domains[row['target']]['struct_list'].count('X')
                intra_distances_domain=calc_loop_lengths(struc,'L')[2]
                intra_distances_domain_whole=calc_loop_lengths(summary_domains[row['target']]['struct'],'L')[2]
                max_length=0
                for item in intra_distances_domain:
                    if item > max_length:
                        max_length=item
                print>>f3, row['target']+'\tpfam\t'+format(max_length)
                max_length=0
                for item in intra_distances_domain_whole:
                    if item > max_length:
                        max_length=item
                print>>f3, row['target']+'\tpfam_whole\t'+format(max_length)
            else:
                struc='NA'
                group='?'
                i+=1
        else:
            struc='NA'
            group='?'
            i+=1
        word=create_word(struc, True)[0]
        print>>f1, row['target']+'\t'+row['q_accession']+'\t'+format(row['qlen'])+'\t'+group+'\t'+format(row['alncoord_from'])+'\t'+format(row['alncoord_to'])+'\t'+struc+'\t'+format(struc.count('L'))+'\t'+format(struc.strip('L').count('L'))+'\t'+format(struc.count('H'))+'\t'+format(struc.count('E'))+'\t'+row['target_decsription']+'\t'+word[0]

f1.close()


f4=open(main_f+"/results/inter_domains_loop_lengths.txt", 'w')
print>>f4,'protein\tloop_length'
f7=open(main_f+"/results/words_proteins.txt", 'w')
print>>f7,'protein\tword\tword_num\tage\tnum_domains\tprotein_length'

#Extract the length of internal "loops" in all proteins and also in all secondary structures that overlap a pfam domain.
words_dict={}
for prot in summary_domains:
    print>>f0, prot+'\t'+summary_domains[prot]['age']+'\t'+format(summary_domains[prot]['num_domains'])+'\t'+format(summary_domains[prot]['total_domain_length'])+'\t'+format(summary_domains[prot]['protein_length'])
    intra_distances_domain=calc_loop_lengths(summary_domains[prot]['struct'],'L')[2]
    max_length=0
    for item in intra_distances_domain:
         if item > max_length:
             max_length=item
    print>>f3, prot+'\twhole\t'+format(max_length)
    if summary_domains[prot]['num_domains']>=2:
        lengths = inter_domains_loop_length(summary_domains[prot]['struct_list'])
        for length in lengths:
            print>>f4, prot+'\t'+format(length)
    i=0
    summary_domains[prot]['words']=create_word(summary_domains[prot]['struct'], False)[0]
    age = summary_domains[prot]['age']
    if age not in words_dict:
        words_dict[age]={}
    for word in summary_domains[prot]['words']:
        word=word.strip('L') #We remove terminal loops in words.
        if word=='':
            word='L'
        if word not in words_dict[age]:
            words_dict[age][word]=1
        else:
            words_dict[age][word]+=1
        i+=1
        print>>f7,prot+'\t'+word+'\t'+format(i)+'\t'+age+'\t'+format(summary_domains[prot]['num_domains'])+'\t'+format(summary_domains[prot]['protein_length'])

f7.close()
f0.close()
f4.close()
f3.close()

f5=open(main_f+"/results/intra_domain_stretch_lengths.txt", 'w')
f6=open(main_f+"/results/intra_domain_stretch_composition.txt", 'w')

print>>f5,'protein\tstretch_length\tletter\tgroup\tpfam_id'
print>>f6,'protein\tnumber_letter\tpfam_domain\tdomain_length\tgroup\tletter'


#Composition of pfam domains
a=pd.concat(dfs) #Put together the different pandas dataframes into a single one
for index, row in a.iterrows(): #Iterate the single pandas dictionary and generate some variables for the different datasets
    if row['target']!='#': #Ignore commented lines / description of the dataset
        if row['target'] in structure_dictionary:
            if structure_dictionary[row['target']][0] != '':
                struc = structure_dictionary[row['target']][0][int(row['alncoord_from'])-1:int(row['alncoord_to'])-1]
                group = structure_dictionary[row['target']][1]
                for letter in word_parts:
                    b = calc_loop_lengths(struc, letter)[2]
                    for item in b:
                        print>>f5,row['target']+'\t'+format(item)+'\t'+letter+'\t'+group+'\t'+row['q_accession']
                    print>>f6, row['target']+'\t'+format(struc.count(letter))+'\t'+row['q_accession']+'\t'+format(len(struc))+'\t'+group+'\t'+letter
f5.close()
f6.close()

f6b=open(main_f+"/results/all_proteins_stretch_composition.txt", 'w')
print>>f6b,'protein\tnumber_letter\tstretch_num\tlength\tgroup\tletter'
f6c=open(main_f+"/results/all_proteins_composition.txt", 'w')
print>>f6c,'protein\tnumber_letter\tlength\tgroup\tletter'

for prot in structure_dictionary:
    structure = structure_dictionary[prot][0]
    group = structure_dictionary[prot][1]
    i=0
    for letter in word_parts:
        print>>f6c, prot+'\t'+format(structure.count(letter))+'\t'+format(len(structure))+'\t'+group+'\t'+letter
    for struct in create_word(structure, False)[1]:
        i+=1
        for letter in word_parts:
            print>>f6b, prot+'\t'+format(struct.count(letter))+'\t'+format(i)+'\t'+format(len(struct))+'\t'+group+'\t'+letter

f6b.close()

f8=open(main_f+"/results/word_count.txt", 'w')
for age in words_dict:
    for word in words_dict[age]:
        print>>f8, word+'\t'+format(words_dict[age][word])+'\t'+age

f8.close()


#Top 20 more frequent words for each group. The empty ones are loop only ''=L. We extract composition.
###
f9=open(main_f+"/results/top20_word_count.txt", 'w')
print>>f9,'word age value L E H word_length'

for age in words_dict:
    for (word, value) in sorted(words_dict[age].items(), key=lambda x:x[1])[-20:]:
        if word == '':
            word='L'
        print>>f9, word, age, value, word.count('L'), word.count('E'),word.count('H'), len(word)

f9.close()

f9=open(main_f+"/results/top20_word_count.txt", 'w')
print>>f9,'word age value L E H word_length'

for age in words_dict:
    for (word, value) in sorted(words_dict[age].items(), key=lambda x:x[1])[-20:]:
        if word == '':
            word='L'
        print>>f9, word, age, value, word.count('L'), word.count('E'),word.count('H'), len(word)

f9.close()

########################################################
#We are only keeping proteins with word lengths 1 to 3.#
########################################################

#Generating all possible 2 and 3 structure combinations
from itertools import product
def reverse(text):
    return text[::-1]

alphabets = ['E','H','L']
ages=["Hs_ps1-2", "Hs_ps13-17", "Hs_ps20-21", "Hs_ps31", "Hs_igen", "Hs_rand"]

combinations_3 = []
combinations_2 = []
#All 3 letter possible combinations#put palindromes together
n=3
for (a,b,c) in product(alphabets, repeat = n):
    comb=(a+b+c).replace('LLL','L').replace('LL','L')
    if len(comb)==n and reverse(comb) not in combinations_3:
        combinations_3.append(comb)

#All 2 letter possible combinations. #put palindromes together
n=2
for (a,b) in product(alphabets, repeat = n):
    comb=(a+b).replace('LL','L')
    if len(comb)==n and reverse(comb) not in combinations_2:
        combinations_2.append(comb)

#Put them together
set_combinations = combinations_2 + combinations_3 + alphabets
word_dict_short = {}
for age in ages:
    if age not in word_dict_short:
        word_dict_short[age]={}
    for comb in set_combinations:
        word_dict_short[age][comb]=0

for age in words_dict:
    for word in words_dict[age]:
        if len(word)<=3:
            dict_word=word
            if word == '':
                word='L'
                dict_word=''
            if reverse(word) in word_dict_short[age]:
                word_dict_short[age][reverse(word)]+=words_dict[age][dict_word]
            elif word in word_dict_short[age]:
                word_dict_short[age][word]+=words_dict[age][dict_word]

f9=open(main_f+"/results/word_count_length_1_to_3.txt", 'w')
print>>f9, 'word\tcount\tage'
for age in word_dict_short:
    for (word, value) in sorted(word_dict_short[age].items(), key=lambda x:x[1])[-20:]:
        print>>f9, word + '\t' + format(value) + '\t'+age

f9.close()