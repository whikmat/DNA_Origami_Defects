# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 17:29:26 2022

@author: Wisam
"""

import pandas
import numpy as np
import itertools
import matplotlib.pyplot as plt
import os
import random
import math
import time
import pickle
import gc

def foo(l,n):
     yield from itertools.product(*([l] * n)) 

def count_substrings(string, substring):
    string_size = len(string)
    substring_size = len(substring)
    count = 0
    for i in range(0,string_size-substring_size+1):
        if string[i:i+substring_size] == substring:
            count+=1
    return count

def check_debruijn(codons,seq):
    for i in codons:
        if count_substrings(seq,i) == 0:
            print("did not find codons: " + i)
        if count_substrings(seq,i) > 1:
            return 0 #means has it multiple times: NO DEBRUIJN
    return 1 #1 means True Debruijn

#DEFINE WHAT HAPPENS WHEN GOING BACK IF SELECTION FAILS
def db_step_back(current,combs,codons,ignored_codons,order):
    #WHEN GOING BACK REMEMBER THAT THE PATH OF THE LAST CODON SHOULD NOT BE TAKEN
    #LAST NUCLEOTIDES (2,3,...)
    ignore_codon=current[-order:]  
    #GO A STEP BACK AND DELETE THE LAST NUCLEOTIDE        
    current_new=current[:-1] 
    #CREATE A NEW POOL OF CODONS TO CHOOSE FROM 
    combs_new=[]
    #INSERT ALL CODONS (THEY ARE CREATED SEPERATED FOR A LENGTH OF N NUCLEOTIDES)
    for i in codons:
        combs_new.append(i)
    #DELETE USED CODONS IN THE CURRENT SEQUENCE, 
    #AMOUNT OF CODONS IS LEN(SEQ) - LEN(CODON) + 1
    for i in range(len(current_new)-1):  
    #PULL CODONS OUT OF CURRENT SEQUENCE
        if current_new[i:i+order] in combs_new:
            combs_new.remove(current_new[i:i+order])
    #PASS IGNORED CODON TO THE NEXT STEP BACK, START WITH IGNORED CODON OF THIS STEP
    ignored_codons_new=[ignore_codon]
    #DELETE ALL CODONS TO BE IGNORED FROM PREVIOUS STEPS BACK AND ADD THEM TO BE PASSED TO THE NEXT STEP
    for i in ignored_codons:
        combs_new.remove(i)
        ignored_codons_new.append(i)
    #DELETE THE IGNORED CODON FROM THIS STEP BACK
    combs_new.remove(ignore_codon)
    #MAKE A SELECTION OF COMPATIBLE CODONS TO THE END OF CURRENT SEQUENCE
    combs_sel=list(filter(lambda x: x.startswith(current_new[-(order-1):]), combs_new)) 
    #REPEAT WHEN NO COMPATIBLE NEXT CODON IS FOUND
    if len(combs_sel)==0:
        current_new, combs_sel=db_step_back(current_new,combs_new,codons,ignored_codons_new,order)
    #WHEN FINISHED, PASS THE NEW STARTING SEQUENCE TO CONTINUE BUILDING FROM 
    #AND THE SELECTION OF COMPATIBLE CODONS BACK TO MAIN PROCESS
    return current_new, combs_sel

#GENERATE RANDOM DEBRUIJN SEQUENCES
def random_debruijn(seed, codons, n, debruijn,order):
    #SETUP A STARTING SEQUENCE
    w_seq=seed
    current=seed
    combs=[]
    #CREATE A POOL 
    for i in codons:
        combs.append(i)
    #REMOVE USED CODONS
    for i in range(len(current)-(order-1)):
        if current[i:i+order] in combs:
            combs.remove(current[i:i+order])
    #FIND COMPATIBLE POOL TO THE END OF CURRENT SEQUENCE 
    combs_sel=list(filter(lambda x: x.startswith(current[-(order-1):]), combs))
    current_new=""
    #IF NO COMPATIBLE EXTENTION IS FOUND
    if len(combs_sel)==0:
        while len(combs_sel)==0:
            #WHEN NOTHING FOUND TRY TO GO BACK AND RESOLVE THE WHOLE THING WITH ALGO ABOVE
            #THIS WILL CHANGE CURRENT SEQ AND SELECTION ACCORDINGLY
            current_new, combs_sel=db_step_back(current,combs,codons,[],order)
    #IF THE NEW SEQUENCE IS DIFFERENT, INSERT IT IN THE NEW ONE TO EXTEND FROM
    if len(current_new)>0:
        w_seq=current_new    
    #EXTEND SEQUENCE BY RANDOM COMPATIBLE CHOICE
    rng=random.choice(combs_sel)
    w_seqRNG=w_seq+rng[-1]
    #IF THE SEQUENCE DID NOT REACH FINAL END, REPEAT
    if len(w_seqRNG)<=(n+order-2):
        random_debruijn(w_seqRNG, codons, n, debruijn,order)
    #IF FINAL LENGTH IS REACHED, ADD THIS SEQUENCE TO DEBRUIJN SEQUENCES
    if len(w_seqRNG)==(n+order-1):
        debruijn.append(w_seq[:len(w_seq)-(order-2)])
    #RETURN DEBRUIJN SEQUENCE (IN A LIST FORM)
    return(debruijn)

def debruijn_rng_generator(order,amount):
    #CREATE ALL CODONS OF DESIRED LENGTH 
    all_combs=[]
    for x in foo('ACGT',order):
        all_combs.append(''.join(x))
      
    start_time = time.time()  
    rng_db=[]
    #CREATE X AMOUNT OF RANDOM SEQUENCES, INSERT K^N FINAL LENGTH AND AN EMPTY DEBRUIJN LIST
    for i in range(amount):
        db=random_debruijn(random.choice(all_combs), all_combs, 4**order, [],order)
        #LINEARIZE DEBRUIJN SEQUENCE AND CHECK IF DEBRUIJN CONDITION IS FULFILLED
        #print(order)
        dbn=db[0]+db[0][:order-1]
        checker=check_debruijn(all_combs, dbn )
        #if checker == 0:
            #print(check_debruijn(all_combs, dbn))
        if checker == 1:
            rng_db.append(db[0])
        #APPEND GENERATED BEBRUIJN SEQUENCE TO OUR RANDOM DB SEQUENCE LIST
    #print(str(len(rng_db)) + " out of " + str(amount) + " trials are true De-Bruijn sequences." )
        #print(db[0])
    
    print("Generated " + str(amount) + " De-Bruijn Sequences of order " + str(order) + " in:\n--- %s seconds ---" % (time.time() - start_time))
    return rng_db

def gen_debruijn(order, amount):
    if order==4 and amount > 2000:
        done=0
        while amount-done > 0:
            rng_db=debruijn_rng_generator(4, min(amount-done,2000))
            with open('order4.txt', 'a') as f:
                for i in rng_db:
                    f.write(i+"\n")
            gc.collect(2)
            done+=2000
    else:
        rng_db=debruijn_rng_generator(order, amount)
        with open('order'+str(order)+'.txt', 'a') as f:
            for i in rng_db:
                f.write(i+"\n")
            
        
    
#generated_db=debruijn_rng_generator(4,amount=10)
#print(generated_db)
#print("ACGTBCDEF"[:10])