# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 15:42:46 2020

@author: Megha Mathur
"""

import sys
from collections import deque
from collections import defaultdict
import os
import  pandas as pd
import argparse  
import warnings

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input",required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and by default it is set to 2")
parser.add_argument("-s","--split",required=True,type=int, help="Enter the split value")
parser.add_argument("-or","--order",type=int, help="Enter the order value and by default it is set to 1")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args() 
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
sq = args.input 
    
if args.kvalue == None:
    k = int(2)
else:
    k = int(args.kvalue)
    
if args.order == None:
    order = int(1)
else:
    order = int(args.order)  
def allp(k):
  n =['A','T','C','G']
  s=[]
  if (k==1):
    return n
  elif(k==2):
    for i in n:
      for j in n:
        se = i+j
        if (se not in s):
          s.append(se)
    return s
  elif(k==3):
    for i in n:
      for j in n:
        for k in n:
          se = i+j+k
          if (se not in s):
            s.append(se)
    return s

  elif(k==4):
    for i in n:
      for j in n:
        for k in n:
          for l in n:
            se = i+j+k+l
            if (se not in s):
              s.append(se)
    return s
def kmer(k,seq,order):
    s=[]
    for i in range(len(seq)):
        se=""
        if (i+k > len(seq)):
            break
        for j in range(k):
            if i+(j*order)>=len(seq):
                break
            else:
                se = se+seq[i+(j*order)]
        if len(se)<k:
            break
        s.append(se)
    return s
split = args.split
def sp(seq,split,dict_n):
    sq1 = []
    for i in seq:
        sq1.append(i)
    
    if(split>len(seq)):
        sys.exit()
    sqnc=[]
    while(split!=0):
        i = int(len(sq1)/split)
        s=""
        for j in range(i):
            s=s+sq1[j]
            sq1[j]=0
        sq1=deque(sq1)
        for j in range(i):
            sq1.popleft()
      
        sqnc.append(s)
        split=split-1
    a=allp(k)
    count=0
    for sq in sqnc:
        if(len(sq)<k):
            print("Invalid length")
            sys.exit()
        count=count+1
        rs = kmer(k,sq,order)
        
    
    #calculating basic k-mer composition
        for i in a:
            ct =(rs.count(i)/len(rs))*100
            dict_n["CDK_Split_s"+str(count)+"_"+i].append(ct)
    
        
filename, file_extension = os.path.splitext(sq)
sequence1=[]
sequence=[]
f2 = open(out,'r')
dr = f2.readlines()
if len(dr) == 0:
    edk = pd.DataFrame()
else:
    edk = pd.read_csv(out)
dict_n = defaultdict(list)
if(file_extension==""):
    sq=sq.upper()
    alphabet=['A','C','G','T']
    for i in sq:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    sequence.append(sq)
    sp(sq,split,dict_n)
    if 'Sequence' not in edk.columns:
        edk['Sequence'] = sequence
    for i in dict_n.keys():
        edk[i] = dict_n[i]
    edk.to_csv(out,index=False)
else:
    f=open(sq,"r")
    b= f.readlines()
    s_id =[]
    s_id1 =[]
    s=""
    f.close()
    for i in b:
        if i[0] == '>':
            i=i.split("\n")
            s_id.append(i[0])
            if s!= "":
                sequence.append(s)
                s=""
                
            else:
                continue
        else:
            for j in i:
                j=j.capitalize()
                if(j in ['A','G','C','T']):
                    s = s+j
    if s!="":
        sequence.append(s)
    
    for i in sequence:
        sp(i,split,dict_n)
        
        
    if 'Sequence_ID' not in edk.columns:
        edk['Sequence_ID'] = s_id
    
    for i in dict_n.keys():
        edk[i]= dict_n[i]
    edk.to_csv(out,index=False)
    


        