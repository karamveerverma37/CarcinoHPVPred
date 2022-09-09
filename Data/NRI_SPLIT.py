# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 16:34:05 2020

@author: Megha Mathur
"""

import argparse  
import warnings
from collections  import defaultdict
from collections import deque
import  pandas as pd
import os

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-split","--svalue",type=int, help="Enter the split value")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
f1 = args.input 
if args.svalue == None:
    split = int(0)
else:
    split = int(args.svalue)
        
alphabet=['A','C','G','T']
def sp(se,split,NRI):
    sqnc=[]
    while(split!=0):
        i = int(len(se)/split)
        s1=""
        for j in range(i):
            s1=s1+se[j]
        
        se=deque(se)
        for j in range(i):
            se.popleft()
      
        sqnc.append(s1)
        split=split-1
    count=0
    for sq1 in sqnc:
        count=count+1
        for i in alphabet:
            c=0
            s1=0
            s2=0
            j=0
            while(j<(len(sq1))):
                if(sq1[j]==i):
                    c=c+1
                elif(sq1[j]!=i):
                    s1 = (c**2)+s1
                    s2 = s2+c
                    c=0
                j=j+1
            if(c!=0):
                s1 = s1+(c**2)
                s2 = s2+c
                c=0
            if(s2!=0):
                t = s1/s2
            else:
                t=0.0
            NRI["NRI_SPLIT_s"+str(count)+"_"+i].append(t)
seq=[]
filename, file_extension = os.path.splitext(f1)
cdk = pd.DataFrame()
if(file_extension==""):
    f1=f1.upper()
    alphabet=['A','C','G','T']
    for i in f1:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    seq.append(f1)
    cdk['Sequence'] = seq
    
else:
    f=open(f1,"r")
    b= f.readlines()
    s_id =[]
    s=""
    f.close()
    for i in b:
        if i[0] == '>':
            i=i.split("\n")
            s_id.append(i[0])
            if s!= "":
                seq.append(s)
                s=""
                
            else:
                continue
        else:
            for j in i:
                j=j.capitalize()
                if(j in ['A','G','C','T']):
                    s = s+j
    if s!="":
        seq.append(s)
        cdk['Sequence_ID'] =s_id
NRI= defaultdict(list)
for i in seq:
    if(split>len(i)):
        continue
    else:
        sp(i,split,NRI)
for i in NRI.keys():
    cdk[i]=NRI[i]
cdk.to_csv(out,index=False)
