# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:00:33 2020

@author: Megha Mathur
"""

import argparse  
import warnings
from collections  import defaultdict
import  pandas as pd
import os

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
f1 = args.input 

NRI = defaultdict(list)
alphabet=['A','C','G','T']
seq=[]
filename, file_extension = os.path.splitext(f1)

f2 = open(out,'r')
dr = f2.readlines()
if len(dr) == 0:
    cdk = pd.DataFrame()
else:
    cdk = pd.read_csv(out)

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

for s in seq:
    for i in alphabet:
        c=0
        s1=0
        s2=0
        j=0
        while(j<(len(s))):
            if(s[j]==i):
                c=c+1
            elif(s[j]!=i):
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
        NRI["NRI_"+i].append(t)
for i in NRI.keys():
    cdk[i]= NRI[i]
cdk.to_csv(out,index =False) 

