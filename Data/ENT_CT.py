# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 09:17:43 2020

@author: Megha Mathur
"""


import argparse  
import warnings
from collections  import defaultdict
import  pandas as pd
import math
import os

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-c","--cvalue",type=int, help="Enter the n value")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
f1 = args.input 
alphabet=['A','C','G','T']

if args.cvalue == None:
    n = 2
else:
    n = int(args.cvalue)
        
#sq1=""
#if n>len(s):
#    print("Invalid Sequence Length")
#    exit()
#i=len(s)-1
#for i in range(n):
#    sq1=sq1+s[i]
#su=0
#for i in alphabet:
#    c=0
#    j=0
#    while(j<(len(sq1))):
#        if(sq1[j]==i):
#            c=c+1
#        
#        j=j+1
#    p = c/(len(sq1))
#    if p>0:
#        su=su+(-(p*(math.log2(p))))
#print(su) 
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
ENT = []
for s in seq:
    sq1=""
    i=len(s)-1
    for i in range(n):
        sq1=sq1+s[i]
    su=0
    for i in alphabet:
        c=0
        j=0
        while(j<(len(sq1))):
            if(sq1[j]==i):
                c=c+1
            
            j=j+1
        p = c/(len(sq1))
        if p>0:
            su=su+(-(p*(math.log2(p))))
    ENT.append(su)
cdk["Entropy_CT"]=ENT
cdk.to_csv(out,index =False) 
