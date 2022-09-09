# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 10:26:53 2020

@author: Megha Mathur
"""

import argparse  
import warnings
from collections  import defaultdict
import  pandas as pd
import os
from collections import deque
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

alphabet=['A','C','G','T']

if args.svalue == None:
    split = int(0)
else:
    split = int(args.svalue)
   
def sp(s,split,DDON):
    sqnc=[]
    while(split!=0):
        i = int(len(s)/split)
        s1=""
        for j in range(i):
            s1=s1+s[j]
        
        s=deque(s)
        for j in range(i):
            s.popleft()
      
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
            o=0
            while(j<(len(sq1))):
                if(sq1[j]!=i):
                    c=c+1
                elif(sq1[j]==i):
                    o=o+1
                    s1 = (c**2)+s1
                    c=0
                j=j+1
            if(c!=0):
                s1 = s1+(c**2)
                c=0
            s2 = (len(sq1)-o+1)
            t = s1/s2
            if(o==0):
                DDON["DDON_s"+str(count)+"_"+i].append(0.0)
            else:
                DDON["DDON_s"+str(count)+"_"+i].append(t)        
    
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
ENT= defaultdict(list)
for i in seq:
    if(split>len(i)):
        continue
    else:
        sp(i,split,ENT)
for i in ENT.keys():
    cdk[i]=ENT[i]
cdk.to_csv(out,index=False)