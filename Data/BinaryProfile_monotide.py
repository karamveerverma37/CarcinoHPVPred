# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:07:37 2020

@author: Megha Mathur
"""
import argparse  
import warnings
import os
import csv
import  pandas as pd
from collections  import defaultdict

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 
## Read Arguments from command
parser.add_argument("-i", "--input",required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
f1 = args.input 
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
seq = args.input 
nuc=['A','G','C','T']
f2 = open(out,'r')
dr = f2.readlines()
if len(dr) == 0:
    cdk = pd.DataFrame()
else:
    cdk = pd.read_csv(out)
seq=[]
filename, file_extension = os.path.splitext(f1)
f2 = open(out,"w")
max_l=0
if(file_extension==""):
    f1=f1.upper()
    alphabet=['A','C','G','T']
    for i in f1:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    seq.append(f1)
    s_id =[]
    s_id.append(f1)
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
                if len(s)>max_l:
                    max_l=len(s)
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
        if(len(s)>max_l):
            max_l=len(s)
    cdk['Sequence_ID'] =s_id
d1= defaultdict(list)
for p,s in enumerate(seq):
    count=0
    for  i in s:
        count=count+1
        for j in nuc:
            t ='P'+str(count)+'_'+j
            
            if j==i:
                d1[t].append(str(1))
                #d = d+","+str(1)
            else:
                d1[t].append(str(0))
                # d= d+","+str(0)
    if count<max_l:
        while count!=max_l:
            count=count+1
            for j in nuc:
                t = 'P' + str(count) + '_' + j
                d1[t].append(str(0))
for i in d1.keys():
    cdk[i]= d1[i]
cdk.to_csv(out,index =False)


