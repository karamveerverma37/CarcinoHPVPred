# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:44:18 2020

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
parser.add_argument("-i", "--input", type=str,help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and by default it is set to 2")
parser.add_argument("-or","--order",type=int, help="Enter the order value and by default it is set to 1")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
f1 = args.input 
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
#f=open(f1,"r")
#a= f.readlines()
#sequence=[]
#s=""
#f.close()
#for i in a:
#    if s!= "":
#        sequence.append(s)
#        s=""
#    if(i[0]=='>'):
#        continue
#    else:
#        for j in i:
#            j=j.capitalize()
#            if(j in ['A','G','C','T']):
#                s = s+j
#if s!="":
#    sequence.append(s)


if args.kvalue == None:
    k = int(2)
else:
    k = int(args.kvalue)
if args.order == None:
    order = int(1)
else:
    order = int(args.order)
#if(order>=len(f1)):
#    print("Invalid Sequence length")
#    exit()
#if(order>=(len(f1)-k+1)):
#    print("Invalid Sequence length")
#    exit()
    
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

def k_mer_comp(f1,k,order,dict_n):
    a=allp(k)
    rs = kmer(k,f1,order)
    #calculating basic k-mer composition
    for i in a:
        ct =(rs.count(i)/len(rs))*100
        dict_n[i].append(ct)   
dict_n = defaultdict(list)
a=allp(k)
filename, file_extension = os.path.splitext(f1)
f2 = open(out,'r')
check = f2.readlines()
if(len(check)==0):
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
    k_mer_comp(f1,k,order,dict_n)
    s = []
    s.append(f1)
    if 'Sequence'  not in cdk.columns:
        cdk['Sequence']=s
    for i in a:
        cdk["CDK_"+i]= dict_n[i]
    cdk.to_csv(out,index=False)
else:
    f=open(f1,"r")
    b= f.readlines()
    sequence=[]
    s_id =[]
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
        k_mer_comp(i,k,order,dict_n)
    if 'Sequence_ID' not in cdk.columns:
        cdk['Sequence_ID'] = s_id
    
    for i in a:
        cdk["CDK_"+i]= dict_n[i]
    cdk.to_csv(out,index=False)

exit()
