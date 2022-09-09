# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 16:56:47 2020

@author: Megha Mathur
"""

import sys
from collections import defaultdict
import argparse  
import warnings
import os
import pandas as pd
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input",required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and by default it is set to 2")
parser.add_argument("-n","--nvalue",required=True,type=int, help="Enter the n value")
parser.add_argument("-or","--order",type=int, help="Enter the order value and by default it is set to 1")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()

if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output

sq=args.input
if args.kvalue == None:
    k = int(2)
else:
    k = int(args.kvalue)
if args.order == None:
    order = int(1)
else:
    order = int(args.order)    
n= int(args.nvalue)
sq1=""
if n>len(sq):
    print("Invalid Sequence Length")
    sys.exit()
if n<=order:
    print("Invalid Order and N value")
    sys.exit()

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
def end(seq,n,dict_n):
    i=len(seq)-1
    sq1=""
    while(n!=0):
        sq1=seq[i]+sq1
        i=i-1
        n=n-1
 
    a=allp(k)
    if(len(sq1)<k):
        print("Invalid length")
        sys.exit()
    rs = kmer(k,sq1,order)
  
    #calculating basic k-mer composition
    for i in a:
        ct =(rs.count(i)/len(rs))*100
        temp = "CDK_NT_"+i
        dict_n[temp].append(ct)
filename, file_extension = os.path.splitext(sq)
f2 = open(out,'r')
check = f2.readlines()
if(len(check)==0):
    cdk = pd.DataFrame()
else:
    cdk = pd.read_csv(out)
sequence=[]
dict_n = defaultdict(list)
if(file_extension==""):
    sq=sq.upper()
    alphabet=['A','C','G','T']
    for i in sq:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    sequence.append(sq)
    if 'Sequence'  not in cdk.columns:
        cdk['Sequence']=sequence
    end(sq,n,dict_n)
    for i in dict_n.keys():
        cdk[i] = dict_n[i]
    cdk.to_csv(out,index=False)
else:
    f=open(sq,"r")
    b= f.readlines()
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
        end(i,n,dict_n)
    if 'Sequence_ID' not in cdk.columns:
        cdk['Sequence_ID'] = s_id
    for i in dict_n.keys():
        cdk[i]= dict_n[i]
    cdk.to_csv(out,index=False)

    