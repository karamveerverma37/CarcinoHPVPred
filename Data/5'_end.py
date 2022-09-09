# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 16:46:10 2020

@author: Megha Mathur
"""

import sys
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
parser.add_argument("-c","--cvalue",required=True,type=int, help="Enter the c value")
parser.add_argument("-or","--order",type=int, help="Enter the order value and by default it is set to 1")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
f1 = args.input 
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
sq = args.input 

if args.kvalue == None:
    k = int(2)
else:
    k = int(args.kvalue)
    
n= int(args.cvalue)
if args.order == None:
    order = int(1)
else:
    order = int(args.order)   

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
    for i in range(n):
        sq1=sq1+seq[i]
    a=allp(k)
    if(len(sq1)<k):
        print("Invalid length")
        sys.exit()
    rs = kmer(k,sq1,order)
    #calculating basic k-mer composition
    for i in a:
        ct =(rs.count(i)/len(rs))*100
        temp = "CDK_CT_"+i
        dict_n[temp].append(ct)
filename, file_extension = os.path.splitext(sq)
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
    end(sq,n,dict_n)
    if 'Sequence' not in edk.columns:
        edk['Sequence']=sequence
    for i in dict_n.keys():
        edk[i] = dict_n[i]
    edk.to_csv(out,index=False)
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
    if 'Sequence_ID' not in edk.columns:
        edk['Sequence_ID'] = s_id
    for i in dict_n.keys():
        edk[i]= dict_n[i]
    edk.to_csv(out,index=False)


