# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:23:59 2020

@author: Megha Mathur
"""
import argparse
import warnings
import os
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
def kmer(k,seq):
  s=[]
  for i in range(len(seq)):
    se=""
    if (i+k > len(seq)):
      break
    for j in range(k):
      se = se+seq[i+j]
    s.append(se)
  return s
dict_n= defaultdict(list)
seq=[]
filename, file_extension = os.path.splitext(f1)
f2 = open(out,'r')
dr = f2.readlines()
if len(dr) == 0:
    cdk = pd.DataFrame()
else:
    cdk = pd.read_csv(out)
max_l=0
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
        if len(s) > max_l:
            max_l = len(s)
        cdk['Sequence_ID'] =s_id
d1= defaultdict(list)
max_l=max_l-2
for s in seq:
    d=s
    allpairs=allp(3)

    rs = kmer(3,s)
    res=[]
    count=0
    for i in rs:
        count=count+1
        for j in allpairs:
            t = 'P' + str(count) + '_' + j

            if j==i:
                d1[t].append(str(1))
                #d = d+","+str(1)
            else:
                d1[t].append(str(0))
    if count<max_l:
        while count!=max_l:
            count=count+1
            for j in allpairs:
                t = 'P' + str(count) + '_' + j
                d1[t].append(str(0))
for i in d1.keys():
    cdk[i]= d1[i]
cdk.to_csv(out,index =False)