# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 15:40:04 2020

@author: Megha Mathur
"""
from collections import defaultdict
import os
import pandas as pd
import argparse  
import warnings
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed')
parser.add_argument("-i", "--input",required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and by default it is set to 2") 
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
    
def reverse(sq):
    temp=""
    for i in reversed(range(len(sq))):
        temp+=sq[i]
    return temp
         
         
def reverse_mer(sq):
    sq1=reverse(sq)
    temp_sq=""
    for i in sq1:
        if i== 'A':
            temp_sq += 'T'
        elif i== 'T':
            temp_sq += 'A'
        elif i== 'C':
            temp_sq += 'G'
        elif i== 'G':
            temp_sq += 'C'
    if(temp_sq<=sq):
        return temp_sq
    else:
        return sq
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
dict_n = defaultdict(list)
def rev_c(seq,k,dict_n):
    a=allp(k)
    a.sort()
    r = []
    r_dict = defaultdict(list)
    mers= kmer(k,seq,1)
    for i in a:
        rev = reverse_mer(i)
        if rev not in r:
            r.append(rev)
    r.sort()
    for i in r:
        r_dict[i].append(0)
    
    for i in mers:
        rev = reverse_mer(i)
        val = r_dict[rev][0]+1
        r_dict[rev].pop()
        r_dict[rev].append(val)
    r_dict1 = defaultdict(list)
    count =0
    for i in r_dict.keys():
        count = count+1
        st = i
        val = r_dict[i][0]
        v = val*100/(len(mers))
        r_dict1[st].append(v)
    for i in r_dict1.keys():
        dict_n[i].append(r_dict1[i][0])

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
    rev_c(f1,k,dict_n)
    s = []
    s.append(f1)
    if 'Sequence' not in cdk.columns:
        cdk['Sequence']=s
    for i in dict_n.keys():
        cdk["RDK_"+i]= dict_n[i]
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
        rev_c(i,k,dict_n)
    if 'Sequence_ID' not in cdk.columns:
        cdk['Sequence_ID'] = s_id
    
    for i in dict_n.keys():
        cdk["RDK_"+i]= dict_n[i]
    cdk.to_csv(out,index=False)


    