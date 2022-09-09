# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:13:20 2020

@author: Megha Mathur
"""

import argparse  
import warnings
from collections  import defaultdict
import  pandas as pd
import os
import sys

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and its default value is 3")
parser.add_argument("-w","--wvalue",type=float, help="Enter the w value and its default value is 0.5")
parser.add_argument("-lm","--lmvalue",type=int, help="Enter the lamada value and its default value is 1")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
f1 = args.input 
     
    
if args.kvalue == None:
    k = int(3)
else:
    k = int(args.kvalue)
    
if args.wvalue == None:
    w = float(0.5)
else:
    w = float(args.wvalue)
if  w<0.0 or w>1.0:
    print("w value should vary from 0 to 1")
    sys.exit()
if args.lmvalue == None:
    lm = int(1)
else:
    lm = int(args.lmvalue)    
    
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


rs = allp(3) 
rs.sort()
ph_v = {
      "AA": [0.06, 0.5, 0.09, 1.59, 0.11, -0.11],
        "AC": [1.5, 0.5, 1.19, 0.13, 1.29, 1.04],
        "GT": [1.5, 0.5, 1.19, 0.13, 1.29, 1.04],
        "AG": [0.78, 0.36, -0.28, 0.68, -0.24, -0.62],
        "CC": [0.06, 1.08, -0.28, 0.56, -0.82, 0.24],
        "CA": [-1.38, -1.36, -1.01, -0.86, -0.62, -1.25],
        "CG": [-1.66, -1.22, -1.38, -0.82, -0.29, -1.39],
        "TT": [0.06, 0.5, 0.09, 1.59, 0.11, -0.11],
        "GG": [0.06, 1.08, -0.28, 0.56, -0.82, 0.24],
        "GC": [-0.08, 0.22, 2.3, -0.35, 0.65, 1.59],
        "AT": [1.07, 0.22, 0.83, -1.02, 2.51, 1.17],
        "GA": [-0.08, 0.5, 0.09, 0.13, -0.39, 0.71],
        "TG": [-1.38, -1.36, -1.01, -0.86, -0.62, -1.25],
        "TA": [-1.23, -2.37, -1.38, -2.24, -1.51, -1.39],
        "TC": [-0.08, 0.5, 0.09, 0.13, -0.39, 0.71],
        "CT": [0.78, 0.36, -0.28, 0.68, -0.24, -0.62],
        }
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
        
res = defaultdict(list)
for s in seq:
    if len(s) < k or lm + k > len(s):
        continue
    mer = kmer(2,s)
    mer2 = kmer(3,s)
    fre=[]
    done=[]
    for i in rs:
        fre.append(mer2.count(i))
    
    fre_sum = sum(fre)
    fre = [(f/fre_sum) for f in fre]
    theta =[]
    for i in range(1,lm+1):
        temp_sum =0.0
        for j in range(len(s)-1-lm):
            n1 = mer[j]
            n2 = mer[j+i]
            temp = 0.0
            for y in range(6):
                temp += ((ph_v[n1][y]-ph_v[n2][y])**2)
            temp = temp/6
            temp_sum += temp
        temp_sum = (temp_sum/(len(s)-i-1))
        theta.append(temp_sum)
    t_sum = sum(theta)
    dm = 1 + w*t_sum
    temp_vec =[f/dm for f in fre]
    for i in theta:
        temp_vec.append(w*i/dm)
    for i in range(0,len(temp_vec)):
        if i<64:
            st="PKNC_"+str(rs[i])
        else:
            st = "PKNC_lm_"+str(i-63)
        res[st].append(temp_vec[i-1])
for i in res.keys():
    cdk[i]= res[i]
cdk.to_csv(out,index =False) 
