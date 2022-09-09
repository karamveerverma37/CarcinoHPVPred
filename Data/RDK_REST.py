# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:30:13 2020

@author: Megha Mathur
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:04:39 2020

@author: Megha Mathur
"""


import sys
from collections import defaultdict
import argparse  
import warnings
import os
import  pandas as pd
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input",required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and by default it is set to 2")
parser.add_argument("-n","--nvalue",required=True,type=int, help="Enter the n value")
parser.add_argument("-c","--cvalue",required=True,type=int, help="Enter the m value")
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
m = int(args.cvalue)
sq1=""
if n+m>len(sq):
    print("Invalid Sequence Length")
    sys.exit()
if n+m<=order:
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

def end(seq,n,dict_n,k):
    i=len(seq)-1
    sq1=""
    for i in range(n,len(seq)-m):
        sq1=sq1+seq[i]
    if(len(sq1)<k):
        print("Invalid length")
        sys.exit()
    rev_c(sq1,k,dict_n)
  
    #calculating basic k-mer composition
    
sequence=[]

dict_n=defaultdict(list)
filename, file_extension = os.path.splitext(sq)
f2 = open(out,'r')
dr = f2.readlines()
if len(dr) == 0:
    cdk = pd.DataFrame()
else:
    cdk = pd.read_csv(out)
if(file_extension==""):
    f1=sq.upper()
    alphabet=['A','C','G','T']
    for i in f1:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    end(f1,n,dict_n,k)
    s = []
    s.append(f1)
    if 'Sequence' not in cdk.columns:
        cdk['Sequence']=s
    for i in dict_n.keys():
        cdk["RDK_REST_"+i]= dict_n[i]
    cdk.to_csv(out,index=False)
else:
    f=open(sq,"r")
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
        end(i,n,dict_n,k)
    if 'Sequence_ID' not in cdk.columns:
        cdk['Sequence_ID'] = s_id
    
    for i in dict_n.keys():
        cdk["RDK_REST_"+i]= dict_n[i]
    cdk.to_csv(out,index=False)    