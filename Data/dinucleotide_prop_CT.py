# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 23:27:42 2020

@author: Megha Mathur
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 23:27:43 2020

@author: Megha Mathur
"""

import  pandas as pd
import argparse  
import warnings
from collections  import defaultdict
import os
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input",required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-p", "--property",type=str,nargs='+',required=True,help=" Refer the property list of dipeptides and enter property names having space in between")
parser.add_argument("-or", "--order",type=int,help=" Enter the order of nucleotide")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
parser.add_argument("-c","--cvalue",required=True,type=int, help="Enter the c value")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
if args.order == None:
    order = int(1)
else:
    order = int(args.order)
f1= args.input 
n= int(args.cvalue)
py= args.property
property_list=[]

data = pd.DataFrame()
data['Physicochemical properties']=['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28','p29','p30','p31','p32','p33','p34','p35','p36','p37','p38']
data['GA']=[-0.654, -0.14400000000000002, 1.112, 0.0, 1.023, -0.27399999999999997, 0.27, -0.5, -0.21899999999999997, -0.11800000000000001, 0.8270000000000001, 0.15, -0.19, -0.5, -1.105, -1.251, 0.0, -0.413, -0.495, 1.023, -0.402, 0.0, 0.17, -0.026000000000000002, -0.036000000000000004, -0.434, 0.516, 0.6509999999999999, -0.23600000000000002, 0.08800000000000001, 0.191, 0.0, -0.081, 0.502, 0.266, 0.126, -0.39399999999999996, 0.711]
data['GC']=[-2.455, -0.301, 0.7859999999999999, 1.369, 0.322, 0.47200000000000003, -1.232, 1.4040000000000001, 2.353, 0.6659999999999999, -0.22399999999999998, 0.588, 1.5190000000000001, 1.4040000000000001, 1.35, 1.162, 1.369, -1.7069999999999999, 1.4080000000000001, 0.322, -1.7069999999999999, 1.369, 1.9969999999999999, 2.354, 2.1, 0.076, 2.517, 2.45, 1.3530000000000002, 1.04, 0.8440000000000001, 1.369, -0.081, 0.215, 1.331, -0.348, 0.6459999999999999, 1.585]
data['GG']=[-0.07, 0.355, -0.055999999999999994, 1.369, -1.1909999999999998, 1.3969999999999998, -1.232, 1.4040000000000001, 0.67, 2.076, -0.496, -0.579, -0.498, 1.4040000000000001, 1.306, 1.1440000000000001, 1.369, -0.485, 1.4080000000000001, -1.1909999999999998, -0.488, 1.369, 0.86, -0.726, -1.276, -0.789, 0.49700000000000005, 0.068, 0.956, 1.264, 1.2930000000000001, 1.369, 0.063, 1.077, 0.08900000000000001, 0.56, -0.8220000000000001, 0.242]
data['GT']=[-0.9179999999999999, -0.831, -0.653, 0.0, -1.359, -0.156, 0.27, -0.8809999999999999, 1.107, -0.11800000000000001, -1.041, 1.025, 0.259, -0.8809999999999999, -0.703, -0.556, 0.0, -1.276, -0.887, -1.359, -1.278, 0.0, 0.10300000000000001, 0.604, 0.6759999999999999, 0.852, 0.971, 0.915, -0.23600000000000002, 0.424, 0.6409999999999999, 0.0, 1.5019999999999998, 0.502, 0.799, 0.126, 1.2890000000000001, 1.044]
data['AA']=[1.0190000000000001, -0.644, -0.002, -1.369, 0.995, -1.8869999999999998, 0.833, -0.11900000000000001, -0.8420000000000001, -0.9009999999999999, 0.36, 0.515, -0.9329999999999999, -0.11900000000000001, 0.45799999999999996, 0.6679999999999999, -1.369, 0.593, -0.132, 0.995, 0.593, -1.369, -0.81, 0.46399999999999997, 0.8340000000000001, -0.7, -0.77, -1.02, -0.831, -0.361, -0.16899999999999998, -1.369, 0.063, 0.502, 0.266, 1.587, 0.111, -0.109]
data['AC']=[-0.9179999999999999, -0.831, -0.653, 0.0, -1.359, -0.156, 0.27, -0.8809999999999999, 1.107, -0.11800000000000001, -1.041, 1.025, 0.259, -0.8809999999999999, -0.703, -0.556, 0.0, -1.276, -0.887, -1.359, -1.278, 0.0, 0.10300000000000001, 0.604, 0.6759999999999999, 0.852, 0.971, 0.915, -0.23600000000000002, 0.424, 0.6409999999999999, 0.0, 1.5019999999999998, 0.502, 0.799, 0.126, 1.2890000000000001, 1.044]
data['AG']=[0.488, -0.894, -1.3319999999999999, 0.0, -0.799, -0.436, 0.27, -0.5, 0.016, -0.11800000000000001, -0.885, 0.15, -0.99, -0.5, -0.12300000000000001, 0.083, 0.0, 0.23399999999999999, -0.495, -0.799, 0.233, 0.0, -0.498, -1.147, -1.1440000000000001, -0.5670000000000001, -0.612, -0.489, -0.23600000000000002, -1.145, -1.406, 0.0, 0.7829999999999999, 0.359, 0.08900000000000001, 0.679, -0.24100000000000002, -0.623]
data['AT']=[0.5670000000000001, -1.05, 2.089, -1.369, -0.098, -0.7509999999999999, 1.396, -1.3880000000000001, -0.5760000000000001, -1.371, -1.896, 1.973, 1.03, -0.627, 0.23399999999999999, 0.65, -1.369, -0.485, -0.615, -0.098, -0.488, -1.369, -1.456, -0.866, -0.43200000000000005, 3.159, -0.669, -0.568, -1.4269999999999998, -1.705, -1.676, -1.369, 1.071, 0.215, 0.621, -1.0190000000000001, 2.513, 1.171]
data['CA']=[0.5670000000000001, 1.51, 0.596, 0.0, 1.1909999999999998, 0.98, -0.106, -0.11900000000000001, -0.915, -0.11800000000000001, 1.216, -1.38, 0.45399999999999996, -0.11900000000000001, -1.015, -1.361, 0.0, 1.0959999999999999, -0.132, 1.1909999999999998, 1.091, 0.0, -0.008, -0.23600000000000002, -0.3, 0.032, -0.043, -0.568, 0.161, -0.249, -0.371, 0.0, -1.376, -1.364, -0.266, -0.861, -0.623, -1.254]
data['CC']=[-0.07, 0.355, -0.055999999999999994, 1.369, -1.1909999999999998, 1.3969999999999998, -1.232, 1.4040000000000001, 0.67, 2.076, -0.496, -0.579, -0.498, 1.4040000000000001, 1.306, 1.1440000000000001, 1.369, -0.485, 1.4080000000000001, -1.1909999999999998, -0.488, 1.369, 0.86, -0.726, -1.276, -0.789, -0.762, 0.068, 0.956, 1.264, 1.2930000000000001, 1.369, 0.063, 1.077, 0.08900000000000001, 0.56, -0.8220000000000001, 0.242]
data['CG']=[-0.579, 2.229, -1.1420000000000001, 1.369, -0.266, 0.799, -2.17, 2.039, 0.187, 0.6659999999999999, 0.7490000000000001, -1.818, 2.36, 2.039, 1.7069999999999999, 1.3630000000000002, 1.369, 0.665, 2.042, -0.266, 0.662, 1.369, 1.5730000000000002, 1.6540000000000001, 1.335, -0.41200000000000003, -0.762, 0.606, 2.346, 1.768, 1.4280000000000002, 1.369, -1.6640000000000001, -1.22, -0.444, -0.8220000000000001, -0.287, -1.389]
data['CT']=[0.488, -0.894, -1.3319999999999999, 0.0, -0.799, -0.436, 0.27, -0.5, 0.016, -0.11800000000000001, -0.885, 0.15, -0.99, -0.5, -0.12300000000000001, 0.083, 0.0, 0.23399999999999999, -0.495, -0.799, 0.233, 0.0, -0.498, -1.147, -1.1440000000000001, -0.5670000000000001, 0.49700000000000005, -0.489, -0.23600000000000002, -1.145, -1.406, 0.0, 0.7829999999999999, 0.359, 0.08900000000000001, 0.679, -0.24100000000000002, -0.623]
data['TA']=[1.6030000000000002, 0.418, -1.061, -1.369, 0.322, 0.233, 1.396, -0.627, -1.598, -1.371, 1.41, -0.506, -1.114, -1.3880000000000001, -0.9259999999999999, -0.629, -1.369, 2.031, -1.37, 0.322, 2.036, -1.369, -1.746, -1.006, -0.511, 0.387, -1.486, -1.6030000000000002, -1.4269999999999998, -1.145, -0.956, -1.369, -1.2329999999999999, -2.3680000000000003, -0.444, -2.2430000000000003, -1.511, -1.389]
data['TC']=[-0.654, -0.14400000000000002, 1.112, 0.0, 1.023, -0.27399999999999997, 0.27, -0.5, -0.21899999999999997, -0.11800000000000001, 0.8270000000000001, 0.15, -0.19, -0.5, -1.105, -1.251, 0.0, -0.413, -0.495, 1.023, -0.402, 0.0, 0.17, -0.026000000000000002, -0.036000000000000004, -0.434, 0.516, 0.6509999999999999, -0.23600000000000002, 0.08800000000000001, 0.191, 0.0, -0.081, 0.502, 0.266, 0.126, -0.39399999999999996, 0.711]
data['TG']=[0.5670000000000001, 1.51, 0.596, 0.0, 1.1909999999999998, 0.98, -0.106, -0.11900000000000001, -0.915, -0.11800000000000001, 1.216, -1.38, 0.45399999999999996, -0.11900000000000001, -1.015, -1.361, 0.0, 1.0959999999999999, -0.132, 1.1909999999999998, 1.091, 0.0, -0.008, -0.23600000000000002, -0.3, 0.032, -0.612, -0.568, 0.161, -0.249, -0.371, 0.0, -1.376, -1.364, -0.266, -0.861, -0.623, -1.254]
data['TT']=[1.0190000000000001, -0.644, -0.002, -1.369, 0.995, -1.8869999999999998, 0.833, -0.11900000000000001, -0.8420000000000001, -0.9009999999999999, 0.36, 0.515, -0.9329999999999999, -0.11900000000000001, 0.45799999999999996, 0.6679999999999999, -1.369, 0.593, -0.132, 0.995, 0.593, -1.369, -0.81, 0.46399999999999997, 0.8340000000000001, -0.7, -0.77, -1.02, -0.831, -0.361, -0.16899999999999998, -1.369, 0.063, 0.502, -3.284, 1.587, 0.111, -0.109]

for i in data['Physicochemical properties']:
    property_list.append(i)
al=['all']    
for i in py:
    if i not in property_list:
        if i not in al:
            print("No such property found")
            exit()
    
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
def end(seq,n,k):
    i=len(seq)-1
    sq1=""
    for i in range(n):
        sq1=sq1+seq[i]
    if(len(sq1)<k):
        print("Invalid length")
        sys.exit()
    return sq1
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
seq=[]
max_len=0
filename, file_extension = os.path.splitext(f1)
cdk = pd.DataFrame()
if(file_extension==""):
    f1=f1.upper()
    alphabet=['A','C','G','T']
    for i in f1:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    seq.append(end(f1,n,2))
    if max_len < len(end(f1,n,2)):
        max_len = len(end(f1,n,2))
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
                seq.append(end(s,n,2))
                if max_len < len(end(s, n, 2)):
                    max_len = len(end(s, n, 2))
                s=""
                
            else:
                continue
        else:
            for j in i:
                j=j.capitalize()
                if(j in ['A','G','C','T']):
                    s = s+j
    if s!="":
        seq.append(end(s,n,2))
        if max_len < len(end(s, n, 2)):
            max_len = len(end(s, n, 2))
        cdk['Sequence_ID'] =s_id
        
dict_n=defaultdict(list)
for sq in seq:
    pl = allp(2)
    rs = kmer(2,sq,order)
    a=[]
    for p in py:
        if p=="all":
            for q in property_list:
                count=0
                for i in rs:
                    count=count+1
                    for j in pl:
                        t='P'+str(count)+'_'+j
                        if j == i:
                            c=-1
                            for l in data['Physicochemical properties']:
                                c=c+1
                                if(l==q):
                                    w = q.split('p')
                                    s = "PR"+str(w[1])
                                    t1 = t+"_"+s+"_CT"
                                    dict_n[t1].append(data[i][c])
                        else:
                            w = q.split('p')
                            s = "PR"+str(w[1])
                            t1=t+"_"+s+"_CT"
                            dict_n[t1].append(0)
                if count != max_len - 1:
                    while count != max_len - 1:
                        count = count + 1
                        for j in pl:
                            t = 'P' + str(count) + '_' + j
                            w = q.split('p')
                            s = "PR" + str(w[1])
                            t1 = t + "_" + s+"_CT"
                            dict_n[t1].append(0)


        else:
            count=0
            for i in rs:
                count=count+1
                for j in pl:
                    t='P'+str(count)+'_'+j
                    if j == i:
                        c=-1
                        for l in data['Physicochemical properties']:
                            c=c+1
                            if(l==p):
                                w = p.split('p')
                                s = "PR"+str(w[1])
                                t1 = t+"_"+s+"_CT"
                                dict_n[t1].append(data[i][c])
                    else:
                        w = p.split('p')
                        s = "PR"+str(w[1])
                        t1=t+"_"+s+"_CT"
                        dict_n[t1].append(0)
            if count != max_len - 1:
                while count != max_len - 1:
                    count=count+1
                    for j in pl:
                        t='P'+str(count)+'_'+j
                        w = p.split('p')
                        s = "PR" + str(w[1])
                        t1 = t + "_" + s+"_CT"
                        dict_n[t1].append(0)


for i in dict_n.keys():
    cdk[i]= dict_n[i]
cdk.to_csv(out,index =False)
                             
    
