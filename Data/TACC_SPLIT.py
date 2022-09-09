# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:37:36 2020

@author: Megha Mathur
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:37:33 2020

@author: Megha Mathur
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:31:14 2020

@author: Megha Mathur
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 15:30:43 2020

@author: Megha Mathur
"""

import argparse  
import warnings
from collections  import defaultdict
import  pandas as pd
from collections import deque
import sys
import os

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-p", "--property",type=str,nargs='+',required=True,help=" Refer the property list of dipeptides and enter property names having space in between")
parser.add_argument("-l","--lagvalue",type=int, help="Enter the lag value and its default value is 2")
parser.add_argument("-s","--split",required=True,type=int, help="Enter the split value")
parser.add_argument("-o","--output",type=str, help="Enter the output file name")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
if args.output == None:
    out= "outfile.csv" 
else:
    out= args.output
f1 = args.input 
split = args.split
     
if args.lagvalue == None:
    lag = int(2)
else:
    lag = int(args.lagvalue)
prop= args.property
property_list=[]
# Creating a dataframe for all the dinucleotides properties
data = pd.DataFrame()
data['Physicochemical properties']=['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12']
data['AAA']=[-2.0869999999999997, -2.745, -1.732, -2.349, -2.7439999999999998, -2.7439999999999998, 2.274, 2.1180000000000003, -1.0, -1.0, -2.342, 2.386]
data['AAC']=[-1.5090000000000001, -1.354, -0.5770000000000001, -0.561, -1.3630000000000002, -1.3630000000000002, 1.105, 1.516, -1.0, -1.0, -0.555, 0.5479999999999999]
data['AAG']=[-0.506, -0.257, -0.5770000000000001, 0.155, -0.26, -0.26, 0.193, 0.493, -1.0, -1.0, 0.16899999999999998, -0.179]
data['AAT']=[-2.126, -2.585, -1.732, -1.9909999999999999, -2.591, -2.591, 2.141, 2.158, -1.0, -1.0, -2.004, 2.032]
data['ACA']=[0.111, 0.171, -0.5770000000000001, 0.155, 0.16399999999999998, 0.16399999999999998, -0.153, -0.12300000000000001, 1.0, 1.0, 0.16899999999999998, -0.179]
data['ACC']=[-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.071, 0.071, -0.078, 0.107, 1.0, 1.0, 0.266, -0.275]
data['ACG']=[-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.065, 0.065, -0.07400000000000001, 0.107, 1.0, 1.0, 0.266, -0.275]
data['ACT']=[-1.354, -0.685, -0.5770000000000001, 0.45299999999999996, -0.6759999999999999, -0.6759999999999999, 0.536, 1.357, 1.0, 1.0, 0.45899999999999996, -0.466]
data['AGA']=[0.381, -0.15, -0.5770000000000001, -0.74, -0.158, -0.158, 0.109, -0.389, 1.0, 1.0, -0.748, 0.743]
data['AGC']=[0.304, 0.92, 0.5770000000000001, 1.287, 0.9109999999999999, 0.9109999999999999, -0.753, -0.313, 1.0, 1.0, 1.28, -1.272]
data['AGG']=[-0.313, -0.07, 0.5770000000000001, 0.27399999999999997, -0.07, -0.07, 0.039, 0.3, 1.0, 1.0, 0.266, -0.275]
data['AGT']=[-1.354, -0.685, -0.5770000000000001, 0.45299999999999996, -0.6759999999999999, -0.6759999999999999, 0.536, 1.357, 1.0, 1.0, 0.45899999999999996, -0.466]
data['ATA']=[1.615, 0.5720000000000001, -1.732, -0.978, 0.584, 0.584, -0.491, -1.585, -1.0, -1.0, -0.99, 0.988]
data['ATC']=[-0.737, -0.391, -0.5770000000000001, 0.214, -0.397, -0.397, 0.307, 0.727, -1.0, -1.0, 0.217, -0.22699999999999998]
data['ATG']=[1.229, 1.348, -0.5770000000000001, 0.87, 1.358, 1.358, -1.112, -1.215, -1.0, -1.0, 0.893, -0.894]
data['ATT']=[-2.126, -2.585, -1.732, -1.9909999999999999, -2.591, -2.591, 2.141, 2.158, -1.0, -1.0, -2.004, 2.032]
data['CAA']=[0.265, -0.231, -0.5770000000000001, -0.74, -0.226, -0.226, 0.166, -0.275, -1.0, -1.0, -0.748, 0.743]
data['CAC']=[0.496, 0.7859999999999999, 0.5770000000000001, 0.81, 0.773, 0.773, -0.6459999999999999, -0.503, -1.0, -1.0, 0.797, -0.8]
data['CAG']=[1.5759999999999998, 0.92, 0.5770000000000001, -0.322, 0.92, 0.92, -0.762, -1.5490000000000002, -1.0, -1.0, -0.314, 0.304]
data['CAT']=[1.229, 1.348, -0.5770000000000001, 0.87, 1.358, 1.358, -1.112, -1.215, -1.0, -1.0, 0.893, -0.894]
data['CCA']=[-1.8559999999999999, -1.14, 0.5770000000000001, 0.27399999999999997, -1.139, -1.139, 0.917, 1.876, 1.0, 1.0, 0.266, -0.275]
data['CCC']=[0.07200000000000001, 0.358, 1.732, 0.5720000000000001, 0.345, 0.345, -0.3, -0.084, 1.0, 1.0, 0.555, -0.562]
data['CCG']=[-0.9690000000000001, -0.7120000000000001, 1.732, -0.084, -0.705, -0.705, 0.5579999999999999, 0.9620000000000001, 1.0, 1.0, -0.07200000000000001, 0.062]
data['CCT']=[-0.313, -0.07, 0.5770000000000001, 0.27399999999999997, -0.07, -0.07, 0.039, 0.3, 1.0, 1.0, 0.266, -0.275]
data['CGA']=[0.111, 1.0, 0.5770000000000001, 1.645, 1.012, 1.012, -0.8340000000000001, -0.12300000000000001, 1.0, 1.0, 1.666, -1.646]
data['CGC']=[-0.46799999999999997, 0.385, 1.732, 1.287, 0.379, 0.379, -0.326, 0.455, 1.0, 1.0, 1.28, -1.272]
data['CGG']=[-0.9690000000000001, -0.7120000000000001, 1.732, -0.084, -0.705, -0.705, 0.5579999999999999, 0.9620000000000001, 1.0, 1.0, -0.07200000000000001, 0.062]
data['CGT']=[-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.065, 0.065, -0.07400000000000001, 0.107, 1.0, 1.0, 0.266, -0.275]
data['CTA']=[0.882, -0.09699999999999999, -0.5770000000000001, -1.276, -0.09699999999999999, -0.09699999999999999, 0.062, -0.88, -1.0, -1.0, -1.28, 1.285]
data['CTC']=[0.419, 0.43799999999999994, 0.5770000000000001, 0.27399999999999997, 0.42700000000000005, 0.42700000000000005, -0.365, -0.42700000000000005, -1.0, -1.0, 0.266, -0.275]
data['CTG']=[1.5759999999999998, 0.92, 0.5770000000000001, -0.322, 0.92, 0.92, -0.762, -1.5490000000000002, -1.0, -1.0, -0.314, 0.304]
data['CTT']=[-0.506, -0.257, -0.5770000000000001, 0.155, -0.26, -0.26, 0.193, 0.493, -1.0, -1.0, 0.16899999999999998, -0.179]
data['GAA']=[-0.159, -0.605, -0.5770000000000001, -0.9179999999999999, -0.6, -0.6, 0.474, 0.146, -1.0, -1.0, -0.893, 0.89]
data['GAC']=[0.034, 0.171, 0.5770000000000001, 0.27399999999999997, 0.17800000000000002, 0.17800000000000002, -0.165, -0.046, -1.0, -1.0, 0.266, -0.275]
data['GAG']=[0.419, 0.43799999999999994, 0.5770000000000001, 0.27399999999999997, 0.42700000000000005, 0.42700000000000005, -0.365, -0.42700000000000005, -1.0, -1.0, 0.266, -0.275]
data['GAT']=[-0.737, -0.391, -0.5770000000000001, 0.214, -0.397, -0.397, 0.307, 0.727, -1.0, -1.0, 0.217, -0.22699999999999998]
data['GCA']=[0.7659999999999999, 0.8390000000000001, 0.5770000000000001, 0.5720000000000001, 0.8420000000000001, 0.8420000000000001, -0.7020000000000001, -0.767, 1.0, 1.0, 0.555, -0.562]
data['GCC']=[1.036, 2.097, 1.732, 2.479, 2.089, 2.089, -1.6869999999999998, -1.0290000000000001, 1.0, 1.0, 2.487, -2.4330000000000003]
data['GCG']=[-0.46799999999999997, 0.385, 1.732, 1.287, 0.379, 0.379, -0.326, 0.455, 1.0, 1.0, 1.28, -1.272]
data['GCT']=[0.304, 0.92, 0.5770000000000001, 1.287, 0.9109999999999999, 0.9109999999999999, -0.753, -0.313, 1.0, 1.0, 1.28, -1.272]
data['GGA']=[0.265, -0.09699999999999999, 0.5770000000000001, -0.501, -0.10300000000000001, -0.10300000000000001, 0.066, -0.275, 1.0, 1.0, -0.507, 0.499]
data['GGC']=[1.036, 2.097, 1.732, 2.479, 2.089, 2.089, -1.6869999999999998, -1.0290000000000001, 1.0, 1.0, 2.487, -2.4330000000000003]
data['GGG']=[0.07200000000000001, 0.358, 1.732, 0.5720000000000001, 0.345, 0.345, -0.3, -0.084, 1.0, 1.0, 0.555, -0.562]
data['GGT']=[-0.121, 0.064, 0.5770000000000001, 0.27399999999999997, 0.071, 0.071, -0.078, 0.107, 1.0, 1.0, 0.266, -0.275]
data['GTA']=[0.342, -0.07, -0.5770000000000001, -0.561, -0.062, -0.062, 0.031, -0.35100000000000003, -1.0, -1.0, -0.555, 0.5479999999999999]
data['GTC']=[0.034, 0.171, 0.5770000000000001, 0.27399999999999997, 0.17800000000000002, 0.17800000000000002, -0.165, -0.046, -1.0, -1.0, 0.266, -0.275]
data['GTG']=[0.496, 0.7859999999999999, 0.5770000000000001, 0.81, 0.773, 0.773, -0.6459999999999999, -0.503, -1.0, -1.0, 0.797, -0.8]
data['GTT']=[-1.5090000000000001, -1.354, -0.5770000000000001, -0.561, -1.3630000000000002, -1.3630000000000002, 1.105, 1.516, -1.0, -1.0, -0.555, 0.5479999999999999]
data['TAA']=[0.6890000000000001, -0.284, -1.732, -1.395, -0.275, -0.275, 0.20600000000000002, -0.6920000000000001, -1.0, -1.0, -1.376, 1.3840000000000001]
data['TAC']=[0.342, -0.07, -0.5770000000000001, -0.561, -0.062, -0.062, 0.031, -0.35100000000000003, -1.0, -1.0, -0.555, 0.5479999999999999]
data['TAG']=[0.882, -0.09699999999999999, -0.5770000000000001, -1.276, -0.09699999999999999, -0.09699999999999999, 0.062, -0.88, -1.0, -1.0, -1.28, 1.285]
data['TAT']=[1.615, 0.5720000000000001, -1.732, -0.978, 0.584, 0.584, -0.491, -1.585, -1.0, -1.0, -0.99, 0.988]
data['TCA']=[1.73, 1.348, -0.5770000000000001, 0.27399999999999997, 1.348, 1.348, -1.103, -1.696, 1.0, 1.0, 0.266, -0.275]
data['TCC']=[0.265, -0.09699999999999999, 0.5770000000000001, -0.501, -0.10300000000000001, -0.10300000000000001, 0.066, -0.275, 1.0, 1.0, -0.507, 0.499]
data['TCG']=[0.111, 1.0, 0.5770000000000001, 1.645, 1.012, 1.012, -0.8340000000000001, -0.12300000000000001, 1.0, 1.0, 1.666, -1.646]
data['TCT']=[0.381, -0.15, -0.5770000000000001, -0.74, -0.158, -0.158, 0.109, -0.389, 1.0, 1.0, -0.748, 0.743]
data['TGA']=[1.73, 1.348, -0.5770000000000001, 0.27399999999999997, 1.348, 1.348, 4.522, -1.696, 1.0, 1.0, 0.266, -0.275]
data['TGC']=[0.7659999999999999, 0.8390000000000001, 0.5770000000000001, 0.5720000000000001, 0.8420000000000001, 0.8420000000000001, -0.7020000000000001, -0.767, 1.0, 1.0, 0.555, -0.562]
data['TGG']=[-1.8559999999999999, -1.14, 0.5770000000000001, 0.27399999999999997, -1.139, -1.139, 0.917, 1.876, 1.0, 1.0, 0.266, -0.275]
data['TGT']=[0.111, 0.171, -0.5770000000000001, 0.155, 0.16399999999999998, 0.16399999999999998, -0.153, -0.12300000000000001, 1.0, 1.0, 0.16899999999999998, -0.179]
data['TTA']=[0.6890000000000001, -0.284, -1.732, -1.395, -0.275, -0.275, 0.20600000000000002, -0.6920000000000001, -1.0, -1.0, -1.376, 1.3840000000000001]
data['TTC']=[-0.159, -0.605, -0.5770000000000001, -0.9179999999999999, -0.6, -0.6, 0.474, 0.146, -1.0, -1.0, -0.893, 0.89]
data['TTG']=[0.265, -0.231, -0.5770000000000001, -0.74, -0.226, -0.226, 0.166, -0.275, -1.0, -1.0, -0.748, 0.743]
data['TTT']=[-2.0869999999999997, -2.745, -1.732, -2.349, -2.7439999999999998, -2.7439999999999998, -2.615, 2.1180000000000003, -1.0, -1.0, -2.342, 2.386]
for i in data['Physicochemical properties']:
    property_list.append(i)
al=['all']    
for i in prop:
    if i not in property_list:
        if i not in al:
            print("No such property found")
            exit()    
# function to calculate all possible kmers
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
# function to create kmers for the given input sequence
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
def sp(seq,split,k):
    sq1 = []
    for i in seq:
        sq1.append(i)
    if(split>len(seq)):
        sys.exit()
    sqnc=[]
    while(split!=0):
        i = int(len(sq1)/split)
        s=""
        for j in range(i):
            s=s+sq1[j]
            sq1[j]=0
        sq1=deque(sq1)
        for j in range(i):
            sq1.popleft()
      
        sqnc.append(s)
        split=split-1
    count=0
    for sq in sqnc:
        if(len(sq)<k):
            print("Invalid length")
            sys.exit()
        count=count+1
    return sqnc

rs = allp(3)  
sq=[]
dac_c=[]
phy=[]
for i in prop:
    if i=='all':
        for i in property_list:
            phy.append(i)
    else:
        phy.append(i)
temp = defaultdict(list)
for p in phy:
    for k in rs:
        count=-1
        for i in data['Physicochemical properties']:
            count=count+1
            if(i==p):
                temp[k].append(data[k][count])
            else:
                continue
            
# For taking input 
filename, file_extension = os.path.splitext(f1)
cdk = pd.DataFrame()
# if input is provided using cmd itself
if(file_extension==""):
    f1=f1.upper()
    alphabet=['A','C','G','T']
    for i in f1:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    temp1 = sp(f1,split,3)
    for k_c,k in enumerate(temp1):
        sq.append(k)
        dac_c.append(k_c+1)
    cdk['Sequence'] = [f1]
# if input sequences are to be read by input file   
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
                temp1 = sp(s,split,3)
                for k_c,k in enumerate(temp1):
                    sq.append(k)
                    dac_c.append(k_c)
                s=""
                
            else:
                continue
        else:
            for j in i:
                j=j.capitalize()
                if(j in ['A','G','C','T']):
                    s = s+j
    if s!="":
        temp1 = sp(s,split,3)
        for k_c,k in enumerate(temp1):
            sq.append(k)
            dac_c.append(k_c)
        
    cdk['Sequence_ID'] =s_id

# calculating the DAC value 
tac= defaultdict(list)
for k_c,seq in enumerate(sq):
    pu= []
    for i in range(len(phy)):
        mean=0.0
        for j in temp.keys():
            mean += temp[j][i]
        mean = mean/len(seq)
        pu.append(mean)  
    mer=kmer(3,seq)   
    count=0
    for l in range(1,lag+1):
        for p in range(len(phy)):
            if((len(seq)-l-2)==0):
                print("Invalid Sequence length")
                exit()
            av=0.0
            for i in range(len(seq)-l-2):
                av += temp[mer[i]][p]
            av= av/len(seq)
            count = count+1
            sum=0.0
            st = "TACC_SPLIT_s"+str(dac_c[k_c]+1)+"_"+ phy[p]+"_lag"+str(l)
            for i in range(len(seq)-l-2):
                pu1 = temp[mer[i]][p]-av
                pu2 = temp[mer[i+l]][p]
                sum += (pu1*pu2)
            sum = sum/(len(seq)-l-2)
            tac[st].append(sum)
for i in tac.keys():
    cdk[i]= tac[i]
tcc= defaultdict(list)
for k_c,seq in enumerate(sq):
    pu= []
    for i in range(len(phy)):
        mean=0.0
        for j in temp.keys():
            mean += temp[j][i]
        mean = mean/len(seq)
        pu.append(mean)  
    mer=kmer(3,seq)   
    count1=count
    for l in range(1,lag+1):
        for p1 in range(len(phy)):
            for p2 in range(len(phy)):
                if((len(seq)-l-2)==0):
                    print("Invalid sequence length")
                    exit()
                if(p1!=p2):
                    av1=0.0
                    av2=0.0
                    for i in range(len(seq)-l-2):
                        av1 += temp[mer[i]][p1]
                        av2 += temp[mer[i]][p2]
                    av1= av1/len(seq)
                    av2=av2/len(seq)
                    count1 = count1+1
                    sum=0.0
                    st = "TACC_SPLIT_s"+str(dac_c[k_c]+1)+"_"+ phy[p1]+"_"+phy[p2]+"_lag"+str(l)
                    for i in range(len(seq)-l-2):
                        pu1 = temp[mer[i]][p1]-av1
                        pu2 = temp[mer[i+l]][p2]-av2
                        sum += (pu1*pu2)
                    sum = sum/(len(seq)-l-2)
                    tcc[st].append(sum)
for i in tcc.keys():
    cdk[i]= tcc[i]
cdk.to_csv(out,index =False)

