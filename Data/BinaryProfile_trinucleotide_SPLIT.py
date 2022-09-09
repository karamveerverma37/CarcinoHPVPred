# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 23:19:26 2020

@author: Megha Mathur
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 11:23:59 2020

@author: Megha Mathur
"""
import argparse
import warnings
import os
import pandas as pd
from collections import defaultdict
from collections import deque
import sys

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments to proceed')
## Read Arguments from command
parser.add_argument("-i", "--input", required=True, type=str, help="Input: protein or peptide sequence in FASTA format")
parser.add_argument("-o", "--output", type=str, help="Enter the output file name")
parser.add_argument("-s", "--split", required=True, type=int, help="Enter the split value")
# Parameter initialization or assigning variable for command level arguments
args = parser.parse_args()
f1 = args.input
if args.output == None:
    out = "outfile.csv"
else:
    out = args.output
seq = args.input
split = args.split
l_ar = [0 for i in range(split)]


def allp(k):
    n = ['A', 'T', 'C', 'G']
    s = []
    if (k == 1):
        return n
    elif (k == 2):
        for i in n:
            for j in n:
                se = i + j
                if (se not in s):
                    s.append(se)
        return s
    elif (k == 3):
        for i in n:
            for j in n:
                for k in n:
                    se = i + j + k
                    if (se not in s):
                        s.append(se)
        return s


def kmer(k, seq):
    s = []
    for i in range(len(seq)):
        se = ""
        if (i + k > len(seq)):
            break
        for j in range(k):
            se = se + seq[i + j]
        s.append(se)
    return s


def sp(seq, split, k, s_t, t, l_arr):
    sq1 = []
    for i in seq:
        sq1.append(i)
    if (split > len(seq)):
        sys.exit()
    sqnc = []
    while (split != 0):
        i = int(len(sq1) / split)
        s = ""
        for j in range(i):
            s = s + sq1[j]
            sq1[j] = 0
        sq1 = deque(sq1)
        for j in range(i):
            sq1.popleft()

        sqnc.append(s)
        split = split - 1
    count = 0
    for k, sq in enumerate(sqnc):
        s_t.append(t)
        if len(sq) > l_arr[k]:
            l_arr[k] = len(sq)
        if (len(sq) < k):
            print("Invalid length")
            sys.exit()
        count = count + 1
    return sqnc


dict_n = defaultdict(list)
seq = []
dac_c = []
s_temp = []
filename, file_extension = os.path.splitext(f1)
cdk = pd.DataFrame()
if (file_extension == ""):
    f1 = f1.upper()
    alphabet = ['A', 'C', 'G', 'T']
    for i in f1:
        if i not in alphabet:
            print("Invalid Character found in the given sequence")
            exit()
    temp1 = sp(f1, split, 3, s_temp, f1, l_ar)
    seq.append(temp1)
    #    for k_c,k in enumerate(temp1):
    #        seq.append(k)
    #        dac_c.append(k_c+1)
    cdk['Sequence'] = seq

else:
    f = open(f1, "r")
    b = f.readlines()
    s_id = []
    s = ""
    f.close()
    s_d = ''
    for i in b:
        if i[0] == '>':
            i = i.split("\n")
            s_d = i[0]
            s_id.append(i[0])
            if s != "":
                temp1 = sp(s, split, 3, s_temp, s_d, l_ar)
                seq.append(temp1)
                #                temp1 = sp(s,split,2,s_temp,s_d)
                #                for k_c,k in enumerate(temp1):
                #                    seq.append(k)
                #                    dac_c.append(k_c)
                s = ""

            else:
                continue
        else:
            for j in i:
                j = j.capitalize()
                if (j in ['A', 'G', 'C', 'T']):
                    s = s + j
    if s != "":
        temp1 = sp(s, split, 3, s_temp, s_d, l_ar)
        seq.append(temp1)
        cdk['Sequence'] = s_id
d1 = defaultdict(list)
for k, sqnc in enumerate(seq):
    for k_c, s in enumerate(sqnc):
        allpairs = allp(3)
        rs = kmer(3, s)
        res = []
        count = 0
        for i in rs:
            count = count + 1
            for j in allpairs:
                t = 'P' + str(count) + '_SPLIT_s' + str(k_c + 1) + '_' + j
                if j == i:
                    d1[t].append(str(1))
                else:
                    d1[t].append(str(0))
        if count != l_ar[k_c]-2:
            while count != l_ar[k_c]-2:
                count = count + 1
                for j in allpairs:
                    t = 'P' + str(count) + '_SPLIT_s' + str(k_c + 1) + '_' + j
                    d1[t].append(str(0))

for i in d1.keys():
    cdk[i] = d1[i]
cdk.to_csv(out, index=False)
