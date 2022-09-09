# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 12:42:37 2020

@author: Megha Mathur
"""
import argparse  
import warnings
import os
from argparse import RawTextHelpFormatter
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='NFDNA....Please provide following arguments to proceed',formatter_class=RawTextHelpFormatter) 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str,required=True,help="Input: Enter DNA sequence in FASTA format")
parser.add_argument("-o", "--output", type=str,help="Output File name")
parser.add_argument("-ft", "--feature",required=True,type=str,
                     help="Select among various features and include _NT,_CT,_REST,_SPLIT after feature name for N-Terminal,C-Terminal,Reamining of N-terminal and C-Terminal and for splitted sequence respectively.\n"
                    "\n"
                    "ALL_COMP : all composition based features\n"
	  "ALL_CORR : all correlation based features\n"
	  "ALL_BIN: all binary profile based features\n"
                    "CDK : kmer composition\n"
                    "RDK : Reverse Compliment kmer composition\n"                   
                    "DAC : Dinucleotide based auto correlation\n"
                    "DCC : Dinucleotide based cross correlation\n"
                    "DACC : Dinucleotide based auto cross correlation\n"
                    "TAC : Trinucleotide based auto correlation\n"
                    "TCC : Trinucleotide based cross correlation\n"
                    "TACC : Trinucleotide based auto cross correlation\n"
                    "PDNC :  Pseudo dinucleotide composition\n"
                    "PKNC :  Pseudo k-tuple nucleotide composition\n"
                    "PC_PDNC : parallel correlation pseudo dinucleotide composition\n"
                    "PC_PTNC : parallel correlation pseudo trinucleotide composition\n"
                    "SC_PDNC : Serial correlation pseudo dinucleotide composition\n"
                    "SC_PTNC : Serial correlation pseudo trinucleotide composition\n"
                    "NMBAC : Normalized Moreau–Broto autocorrelation\n"
                    "MAC :  Moran autocorrelation\n"
                    "GAC : Geary autocorrelation\n"
                    "NRI : Nucleotide repeat index\n"
                    "ES : Sequence level entropy of whole sequence\n"
                    "EN : Nucleotide level entropy of whole sequence\n"
                    "DDN : Distance Distribution of whole sequence\n"
                    "BPD : Binary profile dinucleotide\n"
                    "BPM : Binary profile Monotide\n"
                    "BPT : binary profile trinucleotide\n"
                    "BP_TP : Trinucleotide properties\n"
                    "BP_DP : Dinucleotide properties")
parser.add_argument("-k","--kvalue",type=int, help="Enter the k value and by default it is set to 2")
parser.add_argument("-n","--nvalue",type=int, help="Enter the n value")
parser.add_argument("-c","--cvalue",type=int, help="Enter the c value")
parser.add_argument("-s","--split",type=int, help="Enter the split value")
parser.add_argument("-or","--order",type=int, help="Enter the order value and by default it is set to 1")
parser.add_argument("-p", "--property",type=str,nargs='+',help=" Refer the property list of dipeptides and enter property names having space in between")
parser.add_argument("-l","--lagvalue",type=int, help="Enter the lag value and its default value is 2")
parser.add_argument("-w","--wvalue",type=float, help="Enter the w value and its default value is 0.05")
parser.add_argument("-lm","--lmvalue",type=int, help="Enter the lamada value and its default value is 1")
parser.add_argument("-path","--pythonpath",type=str, help="Enter the python environment path")


args = parser.parse_args()
sequence = args.input
feature = args.feature
feature=feature.upper()
f_list=['ALL_COMP','ALL_CORR','ALL_BIN','CDK','CDK_NT','RDK_NT','RDK_CT','RDK_REST','RDK_SPLIT','CDK_CT','CDK_REST','CDK_SPLIT','RDK','NRI','NRI_NT','NRI_CT','NRI_REST','NRI_SPLIT','DDN','DDN_NT','DDN_CT','DDN_REST','DDN_SPLIT','ES','ES_NT','ES_CT','ES_REST','ES_SPLIT','EN','EN_NT','EN_CT','EN_REST','EN_SPLIT','BPM','BPM_NT','BPM_CT','BPM_REST','BPM_SPLIT','BPD','BPD_NT','BPD_CT','BPD_REST','BPD_SPLIT','BPT','BPT_NT','BPT_CT','BPT_REST','BPT_SPLIT','BP_DP','BP_DP_NT','BP_DP_CT','BP_DP_REST','BP_DP_SPLIT','BP_TP_NT','BP_TP_CT','BP_TP_REST','BP_TP','BP_TP_SPLIT','DAC','DAC_NT','DAC_CT','DAC_REST','DAC_SPLIT','DACC','DACC_NT','DACC_CT','DACC_REST','DACC_SPLIT','DCC','DCC_NT','DCC_CT','DCC_REST','DCC_SPLIT','TAC','TAC_NT','TAC_CT','TAC_REST','TAC_SPLIT','TACC','TACC_NT','TACC_CT','TACC_REST','TACC_SPLIT','TCC','TCC_NT','TCC_CT','TCC_REST','TCC_SPLIT','MAC','MAC_NT','MAC_CT','MAC_REST','MAC_SPLIT','GAC','GAC_NT','GAC_CT','GAC_REST','GAC_SPLIT','NMBAC','NMBAC_NT','NMBAC_CT','NMBAC_REST','NMBAC_SPLIT','PDNC','PDNC_NT','PDNC_CT','PDNC_REST','PDNC_SPLIT','PKNC','PKNC_NT','PKNC_CT','PKNC_REST','PKNC_SPLIT','PC_PDNC','PC_PDNC_NT','PC_PDNC_CT','PC_PDNC_REST','PC_PDNC_SPLIT','PC_PTNC','PC_PTNC_NT','PC_PTNC_CT','PC_PTNC_REST','PC_PTNC_SPLIT','SC_PDNC','SC_PDNC_NT','SC_PDNC_CT','SC_PDNC_REST','SC_PDNC_SPLIT','SC_PTNC','SC_PTNC_NT','SC_PTNC_CT','SC_PTNC_REST','SC_PTNC_SPLIT']
if feature not in f_list:
    print("No Such Feature Exists")
    exit()   
if(args.output==None):
    out="Output.csv"
else:
    out=args.output
if(args.pythonpath==None):
    path='python'
else:
    path=args.pythonpath
# Calling each features
if feature == 'CDK':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    a = path+" Data/K_Mer.py "+"-i "+sequence+" -k "+k+" -or "+order+" -o "+out
    os.system(a)
    
if feature == 'RDK':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    a = path+" Data/RDK.py "+"-i "+sequence+" -k "+k+" -o "+out
    os.system(a)
    
if feature == 'ALL_COMP':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    file = open(out,'w')
    file.close()
    a = path+" Data/K_Mer.py "+"-i "+sequence+" -k "+k+" -or "+order+" -o "+out
    b= path+" Data/NRI.py "+"-i "+sequence+" -o "+out
    c = path+" Data/DDON.py "+"-i "+sequence+" -o "+out
    d= path+" Data/RDK.py "+"-i "+sequence+" -k "+k+" -o "+out
    e = path+" Data/ENT.py "+"-i "+sequence+" -o "+out
    f = path+" Data/ENT_NL.py "+"-i "+sequence+" -o "+out
    k1 = path+" Data/psednc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out
    l = path+" Data/pseknc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
    os.system(b)
    os.system(c)
    os.system(d)
    os.system(e)
    os.system(f)
    os.system(k1)
    os.system(l)
if feature == 'ALL_BIN':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    file = open(out,'w')
    file.close()
    a = path+" Data/BinaryProfile_monotide.py "+"-i "+sequence+" -o "+out
    b= path+" Data/BinaryProfile_dinucleotide.py "+"-i "+sequence+" -o "+out
    c = path+" Data/BinaryProfile_trinucleotide.py "+"-i "+sequence+" -o "+out
    d = path+" Data/dinucleotide_prop.py "+"-i "+sequence+" -or "+order+" -o "+out+" -p"   
    e = path+" Data/trinucleotide_prop.py "+"-i "+sequence+" -or "+order+" -o "+out+" -p"
    for i in prop:
        j='"'+i+'"'
        d=d+" "+j
        e=e+" "+j
    os.system(a)
    os.system(b)
    os.system(c)
    os.system(d)
    os.system(e)
if feature == 'ALL_CORR':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    file = open(out,'w')
    file.close()
    a = path+" Data/DAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    b = path+" Data/TAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    c = path+" Data/MAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -k "+k+" -p"
    d = path+" Data/GAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -k "+k+" -p"
    e = path+" Data/NMBAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -k "+k+" -p"
    f = path+" Data/DCC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    g = path+" Data/TCC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    h = path+" Data/DACC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    j1 = path+" Data/TACC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    m = path+" Data/pcpsednc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -o "+out+" -lm "+lm+" -p"
    n = path+" Data/pcpsetnc.py "+"-i "+sequence+" -w "+w+" -k "+str(3)+" -lm "+lm+" -o "+out+" -p" 
    q = path+" Data/scpsednc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    r = path+" Data/scpsetnc.py "+"-i "+sequence+" -w "+w+" -k "+str(3)+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
        b=b+" "+j
        c=c+" "+j
        d=d+" "+j
        e=e+" "+j
        f=f+" "+j
        g=g+" "+j
        h=h+" "+j
        j1=j1+" "+j
        m=m+" "+j
        n=n+" "+j
        q=q+" "+j
        r=r+" "+j 
    os.system(a)
    os.system(b)
    os.system(c)
    os.system(d)
    os.system(e)
    os.system(f)
    os.system(g)
    os.system(h)
    os.system(j1)
    os.system(m)
    os.system(n)
    os.system(q)
    os.system(r)
  
if feature == 'CDK_NT':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.nvalue==None:
        n = str(len(sequence))  
    else:
        n=str(args.nvalue)
    b = path+" Data/3_end.py "+"-i "+sequence+" -k "+k+" -n "+n+" -or "+order+" -o "+out
    os.system(b)
if feature == 'RDK_NT':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.nvalue==None:
        n = str(len(sequence))  
    else:
        n=str(args.nvalue)
    b = path+" Data/RDK_NT.py "+"-i "+sequence+" -k "+k+" -n "+n+" -or "+order+" -o "+out
    os.system(b)
if feature == 'RDK_CT':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.cvalue==None:
        n = str(len(sequence))  
    else:
        n=str(args.cvalue)
    b = path+" Data/RDK_CT.py "+"-i "+sequence+" -k "+k+" -c "+n+" -or "+order+" -o "+out
    os.system(b)
if feature == 'RDK_REST':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    b = path+" Data/RDK_REST.py "+"-i "+sequence+" -k "+k+" -n "+n+" -c "+m+" -or "+order+" -o "+out
    os.system(b)
if feature == 'RDK_SPLIT':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)    
    e = path+" Data/RDK_SPLIT.py "+"-i "+sequence+" -k "+k+" -s "+s+" -or "+order+" -o "+out
    os.system(e)
if feature == 'CDK_CT':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.cvalue==None:
        n = str(len(sequence))  
    else:
        n=str(args.cvalue)
    b = path+" Data/5_end.py "+"-i "+sequence+" -k "+k+" -c "+n+" -or "+order+" -o "+out
    os.system(b)
if feature == 'CDK_REST':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    b = path+" Data/rest.py "+"-i "+sequence+" -k "+k+" -n "+n+" -c "+m+" -or "+order+" -o "+out
    os.system(b)
if feature == 'CDK_SPLIT':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if args.order == None:
        order = str(1)
    else:
        order = str(args.order)    
    e = path+" Data/split.py "+"-i "+sequence+" -k "+k+" -s "+s+" -or "+order+" -o "+out
    os.system(e)
    
if feature == 'DAC':
    file = open(out,'w')
    file.close()
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DAC_NT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DAC_NT.py "+"-i "+sequence+" -n "+n+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DAC_CT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DAC_CT.py "+"-i "+sequence+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DAC_REST':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DAC_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DAC_SPLIT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DAC_SPLIT.py "+"-i "+sequence+" -s "+s+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DCC':
    file = open(out,'w')
    file.close()
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DCC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DCC_NT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DCC_NT.py "+"-i "+sequence+" -n "+n+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DCC_CT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.cvalue==None:
        m = str(2) 
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DCC_CT.py "+"-i "+sequence+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DCC_REST':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DCC_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DCC_SPLIT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DCC_SPLIT.py "+"-i "+sequence+" -s "+s+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DACC':
    file = open(out,'w')
    file.close()
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/DACC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DACC_NT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DACC_NT.py "+"-i "+sequence+" -n "+n+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DACC_CT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DACC_CT.py "+"-i "+sequence+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DACC_REST':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DACC_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'DACC_SPLIT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/DACC_SPLIT.py "+"-i "+sequence+" -s "+s+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TAC':
    file = open(out,'w')
    file.close()
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/TAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TAC_NT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TAC_NT.py "+"-i "+sequence+" -n "+n+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TAC_CT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TAC_CT.py "+"-i "+sequence+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TAC_REST':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TAC_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TAC_SPLIT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TAC_SPLIT.py "+"-i "+sequence+" -s "+s+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TACC':
    file = open(out,'w')
    file.close()
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/TACC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TACC_NT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TACC_NT.py "+"-i "+sequence+" -n "+n+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TACC_CT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TACC_CT.py "+"-i "+sequence+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TACC_REST':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TACC_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TACC_SPLIT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TACC_SPLIT.py "+"-i "+sequence+" -s "+s+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TCC':
    file = open(out,'w')
    file.close()
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/TCC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TCC_NT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TCC_NT.py "+"-i "+sequence+" -n "+n+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TCC_CT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TCC_CT.py "+"-i "+sequence+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TCC_REST':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TCC_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'TCC_SPLIT':
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop= args.property
    a = path+" Data/TCC_SPLIT.py "+"-i "+sequence+" -s "+s+" -l "+lag+" -o "+out+" -p "
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'GAC':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/GAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'GAC_NT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/GAC_NT.py "+"-i "+sequence+" -l "+lag+" -n "+n+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'GAC_CT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/GAC_CT.py "+"-i "+sequence+" -l "+lag+" -c "+m+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'GAC_REST':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/GAC_REST.py "+"-i "+sequence+" -l "+lag+" -n "+n+" -c "+m+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'GAC_SPLIT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/GAC_SPLIT.py "+"-i "+sequence+" -l "+lag+" -s "+s+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'MAC':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/MAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'MAC_NT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/MAC_NT.py "+"-i "+sequence+" -l "+lag+" -n "+n+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'MAC_CT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/MAC_CT.py "+"-i "+sequence+" -l "+lag+" -c "+m+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'MAC_REST':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/MAC_REST.py "+"-i "+sequence+" -l "+lag+" -n "+n+" -c "+m+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'MAC_SPLIT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/MAC_SPLIT.py "+"-i "+sequence+" -l "+lag+" -s "+s+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'NMBAC':
    file = open(out,'w')
    file.close()
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/NMBAC.py "+"-i "+sequence+" -l "+lag+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'NMBAC_NT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/NMBAC_NT.py "+"-i "+sequence+" -l "+lag+" -n "+n+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'NMBAC_CT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/NMBAC_CT.py "+"-i "+sequence+" -l "+lag+" -c "+m+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'NMBAC_REST':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/NMBAC_REST.py "+"-i "+sequence+" -l "+lag+" -n "+n+" -c "+m+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'NMBAC_SPLIT':
    if args.kvalue==None:
        k=str(2)
    else:
        k=str(args.kvalue)
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    if args.lagvalue == None:
        lag = str(2)
    else:
        lag = str(args.lagvalue)
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/NMBAC_SPLIT.py "+"-i "+sequence+" -l "+lag+" -s "+s+" -o "+out+" -k "+k+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'BPD':
    file = open(out,'w')
    file.close()
    a = path+" Data/BinaryProfile_dinucleotide.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature == 'BPD_NT':
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    a = path+" Data/BinaryProfile_dinucleotide_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature == 'BPD_CT':
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    a = path+" Data/BinaryProfile_dinucleotide_CT.py "+"-i "+sequence+" -c "+m+" -o "+out
    os.system(a)
if feature == 'BPD_REST':
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    a = path+" Data/BinaryProfile_dinucleotide_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature == 'BPD_SPLIT':
    if args.split==None:
        s = str(1)  
    else:
        s=str(args.split)
    a = path+" Data/BinaryProfile_dinucleotide_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
if feature=='BPM':
    file = open(out,'w')
    file.close()
    a=path+" Data/BinaryProfile_monotide.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature == 'BPM_NT':
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    a = path+" Data/BinaryProfile_mononucleotide_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature == 'BPM_CT':
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    a = path+" Data/BinaryProfile_mononucleotide_CT.py "+"-i "+sequence+" -c "+m+" -o "+out
    os.system(a)
if feature == 'BPM_REST':
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    a = path+" Data/BinaryProfile_mononucleotide_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature == 'BPM_SPLIT':
    if args.split==None:
        s = str(1)  
    else:
        s=str(args.split)
    a = path+" Data/BinaryProfile_mononucleotide_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
if feature == 'BPT':
    file = open(out,'w')
    file.close()
    a=path+" Data/BinaryProfile_trinucleotide.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature == 'BPT_NT':
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    a = path+" Data/BinaryProfile_trinucleotide_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature == 'BPT_CT':
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    a = path+" Data/BinaryProfile_trinucleotide_CT.py "+"-i "+sequence+" -c "+m+" -o "+out
    os.system(a)
if feature == 'BPT_REST':
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    a = path+" Data/BinaryProfile_trinucleotide_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature == 'BPT_SPLIT':
    if args.split==None:
        s = str(1)  
    else:
        s=str(args.split)
    a = path+" Data/BinaryProfile_trinucleotide_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
if feature == 'PDNC':
    file = open(out,'w')
    file.close()
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(3)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/psednc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PDNC_NT':
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(3)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/psednc_NT.py "+"-i "+sequence+" -w "+w+" -n "+n+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PDNC_CT':
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(3)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/psednc_CT.py "+"-i "+sequence+" -w "+w+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PDNC_REST':
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(3)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/psednc_REST.py "+"-i "+sequence+" -w "+w+" -n "+n+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PDNC_SPLIT':
    if args.split==None:
        s = str(1)  
    else:
        s=str(args.split)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(3)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/psednc_SPLIT.py "+"-i "+sequence+" -w "+w+" -s "+s+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PKNC':
    file = open(out,'w')
    file.close()
    if args.wvalue == None:
        w = str(0.5)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/pseknc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PKNC_NT':
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if args.wvalue == None:
        w = str(0.5)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/pseknc_NT.py "+"-i "+sequence+" -w "+w+" -n "+n+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PKNC_CT':
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if args.wvalue == None:
        w = str(0.5)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/pseknc_CT.py "+"-i "+sequence+" -w "+w+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PKNC_REST':
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if args.wvalue == None:
        w = str(0.5)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/pseknc_REST.py "+"-i "+sequence+" -w "+w+" -n "+n+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PKNC_SPLIT':
    if args.split==None:
        s = str(1)  
    else:
        s=str(args.split)
    if args.wvalue == None:
        w = str(0.5)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)        
    a = path+" Data/pseknc_SPLIT.py "+"-i "+sequence+" -w "+w+" -s "+s+" -k "+k+" -lm "+lm+" -o "+out
    os.system(a)
if feature == 'PC_PDNC':
    file = open(out,'w')
    file.close()
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsednc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -o "+out+" -lm "+lm+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'PC_PDNC_NT':
    if args.nvalue==None:
        n = str(2)  
    else:
        n=str(args.nvalue)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsednc_NT.py "+"-i "+sequence+" -w "+w+" -n "+n+" -k "+k+" -o "+out+" -lm "+lm+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'PC_PDNC_CT':
    if args.cvalue==None:
        m = str(2)  
    else:
        m=str(args.cvalue)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsednc_CT.py "+"-i "+sequence+" -w "+w+" -c "+m+" -k "+k+" -o "+out+" -lm "+lm+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature == 'PC_PDNC_REST':
    if args.nvalue==None:
        n = str(0)  
    else:
        n=str(args.nvalue)
    if args.cvalue==None:
        m = str(0)  
    else:
        m=str(args.cvalue)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsednc_REST.py "+"-i "+sequence+" -w "+w+" -n "+n+" -c "+m+" -k "+k+" -o "+out+" -lm "+lm+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature == 'PC_PDNC_SPLIT':
    if args.split==None:
        s = str(1)  
    else:
        s=str(args.split)
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsednc_SPLIT.py "+"-i "+sequence+" -w "+w+" -s "+s+" -k "+k+" -o "+out+" -lm "+lm+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'PC_PTNC':
    file = open(out,'w')
    file.close()
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsetnc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'PC_PTNC_NT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.nvalue == None:
        n = str(2)
    else:
        n = str(args.nvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsetnc_NT.py "+"-i "+sequence+" -w "+w+" -n "+n+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'PC_PTNC_CT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.cvalue == None:
        m = str(2)
    else:
        m = str(args.cvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsetnc_CT.py "+"-i "+sequence+" -w "+w+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'PC_PTNC_REST':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsetnc_REST.py "+"-i "+sequence+" -w "+w+" -n "+n+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature == 'PC_PTNC_SPLIT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.split == None:
        s = str(1)
    else:
        s = str(args.split)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/pcpsetnc_SPLIT.py "+"-i "+sequence+" -w "+w+" -s "+s+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PDNC':
    file = open(out,'w')
    file.close()
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsednc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PDNC_NT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.nvalue == None:
        n = str(2)
    else:
        n = str(args.nvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsednc_NT.py "+"-i "+sequence+" -w "+w+" -n "+n+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PDNC_CT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.cvalue == None:
        m = str(2)
    else:
        m = str(args.cvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsednc.py_CT "+"-i "+sequence+" -w "+w+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature == 'SC_PDNC_REST':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsednc_REST.py "+"-i "+sequence+" -w "+w+" -n "+n+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PDNC_SPLIT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.split == None:
        s = str(1)
    else:
        s = str(args.split)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(2)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsednc_SPLIT.py "+"-i "+sequence+" -w "+w+" -s "+s+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PTNC':
    file = open(out,'w')
    file.close()
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsetnc.py "+"-i "+sequence+" -w "+w+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PTNC_NT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.nvalue == None:
        n = str(2)
    else:
        n = str(args.nvalue)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsetnc_NT.py "+"-i "+sequence+" -w "+w+" -n "+n+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature == 'SC_PTNC_CT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.cvalue == None:
        m = str(2)
    else:
        m = str(args.cvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsetnc_CT.py "+"-i "+sequence+" -w "+w+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature == 'SC_PTNC_REST':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsetnc_REST.py "+"-i "+sequence+" -w "+w+" -n "+n+" -c "+m+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature == 'SC_PTNC_SPLIT':
    if args.wvalue == None:
        w = str(0.05)
    else:
        w = float(args.wvalue)
        w=str(w)
    if args.split == None:
        s = str(1)
    else:
        s = str(args.split)

    if args.lmvalue == None:
        lm = str(1)
    else:
        lm = int(args.lmvalue) 
        lm=str(lm)
    if args.kvalue == None:
        k = str(3)
    else:
        k = int(args.kvalue)
        k=str(k)  
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    a = path+" Data/scpsetnc_SPLIT.py "+"-i "+sequence+" -w "+w+" -s "+s+" -k "+k+" -lm "+lm+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature=='BP_DP':
    file = open(out,'w')
    file.close()
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/dinucleotide_prop.py "+"-i "+sequence+" -or "+order+" -o "+out+" -p"
    for i in prop:
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature=='BP_DP_NT':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.nvalue == None:
        n = str(2)
    else:
        n = str(args.nvalue)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/dinucleotide_prop_NT.py "+"-i "+sequence+" -n "+n+" -or "+order+" -o "+out+" -p"
    for i in prop:
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_DP_CT':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.cvalue == None:
        m = str(2)
    else:
        m = str(args.cvalue)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/dinucleotide_prop_CT.py "+"-i "+sequence+" -c "+m+" -or "+order+" -o "+out+" -p"
    for i in prop:
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_DP_REST':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/dinucleotide_prop_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -or "+order+" -o "+out+" -p"
    for i in prop:
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_DP_SPLIT':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.split == None:
        s = str(1)
    else:
        s = str(args.split)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/dinucleotide_prop_SPLIT.py "+"-i "+sequence+" -s "+s+" -or "+order+" -o "+out+" -p"
    for i in prop:
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)
if feature=='BP_TP':
    file = open(out,'w')
    file.close()
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/trinucleotide_prop.py "+"-i "+sequence+" -or "+order+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_TP_NT':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.nvalue == None:
        n = str(2)
    else:
        n = str(args.nvalue)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/trinucleotide_prop_NT.py "+"-i "+sequence+" -n "+n+" -or "+order+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_TP_CT':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.cvalue == None:
        m = str(2)
    else:
        m = str(args.cvalue)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/trinucleotide_prop_CT.py "+"-i "+sequence+" -c "+m+" -or "+order+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_TP_REST':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/trinucleotide_prop_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -or "+order+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a) 
if feature=='BP_TP_SPLIT':
    if(args.property==None):
        print("Invalid number of Arguments.Enter property name in the format -p propertyname")
        exit()
    else:
        prop=args.property
    if args.split == None:
        s = str(1)
    else:
        s = str(args.split)
    if args.order == None:
        order = str(1)
    else:
        order = int(args.order)
        order=str(order)
    a = path+" Data/trinucleotide_prop_SPLIT.py "+"-i "+sequence+" -s "+s+" -or "+order+" -o "+out+" -p"
    for i in prop:
        i=i.lower()
        j='"'+i+'"'
        a=a+" "+j
    os.system(a)

if feature=='NRI':
    file = open(out,'w')
    file.close()
    a = path+" Data/NRI.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature=='NRI_NT':
    if args.nvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.nvalue)
    a = path+" Data/NRI_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature=='NRI_CT':
    if args.cvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.cvalue)
    a = path+" Data/NRI_CT.py "+"-i "+sequence+" -c "+n+" -o "+out
    os.system(a)
if feature=='NRI_REST':
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    a = path+" Data/NRI_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature=='NRI_SPLIT':
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    a = path+" Data/NRI_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
if feature=='DDN':
    file = open(out,'w')
    file.close()
    a = path+" Data/DDON.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature=='DDN_NT':
    if args.nvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.nvalue)
    a = path+" Data/DDON_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature=='DDN_CT':
    if args.cvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.cvalue)
    a = path+" Data/DDON_CT.py "+"-i "+sequence+" -c "+n+" -o "+out
    os.system(a)
if feature=='DDN_REST':
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    a = path+" Data/DDON_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature=='DDN_SPLIT':
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    a = path+" Data/DDON_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
if feature=='ES':
    file = open(out,'w')
    file.close()
    a = path+" Data/ENT.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature=='ES_NT':
    if args.nvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.nvalue)
    a = path+" Data/ENT_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature=='ES_CT':
    if args.cvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.cvalue)
    a = path+" Data/ENT_CT.py "+"-i "+sequence+" -c "+n+" -o "+out
    os.system(a)
if feature=='ES_REST':
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    a = path+" Data/ENT_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature=='ES_SPLIT':
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    a = path+" Data/ENT_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
if feature=='EN':
    file = open(out,'w')
    file.close()
    a = path+" Data/ENT_NL.py "+"-i "+sequence+" -o "+out
    os.system(a)
if feature=='EN_NT':
    if args.nvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.nvalue)
    a = path+" Data/ENT_NL_NT.py "+"-i "+sequence+" -n "+n+" -o "+out
    os.system(a)
if feature=='EN_CT':
    if args.cvalue == None:
        n = str(len(sequence))
    else:
        n = str(args.cvalue)
    a = path+" Data/ENT_NL_CT.py "+"-i "+sequence+" -c "+n+" -o "+out
    os.system(a)
if feature=='EN_REST':
    if args.nvalue == None:
        n = str(0)
    else:
        n = str(args.nvalue)
    if args.cvalue == None:
        m = str(0)
    else:
        m = str(args.cvalue)
    a = path+" Data/ENT_NL_REST.py "+"-i "+sequence+" -n "+n+" -c "+m+" -o "+out
    os.system(a)
if feature=='EN_SPLIT':
    if args.split==None:
        s=str(1)
    else:
        s=str(args.split)
    a = path+" Data/ENT_NL_SPLIT.py "+"-i "+sequence+" -s "+s+" -o "+out
    os.system(a)
    
    
