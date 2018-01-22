import os
import subprocess
import time, math
import re


def G4HScore(seq):
    i=0
    # baseScore=[0]*len(seq)
    baseScore=[]
    while i<len(seq):
        tractScore=[0]
        k=1
        GTract=False
        # print(i)
        while seq[i]=="G":
            # region derivation from original algorithm: if prev is "C" apply bigger penalty. penalizes GCs
            try:
                if seq[i-k]=="C":baseScore[-1]=-2
            except:
                pass
            # endregion
            tractScore=[(min(k-1,4))]*k #derivation from original algorithm: tractScore=[min(k-1,16)]*k
            k+=1
            i+=1
            GTract=True
            if i==len(seq): break
        if not GTract:
            while seq[i]=="C":
                # region derivation from original algorithm: if prev is "G" apply bigger penalty. penalizes GCs
                try:
                    if seq[i - k] == "G": baseScore[-1] = 2
                except:
                    pass
                # endregion
                tractScore=[max(-k+1,-4)]*k #derivation from original algorithm: tractScore=[max(-k,-16)]*k
                k+=1
                i+=1
                GTract=True
                if i == len(seq): break
        baseScore=baseScore.__add__(tractScore)
        if not GTract: i += 1
        # print baseScore
    Score=0
    for value in baseScore:
        Score+=value
    return float(Score)/len(seq)

print G4HScore("GGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG")


file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
file="testedG4s_4.fa"
# file="test 12fa"
file="testedG4s3.fa"
# file="empty.fa"

quadparserCommand = 'python quadparserModified.v3.py ' #shell command for quadparser command. Can add new regex pAUTtern: -r "\w{1}([gG]{3}\w{1,7}){3,}[gG]{3}\w{1}" # status for reference dataset: MCC:0.249, precision 98.6%
quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{2,}\w{1,7}){3,}[gG]{2,}"'# status for reference dataset: MCC:0.667, precision 94.2%



quadparserCommand = 'python G4HunterModified.v2.py -w 20 -s 1 '  # status for reference dataset: MCC:0.764, precision 94.6%

starttime=time.time()
output= subprocess.check_output(quadparserCommand + ' -f "' + file+'"', shell=True)

from numpy import arange

print(time.time()-starttime)
for parameter in arange(1.0,2.1,0.1):
    print(parameter)
    quadparserCommand = 'python G4HunterModified.v2.py -w 20 -s '+str(parameter)  # status for reference dataset: MCC:0.764, precision 94.6%
    output = subprocess.check_output(quadparserCommand + ' -f "' + file + '"', shell=True)

    G4HScoreTreshold=0.#473
    TP=0
    FP=0
    G4List=[]
    nonG4List=[]
    for line in output.splitlines():
        if line.__contains__("nonGQ_"):
            G4no=int(re.search(r"[0-9]+",line).group(0))
            if G4no not in G4List:
                score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
                if abs(score)>G4HScoreTreshold:
                    G4List.append(G4no)
                    FP+=1
                    # print line,score
        elif line.__contains__("GQ_"):
            G4no=int(re.search(r"[0-9]+",line).group(0))
            if G4no not in nonG4List:
                score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
                if abs(score)>G4HScoreTreshold:
                    nonG4List.append(G4no)
                    TP+=1
                    # print line,score

    if TP==0 and FP==0: exit() #if nothing exists then exit without error.
    # FN=71-TP
    # TN=138-FP
    FN = 298 - TP
    TN = 94 - FP
    print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
    MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    print "MCC:",MCC
    precision=float(TP)/(TP+FP)
    print "precision:",precision*100,"%"
    print  "TPR:", float(TP) / 298, "FPR:", float(FP) / 94
