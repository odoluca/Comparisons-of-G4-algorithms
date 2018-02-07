import os
import subprocess
import time, math
import re


def G4HScoreMod(seq,minRepeat=2,penalizeGC=True):
    i=0
    # baseScore=[0]*len(seq)
    baseScore=[]
    while i<len(seq):
        tractScore=[0]
        k=1
        GTract=False
        # print(i)
        while seq[i]=="G":
            tractScore=[(min(k-minRepeat+1,4))]*k #derivation from original algorithm: tractScore=[min(k-1,16)]*k
            # region derivation from original algorithm: if prev is "C" apply bigger penalty. penalizes GCs
            if penalizeGC:
                try:
                    if seq[i-k]=="C":baseScore[-1]=-2
                except:
                    pass
            # endregion
            k+=1
            i+=1
            GTract=True
            if i==len(seq): break
        if not GTract:
            while seq[i]=="C":
                tractScore=[max(-k+minRepeat-1,-4)]*k #derivation from original algorithm: tractScore=[max(-k,-16)]*k
                # region derivation from original algorithm: if prev is "G" apply bigger penalty. penalizes GCs
                if penalizeGC:
                    try:
                        if seq[i - k] == "G": baseScore[-1] = 2
                    except:
                        pass
                # endregion
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

# print G4HScore("GGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG")

def G4HScore(seq):
    i = 0
    # baseScore=[0]*len(seq)
    baseScore = []
    while i < len(seq):
        tractScore = [0]
        k = 1
        GTract = False
        # print(i)
        while seq[i] == "G":
            tractScore = [(min(k ,4))] * k  # derivation from original algorithm: tractScore=[min(k-1,16)]*k
            k += 1
            i += 1
            GTract = True
            if i == len(seq): break
        if not GTract:
            while seq[i] == "C":
                tractScore = [max(-k , -4)] * k  # derivation from original algorithm: tractScore=[max(-k,-16)]*k
                k += 1
                i += 1
                GTract = True
                if i == len(seq): break
        baseScore = baseScore.__add__(tractScore)
        if not GTract: i += 1
        # print baseScore
    Score = 0
    for value in baseScore:
        Score += value
    return float(Score) / len(seq)


file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
file="testedG4s3.fa"
# file="empty.fa"
# file="testedG4s_UnNvsG.fa"


#region quadparser tests
quadparserCommand = 'python quadparserModified.v3.py ' #shell command for quadparser command. Can add new regex pAUTtern: -r "\w{1}([gG]{3}\w{1,7}){3,}[gG]{3}\w{1}" # status for reference dataset: MCC:0.249, precision 98.6%
quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{2,}\w{1,7}){3,}[gG]{2,}"'# status for reference dataset: MCC:0.667, precision 94.2%
# quadparserCommand = 'python ImGQfinder.v2.py ' #shell command for perfect, buldged or mismAUTched sequences
#
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+ "'   #status for sentetic dataset: Finds all bulged: OK
#                                                                                                                                                                                                                                                                                                             # status for reference dataset: MCC:0.670, precision 99.2%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+ "' # status for reference dataset: MCC:0.670, precision 99.2%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "\w( [G]{3,} | (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))(\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "\w( [G]{3,} | (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))(\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})) )) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})) ))(\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUT][G]{2,}|[G]{2,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.693, precision 99.6%

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "\w( [G]{3,} | (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))(\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+" ' #status for reference dataset: MCC:0.734, precision 93.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.734, precision 93.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.750, precision 94.3%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.780, precision 94.1%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,12}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.780, precision 94.1%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,})) ))+" ' #status for reference dataset: MCC:0.787, precision 94.4%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,11}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,}))))  (\w{1,12}  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,})) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))  ( \w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))  (\w{1,12}  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))+" ' #status for reference dataset: MCC:0.817, precision 95.3%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))  ( \w{1,20}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))  (\w{1,20}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))  (\w{1,20}  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))+" ' #status for reference dataset: MCC:0.821, precision 94.8%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( (?P<edge>[G]{3,})|[G]{2,})  ( \w{1,9}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,9}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))+  (\w{1,9}(?(edge)[G]{2,}|[G]{3,}))" ' #status for reference dataset: MCC:0.809, precision 95.0%


# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  (\w{1,9}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,9}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,9}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.741, precision 92.5%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.747, precision 92.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.747, precision 92.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+" ' #status for reference dataset: MCC:0.747, precision 92.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+" ' #status for reference dataset: MCC:0.747, precision 92.6%
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))(\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+" ' #status for reference dataset: MCC:0.531, precision 92.3%
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,})) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) )) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))(\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.531, precision 92.3%
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w([G]{3,}?| (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,7}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )){3,}\w" '  # only perfect and buldged BEING TESTED
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w([G]{3,}?| (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,5}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )){3,}\w" '  # only perfect and buldged BEING TESTED
# quadparserCommand = 'python ImGQfinder.v2.py -r "([G]{3,}?| (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,5}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )){3,}" '  # status for reference dataset: MCC:0.462, precision 98.8%
quadparserCommand = 'python G4-predictorModified.py -w 30 -s 1.4 '  # status for reference dataset: MCC:0.076, precision 100%
quadparserCommand = 'python G4-predictorModified.py -w 25 -s 1.4 '  # status for reference dataset: MCC:0.116, precision 100%
quadparserCommand = 'python G4-predictorModified.py -w 15 -s 1.5 '  # status for reference dataset: MCC:0.349, precision 98.3%
quadparserCommand = 'python G4-predictorModified.py  -w 20 -s 1.7 '  # status for reference dataset: MCC:0.123, precision 92.1%
quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.4 '  # status for reference dataset: MCC:0.179, precision 100%
quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.2 '  # status for reference dataset: MCC:0.143, precision 100%

quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.0 '  # status for reference dataset: MCC:0.813, precision 96.875%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.2 '  # status for reference dataset: MCC:0.744, precision 97.8%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.4 '  # status for reference dataset: MCC:0.552, precision 98.1%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.5 '  # status for reference dataset: MCC:0.305, precision 96.5%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.6 '  # status for reference dataset: MCC:0.252, precision 95.6%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.7 '  # status for reference dataset: MCC:0.242, precision 98.5%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.0 '  # status for reference dataset: MCC:0.764, precision 94.6%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.2 '  # status for reference dataset: MCC:0.798, precision 97.2%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.4 '  # status for reference dataset: MCC:0.651, precision 97.9%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.5 '  # status for reference dataset: MCC:0.575, precision 97.7%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.6 '  # status for reference dataset: MCC:0.494, precision 97.9%
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.7 '  # status for reference dataset: MCC:0.333, precision 96.8%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1 '  # status for reference dataset: MCC:

# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.4 '  # status for reference dataset: MCC:0.610, precision 96.3%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.5 '  # status for reference dataset: MCC:0.510, precision 96.6%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.6 '  # status for reference dataset: MCC:0.439, precision 96.2%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.7 '  # status for reference dataset: MCC:0.317, precision 94.8%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.4 '  # status for reference dataset: MCC:0.162, precision 100%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.5 '  # status for reference dataset: MCC:0.140, precision 100%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.6 '  # status for reference dataset: MCC:0.130, precision 100%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.7 '  # status for reference dataset: MCC:0.120, precision 100%

quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.7 '  # status for reference dataset: MCC:0.764, precision 94.6%

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))" ' #status for reference dataset: MCC:0.864
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,}))) ))" ' #status for reference dataset: MCC:0.864
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))" ' #status for reference dataset: MCC:0.827
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,}))) ))" ' #status for reference dataset: MCC:0.827
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,}))) ))" ' #status for reference dataset: MCC:0.724
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{1,}[AUTC][G]{1,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{1,}[AUTC][G]{1,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{1,}[AUTC][G]{1,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{1,}[AUTC][G]{1,}))) ))" ' #status for reference dataset: MCC:0.820
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{1,}[AUT][G]{1,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,}))) ))" ' #status for reference dataset: MCC:0.820
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{2,}|[G]{2,}[AUT][G]{1,}))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (( [G]{2,}|(?P<mis>([G]{1,}[AUT][G]{2,}|[G]{2,}[AUT][G]{1,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{2,}|[G]{2,}[AUT][G]{1,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{2,}|[G]{2,}[AUT][G]{1,}))) ))" ' #status for reference dataset: MCC:0.864 DOESNT QUESTION <mis>!!!
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))+" ' #status for reference dataset: MCC:0.864 precision 95.1%, TPR:0.99, FPR:0.16
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{1,}[AUT][G]{1,}|[G]{1,}[AUT][G]{1,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,}|[G]{1,}[AUT][G]{1,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,}|[G]{1,}[AUT][G]{1,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,}|[G]{1,}[AUT][G]{1,}))) ))+" ' #status for reference dataset: MCC:0.842 precision 94.2%, TPR:0.99, FPR:0.19
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{2,7})|(?P<extreme>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(extreme)\w{2,7}|(?P<extreme>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(extreme)\w{2,7}|(?P<extreme>\w{1,30}) )   (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))" ' #status for reference dataset: MCC:0.829 precision 95.1% #I LIKE IT AND ITS LOGIC
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( ((?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  (((?P<lloop>\w{1,30}) )   (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))" ' #status for reference dataset: MCC:0.857

#
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,20}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,20})) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,20})) )   (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))" ' #status for reference dataset: MCC:0.837
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{1,5})|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{3,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(lloop)\w{1,5}|(\w{1,5}|(?P<lloop>\w{1,30})) )  (?(mis)([G]{3,})| ([G]{3,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(lloop)\w{1,5}|(\w{1,5}|(?P<lloop>\w{1,30})) )   (?(mis)([G]{3,})| ([G]{3,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))" ' #status for reference dataset: MCC:0.721 precision: 99.6% #VERY HIGH PRECISION

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{2,7})|(?P<extreme>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(extreme)\w{2,7}|(?P<extreme>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))))) ( (?(extreme)\w{2,7}|(?P<extreme>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))" ' #status for reference dataset: MCC:0.829 precision 95.1% #I LIKE IT AND ITS LOGIC

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))(((\w{1,7})|(?P<lloop>\w{1,20}))(?(mis)([G]{3,})|([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))((?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,20})))(?(mis)([G]{3,})|([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))((?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,20})))(?(mis)([G]{3,})|([G]{2,}|(?P<mis>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))"' #status for reference dataset: MCC:0.837 precision 95.4% #I LIKE IT AND ITS LOGIC
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))(((\w{1,7})|(?P<ExtremeLoop>\w{1,20}))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,20})))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,20})))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))"' #status for reference dataset: MCC:0.837 precision 95.4% #I LIKE IT AND ITS LOGIC

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,}|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))   ((\w{1,7})|(?P<ExtremeLoop>\w{1,20}))    (?(imperfectTract)([G]{3,})|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))   (?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))  (?(imperfectTract)([G]{3,})|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))  (?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))  (?(imperfectTract)([G]{2,})|(?P<imperfectTract>([G]{3,}|[G]+[AUTC][G]+)))"'

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))(((\w{2,7})|(?P<ExtremeLoop>\w{1,20}))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))((?(ExtremeLoop)\w{2,7}|(\w{2,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))((?(ExtremeLoop)\w{2,7}|(\w{2,7}|(?P<ExtremeLoop>\w{1,20})))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))"' #status for reference dataset: MCC:0.849 precision 95.1% #has too short limitation as well for loops
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]{3,}|[G]{3,}[AUTC][G]{2,})))(((\w{0,5})|(?P<ExtremeLoop>\w{1,35}))(?(imperfectTract)([G]{2,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]{3,}|[G]{3,}[AUTC][G]{2,})))))((?(ExtremeLoop)\w{2,7}|(\w{0,5}|(?P<ExtremeLoop>\w{1,35})))(?(imperfectTract)([G]{2,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]{3,}|[G]{3,}[AUTC][G]{2,})))))((?(ExtremeLoop)\w{0,5}|(\w{0,5}|(?P<ExtremeLoop>\w{1,35})))(?(imperfectTract)([G]{2,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]{3,}|[G]{3,}[AUTC][G]{2,})))))"' #status for reference dataset: MCC:0.850 #FOUND BY genetic algorithm

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))(((\w{1,7})|(?P<ExtremeLoop>\w{1,30}))(?(ImperfectTract)[G]{2,}|([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))))((?(ExtremeLoop)\w{1,7}|(?P<ExtremeLoop>\w{1,30}))(?(ImperfectTract)[G]{2,}|([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))))((?(ExtremeLoop)\w{1,7}|(?P<ExtremeLoop>\w{1,30}))(?(ImperfectTract)[G]{2,}|([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))))"' #status for reference dataset: MCC:0.820
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))(((\w{1,7})|(?P<ExtremeLoop>\w{1,30}))(?(ImperfectTract)[G]{3,}|([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))))((?(ExtremeLoop)\w{1,7}|(?P<ExtremeLoop>\w{1,30}))(?(ImperfectTract)[G]{3,}|([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))))((?(ExtremeLoop)\w{1,7}|(?P<ExtremeLoop>\w{1,30}))(?(ImperfectTract)[G]{3,}|([G]{2,}|(?P<ImperfectTract>([G]{1,}[AUTC][G]{1,})))))"' #status for reference dataset: MCC:0.842 precision 94.2% >G4Hs>0.473> MCC:0.888, precision 97.3%, TPR:0.97, FPR: 0.09

#Double allowance to imperfect tracts single allowance to extereme loop:
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis1>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))  ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis1)([G]{3,}|(?P<mis2>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,}))| ([G]{3,}|(?P<mis1>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis1)(?(mis2)[G]{3,}|([G]{3,} |(?P<mis2>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))| ([G]{3,}|(?P<mis1>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis1)(?(mis2)[G]{3,}|(?P<mis2>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))| ([G]{2,}|(?P<mis1>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,}))) ))" ' #status for reference dataset: MCC:0.811 precision: 98.2%

# quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{2,}\w{1,7}){3,}[gG]{2,}"'# status for reference dataset: MCC:0.667, precision 94.2% >G4Hs>0.473> MCC:0.735, precision 97.7%, TPR:0.87, FPR: 0.06
# quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{3,}\w{1,12}){3,}[gG]{3,}"'# status for reference dataset: MCC:0.667, precision 94.2%
# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.2 '  #

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<imperfectTract>([G]{3,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{3,})))(((\w{0,4})|(?P<ExtremeLoop>\w{1,34}))(?(imperfectTract)([G]{1,})|([G]{1,}|(?P<imperfectTract>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))((?(ExtremeLoop)\w{0,4}|(\w{0,4}|(?P<ExtremeLoop>\w{1,34})))(?(imperfectTract)([G]{1,})|([G]{2,}|(?P<imperfectTract>([G]{1,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))((?(ExtremeLoop)\w{0,4}|(\w{0,4}|(?P<ExtremeLoop>\w{1,34})))(?(imperfectTract)([G]{3,})|([G]{2,}|(?P<imperfectTract>([G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})))))"'

# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1 '  # status for reference dataset: MCC:0.764, precision 94.6%, TPR:0.94, FPR: 0.17
#
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1 '  # status for reference dataset: MCC:0.812, precision 96.9%, TPR:0.94, FPR: 0.10

#minG:2, bulged, max imperfect:1
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<imperfectTract>([G]+[AUTC][G]+)))(((\w{1,7})|(?P<ExtremeLoop>\w{1,30}))(?(imperfectTract)([G]{2,})|([G]{2,}|(?P<imperfectTract>([G]+[AUTC][G]+)))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{2,})|([G]{2,}|(?P<imperfectTract>([G]+[AUTC][G]+)))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{2,})|([G]{2,}|(?P<imperfectTract>([G]+[AUTC][G]+)))))"' #status for reference dataset: MCC:0.820, precision 93.3%, TPR:0.99, FPR: 0.22 >G4Hs>0.473> MCC:0.867, precision 96.7%, TPR:0.97, FPR: 0.11 >G4Hs>0.6> MCC:0.34, precision 97.6%, TPR:0.94, FPR: 0.07
#
# #minG: 3, mismatched, max imperfect:1
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,}|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))(((\w{1,7})|(?P<ExtremeLoop>\w{1,30}))(?(imperfectTract)([G]{3,})|([G]{3,}|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{3,})|([G]{3,}|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{3,})|([G]{3,}|(?P<imperfectTract>([G]{2,}|[G]+[AUTC][G]+)))))"' #status for reference dataset: MCC:0.787, precision 98.9%, TPR:0.89, FPR: 0.03 >G4Hs>0.473> no change
#
# #minG: 3, buldged, max imperfect:1, lowest FPR
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))(((\w{1,7})|(?P<ExtremeLoop>\w{1,30}))(?(imperfectTract)([G]{3,})|([G]{3,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{3,})|([G]{3,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))((?(ExtremeLoop)\w{1,7}|(\w{1,7}|(?P<ExtremeLoop>\w{1,30})))(?(imperfectTract)([G]{3,})|([G]{3,}|(?P<imperfectTract>([G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))))"' #status for reference dataset: MCC:0.736, precision 99.6%, TPR:0.84, FPR: 0.01 >G4Hs>0.473> MCC:0.736, precision 99.6%, TPR:0.84, FPR: 0.01

# minG: 2, bulged only, max imperfect:2
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<im1>[G]+[AUTC][G]+)) (\w{1,7}|(?P<extL>\w{1,30})) (?(im1)([G]{2,}|(?P<im2>[G]+[AUTC][G]+))|([G]{2,}|(?P<im1>[G]+[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{2,}|([G]{2,}|(?P<im2>[G]+[AUTC][G]+)))|([G]{2,}|(?P<im1>[G]+[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{2,}|([G]{2,}|(?P<im2>[G]+[AUTC][G]+)))|([G]{2,}|(?P<im1>[G]+[AUTC][G]+)))  "' #status for reference dataset: MCC:0.775, precision 91.6%, TPR:0.99, FPR: 0.29 >G4Hs>0.473> MCC:0.852, precision 96.0%, TPR:0.97, FPR: 0.13 >G4Hs>0.53> MCC:0.856, precision 97.3%, TPR:0.96, FPR: 0.09

# minG: 2, mismatched, max imperfect:2
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<im1>[G]+[AUTC]|[AUTC][G]+)) (\w{1,7}|(?P<extL>\w{1,30})) (?(im1)([G]{2,}|(?P<im2>[G]+[AUTC]|[AUTC][G]+))|([G]{2,}|(?P<im1>[G]+[AUTC]|[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{2,}|([G]{2,}|(?P<im2>[G]+[AUTC]|[AUTC][G]+)))|([G]{2,}|(?P<im1>[G]+[AUTC]|[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{2,}|([G]{2,}|(?P<im2>[G]+[AUTC]|[AUTC][G]+)))|([G]{2,}|(?P<im1>[G]+[AUTC]|[AUTC][G]+)))  "' #status for reference dataset: MCC:0.764, precision 90.0%, TPR:1.0, FPR: 0.35 >G4Hs>0.473> MCC:0.830, precision 95.4%, TPR:0.97, FPR: 0.14

# # minG: 3, bulged only, max imperfect:2
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,}|(?P<im1>[G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})) (\w{1,7}|(?P<extL>\w{1,30})) (?(im1)([G]{3,}|(?P<im2>[G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,}))|([G]{3,}|(?P<im1>[G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,}))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{3,}|([G]{3,}|(?P<im2>[G]{3,}[AUTC][G]+|[G]+[AUTC][G]{2,})))|([G]{3,}|(?P<im1>[G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,}))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{3,}|([G]{3,}|(?P<im2>[G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))|([G]{3,}|(?P<im1>[G]{2,}[AUTC][G]+|[G]+[AUTC][G]{2,})))  "' # status for reference dataset: MCC:0.754, precision 99.2%, TPR:0.86, FPR: 0.02 >G4Hs>0.473> MCC:0.754, precision 99.2%, TPR:0.86, FPR: 0.02

# minG: 3, mismatched, max imperfect:2
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,}|(?P<im1>[G]{2,}|[G]+[AUTC]|[AUTC][G]+)) (\w{1,7}|(?P<extL>\w{1,30})) (?(im1)([G]{3,}|(?P<im2>[G]{2,}|[G]+[AUTC]|[AUTC][G]+))|([G]{3,}|(?P<im1>[G]{2,}|[G]+[AUTC]|[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{3,}|([G]{3,}|(?P<im2>[G]{2,}|[G]+[AUTC]|[AUTC][G]+)))|([G]{3,}|(?P<im1>[G]{2,}|[G]+[AUTC]|[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,30}))) (?(im1)(?(im2)[G]{3,}|([G]{3,}|(?P<im2>[G]{2,}|[G]+[AUTC]|[AUTC][G]+)))|([G]{3,}|(?P<im1>[G]{2,}|[G]+[AUTC]|[AUTC][G]+)))  "'# status for reference dataset: MCC:0.809, precision 96.2%, TPR:0.94, FPR: 0.12 >G4Hs>0.473> MCC:0.819, precision 96.9%, TPR:0.94, FPR: 0.10 >G4Hs>0.473>

# TESTING
# minG: 2, bulged only, max imperfect:2, penalize imperfect
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{2,}|(?P<im1>[G]+[AUTC][G]+)) (\w{1,7}|(?P<extL>\w{1,20})) (?(im1)([G]{3,}|(?P<im2>[G]+[AUTC][G]+))|([G]{2,}|(?P<im1>[G]+[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,20}))) (?(im1)(?(im2)[G]{3,}|([G]{3,}|(?P<im2>[G]+[AUTC][G]+)))|([G]{2,}|(?P<im1>[G]+[AUTC][G]+))) (?(extL)\w{1,7}|(\w{1,7}|(?P<extL>\w{1,20}))) (?(im1)(?(im2)[G]{3,}|([G]{3,}|(?P<im2>[G]+[AUTC][G]+)))|([G]{2,}|(?P<im1>[G]+[AUTC][G]+)))  "'

# TESTING
# minG:2, bulged only, max imperfect:1, penalize imperfect in all successive Gtracts
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{1,}[AUT][G]{1,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)([G]{3,})| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,}))) ))" ' #status for reference dataset: MCC:0.842 precision 94.2%, TPR:0.99, FPR:0.19 >G4Hs>0.473> MCC:0.89, precision:97.3%, TPR:0.98, FPR:0.09

# TESTING
# minG:2, bulged only, max imperfect:1, penalize imperfect once
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{1,}[AUT][G]{1,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)((?P<pen>[G]{3,}))| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)((?(pen)[G]{2,}|(?P<pen>[G]{3,})))| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)((?(pen)[G]{2,}|(?P<pen>[G]{3,})))| ([G]{2,}|(?P<mis>([G]{1,}[AUT][G]{1,}))) ))" ' #status for reference dataset: MCC:0.842 precision 94.2%, TPR:0.99, FPR:0.19 >G4Hs>0.473> MCC:0.89, precision:97.3%, TPR:0.98, FPR:0.09

# TESTING - BEST RESULT
# if Gtract==2, no bulge or mismatch allowed. if Gtract>=3, single bulge or mismatch allowed
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]+[AUTC][G]+))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]+[AUTC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]+[AUTC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]+[AUTC][G]+))))"' #status for reference dataset: MCC:0.864 precision 95.1%, TPR:0.99, FPR:0.16 >G4Hs>0.473> MCC:0.909, precision:98.0%, TPR:0.98, FPR:0.06
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G][AUTC][G]))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G][AUTC][G]))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G][AUTC][G]))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G][AUTC][G]))))"' #status for reference dataset: MCC:0.864 precision 95.1%, TPR:0.99, FPR:0.16 >G4Hs>0.473> MCC:0.909, precision:98.0%, TPR:0.98, FPR:0.06

# TESTING
# if G tract is short, then shorten the loop allowance. if Gtract is long, allow bulge or mismatch
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G][AUTC][G]))  ((?(shrt)\w{1,4}|\w{1,9})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G][AUTC][G]))))  ((?(shrt)\w{1,4}|\w{1,9})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G][AUTC][G]))))  ((?(shrt)\w{1,4}|\w{1,9})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G][AUTC][G]))))"' #status for reference dataset: MCC:0.856 precision 94.8%, TPR:0.99, FPR:0.17 >G4Hs>0.473> MCC:0.909, precision:98.0%, TPR:0.98, FPR:0.06


# TESTING
# min Gtract=2, no bulge or mismatch allowed.
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " [G]{2,} ((\w{1,7})|(?P<lloop>\w{1,30}))  [G]{2,}  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  [G]{2,}  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  [G]{2,}"' #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.97, FPR:0.06
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " [G]{2,} ((\w{1,9})|(?P<lloop>\w{1,30}))  [G]{2,}  (?(lloop)\w{1,9}|(?P<lloop>\w{1,30}))  [G]{2,}  (?(lloop)\w{1,9}|(?P<lloop>\w{1,30}))  [G]{2,}"' #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.97, FPR:0.06

# TESTING
# if Gtract==2, no bulge or mismatch allowed. if Gtract>=3, single bulge allowed
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"' #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.98, FPR:0.06
#
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"' #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.98, FPR:0.06
#
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"' #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.98, FPR:0.06
#
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"' #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.98, FPR:0.06


#G2s are strict but loops typical, G3 one bulge
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
#status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.98, FPR:0.06











#G2s: YES but typical ONE extreme loop, G3s: ONE extreme loop TWO bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+)|(?P<shrt>[G]{2}))  ((\w{1,7})|(?P<lloop>\w{1,30}))                     (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))               (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))               (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))))"'
#status for reference dataset: MCC:0.864 precision 95.1%, TPR:0.99, FPR:0.16 >G4Hs>0.473> MCC:0.909, precision:98.0%, TPR:0.98, FPR:0.06
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+)|(?P<shrt>[G]{2}))   ((\w{1,7})|(?P<lloop>\w{1,30}))                     (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2,}|[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]{2,}|[G]+[ATUC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))               (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2,}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2,}|[G]+[ATUC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))               (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2,}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2,}|[G]+[ATUC][G]+))))"'
# MCC:0.857 precision:94.8 TPR:0.987 FPR:0.170
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  ((\w{1,7})|(?P<lloop>\w{1,30}))                     (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]+[ATUC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))               (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]+))))  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))               (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]+))))"'
# MCC:0.850 precision:94.5 TPR:0.987 FPR:0.181


# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))   (?(shrt)\w{1,7}| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))))   (?(shrt)\w{1,7}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))              (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))))  (?(shrt)\w{1,7}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))          (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2,}|[G]{2,}[ATUC][G]+))))"'



# #G2s are strict, G3 no bulge/mismatch (logic is similar to quadparser classic)
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2}))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|[G]{3,})  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30})) (?(shrt)[G]{2,}|[G]{3,})  (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}))  (?(shrt)[G]{2,}|[G]{3,})"'
# #status for reference dataset: MCC:0.857 precision 95.1%, TPR:0.98, FPR:0.16 >G4Hs>0.473> MCC:0.903, precision:98.0%, TPR:0.98, FPR:0.06

#G2s has short loops only, G3 one bulge
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2})|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))  (?(shrt)\w{1,4}|(\w{1,9}|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}|(?(lloop)\w{1,9}|(\w{1,9}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}|(?(lloop)\w{1,9}|(\w{1,9}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
#status for reference dataset: MCC:0.560 precision 96.1%, TPR:0.74, FPR:0.10

# quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{3,}\w{1,7}){3,}[gG]{3,}"'# status for reference dataset: MCC:0.667, precision 94.2%





#FINAL SETS ARE ADDED HERE

# G2s: typical (adjustable) ONE extreme loop, G3s: ONE extreme loop ONE bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})|(?P<shrt>[G]{2}))  (?(shrt)((\w{1,4})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)(?(lloop)(\w{1,4})|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)(?(lloop)(\w{1,4})|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.872 precision:96.1 TPR:0.980 FPR:0.128

# G2s: typical (adjustable) ONE extreme loop, G3s: ONE extreme loop NO bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2}))  (?(shrt)((\w{1,4})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|([G]{3,}))  (?(shrt)(?(lloop)(\w{1,4})|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|([G]{3,}))   (?(shrt)(?(lloop)(\w{1,4})|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|([G]{3,})) "'
# MCC:0.852 precision:96.0 TPR:0.970 FPR:0.128

# G2s: typical (adjustable) ONE extreme loop, G3s: ONE extreme loop TWO bulge only
# COMMENT: technically bulge only should be unnecessary when G2s are allowed, however, this version finds different scores so better to keep it.
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)|(?P<shrt>[G]{2}))   (?(shrt)((\w{1,7})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)))) (?(shrt) (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+))))  (?(shrt) (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+))))"'
# MCC:0.864 precision:95.1 TPR:0.987 FPR:0.160
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)|(?P<shrt>[G]{2}))   (?(shrt)((\w{1,4})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)))) (?(shrt) (?(lloop)\w{1,4}|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+))))  (?(shrt) (?(lloop)\w{1,4}|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]{2}|[G]{2}[ATUC][G]+))))"'
# MCC:0.879 precision:96.1 TPR:0.983 FPR:0.128


# G2s: typical (adjustable) ONE extreme loop, G3s: ONE extreme loop TWO bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))   (?(shrt)((\w{1,4})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))) (?(shrt) (?(lloop)\w{1,4}|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt) (?(lloop)\w{1,4}|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.850 precision:94.8 TPR:0.983 FPR:0.170

# G2s: typical (adjustable) ONE extreme loop, G3s: ONE extreme loop ONE bulge/mismatch
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)((\w{1,7})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)(?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)(?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.864 precision:95.1 TPR:0.987 FPR:0.160
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)((\w{1,4})|(?P<lloop>\w{1,30}))| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)(?(lloop)(\w{1,4})|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)(?(lloop)(\w{1,4})|(\w{1,4}|(?P<lloop>\w{1,30})))| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.879 precision:96.1 TPR:0.983 FPR:0.128


# G2s: NOT ALLOWED, G3s: NO extreme loop TWO bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))   (\w{1,7})   (?(imp1)([G]{3,}|(?P<imp2>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))   (\w{1,7})  (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))) )|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))   (\w{1,7})  (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))) )|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.659 precision:99.1 TPR:0.779 FPR:0.021

# G2s: NOT ALLOWED, G3s: NO extreme loop ONE bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>(([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (\w{1,7})      (?(imp)[G]{3,}  |([G]{3,}|(?P<imp>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))   (\w{1,7})   (?(imp)[G]{3,}  |([G]{3,}|(?P<imp>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))   (\w{1,7})  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.650 precision:99.6 TPR:0.762 FPR:0.011

# G2s: NOT ALLOWED, G3s: NO extreme loop NO bulge/mismatch
# COMMENT: technically identical to QP37
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,})   (\w{1,7})   ([G]{3,})   (\w{1,7})    ([G]{3,})   (\w{1,7})   ([G]{3,})"'
# MCC:0.553 precision:99.5 TPR:0.658 FPR:0.011

# G2s: NOT ALLOWED, G3s: NO extreme loop TWO bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))  (\w{1,7})    (?(imp1)([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))  (\w{1,7})  (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))  (\w{1,7})  (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))"'
# MCC:0.729 precision:96.7 TPR:0.883 FPR:0.096

# G2s: NOT ALLOWED, G3s: NO extreme loop ONE bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))  (\w{1,7})  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))) (\w{1,7})  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)))   (\w{1,7})  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)))"'
# MCC:0.720 precision:98.8 TPR:0.839 FPR:0.032

# G2s: NOT ALLOWED, G3s: ONE extreme loop TWO bulge only
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))   ((\w{1,7})|(?P<lloop>\w{1,30}))   (?(imp1)([G]{3,}|(?P<imp2>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))  (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))) )|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))) (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))) )|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.754 precision:99.2 TPR:0.859 FPR:0.021

# G2s: NOT ALLOWED, G3s: ONE extreme loop ONE bulge only
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>(([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))    ((\w{1,7})|(?P<lloop>\w{1,30}))      (?(imp)[G]{3,}  |([G]{3,}|(?P<imp>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))   (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))           (?(imp)[G]{3,}  |([G]{3,}|(?P<imp>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))   (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.737 precision:99.6 TPR:0.839 FPR:0.011

# G2s: NOT ALLOWED, G3s: ONE extreme loop NO imperfect tract
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,})   ((\w{1,7})|(?P<lloop>\w{1,30}))  ([G]{3,})   (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))  ([G]{3,})   (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))) ([G]{3,})"'
# MCC:0.671 precision:99.6 TPR:0.782 FPR:0.011
# COMMENT: very low FPR

# G2s: NOT ALLOWED, G3s: ONE extreme loop TWO bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))   ((\w{1,7})|(?P<lloop>\w{1,30}))   (?(imp1)([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))  (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))  (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))  (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))) (?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))"'
# MCC:0.805 precision:96.5 TPR:0.936 FPR:0.106


# G2s: NOT ALLOWED, G3s: ONE extreme loop ONE bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))  ((\w{1,7})|(?P<lloop>\w{1,30}))  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))) (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)))   (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))  (?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)))"'
# MCC:0.787 precision:98.9 TPR:0.889 FPR:0.032

# G2s: Typical (adjustable) loop, G3s: NO extreme loop NO bulge/mismatch
# COMMENT: technically with typical loop for G2 should be identical to original QP27 algorithm, and it is:
# quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{2,}\w{1,7}){3,}[gG]{2,}"'# status for reference dataset: MCC:0.667, precision 94.2%
# MCC:0.691 precision:95.3 TPR:0.879 FPR:0.138
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}|\w{1,7})   (?(shrt)[G]{2,}|[G]{3,})  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|[G]{3,})  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|[G]{3,})"'
# MCC:0.666 precision:96.5 TPR:0.836 FPR:0.096


# G2s: Typical (adjustable) loop, G3s: NO extreme loop TWO bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}|\w{1,7})   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))) (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}|\w{1,7}) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.776 precision:95.5 TPR:0.933 FPR:0.138


# G2s: Typical (adjustable) loop, G3s: NO extreme loop ONE bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.782 precision:96.5 TPR:0.923 FPR:0.106

# G2s: Typical (adjustable) loop, G3s: ONE extreme loop NO bulge/mismatch
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|[G]{3,})  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|[G]{3,})  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|[G]{3,})"'
# MCC:0.813 precision:96.9 TPR:0.936 FPR:0.096

# G2s: Typical (adjustable) loop, G3s: ONE extreme loop ONE bulge/mismatch
# COMMENT: Untested final version
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.860 precision:96.6 TPR:0.966 FPR:0.106

# G2s: Typical (adjustable) loop, G3s: ONE extreme loop TWO bulge/mismatch
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}| ((\w{1,7})|(?P<lloop>\w{1,30})))   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]+[ATUC][G]+))))               (?(shrt)\w{1,4}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]+[ATUC][G]+))))"'
# MCC:0.866 precision:96.3 TPR:0.973 FPR:0.117 >G4Hs>0.473> MCC:0.895, precision:97.3%, TPR:0.98, FPR:0.09
# COMMENT: Untested final version: mismatch is allowed at the edge aswell
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}| ((\w{1,7})|(?P<lloop>\w{1,30})))   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+))|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+)))) (?(shrt)\w{1,4}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))))  (?(shrt)\w{1,4}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}|[G]+[ATUC][G]+)) )|([G]{3,}|(?P<imp1>[G]{2}|[G]+[ATUC][G]+))))"'
# MCC:0.864 precision:95.4 TPR:0.983 FPR:0.149


# G2s: adjustable loop, G3s: ONE extreme loop ONE bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}| ((\w{1,7})|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}| (?(lloop)(\w{1,7})|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.842 precision:96.9 TPR:0.953 FPR:0.096

# G2s: adjustable loop, G3s: ONE extreme loop TWO bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}| ((\w{1,7})|(?P<lloop>\w{1,30})))   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))|([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))) (?(shrt)\w{1,4}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30}))))  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})) )|([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}| (?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})) )|([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.855 precision:96.9 TPR:0.960 FPR:0.096


# G2s: adjustable, G3s: NO extreme loop, ONE bulge
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp)[G]{3,}|([G]{3,}|(?P<imp>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
#                                                structure='([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))|(?P<shrt>[G]{2}))  (?(shrt)(\w{1,30})|(\w{1,7}))  (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))) (?(shrt)(\w{1,30})|(\w{1,7})) (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))) (?(shrt)(\w{1,30})|(\w{1,7})) (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))))'

# MCC:0.734 precision:96.7 TPR:0.886 FPR:0.096

# G2s: adjustable, G3s: NO extreme loop, TWO bulge
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}|\w{1,7})   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))|([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))) (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})) )|([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))  (?(shrt)\w{1,4}|\w{1,7}) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})) )|([G]{3,}|(?P<imp1>[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))"'
# MCC:0.749 precision:96.7 TPR:0.896 FPR:0.096
#endregion

#region set regex
G2sAllowed=False
ExtremeAllowed=True
ExtremeAllowedForG2s=False
ImperfectTractsAllowed=2
BulgedTractsOnly=True
if not G2sAllowed or not ExtremeAllowed:
    ExtremeAllowedForG2s = False
typLoopMax=str(2)
extLoopMax=str(31)
shrtLoopMax=str(2)

#region construct regex
bulgeOnly="[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}"
mismatch="[G]{2,}|[G]+[ATUC][G]+"
Dimp1='?P<imp1>'
Dimp2='?P<imp2>'
imp='('+mismatch+')'
if BulgedTractsOnly:
    imp='('+bulgeOnly+')'
Timp1='?(imp1)'
Timp2='?(imp2)'
shrt='\w{1,'+shrtLoopMax+'}'
Tshrt='?(shrt)'
Dshrt='?P<shrt>[G]{2}'
typ='\w{1,'+typLoopMax+'}'
ext='\w{1,'+extLoopMax+'}'
Dext='?P<lloop>'
Text='?(lloop)'

# Construct Tract 1:
Tract1='[G]{3,}'
if ImperfectTractsAllowed>0: Tract1=Tract1+'|('+Dimp1+imp+')'
if G2sAllowed: Tract1=Tract1+'|('+Dshrt+')'

# Construct Loop 1:
Loop1=typ
if ExtremeAllowed:Loop1=Loop1+'|('+Dext+ext+')'
if G2sAllowed:
    shrtAdd = shrt
    if ExtremeAllowedForG2s:
        shrtAdd = shrtAdd + '|(' + Dext + ext+')'
    Loop1 = Tshrt + '(' + shrtAdd + ')|('+Loop1+')'

# Construct Tract 2:
Tract2='[G]{3,}'
if ImperfectTractsAllowed>1: Tract2=Tract2+'|('+Dimp2+imp+')'
if ImperfectTractsAllowed>0: Tract2=Timp1+'('+Tract2+')|([G]{3,}|('+Dimp1+imp+'))'
if G2sAllowed:
    Tract2=Tshrt+'[G]{2,}|('+Tract2+')'

# Construct Loop 2:
Loop2=typ
if ExtremeAllowed:
    Loop2=Text+Loop2+'|('+typ+'|('+Dext+ext+'))'
if G2sAllowed:
    shrtAdd=shrt
    if ExtremeAllowedForG2s:
        shrtAdd = Text+ shrtAdd + '|('+shrt+'|(' + Dext + ext+'))'
    Loop2=Tshrt+'('+shrtAdd+')|('+Loop2+')'

# Construct Tract 3:
Tract3='[G]{3,}'
if ImperfectTractsAllowed>1:
    Tract3=Timp2+Tract3+'|([G]{3,}|('+Dimp2+imp+'))'
if ImperfectTractsAllowed>0:
    Tract3=Timp1+'('+Tract3+')|([G]{3,}|('+Dimp1+imp+'))'
if G2sAllowed:
    Tract3=Tshrt+'[G]{2,}|('+Tract3+')'

# Combine all regions:
structure='('+Tract1+')  ('+Loop1+')  ('+ Tract2+') ('+Loop2+') ('+Tract3+') ('+Loop2+') ('+Tract3+')'
#endregion
#endregion
#structure='([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))|(?P<shrt>[G]{2}))  (?(shrt)(\w{1,30})|(\w{1,7}))  (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))) (?(shrt)(\w{1,30})|(\w{1,7})) (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))) (?(shrt)(\w{1,30})|(\w{1,7})) (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))))'
quadparserCommand = 'python ImGQfinder.v2.py  -r "'
quadparserCommand=quadparserCommand[:40]+structure+'"'


from numpy import arange
parameter=None
print "parameter\tMCC\tTPR\tFPR"
# print(parameter)
# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1.2 '  # status for reference dataset: MCC:0.764, precision 94.6%

output= subprocess.check_output(quadparserCommand + ' -f "' + file+'"', shell=True)
# print output

G4HScoreTreshold=0#1.2#.473
TP=0
FP=0
G4List=[]
nonG4List=[]
for line in output.splitlines():
    if line.__contains__("not GQ_"):
        G4no=int(re.search(r"[0-9]+",line).group(0))
        if G4no not in nonG4List:
            score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
            if abs(score)>G4HScoreTreshold:
                nonG4List.append(G4no)
                FP+=1
                print line,score
    elif line.__contains__("GQ_"):
        G4no=int(re.search(r"[0-9]+",line).group(0))
        if G4no not in G4List:
            score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
            if abs(score)>G4HScoreTreshold:
                G4List.append(G4no)
                TP+=1
                print line,score

if TP==0 and FP==0: exit() #if nothing exists then exit without error.
FN=298-TP
TN=94-FP
# print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
# print "MCC:",MCC
precision=float(TP)/(TP+FP)
# print "precision:",precision*100,"%"
# print  "TPR:",float(TP)/298,"FPR:",float(FP)/94
print str(parameter)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94)
print "MCC:%.3f precision:%.1f TPR:%.3f FPR:%.3f" % (MCC, precision*100,float(TP)/298,float(FP)/94)