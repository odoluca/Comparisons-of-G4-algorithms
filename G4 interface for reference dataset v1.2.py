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
            tractScore=[(min(k,4))]*k
            k+=1
            i+=1
            GTract=True
            if i==len(seq): break
        if not GTract:
            while seq[i]=="C":
                tractScore=[max(-k,-4)]*k
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

# print G4HScore("GGGGGAGGGGGGGUUUGG")


file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
file="testedG4s3.fa"
# file="empty.fa"

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
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,}))))  (\w{1,12}  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,})) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}  (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%
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
# quadparserCommand = 'python G4-predictorModified.py -w 30 -s 1.4 '  # status for reference dataset: MCC:0.076, precision 100%
# quadparserCommand = 'python G4-predictorModified.py -w 25 -s 1.4 '  # status for reference dataset: MCC:0.116, precision 100%
# quadparserCommand = 'python G4-predictorModified.py -w 15 -s 1.5 '  # status for reference dataset: MCC:0.349, precision 98.3%
# quadparserCommand = 'python G4-predictorModified.py  -w 20 -s 1.7 '  # status for reference dataset: MCC:0.123, precision 92.1%
# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.4 '  # status for reference dataset: MCC:0.179, precision 100%
# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.2 '  # status for reference dataset: MCC:0.143, precision 100%

# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 25 -s 1.0 '  # status for reference dataset: MCC:0.813, precision 96.875%
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

quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{1,} | (?P<mis>[G]{2,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{2,})) (\w{1,6}  (?(mis)[G]{1,}| ([G]{1,}|(?P<mis>[G]{2,}[AUTC][G]{4,}|[G]{4,}[AUTC][G]{2,})) )) (\w{1,11}  (?(mis)[G]{1,}| ([G]{1,}|(?P<mis>[G]{2,}[AUTC][G]{2,}|[G]{4,}[AUTC][G]{2,})) ))(\w{1,11}  (?(mis)[G]{1,}| ([G]{1,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))"  -f "testedG4s3.fa"'  # status for reference dataset: MCC:0.764, precision 94.6%


# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.1 '  #
starttime=time.time()
output= subprocess.check_output(quadparserCommand + ' -f "' + file+'"', shell=True)


print(time.time()-starttime)

G4HScoreTreshold=0.9
TP=0
FP=0
G4List=[]
nonG4List=[]
for line in output.splitlines():
    if line.__contains__("not GQ_"):
        G4no=int(re.search(r"[0-9]+",line).group(0))
        if G4no not in G4List:
            score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
            if score>G4HScoreTreshold:
                G4List.append(G4no)
                FP+=1
                print line,score
    elif line.__contains__("GQ_"):
        G4no=int(re.search(r"[0-9]+",line).group(0))
        if G4no not in nonG4List:
            score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
            if score>G4HScoreTreshold:
                nonG4List.append(G4no)
                TP+=1
                print line,score

if TP==0 and FP==0: exit() #if nothing exists then exit without error.
FN=298-TP
TN=94-FP
print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
print "MCC:",MCC
precision=float(TP)/(TP+FP)
print "precision:",precision*100,"%"