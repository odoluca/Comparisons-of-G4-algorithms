import os
import subprocess
import time, math
import re
file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
file="testedG4s3.fa"



quadparserCommand = 'python quadparserModified.v3.py ' #shell command for quadparser command. Can add new regex pAUTtern: -r "\w{1}([gG]{3}\w{1,7}){3,}[gG]{3}\w{1}" # status for reference dataset: MCC:0.249, precision 98.6%
quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{3,}\w{1,7}){3,}[gG]{3,}"'# status for reference dataset: MCC:0.667, precision 94.2%
# quadparserCommand = 'python ImGQfinder.v2.py ' #shell command for perfect, buldged or mismAUTched sequences
#
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+ "'   #status for sentetic dataset: Finds all bulged: OK
                                                                                                                                                                                                                                                                                                            # status for reference dataset: MCC:0.670, precision 99.2%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "([G]{3,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+ "' # status for reference dataset: MCC:0.670, precision 99.2%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "\w( [G]{3,} | (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))(\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "\w( [G]{3,} | (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))(\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{3,} | (?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})) )) (\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})) ))(\w{1,12}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUT][G]{2,}|[G]{2,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.693, precision 99.6%

# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "\w( [G]{3,} | (?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) )) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[AUTC][G]{1,}|[G]{1,}[AUTC][G]{2,})) ))(\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[AUTC][G]{2,}|[G]{2,}[AUTC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+" ' #status for reference dataset: MCC:0.734, precision 93.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.734, precision 93.6%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.750, precision 94.3%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.780, precision 94.1%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,7}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,})) ))+" ' #status for reference dataset: MCC:0.787, precision 94.4%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,12}   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  (\w{1,12}   (?(mis)[G]{3,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,}))))  (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{2,})) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%


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
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.4 '  # status for reference dataset: MCC:0.610, precision 96.3%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.5 '  # status for reference dataset: MCC:0.510, precision 96.6%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.6 '  # status for reference dataset: MCC:0.439, precision 96.2%
# quadparserCommand = 'python G4HunterModified.v2.py -w 15 -s 1.7 '  # status for reference dataset: MCC:0.317, precision 94.8%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.4 '  # status for reference dataset: MCC:0.162, precision 100%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.5 '  # status for reference dataset: MCC:0.140, precision 100%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.6 '  # status for reference dataset: MCC:0.130, precision 100%
# quadparserCommand = 'python G4HunterModified.v2.py -w 30 -s 1.7 '  # status for reference dataset: MCC:0.120, precision 100%

# quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w 20 -s 1 '  # status for reference dataset: MCC:0.764, precision 94.6%
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( \w{1,12}                          (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  ( \w{1,12}                                (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))   ( \w{1,12}                               (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%
# quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUTC][G]{1,})) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})))) (\w{1,12}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUTC][G]{1,})) ))+" ' #status for reference dataset: MCC:0.747, precision 92.6%
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>[G]{1,}[AUT][G]{1,}))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,}))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AUT][G]{1,})) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%
quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{2,} | (?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))  ( ((\w{1,7})|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ( (?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,})))))  ((?(lloop)\w{1,7}|(?P<lloop>\w{1,30}) )   (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>([G]{2,}[AUT][G]{1,}|[G]{1,}[AUT][G]{2,}))) ))+" ' #status for reference dataset: MCC:0.809, precision 95.0%


# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.2'  #
starttime=time.time()
output= subprocess.check_output(quadparserCommand + ' -f "' + file+'"', shell=True)
# bulged_found=0
# lastline=""
# bulgedList=[]
# for line in output.splitlines():
#     if line.startswith("bulged_"):
#         number = int(re.search("[0-9]+", line).group(0))
#         if number not in bulgedList:
#             bulgedList.append(number)
#             bulged_found+=1
#             print line
# print "bulged out:",bulged_found


print(time.time()-starttime)

TP=0
FP=0
G4List=[]
nonG4List=[]
for line in output.splitlines():
    if line.__contains__("not GQ_"):
        G4no=int(re.search(r"[0-9]+",line).group(0))
        if G4no not in G4List:
            G4List.append(G4no)
            FP+=1
            print(line)
    elif line.__contains__("GQ_"):
        G4no=int(re.search(r"[0-9]+",line).group(0))
        if G4no not in nonG4List:
            nonG4List.append(G4no)
            TP+=1
            print(line)
if TP==0: exit("no true positive")

FN=298-TP
TN=94-FP
print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
print "MCC:",MCC
precision=float(TP)/(TP+FP)
print "precision:",precision*100,"%"