import os
import subprocess
import time, math
import re
file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
#file="testedG4s.fa"


quadparserCommand = 'python quadparserModified.v3.py ' #shell command for quadparser command. Can add new regex pattern: -r "\w{1}([gG]{3}\w{1,7}){3,}[gG]{3}\w{1}"
# quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{2,}\w{1,7}){3,}[gG]{2,}"'
# quadparserCommand = 'python ImGQfinder.v2.py ' #shell command for perfect, buldged or mismatched sequences

quadparserCommand = 'python ImGQfinder.v2.py -r "([G]{3,} | (?P<mis>[G]{1,}[ATC][G]{1,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))+ "' #status for sentetic dataset: bulged: 100/100 mismatched: 100/100

quadparserCommand = 'python ImGQfinder.v2.py -r "([G]{3,} | (?P<mis>[G]{1,}[ATC][G]{1,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))+ "' #
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w( [G]{3,} | (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) ))(\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{2,}|[G]{2,}[ATC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w( [G]{3,} | (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) ))(\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{2,}|[G]{2,}[ATC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[ATC][G]{1,})) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))(\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[AT][G]{1,})) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) )) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))(\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[ATC][G]{1,})) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))(\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[AT][G]{1,})) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) )) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))(\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w([G]{3,}?| (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,7}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )){3,}\w" '  # only perfect and buldged BEING TESTED
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w([G]{3,}?| (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,5}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )){3,}\w" '  # only perfect and buldged BEING TESTED
quadparserCommand = 'python ImGQfinder.v2.py -r "([G]{3,}?| (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,7}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )){3,}" ' # Target:perfect and bulged.
                                                                                                                                                                                                                # status for synthetic dataset: Bulged:100/100 Mismatched:66/100

# quadparserCommand = 'python G4-predictorModified.py -w 30 -s 1.4 '  # shell command for quadparser command
# quadparserCommand = 'python G4-predictorModified.py -w 25 -s 1.2 '  # shell command for quadparser command
# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.2 '  # shell command for quadparser command
# quadparserCommand = 'python G4-predictorModified.v2b.py -w 30 -s 1.4 '  # shell command for quadparser command
# quadparserCommand = 'python G4HunterModified.v2.py -w 25 -s 1.4 '  # shell command for quadparser command
starttime=time.time()
output= subprocess.check_output(quadparserCommand + ' -f "' + file+'"', shell=True)
# print  output

lastline=""
bulgedList=[]
bulged_found=0
mismatchedList=[]
mismatched_found=0
nonGQList=[]
nonGQ_found=0
for line in output.splitlines():
    if line.startswith("bulged_"):
        number = int(re.search("[0-9]+", line).group(0))
        if number not in bulgedList:
            bulgedList.append(number)
            bulged_found+=1
            print line
    if line.startswith("mismatched_"):
        number = int(re.search("[0-9]+", line).group(0))
        if number not in mismatchedList:
            mismatchedList.append(number)
            mismatched_found+=1
            print line
    if line.startswith("nonGQ_"):
        number = int(re.search("[0-9]+", line).group(0))
        if number not in nonGQList:
            nonGQList.append(number)
            nonGQ_found+=1
            print line
print "bulged out:",bulged_found
print "mismatched out:",mismatched_found
print "nonGQ out:",nonGQ_found


print(time.time()-starttime)

# TP=0
# FP=0
# for line in output.splitlines():
#     if line.__contains__("not G4_"):
#         FP+=1
#         print(line)
#     elif line.__contains__("G4_"):
#         TP+=1
#         print(line)
#
# FN=298-TP
# TN=94-FP
# print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
# MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
# print "MCC:",MCC
# precision=float(TP)/(TP+FP)
# print "precision:",precision*100,"%"
