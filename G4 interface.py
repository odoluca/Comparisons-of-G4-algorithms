import os
import subprocess
import time, math

file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
file="testedG4s.fa"


quadparserCommand = 'python quadparserModified.v3.py ' #shell command for quadparser command. Can add new regex pattern: -r "\w{1}([gG]{3}\w{1,7}){3,}[gG]{3}\w{1}"
# quadparserCommand = 'python quadparserModified.v3.py -r "([gG]{2,}\w{1,7}){3,}[gG]{2,}"'
# quadparserCommand = 'python ImGQfinder.v2.py ' #shell command for perfect, buldged or mismatched sequences
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w( [G]{3,} | (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )) (\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) ))(\w{1,7}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{2,}|[G]{2,}[ATC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w( [G]{3,} | (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )) (\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) ))(\w{1,5}  (?(mis)[G]{3,}| ([G]{3,}|(?P<mis>[G]{1,}[ATC][G]{2,}|[G]{2,}[ATC][G]{1,})) ))+\w" ' #only perfect and buldged
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[ATC][G]{1,})) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))(\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[AT][G]{1,})) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) )) (\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))(\w{1,7}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[ATC][G]{1,})) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) )) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))(\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[ATC][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "( [G]{2,} | (?P<mis>[G]{1,}[AT][G]{1,})) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) )) (\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))(\w{1,5}  (?(mis)[G]{2,}| ([G]{2,}|(?P<mis>[G]{1,}[AT][G]{1,})) ))+" '
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w([G]{3,}?| (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,7}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )){3,}\w" '  # only perfect and buldged BEING TESTED
# quadparserCommand = 'python ImGQfinder.v2.py -r "\w([G]{3,}?| (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,5}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )){3,}\w" '  # only perfect and buldged BEING TESTED
# quadparserCommand = 'python ImGQfinder.v2.py -r "([G]{3,}?| (?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) (\w{1,5}?  (?(mis)[G]{3,}?| ([G]{3,}?|(?P<mis>[G]{2,}[ATC][G]{1,}|[G]{1,}[ATC][G]{2,})) )){3,}" '  # only perfect and buldged BEING TESTED
# quadparserCommand = 'python G4-predictorModified.py -w 30 -s 1.4 '  # shell command for quadparser command
# quadparserCommand = 'python G4-predictorModified.py -w 25 -s 1.6 '  # shell command for quadparser command
# quadparserCommand = 'python G4-predictorModified.py -w 20 -s 1.7 '  # shell command for quadparser command
# quadparserCommand = 'python G4-predictorModified.v2b.py -w 30 -s 1.4 '  # shell command for quadparser command
# quadparserCommand = 'python G4HunterModified.v2.py -w 25 -s 1.4 '  # shell command for quadparser command
starttime=time.time()
output= subprocess.check_output(quadparserCommand + ' -f "' + file+'"', shell=True)

print(time.time()-starttime)

TP=0
FP=0
for line in output.splitlines():
    if line.__contains__("not G4_"):
        FP+=1
        print(line)
    elif line.__contains__("G4_"):
        TP+=1
        print(line)

FN=298-TP
TN=94-FP
print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN*FN))
print "MCC:",MCC
precision=float(TP)/(TP+FP)
print "precision:",precision*100,"%"
