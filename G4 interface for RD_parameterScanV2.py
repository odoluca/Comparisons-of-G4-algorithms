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

# print G4HScore("GGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG")


file="Mitochondria_NC_012920_1.fasta"
file="test 1.fa"
file="testedG4s_4.fa"
# file="test 12fa"
file="testedG4s3.fa"
# file="empty.fa"

def ConstructRegex(typLoopMax=7,extLoopMax=30,shrtLoopMax=4):
    G2sAllowed=True
    ExtremeAllowed=False
    ExtremeAllowedForG2s=False
    if not G2sAllowed or not ExtremeAllowed:
        ExtremeAllowedForG2s=False
    ImperfectTractsAllowed=1
    BulgedTractsOnly=True
    typLoopMax=str(typLoopMax)
    extLoopMax=str(extLoopMax)
    shrtLoopMax=str(shrtLoopMax)

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
    structure=r'('+Tract1+')  ('+Loop1+')  ('+ Tract2+') ('+Loop2+') ('+Tract3+') ('+Loop2+') ('+Tract3+')'
    return structure



quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r " ([G]{3,}|(?P<imp1>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2})|(?P<shrt>[G]{2}))  (?(shrt)\w{1,4}|\w{1,7})   (?(shrt)[G]{2,}|(?(imp1)([G]{3,}|(?P<imp2>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2}))|([G]{3,}|(?P<imp1>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2})))) (?(shrt)\w{1,4}|\w{1,7})  (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2})) )|([G]{3,}|(?P<imp1>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2}))))  (?(shrt)\w{1,4}|\w{1,7}) (?(shrt)[G]{2,}|(?(imp1)(?(imp2)[G]{3,}|([G]{3,}|(?P<imp2>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2})) )|([G]{3,}|(?P<imp1>[G]{2}[ATUC][G]+|[G]+[ATUC][G]{2}))))"'
                     # python ImGQfinder.v2.py --noreverse -r "([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))|(?P<shrt>[G]{2}))  (?(shrt)(\w{1,4}|(?P<lloop>\w{1,30}))|(\w{1,7}|(?P<lloop>\w{1,30})))  (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))) (?(shrt)(?(lloop)\w{1,4}|(\w{1,4}|(?P<lloop>\w{1,30})))|(?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}))))) (?(shrt)(?(lloop)\w{1,4}|(\w{1,4}|(?P<lloop>\w{1,30})))|(?(lloop)\w{1,7}|(\w{1,7}|(?P<lloop>\w{1,30})))) (?(shrt)[G]{2,}|(?(imp1)([G]{3,})|([G]{3,}|(?P<imp1>([G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,})))))"


from numpy import arange, random

# for iteration in range(1000):
def iterate(queue):

    for iteration in range(100):
        # print(iteration)
        # print(parameter)
        # typLoopMax=random.randint(1,15)
        # shrtLoopMax=max(2,typLoopMax-random.randint(1,12))
        # extLoopMax=max(typLoopMax,typLoopMax+random.randint(1,40))
        window='25'#str(random.randint(15,30))
        treshold=str(random.uniform(1.0,2))
        quadparserCommand = 'python G4HunterModified.v2.py --noreverse -w '+window+' -s '+treshold  # status for reference dataset: MCC:0.812, precision 96.9%, TPR:0.94, FPR: 0.10

        # quadparserCommand = r'python ImGQfinder.v2.py --noreverse -r "' + ConstructRegex(typLoopMax,shrtLoopMax,extLoopMax) + '"'
        output = subprocess.check_output(quadparserCommand + ' -f "' + file + '"', shell=True)
        G4HScoreTreshold=0.#473
        TP=0
        FP=0
        G4List=[]
        nonG4List=[]
        for line in output.splitlines():
            if line.__contains__("not GQ_"):
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
        # print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
        MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        # print "MCC:",MCC
        precision=float(TP)/(TP+FP)
        # print "precision:",precision*100,"%"
        # print  "TPR:", float(TP) / 298, "FPR:", float(FP) / 94
        # if float(TP)/298>0.8 and float(FP)/94<0.14:
        # print str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94) #,quadparserCommand
        # print "MCC:%.3f precision:%.1f TPR:%.3f FPR:%.3f" % (MCC, precision*100,float(TP)/298,float(FP)/94)
        # queue.put(str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94))
        queue.put(window+"\t"+treshold+"\t"+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94))

import multiprocessing

if __name__=="__main__":
    print(ConstructRegex(7,4,30))
    jobs=[]
    queue=multiprocessing.Queue()
    for i in range(2):
        p=multiprocessing.Process(target=iterate,args=(queue,))
        jobs.append(p)
        p.start()

    while True:
        print(queue.get())
        time.sleep(0.01)
