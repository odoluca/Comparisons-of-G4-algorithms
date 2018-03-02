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
# file="testedG4s_5.fa"

def ConstructRegex(typLoopMax=7,shrtLoopMax=2,extLoopMax=30,typLoopMin=1,shrtLoopMin=1,extLoopMin=1):
    G2sAllowed=False
    ExtremeAllowed=False
    ExtremeAllowedForG2s=False
    ImperfectTractsAllowed=1
    BulgedTractsOnly=False
    typLoopMax=str(typLoopMax)
    extLoopMax=str(extLoopMax)
    shrtLoopMax=str(shrtLoopMax)
    typLoopMin = str(typLoopMin)
    extLoopMin = str(extLoopMin)
    shrtLoopMin = str(shrtLoopMin)
    if not G2sAllowed or not ExtremeAllowed:
        ExtremeAllowedForG2s=False
    #region regex construct
    bulgeOnly="[G]{2,}[ATUC][G]+|[G]+[ATUC][G]{2,}"
    mismatch="[G]{2,}|[G]+[ATUC][G]+"
    Dimp1='?P<imp1>'
    Dimp2='?P<imp2>'
    imp='('+mismatch+')'
    if BulgedTractsOnly:
        imp='('+bulgeOnly+')'
    Timp1='?(imp1)'
    Timp2='?(imp2)'
    shrt='\w{'+shrtLoopMin+','+shrtLoopMax+'}'
    Tshrt='?(shrt)'
    Dshrt='?P<shrt>[G]{2}'
    typ='\w{'+typLoopMin+','+typLoopMax+'}'
    ext='\w{'+extLoopMin+','+extLoopMax+'}'
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
    # print(structure)
    return structure
    #endregion



# for iteration in range(1000):
def iterate(args):

    typLoopMax=7
    shrtLoopMax=2
    extLoopMax=30
    typLoopMin=1
    shrtLoopMin=1
    extLoopMin=1

    if len(args)==3:
        typLoopMax, shrtLoopMax, extLoopMax=args
    elif len(args)==6:
        typLoopMax, shrtLoopMax, extLoopMax, typLoopMin, shrtLoopMin, extLoopMin=args

    quadparserCommand = r'python ImGQfinder.v2.py --noreverse -r "' + ConstructRegex(typLoopMax,shrtLoopMax,extLoopMax,typLoopMin,shrtLoopMin,extLoopMin) + '"'

    output = subprocess.check_output(quadparserCommand + ' -f "' + file + '"', shell=True)

    # print(output)
    G4HScoreTreshold=0#.473#473
    TP=0
    FP=0
    G4List=[]
    nonG4List=[]
    for line in output.splitlines():
        if line.__contains__("not GQ_"):
            G4no=int(re.search(r"[0-9]+",line).group(0))
            if G4no not in G4List:
                score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
                if abs(score)>=G4HScoreTreshold:
                    G4List.append(G4no)
                    FP+=1
                    # print line,score
        elif line.__contains__("GQ_"):
            G4no=int(re.search(r"[0-9]+",line).group(0))
            if G4no not in nonG4List:
                score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
                if abs(score)>=G4HScoreTreshold:
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
    # queue.put(str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94))
    # print "MCC:%.3f precision:%.1f TPR:%.3f FPR:%.3f" % (MCC, precision*100,float(TP)/298,float(FP)/94)
    if len(args)==3:
        # report= str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94)
        report= {"typLoopMax":typLoopMax,"shrtLoopMax":shrtLoopMax,"extLoopMax":extLoopMax,"MCC":MCC,"TPR":float(TP)/298,"FPR":float(FP)/94}


    elif len(args)==6:
        # report= str(typLoopMax) + "\t" + str(shrtLoopMax) + "\t" + str(extLoopMax) + "\t" + str(typLoopMin) + "\t" + str(shrtLoopMin) + "\t" + str(extLoopMin) + "\t" +  str(MCC) + "\t" + str(float(TP) / 298) + "\t" + str(float(FP) / 94)
        report = {"typLoopMax":typLoopMax,"shrtLoopMax":shrtLoopMax,"extLoopMax":extLoopMax,"MCC":MCC,"TPR":float(TP)/298,"FPR":float(FP)/94,"typLoopMin":typLoopMin, "shrtLoopMin":shrtLoopMin,"shrtLoopMin":extLoopMin}

    return report

import multiprocessing
import itertools

if __name__=="__main__":
    print("sample regex:")
    print(ConstructRegex(7,4,30))


    p=multiprocessing.Pool(30)

    allParams=list(itertools.product(*[[t for t in range(1,16)],[s for s in range(1,16)],[e for e in range(15,45)]]))

    results=p.map(iterate,allParams )


    for hit in results:
        print "\t".join([str(hit["typLoopMax"]),str(hit["shrtLoopMax"]),str(hit["extLoopMax"]),str(hit["MCC"]),str(hit["TPR"]),str(hit["FPR"])])

    """THIS CODE FINDS THE HITS THAT HAS THE BEST TPR FOR EACH FPR"""
    #
    # BestTpls = {}
    # for hit in results:
    #     if hit==None: continue
    #     if hit['FPR'] not in BestTpls.keys():
    #         BestTpls.update({hit['FPR']:[hit]})
    #     elif hit['FPR'] in BestTpls.keys():
    #         if BestTpls[hit['FPR']][0]['TPR']==hit['TPR']:
    #             BestTpls[hit['FPR']].append(hit)
    #         elif BestTpls[hit['FPR']][0]['TPR']<hit['TPR']:
    #             BestTpls[hit['FPR']]=[hit]
    #
    #
    # for key, hitList in BestTpls.iteritems():
    #     for hit in hitList:
    #         print "\t".join([str(hit["typLoopMax"]),str(hit["shrtLoopMax"]),str(hit["extLoopMax"]),str(hit["MCC"]),str(hit["TPR"]),str(hit["FPR"])])

    """END"""





