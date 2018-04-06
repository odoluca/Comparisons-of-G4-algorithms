import os
import subprocess
import time, math
import re
import random as rd



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
    G2sAllowed=True
    ExtremeAllowed=True
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


def crateRandomList(limitStart,limitEnd,size):
    result=[]
    i=size
    while i>0:
        p=rd.randint(limitStart,limitEnd)
        if p not in result:
            result.append(p)
            i-=1
    return tuple(result)



# for iteration in range(1000):
def iterate(args,test=False):

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

    useG4List=(52, 282, 67, 71, 133, 114, 90, 155, 63, 32, 124, 196, 26, 176, 184, 5, 287, 252, 147, 154, 128, 2, 1, 85, 293, 277, 88, 197, 172, 113, 77, 212, 14, 123, 258, 236, 265, 270, 246, 37, 9, 152, 86, 187, 170, 220, 243, 129, 189, 11, 24, 60, 193, 110, 87, 262, 80, 255, 36, 232, 8, 156, 198, 45, 263, 249, 218, 75, 93, 16, 217, 102, 141, 27, 120, 4, 104, 12, 219, 92, 190, 41, 74, 58, 195, 158, 138, 175, 18, 259, 247, 201, 62, 53, 54, 191, 43, 188, 185, 205, 20, 226, 70, 165, 275, 231, 108, 57, 202, 178, 145, 199, 222, 294, 122, 118, 61, 94, 139, 149, 146, 33, 95, 192, 79, 257, 285, 166, 292, 82, 210, 279, 107, 281, 66, 241, 125, 31, 163, 148, 106, 10, 251, 127, 143, 235, 213, 101, 103, 299, 264, 186, 223, 164, 271, 261, 204, 266, 267, 34, 40, 171, 121, 221, 157, 227, 6, 215, 142, 83, 200, 291, 216, 162, 240, 42, 237, 272, 283, 168, 268, 136, 91, 278, 173, 161, 30, 245, 48, 286, 203, 68, 132, 297, 3, 151, 225, 159, 22, 98)

    usenonG4List=(330, 313, 328, 314, 388, 315, 375, 346, 387, 333, 308, 371, 302, 344, 342, 393, 351, 319, 372, 321, 369, 312, 359, 376, 367, 365, 368, 341, 317, 383, 310, 331, 301, 361, 377, 347, 364, 350, 378, 352, 385, 336, 307, 323, 316, 374, 366, 382, 326, 337)
    print("--------start--------")
    if  test:
        newG4=[]
        newnonG4=[]
        for i in range(1,299):

            if i not in useG4List:
                newG4.append(i)
            #     print ("test set")
            # else:
            #     print("training set")
        for i in range(299,393):
            if i not in usenonG4List:
                newnonG4.append(i)
            #     print ("test set")
            # else:
            #     print("training set")
        # print("now using test set"),
        # print(newG4)
        # print(newnonG4)
        useG4List=tuple(newG4)
        usenonG4List=tuple(newnonG4)
    print("--------end--------")

    # useG4List=range(298)
    # usenonG4List=range(299,392)

    for line in output.splitlines():
        if line.__contains__("not GQ_"):
            G4no=int(re.search(r"[0-9]+",line).group(0))
            if G4no not in G4List and G4no in usenonG4List :
                score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
                if abs(score)>=G4HScoreTreshold:
                    G4List.append(G4no)
                    FP+=1
                    # print line,G4no,score
        elif line.__contains__("GQ_"):
            G4no=int(re.search(r"[0-9]+",line).group(0))
            if G4no not in nonG4List and G4no in useG4List:
                score=G4HScore(re.search(r"[ATCUG]{5,}",line).group(0))
                if abs(score)>=G4HScoreTreshold:
                    nonG4List.append(G4no)
                    TP+=1
                    # print line,G4no, score

    if TP==0 and FP==0: return #exit() #if nothing exists then exit without error.
    # FN=71-TP
    # TN=138-FP
    FN = len(useG4List)- TP
    TN = len(usenonG4List)- FP
    # print "TP:",TP,"FP:",FP,"FN:",FN,"TN:",TN
    # MCC=(TP*TN-FP*FN)/ math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    Youden=float(TP/float(TP+FN))+float(TN/float(TN+FP))-1 #youden's statistics
    # print((TP/float(TP+FN)),(TN/float(TN+FP)),MCC)
    # print "MCC:",MCC
    precision=float(TP)/(TP+FP)
    # print "precision:",precision*100,"%"
    # print  "TPR:", float(TP) / 298, "FPR:", float(FP) / 94
    # if float(TP)/298>0.8 and float(FP)/94<0.14:
    # print str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94) #,quadparserCommand
    # queue.put(str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94))
    # print "MCC:%.3f precision:%.1f TPR:%.3f FPR:%.3f" % (MCC, precision*100,float(TP)/298,float(FP)/94)
    # print(TP,TN,FP,FN,MCC)
    if len(args)==3:
        # report= str(typLoopMax)+"\t"+str(shrtLoopMax)+"\t"+str(extLoopMax)+"\t"+str(MCC)+"\t"+str(float(TP)/298)+"\t"+str(float(FP)/94)
        report= {"typLoopMax":typLoopMax,"shrtLoopMax":shrtLoopMax,"extLoopMax":extLoopMax,"MCC":Youden,"TPR":float(TP)/len(useG4List),"FPR":float(FP)/len(usenonG4List)}


    elif len(args)==6:
        # report= str(typLoopMax) + "\t" + str(shrtLoopMax) + "\t" + str(extLoopMax) + "\t" + str(typLoopMin) + "\t" + str(shrtLoopMin) + "\t" + str(extLoopMin) + "\t" +  str(MCC) + "\t" + str(float(TP) / 298) + "\t" + str(float(FP) / 94)
        report = {"typLoopMax":typLoopMax,"shrtLoopMax":shrtLoopMax,"extLoopMax":extLoopMax,"MCC":Youden,"TPR":float(TP)/len(useG4List),"FPR":float(FP)/len(usenonG4List),"typLoopMin":typLoopMin, "shrtLoopMin":shrtLoopMin,"shrtLoopMin":extLoopMin}

    return report

import multiprocessing
import itertools

if __name__=="__main__":
    print("sample regex:")
    print(ConstructRegex(7,4,30))
    # useG4List = crateRandomList(1, 298 + 1, 200)  # can be between 1 and 298
    # usenonG4List = crateRandomList(299, 392+1, 50)  # can be between 299 and 392
    # print(useG4List)
    # print(usenonG4List)

    p=multiprocessing.Pool(1)

    allParams=list(itertools.product(*[[t for t in range(1,16)],[s for s in range(1,16)],[e for e in range(15,46)]]))

    results=p.map(iterate,allParams)
    print("--------------------------------------------------------------")
    bestSet=[{"typLoopMax":0,"shrtLoopMax":0,"extLoopMax":0,"MCC":0,"TPR":0,"FPR":0}]
    for hit in results:
        if hit!=None:
            print "\t".join([str(hit["typLoopMax"]),str(hit["shrtLoopMax"]),str(hit["extLoopMax"]),str(hit["MCC"]),str(hit["TPR"]),str(hit["FPR"])])
            if hit["MCC"]>bestSet[0]["MCC"]:
                print(hit["MCC"])
                bestSet=[hit]
            elif hit["MCC"] == bestSet[0]["MCC"]:
                print(hit["MCC"])
                bestSet.append(hit)
    print("--------------------------------------------------------------")

    # allParams=list(itertools.product(*[[t["typLoopMax"] for t in bestSet],[s["shrtLoopMax"] for s in bestSet],[e["extLoopMax"] for e in bestSet]]))
    # print(allParams[0])
    #
    # Tests=p.map(iterate,allParams)
    Tests=[]
    for params in bestSet:
        Tests.append(iterate((params["typLoopMax"],params["shrtLoopMax"],params["extLoopMax"],),True))

    bestTests=[Tests[0]]
    for hit in Tests:
        if hit["MCC"] > bestTests[0]["MCC"]:
            # print(hit["MCC"])
            bestTests = [hit]
        elif hit["MCC"] == bestTests[0]["MCC"]:
            # print(hit["MCC"])
            bestTests.append(hit)
    print("best MCC:",bestTests[0]["MCC"])
    for set in bestTests:
        print(set)

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

"""LATEST CHANGES FINDS BEST PARAMETERS USING A TRAINING SET AND COMPARES TO A TEST SET"""




