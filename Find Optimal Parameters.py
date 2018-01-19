import subprocess
import time, math
import re
from random import randint,choice, uniform

filename="testedG4s3.fa"

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

def constructRegex(parametersList):
    (ifull, ishort, ilong, iloop, mfull, mshort, mlong, mloop, eloop, efull, elong, eshort)=tuple(parametersList[:12])
  # quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{'+str(ifull)+',} | (?P<mis>           ([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))  ( ((\w{1,'+str(iloop)+'})|(?P<lloop>\w{1,'+str(eloop)+'}) )  (?(mis)[G]{'+str(ifull)+',}| ([G]{'+str(ifull)+',}|(?P<mis>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))))  ( (?(lloop)\w{1,'+str(iloop)+'}|(?P<lloop>\w{1,'+str(ilong)+'}) )  (?(mis)[G]{'+str(ifull)+',}| ([G]{'+str(ifull)+',}|(?P<mis>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))))  ((?(lloop)\w{1,'+str(iloop)+'}|(?P<lloop>\w{1,'+str(ilong)+'}) )   (?(mis)[G]{'+str(ifull)+',}| ([G]{'+str(ifull)+',}|(?P<mis>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',}))) ))+" '  # status for reference dataset: MCC:0.809, precision 95.0%
    quadparserCommand = 'python ImGQfinder.v2.py --noreverse -r "( [G]{'+str(ilong)+',} | (?P<imperfectTract>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))(((\w{'+str(mloop-iloop)+','+str(mloop)+'})|(?P<ExtremeLoop>\w{1,'+str(eloop)+'}))(?(imperfectTract)([G]{'+str(mfull)+',})|([G]{'+str(mfull)+',}|(?P<imperfectTract>([G]{'+str(mlong)+',}[AUTC][G]{'+str(mshort)+',}|[G]{'+str(mshort)+',}[AUTC][G]{'+str(mlong)+',})))))((?(ExtremeLoop)\w{'+str(mloop-iloop)+','+str(mloop)+'}|(\w{'+str(mloop-iloop)+','+str(mloop)+'}|(?P<ExtremeLoop>\w{1,'+str(eloop)+'})))(?(imperfectTract)([G]{'+str(mfull)+',})|([G]{'+str(mlong)+',}|(?P<imperfectTract>([G]{'+str(mfull)+',}[AUTC][G]{'+str(mshort)+',}|[G]{'+str(mshort)+',}[AUTC][G]{'+str(mlong)+',})))))((?(ExtremeLoop)\w{'+str(mloop-iloop)+','+str(mloop)+'}|(\w{'+str(mloop-iloop)+','+str(mloop)+'}|(?P<ExtremeLoop>\w{1,'+str(eloop)+'})))(?(imperfectTract)([G]{'+str(efull)+',})|([G]{'+str(elong)+',}|(?P<imperfectTract>([G]{'+str(elong)+',}[AUTC][G]{'+str(eshort)+',}|[G]{'+str(eshort)+',}[AUTC][G]{'+str(elong)+',})))))"'  # status for reference dataset: MCC:0.849 precision 95.1% #has too short limitation as well for loops

    #if you want to search for both strands remove "--noreverse"
    return quadparserCommand

def callCommand(quadparserCommand,filename):
    print quadparserCommand + ' -f "' + filename + '"'
    output = subprocess.check_output(quadparserCommand + ' -f "' + filename + '"', shell=True)
    # for line in output.splitlines():
    #     print(line)
    return output



def MCCCalc(commandOut, G4HScoreTreshold=0.0):
    TP = 0
    FP = 0
    G4List = []
    nonG4List = []
    for line in commandOut.splitlines():
        if line.__contains__("not GQ_"):
            G4no = int(re.search(r"[0-9]+", line).group(0))
            if G4no not in G4List:
                score = G4HScore(re.search(r"[ATCUG]{5,}", line).group(0))
                if score > G4HScoreTreshold:
                    G4List.append(G4no)
                    FP += 1
                    # print line, score
        elif line.__contains__("GQ_"):
            G4no = int(re.search(r"[0-9]+", line).group(0))
            if G4no not in nonG4List:
                score = G4HScore(re.search(r"[ATCUG]{5,}", line).group(0))
                if score > G4HScoreTreshold:
                    nonG4List.append(G4no)
                    TP += 1
                    # print line, score

    if TP == 0 and FP == 0: return 0 # if nothing exists then exit without error.
    FN = 298 - TP
    TN = 94 - FP
    # print "TP:", TP, "FP:", FP, "FN:", FN, "TN:", TN
    MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    # print "MCC:", MCC
    precision = float(TP) / (TP + FP)
    # print "precision:", precision * 100, "%"
    return MCC

def CrossBreed(succesfulParams,childrenCount):
    newParameterList=[]
    for k in range(childrenCount/2):
        parent1=choice(succesfulParams)
        parent2=choice(succesfulParams)
        choices=[parent1,parent2]

        crossoverList=[randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(0,1),randint(5,15)]
        # crossoverList2=[1-x for x in crossoverList1]
        # print(crossoverList1,crossoverList2)
        child1=[]
        child2=[]
        for i in range(13):
            child1.append(choices[1][i]*crossoverList[i]+choices[0][i]*(1-crossoverList[i]))
            child2.append(choices[0][i]*crossoverList[i]+choices[1][i]*(1-crossoverList[i]))
        # print "child1:",child1, "child2:",child2
        newParameterList.append(Mutate(child1,1))
        newParameterList.append(Mutate(child2,1))
    return newParameterList

def Mutate(params,mutation=1):
    for i in range(mutation):
        param=randint(0,len(params)-2)
        params[param]=max(params[param]+randint(-1,1),1)
    return params


# initial_parameters=[3,1,2,12,3,1,2,12,12,3,2,1,0.0] #ifull, ishort, ilong, iloop, mfull, mshort, mlong, mloop, eloop, efull, elong, eshort, GScoreTresh


if __name__=="__main__":
    # print MCCCalc(callCommand(constructRegex(initial_parameters),filename),0.0)

    # print CrossBreed([[1,2,3,4,5,6,7,8,9,10,11,12,0.0],[21,22,23,24,25,26,27,28,29,30,31,32,0.0]])


    ListOfParameters=[]
    for k in range(1000):
        ListOfParameters.append([randint(2,4),randint(1,2),randint(2,4),randint(1,3),randint(1,4),randint(1,2),randint(2,4),randint(2,9),randint(7,40),randint(2,4),randint(2,4),randint(1,2),randint(5,15)])
                                     # ifull, ishort, ilong, iloop, mfull, mshort, mlong, mloop, eloop, efull, elong, eshort, GScoreTresh

    print(ListOfParameters)

    succesfulParams=[]
    succesfulMMCs=[]



    for generation in range(10):
        print "generation:",generation,"length of pool:",len(ListOfParameters)
        for params in ListOfParameters:
            output= callCommand(constructRegex(params),filename)
            MCC=MCCCalc(output,0.0) #If you do not want G4Hunter filtration, set the second parameter to 0.0. The system will still randomize the last parameter but only 0.0 will be used. default=params[-1]/10.0

            if MCC>0.7+generation*0.02:
                succesfulParams.append(params)
                succesfulMMCs.append(MCC)
                print params,MCC
                if MCC>0.78+generation*0.02:
                    succesfulParams.append(params)
                    succesfulMMCs.append(MCC)
                    if MCC>0.85+generation*0.02:
                        succesfulParams.append(params)
                        succesfulMMCs.append(MCC)


        #keep best 20
        bestParams=[]
        bestMMCs=[]
        cutoff=sorted(succesfulMMCs)[min(50,len(succesfulMMCs)-1)]
        print "cutoff:",cutoff
        # print succesfulParams
        for k in range(len(succesfulMMCs)):
            if succesfulMMCs[k]>=cutoff:
                # print "here",succesfulMMCs[k]
                bestMMCs.append(succesfulMMCs[k])
                bestParams.append(succesfulParams[k])
        succesfulMMCs=bestMMCs
        succesfulParams=bestParams
        # print bestParams

        #cross breed to make new children...it can also introduce mutation
        ListOfParameters= CrossBreed(succesfulParams,1000)
    print "RESULTS:"
    maxMCC=succesfulMMCs[0]
    maxMCCno=-1

    for resultNo in range(len(succesfulParams)):
        if succesfulMMCs[resultNo]>maxMCC:
            maxMCC=succesfulMMCs[resultNo]
            maxMCCno=resultNo
        # print succesfulParams[resultNo],succesfulMMCs[resultNo]

    print "most succesful MCC:",maxMCC
    bestParams=[]
    for resultNo in range(len(succesfulParams)):
        if succesfulMMCs[resultNo]==maxMCC and succesfulParams[resultNo] not in bestParams:
            bestParams.append(succesfulParams[resultNo])
            print succesfulParams[resultNo]

