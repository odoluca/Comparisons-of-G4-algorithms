#!/usr/bin/python
########################################################################
"""
    <G4Hunter - a program to search quadruplex-forming regions in DNA.>
    Copyright (C) <2012>  <Bedrat amina  supervised by Dr.Jean Louis Mergny>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
########################################################################
import os, re, sys, getopt
import time
import shutil
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
from Bio import SeqIO

 
def main(argv):

   if not argv:
        sys.stdout.write("Sorry: you must specify at least an argument\n")
        sys.stdout.write("More help avalaible with -h or --help option\n")
        sys.exit(1)

   inputfile = ''
   outputfile = ''
   noreverse=False
   try:
      opts, args = getopt.getopt(argv,"hf:o:w:s:n",["help","ffile=","ofile=","noreverse"])
   except getopt.GetoptError:
      print '\033[1m' +'python G4Hunter.py -f <inputfile> -o <outputrepository> -w <window> -s <score threshold> --noreverse\n'+'\033[0;0m'
      sys.exit(1)
      
   for opt, arg in opts:
      if opt in ('-h',"--help"):
          print  '\033[1m' +'\t ----------------------'+'\033[0;0m'
          print  '\033[1m' +'\n\t  Welcome To G4Hunter :'+'\033[0;0m'
          print  '\033[1m' +'\t ----------------------'+'\033[0;0m'

          print 'G4Hunter takes into account G-richness and G-skewness of a given sequence and gives a quadruplex propensity score as output.'
          print 'To run G4Hunter use the commande line: \n'
          print  '\033[1m' +'python G4Hunter.py -f <inputfile> -o <outputrepository> -w <window> -s <score threshold> --noreverse\n'+'\033[0;0m'
          sys.exit()
      elif opt in ("--noreverse"):
         noreverse=True
      elif opt in ("-f", "--ffile"):
         inputfile= arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-w", "--wind"):
         window = arg
      elif opt in ("-s", "--score"):
         score = arg
             
   return  inputfile, outputfile, int(window), float(score), noreverse

#calcule le score de chaque base dans un fichier qui peut contenir plusieur sequences fasta
#le fichier doit comporter la nom de la seq + la sequence en une seul ligne donc pas de \n ou \r 
class Soft(object):

    def __init__(self):
        pass
     
    def ReadFile(self,Filein):
        ListSeq,LHeader=[],[]
        for record in SeqIO.parse(Filein, "fasta") :
            LHeader.append(record.description)
            ListSeq.append(record.seq)
        return LHeader,ListSeq

    def ReadSeq(self, Seqs):
        ListSeq,LHeader=[],[]
        for record in SeqIO.parse(Seqs, "fasta") :
            LHeader.append(record.description)
            ListSeq.append(record.seq)
        return LHeader,ListSeq

    def GFinder(self,Filein,k):
        LHeader,ListSeq=self.ReadFile(Filein)
        LSeq,LNumber,LScoreSeq,SeqLine=[],[],[],""
        for i in range(len(ListSeq)) :
            Sequence,liste=self.BaseScore(ListSeq[i])
            LSeq.append(Sequence)
            LScoreSeq.append(self.CalScore(liste, k))
            LNumber.append(liste)
        return LScoreSeq, LSeq , LNumber, LHeader
    
    def BaseScore(self,line):
        item, liste=0, []
        #calcule le item de chaque base et la stock dans liste
        while ( item < len(line)):
            #a la fin d une sequence il est possible d avoir des GGG dans se cas
            # on verifie si la secore+1<len(line) car il ya un deuxieme G 
            #et

            if (item < len(line) and (line[item]=="G" or line[item]=="g")):
                liste.append(0)
                #print liste
                if (item + 1 < len(line) and (line[item + 1] == "C" or line[item + 1] == "c")):
                    liste[-1]+=-2
                if(item+1< len(line) and (line[item+1]=="G" or line[item+1]=="g")):
                    liste[item]=1
                    liste.append(1)
                    if (item+2< len(line) and (line[item+2]=="G" or line[item+2]=="g")):
                        liste[item+1]=2
                        liste[item]=2
                        liste.append(2)
                        if (item+3< len(line) and (line[item+3]=="G" or line[item+3]=="g")):
                            liste[item]=3
                            liste[item+1]=3
                            liste[item+2]=3
                            liste.append(3)
                            item=item+1
                        item=item+1
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="G" or line[item]=="g")):
                        liste.append(3)
                        item=item+1
        
            elif (item < len(line) and (line[item]=="T" or line[item]=="A"  or line[item]=="t" or line[item]=="a" or line[item]=="U"or line[item]=="u" or
                                         line[item]=="-" or line[item]=="N" or line[item]=="_" or line[item]=="Y" 
                                         or line[item]=="W" or line[item]=="R" or 
                                        line[item]=="K" or line[item]=="M"or line[item]=="S" or line[item]=="B"
                                        or line[item]=="V"or line[item]=="D"or line[item]=="H"or line[item]=="N")):
                liste.append(0)
                item=item+1
                
            elif(item < len(line) and (line[item]=="C" or line[item]=="c")):
                liste.append(-0)
                if (item + 1 < len(line) and (line[item + 1] == "G" or line[item + 1] == "g")):
                    liste[-1]+=2
                if(item+1< len(line) and (line[item+1]=="C" or line[item+1]=="c" )):
                    liste[item]=-1
                    liste.append(-2)
                    if (item+2< len(line) and (line[item+2]=="C" or line[item+2]=="c" )):
                        liste[item+1]=-2
                        liste[item]=-2
                        liste.append(-2)
                        if (item+3< len(line) and (line[item+3]=="C" or line[item+3]=="c"  )):
                            liste[item]=-3
                            liste[item+1]=-3
                            liste[item+2]=-3
                            liste.append(-3)
                            item=item+1
                        item=item+1   
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="C" or line[item]=="c")):
                    liste.append(-3)
                    item=item+1
            
            else:
                    item=item+1 #la fin du la ligne ou y a des entrers
        return line, liste
 
    def CalScore(self,liste, k):
        Score_Liste=[]
        #calcule de la moynne des scores pour toutes les sequences - les derniers k bases
        for i in range (len (liste)-k):
            j,Sum=0,0
            while (j<k):
                Sum=Sum+liste[i]
                j=j+1
                i=i+1
            Mean=Sum/float(k)
            Score_Liste.append(Mean) 
        return Score_Liste
    
    ###################################################
    ###################################################


    def plot2(self,liste, repert):
        # make a little extra space between the subplots
        plt.subplots_adjust(wspace=1.0)
        dt = 1
        t = np.arange(0, len(liste), dt)
        figure= plt.figure()
        plt.plot(t, liste, 'b-')
        plt.xlim(0,len(liste))
        #figure.suptitle('Score of window nucleotides', fontsize=16)
        plt.xlabel('Position (ntS)')
        plt.ylabel('Score')
        plt.grid(True)
        figure.savefig(repert+'Score_plot.pdf', dpi=figure.dpi)

    # modification: reverse complement function.
    def ReverseComplement(self,seq):
        seq1 = 'ATCGNWSMKRYBDHVatcgnwsmkrybdhv'
        seq2 = 'TAGCNWSKMYRVHDBtagcnwskmyrvhdb'
        seq_dict = {seq1[i]: seq2[i] for i in range(len(seq1))}
        return "".join([seq_dict[base] for base in reversed(seq)])

    """ 
    ###################################################
    ###################################################
    """               
    def GetG4(self,line,liste, Window,k, header, Len ):
        LG4=[]
        SEQ=">"+header+"\n Start \t End \t Sequence\t Length \t Score\n"
        # fileout.write(SEQ)
        for i in range(len(liste)) :
            if (liste[i]>= float(Window) or liste[i]<= - float(Window)):
                seq=line[i:i+k]
                LG4.append(i)
                # self.Write( i, k ,0,0, seq , k, liste[i],"","")
                # fileout.write("\n")
        return LG4
   
    def WriteSeq(self,line,liste, LISTE,header, F, Len , noreverse ):
        i,k,I=0,0,0
        a=b=LISTE[i]
        MSCORE=[]
        #SEQ=header+"\nStart\tEnd\tSequence\tLength\tScore\tNBR\n"
        #fileout.write(SEQ)
        if (len(LISTE)>1):
            c=LISTE[i+1]
            while (i< len(LISTE)-2):
                if(c==b+1):
                    k=k+1
                    i=i+1
                else:
                    I=I+1
                    seq=line[a:a+F+k]
                    sequence,liste2=self.BaseScore(seq)
                    self.Write( a, k ,F,0, seq ,len(seq) , round(np.mean(liste2),2),header,"", noreverse)
                    MSCORE.append(abs(round(np.mean(liste2),2)))
                    # fileout.write("\n")
                    k=0
                    i=i+1
                    a=LISTE[i]
                b=LISTE[i] 
                c=LISTE[i+1] 
            I=I+1
            seq=line[a:a+F+k+1]
            sequence,liste2=self.BaseScore(seq)
            self.Write( a, k ,F,1, seq ,len(seq) , round(np.mean(liste2),2),header,"", noreverse)
            MSCORE.append(abs(round(np.mean(liste2),2)))
            # fileout.write("\t")
            # fileout.write(str(I))
            # fileout.write("\n")
        #dans le cas ou il ya une seul sequence donc pas de chevauchement
            #writing into the file is omitted if there is no result
        # else:
        #     I=I+1
        #     seq=line[a:a+F]
        #     self.Write(fileout, a, 0 ,F,0, seq ,len(seq) , liste[a])
        #     MSCORE.append(abs(liste[a]))
        #     fileout.write("\t")
        #     fileout.write(str(I))
        #     fileout.write("\n")
        return MSCORE
    
    def Write(self, i, k ,F,X, seq ,long, score, source, strand, noreverse):
        if score<0:
            if noreverse:
                return
            strand="-"
            G4sequence=self.ReverseComplement(seq)
        else:
            strand="+"
            G4sequence = seq
        LINE=str(source)+"\t"+str(i)+"\t"+str(i+k+F+X)+"\t"+str(long)+"\t"+str(strand)+"\t"+str(seq) +"\t"+str(G4sequence)+"\t"+str(score)
        # fileout.write(LINE)
        # if fileout==Res2file:
        print(LINE)



       
    #Len dans le cas ou la sequence fini avec des ----- donc il yaura une erreur



if __name__ == "__main__":
    inputfile, outputfile , window, score, noreverse = main(sys.argv[1:])
    # print window
    # OPF= os.listdir(outputfile)
    # flag=False
    # for dir in OPF:
    #     if dir== "Results":
    #         flag=True
    # if flag==True:
    #     shutil.rmtree(outputfile+"/Results/")
    #     os.makedirs(outputfile+"/Results/", mode=0777)        #
    #     print '\033[1m' +"\n \t Re-evaluation of G-quadruplex propensity with G4Hunter " +'\033[0;0m'
    #     print "\n#####################################"
    #     print "#    New Results directory Created  #"
    #     print "#####################################\n"
    # else:
    #     os.makedirs(outputfile+"/Results/", mode=0777)
    #     print "\n########################################################################"
    #     print "#                            Results directory Created                 #"
    #     print "########################################################################\n"

    #================================================================
    # plot=[]
    fname=inputfile.split("/")[-1]
    filefasta=fname.split(".")
    filein=open(inputfile,"r")
    # print "\n Input file:", '\033[1m' + filefasta[0]+'\033[0;0m'
    #repertoire des fichiers de sortie

    # Res1file= open (outputfile +filefasta[0]+"-"+ str(window)+"-nts", "w")
    # Res2file= open (outputfile +filefasta[0]+"-Merged", "w")
    #=========================================
    
    startTime = time.time()
    
    
    soft1=Soft()
    ScoreListe, DNASeq, NumListe, HeaderListe=soft1.GFinder(filein, window)
    for i in range(len(DNASeq)):
        G4Seq=soft1.GetG4(DNASeq[i], ScoreListe[i], float(score), int(window), HeaderListe[i], len(NumListe[i]))
        if (len(G4Seq)>0):
            MSCORE=soft1.WriteSeq(DNASeq[i],ScoreListe[i], G4Seq, HeaderListe[i], int(window), len(NumListe[i]), noreverse)
        # plot.append(MSCORE)
    # soft1.plot2(ScoreListe[0], outputfile +"/Results/")
    filein.close()
    fin=time.time()

    # print "\n Results files and Score Figure are created in:   "#,fin-startTime, "secondes"
    # print '\033[1m' + outputfile,"/Results/ \n "+'\033[0;0m'





    # Res1file.close()
    # Res2file.close()

    
    
