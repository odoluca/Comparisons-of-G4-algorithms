

"""
Original code by JocelynSP at https://github.com/JocelynSP/G4-Hunter/tree/master
modified to report in a standardized format.
"""

"""
Guanine_quadruplex_predictor.py -i <inputFasta> -o <outputDir> -w Window_width -s Threshold_score
Author: J Sietsma Penington
Date: 21 April 2017

Implementation of G4Hunter algorithm from
Re-evaluation of G-quadruplex propensity with G4Hunter
Amina Bedrat  Laurent Lacroix  Jean-Louis Mergny
Nucleic Acids Res (2016) 44 (4): 1746-1759.
DOI:  https://doi.org/10.1093/nar/gkw006

Input Fasta-format file of sequences is scanned for regions with high 'G-richness' +
'G-skewness', implemented as a score based on the number of consecutive G nucleotides.
To simultaneously assess the complementary strand the same score is calculated as a
negative number for consecutive C nucleotides. Skewness is implemented by adding the
positive G scores and negative C scores so that regions with balanced G and C numbers have
small scores.

Two output files are written, one containing all sliding windows with absolute value of
mean score above the given threshold, and one containing extended regions of merged
windows which pass the criterion.

Intervals are written in 0-based, half open format, (]. 'start' is the position
before the interval and 'end' is the last position in the interval.
For example the first nucleotide would be described:
start	end	length
    0	  1	     1

Two decisions that may need changing:
	1. A merged region of adjacent windows has the score for the region calculated, but
	is not filtered by the final score. It is possible a number of windows with
	scores above threshold may combine to a region with mean score below threshold, which
	would still be written to the output file. This is same as Bedrat-Lacroix-Mergny
	output, but may be undesirable.

	2. Score is recalculated for each single-nucleotide advance of the window.
	It would probably be more efficient to calculate the score for the first window, then
	add and subtract scores for just the first and last nucleotides as the window is
	advanced.

"""

import sys
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq


def ScoreSeq(seq, winStart=0, winEnd=-1):
    if winEnd < 0:
        winEnd = len(seq)
    # calculate the total score for a window in the supplied sequence.
    # Default is window=sequence. If seq is larger, leading and trailing nucleotides are
    # used to adjust the score for the window if they extend a run.
    priorGcount = priorCcount = GbeforeWindow = CbeforeWindow = score = 0
    if winStart > 0:  # there are leading nucleotides. Increment length of run
        i = winStart - 1
        if seq[winStart] == "G":
            while i >= 0 and seq[i] == "G" and GbeforeWindow < 4:
                GbeforeWindow += 1
                i -= 1
        elif seq[winStart] == "C":
            while i >= 0 and seq[i] == "C" and CbeforeWindow < 4:
                CbeforeWindow += 1
                i -= 1
    # this corrects for the fact that in the main algorithm the prior counts are used to
    # adjust score for whole run, not just nucleotide under consideration

    for i in range(winStart, winEnd):
        if seq[i] == "G":
            if priorGcount + GbeforeWindow < 4:
                score += priorGcount * 2 + 1 + GbeforeWindow
                priorGcount += 1
            else:  # if run reaches GGGG, each subsequent G adds 4 regardless of number
                score += 4
            priorCcount = 0
            CbeforeWindow = 0  # any C run is terminated
        elif seq[i] == "C":
            if priorCcount + CbeforeWindow < 4:
                score -= priorCcount * 2 + 1 + CbeforeWindow
                priorCcount += 1
            else:  # if run reaches CCCC, each subsequent C subtracts 4 regardless of number
                score -= 4
            priorGcount = 0
            GbeforeWindow = 0  # any G run is terminated
        else:  # nucleotide is A, T or ambiguous: no change to score, any run is terminated
            priorGcount = priorCcount = GbeforeWindow = CbeforeWindow = 0

    if winEnd < len(seq):  # there are trailing nucleotides. Look ahead to adjust score
        i = winEnd
        if seq[winEnd - 1] == "G":
            runlength = priorGcount
            while i < len(seq) and seq[i] == "G" and runlength < 4:
                score += priorGcount
                runlength += 1
                i += 1
        if seq[winEnd - 1] == "C":
            runlength = priorCcount
            while i < len(seq) and seq[i] == "C" and runlength < 4:
                score -= priorCcount
                runlength += 1
                i += 1
    meanScore = float(score) / (winEnd - winStart)
    return meanScore

def ReverseComplement(seq):
    seq1 = 'ATCGNTAGCNatcgntagcn'
    seq_dict = {seq1[i]: seq1[i + 5] for i in range(20) if i < 5 or 10 <= i < 15}
    return "".join([seq_dict[base] for base in reversed(seq)])

def G4predictor(inFasta, win_width, thresh):
    regionCount = 0
    seqRecordCount = 0
    for seq_record in SeqIO.parse(inFasta, "fasta"):
        seqRecordCount += 1
        # out_win_file.write(seq_record.id + '\n')
        # out_win_file.write("Start\tEnd\tSequence\tLength\tScore\n")
        # out_region_file.write(seq_record.id + '\n')
        # print(seq_record.id + '\n')
        # print("Start\tEnd\tSequence\tLength\tScore\n")
        # print(seq_record.description)
        seq_curr = seq_record.upper()
        prevWinPass = False
        G4predCount = predG4start = predG4end = 0
        predG4seq = ""

        # Slide window along sequence, advance 1 nucleotide each loop:
        for i in range(len(seq_curr) - win_width + 1):
            win_curr = seq_curr[i: i + win_width]
            # where possible calculate score using 3 leading and 3 trailing bases, to get full value for bases in runs
            extended_win = seq_curr.seq[max(0, i - 3): min(i + win_width + 3, len(seq_curr))]
            score_curr = ScoreSeq(extended_win, min(3, i), min(win_width + 3, len(extended_win)))
            if abs(score_curr) > thresh:
                if i > predG4end:
                    if G4predCount == 0:
                        pass
                    else:
                        strand="-" if ScoreSeq(predG4seq)<0 else "+"
                        G4FormingStrand=predG4seq if ScoreSeq(predG4seq)<0 else ReverseComplement(predG4seq)
                        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(seq_record.description,
                            predG4start, predG4end, len(predG4seq), strand, predG4seq,
                             G4FormingStrand, ScoreSeq(predG4seq)))
                    G4predCount += 1
                    # Adjust start of new region:
                    if win_curr[0] in ["C", "G"]:
                        #  extend backwards to prior Gs or Cs if they match
                        predG4seq = win_curr.seq
                        predG4start = i
                        while predG4start > 0 and predG4seq[0] == seq_curr.seq[predG4start - 1]:
                            predG4seq = seq_curr.seq[predG4start - 1] + predG4seq
                            predG4start -= 1
                    else:  # strip off leading A and T nucleotides
                        oldseq=predG4seq
                        predG4seq = win_curr.seq.lstrip("AT")

                        # print oldseq[len(predG4seq)]
                        predG4seq = oldseq[len(predG4seq)]+predG4seq

                        predG4start = i + win_width - len(predG4seq)
                else:  # current window extends current G4 region
                    if prevWinPass:  # simple case of adjacent windows
                        predG4seq = predG4seq + win_curr.seq[-1]
                    else:
                        predG4seq = seq_curr[predG4start:i + win_width].seq
                predG4end = i + win_width
                prevWinPass = True
            else:
                if prevWinPass:
                    # a run of 1 or more windows with scores above threshold has just ended
                    if predG4seq[-1] in ["C", "G"]:
                        #  extend forward to next base if it matches
                        if predG4seq[-1] == win_curr.seq[-1]:
                            predG4seq = predG4seq + win_curr.seq[-1]
                    else:  # strip off trailing A and T nucleotides
                        oldseq = predG4seq
                        predG4seq = predG4seq.rstrip("AT")
                        predG4seq = predG4seq+oldseq[len(predG4seq)]
                    predG4end = predG4start + len(predG4seq)

                prevWinPass = False

        # Adjust end of, and write, last region
        if G4predCount == 0:
            pass
            # out_region_file.write("# No predicted quadruplex regions found\n")
        else:
            while (predG4end < len(seq_curr) - 1 and predG4seq[-1] in ["C", "G"] and
                   predG4seq[-1] == seq_curr.seq[predG4end + 1]):
                predG4seq = predG4seq + seq_curr.seq[predG4end + 1]
                predG4end += 1
            predG4seq = predG4seq.rstrip("AT")
            predG4end = predG4start + len(predG4seq)
            # out_region_file.write('{0}\t{1}\t{2}\t{3}\t{4:.2f}\t{5}\n'.format(
            #     predG4start, predG4end, predG4seq, len(predG4seq),
            #     ScoreSeq(predG4seq), G4predCount))
            strand = "-" if ScoreSeq(predG4seq) < 0 else "+"
            G4FormingStrand = predG4seq if ScoreSeq(predG4seq) > 0 else ReverseComplement(predG4seq)
            print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(seq_record.description,
                                                                  predG4start, predG4end, len(predG4seq), strand,
                                                                  predG4seq,
                                                                  G4FormingStrand, ScoreSeq(predG4seq)))
        regionCount += G4predCount
    # print "Number of sequences:", seqRecordCount
    # print "Number of potential quadruplex regions from given parameters:", regionCount


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     'Search a FASTA file for regions likely to form G4 quadruplexes')
    parser.add_argument('-f', dest='fasta_name', metavar='<infile>',
                        help='the input fasta file')
    # parser.add_argument('-o', dest='out_dir', metavar='<outdir>',default="",
    #                     help='directory for output files')
    parser.add_argument('-w', dest='window_width', type=int,
                        help='Width of sliding window. Recommended value 25 nts')
    parser.add_argument('-s', dest='threshold', type=float,
                        help='Threshold for absolute value of score to predict a quadruplex region. Typical value 1.0 to 1.5')

    args = parser.parse_args()
    # basefname = os.path.basename(args.fasta_name).split('.')[0]
    # outDirName = os.curdir
    # outWinName = os.path.join(outDirName,
    #                           basefname + "_" + str(args.window_width) + "nts.tsv")
    # outputWindows = open(outWinName, 'w')
    # outRegionName = os.path.join(outDirName,
    #                              basefname + "_"+str(args.threshold)+"_merged.tsv")
    # outputRegions = open(outRegionName, 'w')

    # print "Input file:", os.path.basename(args.fasta_name)
    G4predictor(args.fasta_name,  args.window_width, args.threshold)