import argparse, sys, re, string, regex



" --------------------------[ Parse arguments ]---------------------- "

parser=argparse.ArgumentParser(
    description="""
    DESCRIPTION
    
        Searches for matches to a G-quadruplex-fitting regex in a fasta file, 
        filters through G4Hunter-like secondary scoring scheme and return a bed file with 
        coordinates of the match, matched sequence, G-quadruplex forming sequence and the score.
        
        the defult regex is 
        
        
        Output bed file has the following columns:
        1. description of the fasta sequence (e.g. NC_00024.11 Y chromosome)
        2. start of the match
        3. end of the match
        4. size of the match
        5. strand of the match (e.g. +)
        6. positive strand sequence of the match (e.g. CCCTTCCCTTTCCCTCCC)
        7. matched G-quadruplex-forming sequence (e.g. GGGAGGGAAAGGGAAGGG)
        8. score of the matched G-quadruplex-forming sequence based on selected scoring scheme
        
    EXAMPLE
        ##Test data:
        echo '>mychr' > /tmp/mychr.fa
        echo 'TTGGGTTGGGACTGGGTACGGGAATA' >> /tmp/mychr.fa
        
        G4Catchall.py -f /tmp/mychr.fa 
            mychr   2   22  20  +   GGGTTGGGACTGGGTACGGG    GGGTTGGGACTGGGTACGGG    2.11

        G4Catchall.py -f /tmp/mychr.fa --min_G-tract 4
            
    DOWNLOAD
        G4Catchall.p is hosted at 

""",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--min_G_tract_length','-t',
type=str,
help="""Minimum length of guanine tracts in the search.
Please not that if bulged or mismatched G-tracts are allowed, minimum length 
may be lower for at most one guanine tract. default=2

""",
default="2")

parser.add_argument('--no_bulge','-b',
action='store_true',
help="""By default the program allows a maximum of one G-tract 
with one bulge. If used, mismatch is 
automatically not allowed as well.

""")

parser.add_argument('--no_mismatch','-m',
action='store_true',
help="""By default the program allows a maximum of one G-tract 
with one mismatch. Does not affect bulge search.

""")

parser.add_argument('--min_typical_loop_length',
type=str,
help="""Allows to set the minimum typical loop length. 
default=1

""", default="1")

parser.add_argument('--max_typical_loop_length',
type=str,
help="""Allows to set the maximum typical loop length. 
default=7

""", default="7")

parser.add_argument('--no_extreme_loop',
action='store_true',
help="""
By default the program allows a single loop to have
 extreme lengths (too short or too long) If used only no 
 extreme loop will be allowed.

""")

parser.add_argument('--min_extreme_loop_length',
type=str,
help="""Allows to set the minimum loop length of an extreme loop. 
Please note that only one extreme loop is allowed by default.
default=1

""",default="1")

parser.add_argument('--max_extreme_loop_length',
type=str,
help="""Allows to set the maximum loop length of an extreme loop. 
Please note that only one extreme loop is allowed by default.
default=30

""",default="30")

parser.add_argument('--no_reverse',
action='store_true',
help="""By default the program searches both strands by 
reversing the regex. If used only + strand is searched for matches.

""")

parser.add_argument('--dont_merge_overlapping',
action='store_true',
help="""Putative G-quadruplex-forming sequences may be found 
overlapping on the same strand. By default the program merges 
these sequences. If used, these are not overlapped. Using may 
result in huge number of matches and cause memory issues.

""")

parser.add_argument('--include_flanks',
action='store_true',
help="""By default the program extracts only matching sequences. 
If used flanking nucleotides are also included in the search. 
Please note if used G-quadruplex-forming sequences at the beginning
or ending of the sequences may be missed. Consider adding "N" to the 
edges of the sequence if G-quadruplex forming sequences are expected 
to be found at the very edge of the target sequence.

""")

args=parser.parse_args()

" --------------------------[ End of Parse arguments ]---------------------- "


" ------------------------------[  Reverse complement ]--------------------------------- "
def ReverseComplement(seq):
    seq1 = 'ATCGNTAGCNatcgntagcn'
    seq_dict = {seq1[i]: seq1[i + 5] for i in range(20) if i < 5 or 10 <= i < 15}
    return "".join([seq_dict[base] for base in reversed(seq)])

" ------------------------------[  End of Reverse complement  ]--------------------------------- "


""" ------------------------------[  Sorter ]--------------------------------- 
Code to sort list of lists
see http://www.saltycrane.com/blog/2007/12/how-to-sort-table-by-columns-in-python/
"""
import operator


def sort_table(table, cols):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols:  a list (or tuple) specifying the column numbers to sort by
               e.g. (1,0) would sort by column 1, then by column 0
    """
    for col in reversed(cols):
        table = sorted(table, key=operator.itemgetter(col))
    return (table)

" --------------------------[ End of Sorter ]---------------------- "


" ------------------------------[  Set regex  ]--------------------------------- "
MinG=args.min_G_tract_length
NoBulge=args.no_bulge
NoMismatch=args.no_mismatch
MinTypLL=args.min_typical_loop_length
MaxTypLL=args.max_typical_loop_length
NoExtLoop=args.no_extreme_loop
MinExtLL=args.min_extreme_loop_length
MaxExtLL=args.max_extreme_loop_length
NoReverse=args.no_reverse
NoMerge=args.dont_merge_overlapping
IncFlank=args.include_flanks




if (not NoBulge):
    ImperfectTract='(?P<imperfectTract>'
    for k in range(1,int(MinG)):
        ImperfectTract+='[G]{'+str(k)+',}[AUTC][G]{'+str(int(MinG)-k)+',}'
        if k!=int(MinG)-1:
            ImperfectTract += '|'
    ImperfectTract += ')'
# print(ImperfectTract)


if (not NoBulge) and (not NoMismatch):
    ImperfectTract='(?P<imperfectTract>'
    for k in range(1,int(MinG)-1):
        ImperfectTract+='[G]{'+str(k)+',}[AUTC][G]{'+str(int(MinG)-k-1)+',}'
        if k!=int(MinG)-2:
            ImperfectTract += '|'
    ImperfectTract += ')'
# print(ImperfectTract)

if not NoExtLoop:
    ExtremeLoop='(?P<ExtremeLoop>\w{'+MinExtLL+','+MaxExtLL+'})'

finalRegex='([G]{'+MinG+',} | )





# chosenRegex='( [G]{'+MinG+',} | (?P<imperfectTract>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))(((\w{'+str(mloop-iloop)+','+str(mloop)+'})|(?P<ExtremeLoop>\w{1,'+str(eloop)+'}))(?(imperfectTract)([G]{'+str(ifull)+',})|([G]{'+str(ifull)+',}|(?P<imperfectTract>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))))((?(ExtremeLoop)\w{'+str(mloop-iloop)+','+str(mloop)+'}|(\w{'+str(mloop-iloop)+','+str(mloop)+'}|(?P<ExtremeLoop>\w{1,'+str(eloop)+'})))(?(imperfectTract)([G]{'+str(ifull)+',})|([G]{'+str(ilong)+',}|(?P<imperfectTract>([G]{'+str(ifull)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))))((?(ExtremeLoop)\w{'+str(mloop-iloop)+','+str(mloop)+'}|(\w{'+str(mloop-iloop)+','+str(mloop)+'}|(?P<ExtremeLoop>\w{1,'+str(eloop)+'})))(?(imperfectTract)([G]{'+str(ifull)+',})|([G]{'+str(ilong)+',}|(?P<imperfectTract>([G]{'+str(ilong)+',}[AUTC][G]{'+str(ishort)+',}|[G]{'+str(ishort)+',}[AUTC][G]{'+str(ilong)+',})))))'  # status for reference dataset: MCC:0.849 precision 95.1% #has too short limitation as well for loops
