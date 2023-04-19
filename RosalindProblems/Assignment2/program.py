#!/usr/bin/env python3

# Parser and main program scaffold taken from program.py on UCSC Canvas website - some commentary and former parts of scaffold have been deleted
########################################################################
# CommandLine
########################################################################

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output' 
                                             )
        self.parser.add_argument('-m', '--minMotif', type=int, default=3, choices=range(3, 9), action='store',
                            help='minimum length of motif')
        self.parser.add_argument('-M', '--maxMotif', type=int, default=8, choices=range(3, 9), action='store',
                            help='maximum length of motif')
        self.parser.add_argument('-c', '--cutoff', type=float, default=-5, action='store', help='Z-Score cutoff')
        self.parser.add_argument('-i', '--input', type=str, default = '', action='store', help='Input File')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
  

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg


########################################################################
# imports and helper functions
#  Here are the imports and helper functions/ generators for the main program
#
#
########################################################################

import math
import fastaReader


# generate all motifs of sizes minMotif to maxMotif
class SeqInit:
    def __init__(self, minsize, maxMotif):
        self.minsize = minsize
        self.maxMotif = maxMotif

    def InitSeq(self):
        sequences = {k: 0 for n in range(self.minsize, self.maxMotif) for k in self.MotifGen(n)}
        return sequences

    def MotifGen(self, n):
        if (n == 0):
            yield ""
        else:
            for c in "ATGC":
                for m in self.MotifGen(n - 1):
                    yield c + m

class Missing:
    def __init__(self, sequences, genomeLen, minMotif, maxMotif):
        self.minMotif = minMotif
        self.maxMotif = maxMotif
        self.genomeLen = genomeLen
        self.sequences = sequences

    # finds expected value and z-score using dictionary s, key k, and genome length glen
    def ScoreFunc(self, s, k, glen):
        denom = s[k[1:-1]]  # to prevent calculating twice
        if (denom > 0):
            expv = (s[k[1:]] * s[k[:-1]]) / denom
            n = glen - len(k) + 1  # of motifs of length k
            if (expv > 0 and expv < n):
                z_score = ((s[k] - expv) / math.sqrt(expv * (1 - (expv / n))))
            else:
                z_score = 0
        else:
            expv = 0
            z_score = 0
        return expv, z_score


    def DictToSortList (self):
        if self.minMotif < 3:
            raise Usage('minMotif must be >3')
        # return error
        else:

            # make a list of tuples of (key, actual count, expected count, Z-value) seqlist
            seqlist = [(k, self.sequences[k], *self.ScoreFunc(self.sequences, k, self.genomeLen)) for k in self.sequences.keys() if
                       len(k) >= self.minMotif]

            # sort seqlist into sortseqlist by reverse length of motif, then reverse count, then reverse z-score
            seqlist.sort(key=lambda x: (-len(x[0]), -x[1], -x[3]))

            return seqlist



########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.


    try:

        fr = fastaReader.FastAreader(myCommandLine.args.input)
        genome = ""
        seguences = list()

        minMotif = myCommandLine.args.minMotif  # must be at least 3
        maxMotif = myCommandLine.args.maxMotif  # up to 9
        cutoff = myCommandLine.args.cutoff
        totalgenomeLen = 0
        if minMotif < 3:
            raise Usage('minMotif must be >3')
        # return error
        else:
            minsize = minMotif - 2
            #find all lengths of sequences ranging from minsize to maxmotif
            seq = SeqInit(minsize,maxMotif)
            sequences = seq.InitSeq()
            for header, line in fr.readFasta(): #for all genomes in all lines of fasta:
                print("Header:", header, len(line))
                genome = line
                genomeLen = len(genome)
                totalgenomeLen += genomeLen
                # for all lengths of sequences ranging from minsize to max Motif, count # of occurences in each line - put in sequences dictionary
                # include non-standard base pairs
                for start in range(0, genomeLen):
                    for length in range(minsize, maxMotif):
                        if start + length <= genomeLen:
                            curr = genome[start:start + length]
                            sequences[curr] = sequences.get(curr, 0) + 1
            # end of input section
        MMotifs = Missing(sequences,totalgenomeLen,minMotif,maxMotif)
        seqlist = MMotifs.DictToSortList()
        # print all values in sortedseqlist with z-score greater than cutoff, in order given using assigned format
        for l in seqlist:
            if (l[3] <= cutoff):
                print('{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2f}'.format(l[0], l[1], l[2], l[3]))


    except Usage as err:
       print (err.msg)

if __name__ == "__main__":
    main()

