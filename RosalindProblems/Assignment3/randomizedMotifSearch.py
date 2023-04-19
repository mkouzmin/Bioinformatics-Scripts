#!/usr/bin/env python3
# scaffold structure taken from canvas website

# Some Notes on storage:
# Profiles are stored as dictionaries of lists:
# Motifs stored as list of strings
# Dna is stored as dictionary of strings to names
########################################################################
# import packages
########################################################################
import math
import fastaReader
import argparse
import random
########################################################################
# CommandLine
########################################################################
class CommandLine() :
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implements a parser to interpret the command line argv string using argparse.
        '''

        self.parser = argparse.ArgumentParser(description = 'Program prolog - Program calculates motifs in each line closest to concensus motif and scores the result',
                                             add_help = True,
                                             prefix_chars = '-',
                                             usage = '%(prog)s [options] <input >output'
                                             )
        self.parser.add_argument('-i', '--iterations', type=int, default=10000, action='store',
                            help='number oftimes to run search')
        self.parser.add_argument('-p', '--pseudo', type=int, default=1, action='store', help='pseudocount to use')
        self.parser.add_argument('-k', '--motifLen', type=int, default=13, action='store', help='length of motif')
        self.parser.add_argument('-f', '--file', type=str, default='', action='store', help='Input File')
        self.parser.add_argument('-r', '--random', type=bool, default=0, action='store', help='Input File')
        self.parser.add_argument('-m', '--motif', type=bool, default=0, action='store', help='print motifs')

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
#  helper functions
#  Here are the helper classes for the main program
#
#
#
########################################################################
class ConcMotifFind: # finds concensus motifs given profile, motifLen
    def __init__(self, Profile, mlen):
        self.Profile = Profile
        self.mlen = mlen

    def OptMotif(self): #generate all concensus motifs based on profile and motif length
        if (self.mlen == 0):
            yield ""
        else:
            for c in self.MaxBase(self.mlen):
                CMFind = ConcMotifFind(self.Profile, self.mlen-1)
                for m in CMFind.OptMotif():
                    yield m + c

    def MaxBase(self, pos): # find maximaly likely bases in profile position
        maxbases = list()
        maxbaseprob = 0
        for c in self.Profile.items():
            prob = c[1][pos-1]
            if prob > maxbaseprob:
                maxbaseprob = prob
                maxbases = [c[0]]
            elif prob == maxbaseprob:
                maxbases.append(c[0])
        return maxbases

class RandomSearch:
    def __init__(self, dna, pseudo, iterations, motiflen): # does main search
        self.dna = dna
        self.pseudo = pseudo
        self.iterations = iterations
        self.motifLen = motiflen
        self.motifCount = len(dna)
        self.dnasize = len(dna[0][1])
        for s in dna:  # check all lengths of dna strands are equal
            if len(s[1]) != self.dnasize:
                raise Usage("All dna strings must be same size")
        self.bestiterscore = 2 * self.motifLen
        self.bestitermotifs = list()
        self.bestiterprof = dict()

    def MotifSearch(self, keepmotifs): #main fuction of class
        for c in range(self.iterations): # repeat finding best motif iterations times
            motifs = list()
            for s in self.dna: # generate random motif for each string in dna
                i = random.randint(0,self.dnasize-self.motifLen)
                motifs.append(s[1][i:i+self.motifLen])
            bestmotifs = list(motifs) # give default values
            bestscore = 2 * self.motifLen
            bestprofile = dict()
            score = 2 * self.motifLen
            while True: # infinite loop
                profile = self.ProfFunc(motifs,self.motifLen,self.motifCount) #find profile of motifs
                for n in range(0, self.motifLen):# change profile to pseudocounts and probabilities
                    for l in "ACGT":
                        profile[l][n] = (profile[l][n] + self.pseudo) / (self.motifCount + self.pseudo * 4)
                motifs = self.MotifFunc(profile) # find motifs for new profile
                score = self.ScoreFunc(profile) # score new profile
                if score < bestscore: # see if new profile better than previous
                    bestscore = score
                    if(keepmotifs):
                        bestmotifs = list(motifs)
                    bestprofile = dict(profile)
                else: # if new profile not better, check best score for iteration against best overall score, end iteration
                    if self.bestiterscore > bestscore:
                        self.bestiterscore = bestscore
                        if (keepmotifs):
                            self.bestitermotifs = list(bestmotifs)
                        self.bestiterprof = dict(bestprofile)
                    break

    def ProfFunc(self, motifs, mlen, mcount): #profile function, stores counts, not probability, given motifs, motif individual length, and motif total count(to avoid repeat counting)
        profile = {l:[0 for i in range(0,mlen)]for l in "ACGT"}
        for c in range(0,mcount):
            for n in range(0, mlen):
                b = motifs[c][n]
                if b in "ACGT":
                    profile[b][n] += 1
        return profile


    def ScoreMotif(self, pseuprofile, dna, motifstart, mcount, mlen, pseudo): # scores motif according to profile with pseudocounts accounted for
        score = 1
        for i in range(0, mlen):
            try:
                score = score * pseuprofile[dna[motifstart+i]][i]
            except KeyError:
                score = score * (1/(mcount+(4*pseudo)))
        return score



    def MotifFunc(self, pseuprofile): #stores motifs in list given profile and dna, mcount - motif total count, mlen - motif individual length, pseudocount, and dnasize - size of dna strands (50 default)
        motifs = list()
        for s in self.dna: # find highest-scoring motif in line
            linescore = 0
            max_motif_pos = 0 # give default value to max_motif_pos in line - possible change required
            for n in range(0, self.dnasize-self.motifLen):
                mscore = self.ScoreMotif(pseuprofile, s[1], n, self.motifCount, self.motifLen, self.pseudo)
                if mscore > linescore:
                    linescore = mscore
                    max_motif_pos = n # moc_motif recorded only as start position
                elif mscore == linescore:
                    if random.getrandbits(1):# implement random number generation if 2 motifs have equal scores
                        linescore = mscore
                        max_motif_pos = n
            motifs.append(s[1][max_motif_pos:max_motif_pos + self.motifLen])# add maximum motif in dna string s to motifs
        return motifs

    def ScoreFunc(self, pseuprofile):# score a profile made of probabilities w/ pseudocounts accounted for
        score = 0
        for n in range (0, self.motifLen):
            for l in "ACGT":
                p = pseuprofile[l][n]
                if p != 0:
                    score = score + (p * math.log(p, 2))
        score = -score
        return score




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

        dna = list()
        fr = fastaReader.FastAreader(myCommandLine.args.file)
        random.seed()
        for header,line in fr.readFasta():
            if myCommandLine.args.random:
                dna.append((header[0: header.index(' ')], random.shuffle(line)))
            else:
                dna.append((header[0: header.index(' ')],line))
        iterations = myCommandLine.args.iterations
        pseudo = myCommandLine.args.pseudo
        motifLen = myCommandLine.args.motifLen
        rsearch = RandomSearch(dna,pseudo,iterations, motifLen)
        rsearch.MotifSearch(myCommandLine.args.motif)
        print("Best score in " + str(iterations) + " repetitions: " + str(rsearch.bestiterscore))
        if (myCommandLine.args.motif):
            print("Best motifs found:")
            for i in range (0, rsearch.motifCount):
                print(dna[i][0] + ": " + rsearch.bestitermotifs[i])
        print("Possible Consensus Motifs Found:")
        CMFind = ConcMotifFind(rsearch.bestiterprof, motifLen)
        for k in CMFind.OptMotif():#output concensus motif according to best profile
            print(k)




    except Usage as err:
       print (err.msg)

if __name__ == "__main__":
    main()
#   main(['-r'])  # this would make this program.py behave as written

