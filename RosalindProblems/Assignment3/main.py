
# Profiles are stored as dictionaries of lists:
# Motifs stored as list of strings
# Dna is stored as list of strings
import math
import fastaReader
import argparse
import sys
import random

def ProfFunc(motifs, mlen, mcount): #profile function, stores counts, not probability, given motifs, motif individual length, and motif total count(to avoid repeat counting)
    profile = {l:[0 for i in range(0,mlen)]for l in "ACGT"}
    for c in range(0,mcount):
        for n in range(0, mlen):
            b = motifs[c][n]
            if b in "ACGT":
                profile[b][n] += 1
    return profile

def OptMotif(Profile, mlen): #generate all concensus motifs based on profile and motif length
    if (mlen == 0):
        yield ""
    else:
        for c in MaxBase(Profile, mlen):
            for m in OptMotif(Profile,mlen - 1):
                yield m + c

def MaxBase(Profile, pos): # find maximaly likely bases in profile position
    maxbases = list()
    maxbaseprob = 0
    for c in Profile.items():
        prob = c[1][pos-1]
        if prob > maxbaseprob:
            maxbaseprob = prob
            maxbases = [c[0]]
        elif prob == maxbaseprob:
            maxbases.append(c[0])
    return maxbases

def ScoreMotif(pseuprofile, dna, motifstart, mcount, mlen): # scores motif according to profile with pseudocounts accounted for
    score = 1
    for i in range(0, mlen):
        try:
            score = score * pseuprofile[dna[motifstart+i]][i]
        except KeyError:
            score = score * 1/mcount # account for non-standard base pairs, by giving them a value of 0: 1/mcount
    return score



def MotifFunc(pseuprofile, dna, mcount, mlen, dnasize): #stores motifs in list given profile and dna, mcount - motif total count, mlen - motif individual length, pseudocount, and dnasize - size of dna strands (50 default)
    motifs = list()
    for s in dna: # find highest-scoring motif in line
        linescore = 0
        max_motif_pos = 0 # give default value to max_motif_pos in line - possible change required
        for n in range(0, dnasize-mlen):
            mscore = ScoreMotif(pseuprofile, s, n, mcount, mlen)
            if mscore > linescore:
                linescore = mscore
                max_motif_pos = n # moc_motif recorded only as start position
            elif mscore == linescore:
                if random.getrandbits(1):# implement random number generation if 2 motifs have equal scores
                    linescore = mscore
                    max_motif_pos = n
        motifs.append(s[max_motif_pos:max_motif_pos + mlen])# add maximum motif in dna string s to motifs
    return motifs

def ScoreFunc(pseuprofile, mlen):# score a profile made of probabilities w/ pseudocounts accounted for
    score = 0
    for n in range (0, mlen):
        for l in "ACGT":
            p = pseuprofile[l][n]
            if p != 0:
                score = score + (p * math.log(p, 2))
    score = -score
    return score


#begin parser
parser = argparse.ArgumentParser(
    description='Program prolog - Program calculates motifs in each line closest to concensus motif and scores the result',
    add_help=True,  # default is True
    prefix_chars='-',
    usage='%(prog)s [options] <input >output'
    )
parser.add_argument('-i', '--iterations', type=int, default = 10000, action='store',help='number oftimes to run search')
parser.add_argument('-p', '--pseudo', type=int, default = 1, action='store', help='pseudocount to use')
parser.add_argument('-k', '--motifLen', type=int, default = 13, action='store',help='length of motif')
parser.add_argument('-f', '--file', type = str, default = '', action = 'store', help='Input File')
parser.add_argument('-r', '--random', type = bool, default = 0, action = 'store', help='Randomize DNA')
parser.add_argument('-m', '--motif', type = bool, default = 0, action = 'store', help='print motifs')


args = parser.parse_args()

# begin main program

dna = []
fr = fastaReader.FastAreader(args.file)
for header,line in fr.readFasta():
    print(header, line)# remove in final
    dna.append(line)
random.seed()
if args.random:
    random.shuffle(dna)
iterations = args.iterations
pseudo = args.pseudo
motifLen = args.motifLen
motifCount = len(dna)
dnasize = len(dna[0])
for s in dna: # check all lengths of dna strands are equal
    if len(s) != dnasize:
        print("All dna strings must be same size", file=sys.stderr)
        quit()

bestiterscore = 2 * motifLen
bestitermotifs = list()
bestiterprof = dict()
for c in range(iterations): # repeat finding best motif iterations times
    motifs = list()
    for s in dna: # generate random mmotif for each string in dna
        i = random.randint(0,dnasize-motifLen)
        motifs.append(s[i:i+motifLen])
    bestmotifs = list(motifs) # give default values
    bestscore = 2 * motifLen
    bestprofile = dict()
    score = 2 * motifLen
    while True: # infinite loop
        profile = ProfFunc(motifs,motifLen,motifCount) #find profile of motifs
        for n in range(0, motifLen):# change profile to pseudocounts and probabilities
            for l in "ACGT":
                profile[l][n] = (profile[l][n] + pseudo) / (motifCount + pseudo * 4)
        motifs = MotifFunc(profile, dna, motifCount, motifLen, dnasize) # find motifs for new profile
        score = ScoreFunc(profile, motifLen) # score new profile
        if score < bestscore: # see if new profile better than previous
            bestscore = score
            bestmotifs = list(motifs)
            bestprofile = dict(profile)
        else: # if new profile not better, check best score for iteration against best overall score, end iteration
            if bestiterscore > bestscore:
                bestiterscore = bestscore
                bestitermotifs = list(bestmotifs)
                bestiterprof = dict(profile)
            break
print("Best score in " + str(iterations) + " repetitions: " + str(bestiterscore))
print("Best motifs found:")
for i in range (0, motifCount):
    print(bestitermotifs[i])
print("Possible Consensus Motifs Found:")
for k in OptMotif(bestiterprof, motifLen):
    print(k)