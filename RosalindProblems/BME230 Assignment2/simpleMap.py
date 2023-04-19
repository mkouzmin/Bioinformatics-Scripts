import array
import sys
import numpy
#import pysam
import argparse
import logging
import time

logger = logging.getLogger()

"""See the comments below to see the code you need to complete.
"""


class MinimizerIndexer(object):
    """ Simple minimizer based substring-indexer.

    Please read: https://doi.org/10.1093/bioinformatics/bth408

    Related to idea of min-hash index and other "sketch" methods.
    """

    def __init__(self, targetString, w, k, t):
        """ The target string is a string/array of form "[ACGT]*".

        Stores the lexicographically smallest k-mer in each window of length w, such that w >= k positions. This
        smallest k-mer is termed a minmer.

        If a minmer occurs more than t times then it is omitted from the index.
        """

        self.targetString = targetString
        self.w = w
        self.k = k
        self.t = t  # If a minmer occurs more than t times then its entry is removed from the index
        # This is a heuristic to remove repetitive minmers that would create many spurious alignments between
        # repeats

        # Hash of minmers to query locations, stored as a map whose keys
        # are minmers and whose values are lists of the start indexes of
        # occurrences of the corresponding minmer in the targetString,
        # sorted in ascending order of index in the targetString.
        #
        # For example if k = 2 and w = 4 and targetString = "GATTACATTT"
        #
        # GATTACATTT
        # GATT (AT)
        #  ATTA (AT)
        #   TTAC (AC)
        #    TACA (AC)
        #     ACAT (AC)
        #      CATT (AT)
        #       ATTT (AT)
        #
        # then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        self.minimizerMap = {}

        todelete = set()
        for a,b in self.iter_minmer(targetString):
            if a not in todelete:
                l=self.minimizerMap.get(a,None)
                if l:
                    if len(l) >= t:
                        todelete.add(a)
                        del self.minimizerMap[a]
                    else:
                        l.append(b)
                else:
                    self.minimizerMap[a] = [b]

    def iter_minmer(self,in_str):
        prevstring = ''
        prev_pos = -1
        for i in range((len(in_str)-self.w)+1):
            if prev_pos<i:
                best_pos = i
                bestnow = in_str[i:i+self.k]
                for j in range(1,(self.w-self.k)+1):
                    curr = in_str[i+j:i+j+self.k]
                    if curr<bestnow:
                        bestnow = curr
                        best_pos = i+j
                prevstring = bestnow
                prev_pos = best_pos
                yield(bestnow,best_pos)
            else:
                curr = in_str[i+(self.w-self.k):i+self.w]
                if curr<prevstring:
                    yield (curr,i+(self.w-self.k))
                    prevstring = curr
                    prev_pos = i+(self.w-self.k)


        # Code to complete to build index - you are free to define additional functions

    def getMatches(self, searchString):
        """ Iterates through search string finding minmers in searchString and
        yields their list of occurrences in targetString, each as a pair of (x, (y,)*N),
        where x is the index searchString and y is an occurrence in targetString.

        For example if k = 2 and w = 4 and targetString = "GATTACATTT" and searchString = "GATTTA"
        then self.minimizerMap = { "AT":(1,6), "TA":(3,), "AC":(4,), "CA":(5,), }
        and getMatches will yield the following sequence:
        (1, (1,6)), (4, (3,))

        You will need to use the "yield" keyword
        """
        # Code to complete - you are free to define additional functions
        for a,b in self.iter_minmer(searchString):
            if a in self.minimizerMap:
                yield(b,tuple(self.minimizerMap[a]))


class SeedCluster:
    """ Represents a set of seeds between two strings.
    """

    def __init__(self, seeds):
        """ Seeds is a list of pairs [ (x_1, y_1), (x_2, y_2), ..., ], each is an instance of a seed
        (see static cluster seeds method below: static methods: https://realpython.com/blog/python/instance-class-and-static-methods-demystified/)
        """
        seeds = list(seeds)
        seeds.sort()
        self.seeds = seeds
        # Gather the minimum and maximum x and y coordinates
        self.minX = seeds[0][0]
        self.maxX = seeds[-1][0]
        ys = map(lambda (x, y): y, seeds)
        self.minY = min(ys)
        self.maxY = max(ys)

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.

        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of shared k-mer in both strings, termed a *seed*, such that the k-mer
        occurrence starts at position x in the first string and starts at position y_i in the second string.

        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.

        Consider a graph in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).

        The return value is a Python set of components, each component is a SeedCluster object.

        (QUESTION 1): The clustering of seeds is very simplistic. Can you suggest alternative strategies by
        which the seeds could be clustered, and what the potential benefits such alternative strategies could
        have? Consider the types of information you could use. Since we need to find points where both xs and ys are close to each other,
         we could use a k-means strategy to cluster seeds if we knew approximately how many groups to expect on average.
         We could also use a density-based clustering to leave out points that are not within a set distnace to a set number of other points.
         This would allow us to exclude some outlying minmersfrom our clustering sets,
          and potentially give us more and smaller clusters if they are only joined by a small amount of points in-between
        """

        # Code to complete - you are free to define other functions as you like
        clusters = set()
        for xcoord,b in seeds:
            for ycoord in b:
                sub_clusters = set(a for a in clusters if a.can_append(xcoord,ycoord,l))
                if len(sub_clusters)>0:
                    elem = sub_clusters.pop()
                    clusters.difference_update(sub_clusters)
                    for a in sub_clusters:
                        elem.add_cluster(a)
                    elem.append_element(xcoord,ycoord)
                else:
                    clusters.add(SeedCluster([(xcoord,ycoord)]))
        return clusters

    def append_element(self,xcoord,ycoord):
        self.seeds.append((xcoord,ycoord))
        self.maxX = xcoord
        self.minY = min(self.minY,ycoord)
        self.maxY = max(self.maxY,ycoord)

    def add_cluster(self,other):
        i = 0
        j = 0
        selflen = len(self.seeds)
        otherlen = len(other.seeds)
        new_seeds = list()
        while i < selflen and j<otherlen:
            if self.seeds[i]<other.seeds[j]:
                new_seeds.append(self.seeds[i])
                i+=1
            else:
                new_seeds.append(other.seeds[j])
                j+=1
        new_seeds+=self.seeds[i:]
        new_seeds+=other.seeds[j:]
        self.seeds = new_seeds
        self.minX = min(self.minX,other.minX)
        self.maxX = max(self.maxX,other.maxX)
        self.minY = min(self.minY,other.minY)
        self.maxY = max(self.maxY,other.maxY)




    def can_append(self, xcoord,ycoord,l):
        for i in reversed(self.seeds):
            if i[0] < xcoord -l:
                return False
            else:
                if i[1]>=ycoord - l and i[1] <= ycoord + l:
                    return True

class SmithWaterman(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3):
        """ Finds an optimal local alignment of two strings.

        Implements the Smith-Waterman algorithm:
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

        (QUESTION 2): The Smith-Waterman algorithm finds the globally optimal local alignment between to
        strings, but requires O(|string1| * |string2|) time. Suggest alternative strategies you could implement
        to accelerate the finding of reasonable local alignments. What drawbacks might such alternatives have?
        You could try splitting the querry string into sections and looking for the smaller subsections of the search string first,
         and extending from there. This would fail more easily, or be less accurate if the segment of the search string
        we are attempting to find never appears in the target string. We could also attempt to do a minmer allignment,
         or somehow trying to find important string sections, such as those that tend to repeat in a given search string
          directly in both strings, then aligning the other parts once we have those.
        """
        # Code to complete to compute the edit matrix
        # string 2 vertical, string 1 horizontal
        self.SWMat = numpy.zeros((len(string2)+1,len(string1)+1))
        self.gapScore = gapScore
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore
        self.string1 = string1
        self.string2 = string2
        self.maxi = 0
        self.maxj = 0
        self.maxScore = 0
        for i in range(len(self.string2)):
            for j in range(len(self.string1)):
                if string1[j] == string2[i]:
                    val = max(self.SWMat[i][j] + self.matchScore,self.SWMat[i][j+1] + self.gapScore,
                                          self.SWMat[i+1][j] + self.gapScore,0)
                    self.SWMat[i+1][j+1] = val
                    if val>self.maxScore:
                        self.maxScore = val
                        self.maxi = i
                        self.maxj = j
                else:
                    val = max(self.SWMat[i][j] + self.mismatchScore, self.SWMat[i][j + 1] + self.gapScore,
                                              self.SWMat[i + 1][j] + self.gapScore, 0)
                    self.SWMat[i+1][j+1] = val
                    if val>self.maxScore:
                        self.maxScore = val
                        self.maxi = i
                        self.maxj = j




    def getAlignment(self):
        """ Returns an optimal local alignment of two strings. Alignment
        is returned as an ordered list of aligned pairs.

        e.g. For the two strings GATTACA and CTACC an optimal local alignment
        is (GAT)TAC(A)
             (C)TAC(C)
        where the characters in brackets are unaligned. This alignment would be returned as
        [ (3, 1), (4, 2), (5, 3) ]
        """
        # Code to complete - generated by traceback through matrix to generate aligned pairs
        i = self.maxi
        j = self.maxj
        rev_list = list()
        value = self.maxScore
        while value != 0:
            if self.string1[j] == self.string2[i] and self.SWMat[i][j] + self.matchScore == value:
                rev_list.append((j, i))
                i = i-1
                j = j-1
            elif self.string1[j] != self.string2[1] and self.SWMat[i][j] + self.mismatchScore == value:
                rev_list.append((j, i))
                i = i-1
                j = j-1
            elif self.SWMat[i][j+1] + self.gapScore == value:
                i = i-1
            elif self.SWMat[i+1][j] + self.gapScore == value:
                j = j-1

            value = self.SWMat[i+1][j+1]
        rev_list.reverse()
        return rev_list
    def getMaxAlignmentScore(self):
        """ Returns the maximum alignment score
        """
        # Code to complete
        return self.maxScore



def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.

    Maps the string in both its forward and reverse complement orientations.

    (QUESTION 3): The code below is functional, but very slow. Suggest ways you could potentially accelerate it,
    and note any drawbacks this might have.
    We could shorten the lengths of large query and target strings by splitting them into sections and running them individually,
     eventually stitching the alignments together. This would require you to implement another algorithm for this stitching,
      and may therefore be even slower. It could also potentially require more data to do so.
      You could also run the Smith-Waterman only on the spaces between minimizers, but this would lead to inaccuracies
       (due to looking in the wrong locations)if we missed a minimizer.
    """
    bestAlignment = [None]

    def mapForwards(queryString):
        """ Maps the query string forwards
        """
        # Find seed matches, aka "aligned kmers"
        seeds = list(minimizerIndex.getMatches(queryString))

        # For each cluster of seeds
        for seedCluster in SeedCluster.clusterSeeds(list(seeds), l=config.l):

            # Get substring of query and target to align
            queryStringStart = max(0, seedCluster.minX - config.c)  # Inclusive coordinate
            queryStringEnd = min(len(queryString), seedCluster.maxX + config.k + config.c)  # Exclusive coordinate
            querySubstring = queryString[queryStringStart:queryStringEnd]

            targetStringStart = max(0, seedCluster.minY - config.c)  # Inclusive coordinate
            targetStringEnd = min(len(targetString), seedCluster.maxY + config.k + config.c)  # Exclusive coordinate
            targetSubstring = targetString[targetStringStart:targetStringEnd]

            # print "target_aligning", targetStringStart, targetStringEnd, targetSubstring
            # print "query_aligning", queryStringStart, queryStringEnd, querySubstring

            # Align the genome and read substring
            alignment = SmithWaterman(targetSubstring, querySubstring,
                                      gapScore=config.gapScore,
                                      matchScore=config.matchScore,
                                      mismatchScore=config.mismatchScore)

            # Update best alignment if needed
            if bestAlignment[0] == None or alignment.getMaxAlignmentScore() > bestAlignment[0].getMaxAlignmentScore():
                bestAlignment[0] = alignment

        return bestAlignment

    def reverseComplement(string):
        """Computes the reverse complement of a string
        """
        rMap = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(rMap[i] for i in string[::-1])

    # Run mapping forwards and reverse
    mapForwards(queryString)
    mapForwards(reverseComplement(queryString))

    return bestAlignment[0]


class Config():
    """ Minimal configuration class for handing around parameters
    """

    def __init__(self):
        self.w = 30
        self.k = 20
        self.t = 10
        self.l = 30
        self.c = 100
        self.gapScore = -2
        self.matchScore = 3
        self.mismatchScore = -3
        self.logLevel = "INFO"


def main():
    # Read parameters
    config = Config()

    # Parse the inputs args/options
    parser = argparse.ArgumentParser(usage="target_fasta query_fastq [options]", version="%prog 0.1")

    parser.add_argument("target_fasta", type=str,
                        help="The target genome fasta file.")
    parser.add_argument("query_fastq", type=str,
                        help="The query sequences.")

    parser.add_argument("--w", dest="w", help="Length of minimizer window. Default=%s" % config.w, default=config.w)
    parser.add_argument("--k", dest="k", help="Length of k-mer. Default=%s" % config.k, default=config.k)
    parser.add_argument("--t", dest="t", help="Discard minmers that occur more frequently "
                                              "in the target than t. Default=%s" % config.w, default=config.w)
    parser.add_argument("--l", dest="l", help="Cluster two minmers into the same cluster if within l bases of"
                                              " each other in both target and query. Default=%s" % config.l,
                        default=config.l)
    parser.add_argument("--c", dest="c", help="Add this many bases to the prefix and suffix of a seed cluster in the"
                                              " target and query sequence. Default=%s" % config.c, default=config.c)
    parser.add_argument("--gapScore", dest="gapScore", help="Smith-Waterman gap-score. Default=%s" %
                                                            config.gapScore, default=config.gapScore)
    parser.add_argument("--matchScore", dest="matchScore", help="Smith-Waterman match-score. Default=%s" %
                                                                config.gapScore, default=config.gapScore)
    parser.add_argument("--mismatchScore", dest="mismatchScore", help="Smith-Waterman mismatch-score. Default=%s" %
                                                                      config.mismatchScore,
                        default=config.mismatchScore)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" %
                                                       config.logLevel, default=config.logLevel)

    options = parser.parse_args()

    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)

    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")

    startTime = time.time()

    # Parse the target sequence and read the first sequence
    with pysam.FastaFile(options.target_fasta) as targetFasta:
        targetString = targetFasta.fetch(targetFasta.references[0])
    logger.info("Parsed target string. Length: %s" % len(targetString))

    # Build minimizer index
    minimizerIndex = MinimizerIndexer(targetString.upper(), w=options.w, k=options.k, t=options.t)
    minmerInstances = sum(map(len, minimizerIndex.minimizerMap.values()))
    logger.info("Built minimizer index in %s seconds. #minmers: %s, #minmer instances: %s" %
                ((time.time() - startTime), len(minimizerIndex.minimizerMap), minmerInstances))

    # Open the query files
    alignmentScores = []  # Array storing the alignment scores found
    with pysam.FastqFile(options.query_fastq) as queryFastq:
        # For each query string build alignment
        for query, queryIndex in zip(queryFastq, xrange(sys.maxint)):
            print(queryIndex)
            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.getMaxAlignmentScore()
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" %
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore))
            # Comment this out to test on a subset
            # if queryIndex > 100:
            #    break

    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" %
                    (time.time() - startTime, float(sum(alignmentScores)) / len(alignmentScores)))


if __name__ == '__main__':
    main()