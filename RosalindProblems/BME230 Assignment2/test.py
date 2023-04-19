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
                bestlist = [i]
                bestnow = in_str[i:i+self.k]
                for j in range(1,(self.w-self.k)+1):
                    curr = in_str[i+j:i+j+self.k]
                    if curr<bestnow:
                        bestnow = curr
                        bestlist=[i+j]
                    elif curr==bestnow:
                        bestlist.append(i+j)
                prevstring = bestnow
                prev_pos = bestlist[-1]
                for j in bestlist:
                    yield(bestnow,j)
            else:
                curr = in_str[i+(self.w-self.k):i+self.w]
                if curr<=prevstring:
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

def main():
    minimizerIndex = MinimizerIndexer("TACCCCTCAGATGCTTAAGC", w=5, k=3, t=1000)
    print(minimizerIndex.minimizerMap)

if __name__ == '__main__':
    main()