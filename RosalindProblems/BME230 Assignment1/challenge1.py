import sys
import numpy

"""The following uses Python to challenge you to create an algorithm for finding
matches between a set of aligned strings. Minimal familiarity with Python is 
necessary, notably list and Numpy array slicing. 
"""

"""Problem 1.

Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length 
N. 

For integer 0<=i<=N we define an ith prefix sort as a lexicographic sort 
(here 0 precedes 1) of the set of ith prefixes: { x[:i] | x in X }.
Similarly an ith reverse prefix sort is a lexicographic sort of the set of
ith prefixes after each prefix is reversed.

Let A be an Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the 
index in X of the ith string ordered by jth reverse prefix. To break ties 
(equal prefixes) the ordering of the strings in X is used. 

Complete code for the following function that computes A for a given X.

Here X is a Python list of Python strings. 
To represent A we use a 2D Numpy integer array.

Example:

>>> X = getRandomX() #This is in the challenge1UnitTest.py file
>>> X
['110', '000', '001', '010', '100', '001', '100'] #Binary strings, M=7 and N=3
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> 

Hint:
Column j (0 < j <= N) of the matrix can be constructed from column j-1 and the 
symbol in each sequence at index j-1.  

Question 1: In terms of M and N what is the asymptotic cost of your algorithm? O(MN)
"""


def constructReversePrefixSortMatrix(X):
    # Creates the Mx(N+1) matrix
    A = numpy.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0]) + 1], dtype=int)
    M = len(X)
    N = 0 if len(X) == 0 else len(X[0])
    for i in range(0,M):#set left row order
        A[i][0] = i
    for j in range(0,N):  # for each column
        list0 = list()
        list1 = list()
        for i in range(0,M):  # for each value in a given column
            if X[int(A[i][j])][j] == '0':  # if jth term of indexed string = 0, add to first list, else add to second
                list0.append(int(A[i][j]))
            else:
                list1.append(int(A[i][j]))
        c = 0  # keep track of place in column
        for ind in list0:  # fill column with the 2 lists
            A[c][j + 1] = ind
            c += 1
        for ind in list1:
            A[c][j + 1] = ind
            c += 1

    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.

    return A


"""Problem 2: 

Following on from the previous problem, let Y be the MxN matrix such that for 
all 0 <= i < M, 0 <= j < N, Y[i,j] = X[A[i,j]][j].

Complete the following to construct Y for X. 

Hint: You can either use your solution to constructReversePrefixSortMatrix() 
or adapt the code from that algorithm to create Y without using 
constructReversePrefixSortMatrix().

Question 2: In terms of M and N what is the asymptotic cost of your algorithm? O(MN)
"""


def constructYFromX(X):
    # Creates the MxN matrix
    Y = numpy.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0])], dtype=int)
    A = numpy.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0]) + 1], dtype=int)
    M = len(X)
    N = 0 if len(X) == 0 else len(X[0])
    for i in range(0,M):
        A[i][0] = i
    for j in range(0,N):  # for each column
        list0 = list()
        list1 = list()
        for i in range(0,M):  # for each value in a given column
            if X[int(A[i][j])][j] == '0':  # if jth term of indexed string = 0, add to first list, else add to second
                list0.append(int(A[i][j]))
            else:
                list1.append(int(A[i][j]))
        c = 0  # keep track of place in column
        for ind in list0:  # fill column with the 2 lists for both matrices
            A[c][j + 1] = ind
            Y[c][j] = X[A[c][j]][j]# 1 of only 2 lines different
            c += 1
        for ind in list1:
            A[c][j + 1] = ind
            Y[c][j] = X[A[c][j]][j]
            c += 1

    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.

    return Y

"""Problem 3.

Y is a transformation of X. Complete the following to construct X from Y, 
returning X as a list of strings as defined in problem 1.
Hint: This is the inverse of X to Y, but the code may look very similar.

Question 3a: In terms of M and N what is the asymptotic cost of your algorithm?O(NM)

Question 3b: What could you use the transformation of Y for? We can immediately use this data for sequence matching, similarly to using the BWT, if we build an FM index for it.
Hint: consider the BWT.

Question 3c: Can you come up with a more efficient data structure for storing Y? Y can be stored with a list, or another matrix row that stores the number of 1s in each column, so that we do not have to count this and can reduce reconstruction time by O(NM). This can also be done by creating an FM index for each column, allowing potentially faster sequence matching. These choices would increase memory usage slighly, but reduce time needed to use this data 
"""


def constructXFromY(Y): #saves memory space - constructs in place, using 2 lists of length M instead of array
    # Creates the MxN matrix
    M = Y.shape[0]
    N = Y.shape[1]
    X= ["" for m in range (M)]# create list of M empty strings
    for j in reversed(range (0,N)):
        num0 = 0 #count number of 1s in each column
        A = ["" for m in range(M)]
        if j == N-1:
            for i in range(0,M):
                X[i] = Y[i][j]
        else:
            for i in range(0,M):#add string of given Y location to start of existing list item, then sort numerically
                if Y[i][j] == 0:
                    num0+=1
            num1s = 0  # count number of 1s in current column while building list items
            for i in range(0,M):
                if Y[i][j] == 0:
                    A[i] = "0" + str(X[i-num1s])
                else:
                    A[i] = "1" + str(X[num0 + num1s])
                    num1s = num1s+1
            X = A
    for i in range(len(X)):
        X[i] = str(X[i])


    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.

    return X


"""Problem 4.

Define the common suffix of two strings to be the maximum length suffix shared 
by both strings, e.g. for "10110" and "10010" the common suffix is "10" because 
both end with "10" but not both "110" or both "010". 

Let D be a Mx(N+1) Numpy integer array such that for all 1<=i<M, 1<=j<=N, 
D[i,j] is the length of the common suffix between the substrings X[A[i,j]][:j] 
and X[A[i-1,j]][:j].  

Complete code for the following function that computes D for a given A.

Example:

>>> X = getRandomX()
>>> X
['110', '000', '001', '010', '100', '001', '100']
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> D = constructCommonSuffixMatrix(A, X)
>>> D
array([[0, 0, 0, 0],
       [0, 1, 2, 2],
       [0, 1, 2, 3],
       [0, 1, 1, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 1, 3]])

Hints: 

As before, column j (0 < j <= N) of the matrix can be constructed from column j-1 
and thesymbol in each sequence at index j-1.

For an efficient algorithm consider that the length of the common suffix 
between X[A[i,j]][:j] and X[A[i-k,j]][:j], for all 0<k<=i is 
min(D[i-k+1,j], D[i-k+2,j], ..., D[i,j]).

Question 4: In terms of M and N what is the asymptotic cost of your algorithm?O(NM^2)
"""


def constructCommonSuffixMatrix(A, X):
    D = numpy.zeros(shape=A.shape, dtype=int)  # Creates the Mx(N+1) D matrix

    # Code to write - you're free to define extra functions
    # (inline or outside of this function) if you like.
    M = len(X)
    N = 0 if len(X) == 0 else len(X[0])
    for j in range(1, N+1):  # for each column
        for i in range(1,M):  # for each value in a given column, except initial 0
            if X[A[i][j]][j-1] == X[A[i-1][j]][j-1]:  # if char in current list item is value in prev list item, find the 2 strings in the previous column, report min match distance between them
                first = 0
                for k in range(0,M):
                    if A[k][j-1] == A[i][j] or A[k][j-1] == A[i-1][j]:
                        first = k
                        break
                second = first
                for l in range(first+1,M):
                    if A[l][j-1] == A[i-1][j]or A[l][j-1] == A[i][j]:
                        second = l
                        break
                minN = int(D[first+1][j-1])
                for l in range(first+2,second+1):#find minimum in interval between the 2 found numbers in previous row
                    if int(D[l][j-1]) < minN:
                        minN = int(D[l][j-1])
                D[i][j] = minN+1#time = M
            else:#else, current item value = 0
                D[i][j] = 0
    return D

"""Problem 5.

For a pair of strings X[x], X[y], a long match ending at j is a common substring
of X[x] and X[y] that ends at j (so that X[x][j] != X[y][j] or j == N) that is longer
than a threshold 'minLength'. E.g. for strings "0010100" and "1110111" and length
threshold 2 (or 3) there is a long match "101" ending at 5.

The following algorithm enumerates for all long matches between all substrings of
X, except for simplicity those long matches that are not terminated at
the end of the strings.

Question 5a: What is the asymptotic cost of the algorithm in terms of M, N and the
number of long matches?O(max(MN*#longmatches))

Question 5b: Can you see any major time efficiencies that could be gained by
refactoring? We could consider long matches in range over j starting only at minLength(i.e. reduce number of columns considered at start of matrix), and reduce time by O(M*minLength)

Question 5c: Can you see any major space efficiencies that could be gained by
refactoring? You could have constructCommonSuffixMatrix take Y as a parameter, and construct Y instead of A, increasing data compression

Question 5d: Can you imagine alternative algorithms to compute such matches?, The Naive approach would be to compute the suffix matches pairwise individually for each pair of strings. This would take a time of O(NM^2) and use O(M^2) space.
if so, what would be the asymptotic cost and space usage?
"""


def getLongMatches(X, minLength):
    assert minLength > 0

    A = constructReversePrefixSortMatrix(X)
    D = constructCommonSuffixMatrix(A, X)

    # For each column, in ascending order of column index
    for j in range(1, 0 if len(X) == 0 else len(X[0])):
        # Working arrays used to store indices of strings containing long matches
        # b is an array of strings that have a '0' at position j
        # c is an array of strings that have a '1' at position j
        # When reporting long matches we'll report all pairs of indices in b X c,
        # as these are the long matches that end at j.
        b, c = [], []

        # Iterate over the aligned symbols in column j in reverse prefix order
        for i in range(len(X)):
            # For each string in the order check if there is a long match between
            # it and the previous string.
            # If there isn't a long match then this implies that there can
            # be no long matches ending at j between sequences indices in A[:i,j]
            # and sequence indices in A[i:,j], thus we report all long matches
            # found so far and empty the arrays storing long matches.
            if D[i, j] < minLength:
                for x in b:
                    for y in c:
                        # The yield keyword converts the function into a
                        # generator - alternatively we could just to append to
                        # a list and return the list

                        # We return the match as tuple of two sequence
                        # indices (ordered by order in X) and coordinate at which
                        # the match ends
                        yield (x, y, j) if x < y else (y, x, j)
                b, c = [], []

            # Partition the sequences by if they have '0' or '1' at position j.
            if X[A[i, j]][j] == '0':
                b.append(A[i, j])
            else:
                c.append(A[i, j])

        # Report any leftover long matches for the column
        for x in b:
            for y in c:
                yield (x, y, j) if x < y else (y, x, j)

