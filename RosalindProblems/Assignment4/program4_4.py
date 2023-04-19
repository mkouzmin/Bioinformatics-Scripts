from collections import defaultdict
import sys
readf = sys.stdin.readline().strip()
read = open (readf, "r")
writef = sys.stdout
k = int(read.readline())
l = k-1
genome = read.readline()[:-1]
genomeLen = len(genome)
dBGraph = defaultdict(list)
for i in range (0, genomeLen-l):
    dBGraph[genome[i:i+l]].append(genome[i+1:i+l+1])
for i in sorted(dBGraph.items()):
    writef.write(i[0] + " -> ")
    writef.write(i[1][0])
    if (len(i[1]) > 1):
        for j in i[1][1:]:
            writef.write("," + j)
    writef.write("\n")