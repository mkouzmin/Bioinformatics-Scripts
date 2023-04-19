from collections import defaultdict
import sys
readf = sys.stdin.readline().strip()
read = open (readf, "r")
writef = sys.stdout
dBGraph = defaultdict(list)
Sequences = list()
for line in read:
    Sequences.append(line[:-1])
for i in Sequences:
    dBGraph[i[:-1]].append(i[1:])
for i in dBGraph.items():
    writef.write(i[0] + " -> ")
    writef.write(i[1][0])
    if (len(i[1]) > 1):
        for j in i[1][1:]:
            writef.write("," + j)
    writef.write("\n")