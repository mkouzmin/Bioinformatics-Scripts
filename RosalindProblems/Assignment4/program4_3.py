import sys
readf = sys.stdin.readline().strip()
read = open (readf, "r")
writef = sys.stdout
Sequences = list()
SeqDict = dict()
for line in read:
    Sequences.append(line[:-1])
for s in Sequences:
    for p in Sequences:
        if s[1:] == p[:-1]:
            SeqDict[s] = p
for i in SeqDict.items():
    writef.write(i[0] + " -> " + i[1] + "\n")