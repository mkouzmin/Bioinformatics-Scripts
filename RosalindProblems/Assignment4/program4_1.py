import sys
readf = sys.stdin.readline().strip()
read = open (readf, "r")
writef = sys.stdout
k = int(read.readline())
genome = read.readline()
genomeLen = len(genome)
for i in range (0, genomeLen-k):
    writef.write(genome[i:i+k]+"\n")
