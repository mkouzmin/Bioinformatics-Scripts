import sys
readf = sys.stdin.readline().strip()
read = open (readf, "r")
writef = sys.stdout
First = True
combS = ""
for line in read:
    if First:
        combS = line[:-1]
        First = False
    else:
        combS = combS + line[-2:-1]
writef.write(combS)
