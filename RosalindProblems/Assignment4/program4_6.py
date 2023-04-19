class CalcPath:
    def __init__(self,dBGraph):
        self.dBGraph = dBGraph

    # takes dBGraph in form dictionary{istring:[[ostring1,count][ostring2,count]], istringcount},
    #  where istringcount is # of times istring is referenced elsewhere

    def findPath(self):
        # find start + end nodes
        start = ""
        end = ""
        startcount = 0
        endcount = 0
        for k, v in self.dBGraph.items():
            osum = sum(j[1] for j in v[0])  # count outgoing edges
            ival = v[1]
            if (osum > ival + 1 or ival > osum + 1):
                print("A path is impossible, difference between number of exits and entrances in node >1 ",
                      file=sys.stderr)
                quit()
            elif (osum == ival + 1):
                start = k
                startcount += 1
            elif (ival == osum + 1):
                end = k
                endcount += 1
        if startcount != endcount or startcount > 1:
            print("A path is impossible, more than one node each must be a start or an end", file=sys.stderr)
            quit()
        if startcount == 0:
            for k in self.dBGraph.keys():
                end = k
                break
            start = end
        # end of finding start, end
        # find path using found start,end
        now = end
        next = start
        path = list()
        pos = 0
        checkpos = 0
        Endpath = False
        while Endpath == False:  # while edges left
            path.insert(pos, next)  # insert last found edge
            deadend = True
            now = next
            if self.dBGraph[now][0]:  # check if list is empty: follow path till deadend, adding nodes to path
                i = self.dBGraph[now][0][-1]
                next = i[0]
                pos += 1
                deadend = False
                if i[1] != 1:
                    i[1] -= 1
                else:
                    self.dBGraph[now][0].pop()
            if (deadend):  # if no paths found from node
                Endpath = True  # end when goes through path without finding loops
                while checkpos < len(
                        path):  # go through graph, find first spot with free path, continue inserting values into path from there
                    i = path[checkpos]
                    if self.dBGraph[i][0]:  # check if list is empty / check if any edges left in node
                        j = self.dBGraph[i][0][-1]
                        next = j[0]
                        pos = checkpos + 1
                        Endpath = False
                        if j[1] != 1:
                            j[1] -= 1
                        else:
                            self.dBGraph[i][0].pop()
                        break
                    checkpos += 1
        return path



import sys
readf = "input.txt"
read = open (readf, "r")
writef = open("output.txt", "w")
dBGraph = dict()
length = 0
for line in read:
    words = line.split("->")
    k = words[0].strip()
    numlist = list(words[1].strip().split(","))
    countlist = list()
    for i in numlist:
        inlist = False
        for j in countlist:
            if i == j[0]:
                j[1] += 1
                inlist = True
                break
        if not inlist:
            countlist.append([i,1])

    for i in countlist:
        if i[0] in dBGraph:

            dBGraph[i[0]][1] += i[1]
        else:
            dBGraph[i[0]] = [list(),i[1]] #counts incoming edges before forming node in graph
    dBGraph[k] = [countlist,dBGraph.get(k, [list(),0])[1]] #form node in graph

#end graph input

cp = CalcPath(dBGraph)
path = cp.findPath()
writef.write(path[0])
for i in path[1:]:
    writef.write("->"+i)