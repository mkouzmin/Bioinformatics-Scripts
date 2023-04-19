import sys
import collections

class Edge:
    def __init__(self, start, end, weight):
        self.start = start
        self.end = end
        self.weight = weight

class Node:
    def __init__(self):
        self.frontedges = list()
        self.bestvalue = float("-inf")
        self.incount = 0
        self.besttonode = str()

class SolveGraph: # finds best path, returns list of nodes in reverse order
    def __init__(self, graph, startnode, endnode):
        self.graph = graph
        self.startnode = startnode
        self.endnode = endnode
        self.bestpathvalue = 0
        self.bestpathrev = list()
    def SolveG(self):
        self.graph[startnode].bestvalue = 0
        waiting = collections.deque()
        for k,v in self.graph.items():
            if v.incount == 0:
                waiting.append(k)
        if startnode not in waiting:
            waiting.appendleft(startnode)
        endwait = False
        while waiting and endwait == False: # iterate while nodes left without incoming edges unacounted for.
            currnode = waiting.popleft()
            if currnode == endnode:
                endwait = True
                break
            for i in self.graph[currnode].frontedges:
                if self.graph[i.end].incount < 0:# if an edge is followrd more than once, graph not acyclic
                    if self.graph[currnode].bestvalue != float("-inf"): # if loopin possible path location, return error
                        print("this is not an acyclic graph", file=sys.stderr)
                        quit()
                    else:
                        continue # if  loop not in path section of graph, continue analysis with other edges from loop
                self.graph[i.end].incount-=1
                valuable = False
                if self.graph[i.end].bestvalue < self.graph[currnode].bestvalue + i.weight:
                    if self.graph[i.end] == startnode:
                        continue #skip evaluation of edges going into startnode
                    else:
                        self.graph[i.end].bestvalue = self.graph[currnode].bestvalue + i.weight
                        valuable = True
                        self.graph[i.end].besttonode = currnode
                if self.graph[i.end].incount == 0:
                    if valuable:
                        waiting.appendleft(i.end)
                    else:
                        waiting.append(i.end)
        self.bestpathvalue = self.graph[endnode].bestvalue
        evalnode = endnode
        while evalnode != startnode:
            self.bestpathrev.append(evalnode)
            evalnode = self.graph[evalnode].besttonode
        self.bestpathrev.append(startnode)



readf = "input.txt"
read = open(readf, "r")
writef = open("output.txt", "w")
startnode = read.readline()[:-1]
endnode = read.readline()[:-1]
DAGraph = dict()
for line in read:
    words = line.split("->")
    fromnode = words[0].strip()
    words2 = words[1].strip().split(":")
    tonode = words2[0].strip()
    edgeweight = float(words2[1].strip())
    newedge = Edge(fromnode, tonode, edgeweight)
    DAGraph[fromnode] = DAGraph.get(fromnode, Node())
    DAGraph[fromnode].frontedges.append(newedge)
    DAGraph[tonode] = DAGraph.get(tonode, Node())
    DAGraph[tonode].incount += 1
toSolve = SolveGraph(DAGraph,startnode,endnode)
toSolve.SolveG()
writef.write(str(int(toSolve.bestpathvalue)) + "\n")
writef.write(toSolve.bestpathrev[-1])
for i in reversed(toSolve.bestpathrev[:-1]):
    writef.write("->" + i)
