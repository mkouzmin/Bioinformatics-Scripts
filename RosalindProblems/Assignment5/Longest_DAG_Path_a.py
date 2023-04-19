# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys
class CommandLine():
    def __init__(self, inOpts=None):

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program reads nodes into Directed Acyclic Graph, and outputs longest path in Graph',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-i', '--input', type=str, default=None, action='store',
                                 help='Input File')
        self.parser.add_argument('-o', '--output', type=str, default=None, action='store',
                                 help='Output File')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''

    def __init__(self, msg):
        self.msg = msg


########################################################################
# My classes
# Here are the original classes
#
#
########################################################################
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
        self.graph[self.startnode].bestvalue = 0
        waiting = collections.deque()
        for k,v in self.graph.items():
            if v.incount == 0:
                waiting.append(k)
        if self.startnode not in waiting:
            waiting.appendleft(self.startnode)
        endwait = False
        while waiting and endwait == False: # iterate while nodes left without incoming edges unacounted for.
            currnode = waiting.popleft()
            if currnode == self.endnode:
                endwait = True
                break
            for i in self.graph[currnode].frontedges:
                if self.graph[i.end].incount < 0:# if an edge is followrd more than once, graph not acyclic
                    if self.graph[currnode].bestvalue != float("-inf"): # if loopin possible path location, return error
                       raise Usage("this is not an acyclic graph")
                    else:
                        continue # if  loop not in path section of graph, continue analysis with other edges from loop
                self.graph[i.end].incount-=1
                valuable = False # use to evaluate edges with positive values first
                endn = False #use to evaluate endnode immediately upon following all edges to it
                if self.graph[i.end].bestvalue < self.graph[currnode].bestvalue + i.weight:
                    if self.graph[i.end] == self.startnode:
                        continue #skip evaluation of edges going into startnode
                    else:
                        self.graph[i.end].bestvalue = self.graph[currnode].bestvalue + i.weight
                        valuable = True
                        self.graph[i.end].besttonode = currnode
                        if self.graph[i.end] == self.endnode:
                            endn = True
                if self.graph[i.end].incount == 0:
                    if valuable:
                        waiting.appendleft(i.end)
                        if endn:
                            break
                    else:
                        waiting.append(i.end)
        self.bestpathvalue = self.graph[self.endnode].bestvalue
        evalnode = self.endnode
        while evalnode != self.startnode:
            self.bestpathrev.append(evalnode)
            evalnode = self.graph[evalnode].besttonode
        self.bestpathrev.append(self.startnode)

########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(
            myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    try:
        readf = myCommandLine.args.input
        if myCommandLine.args.input is None:
            read = sys.stdin
        else:
            read = open(readf, "r")
        if myCommandLine.args.output is None:
            writef = sys.stdout
        else:
            writef = open(myCommandLine.args.output, "w")
        startnode = read.readline()[:-1]
        endnode = read.readline()[:-1]
        DAGraph = dict()
        for line in read:
            if line.strip() is not "":
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
            else:
                break
        toSolve = SolveGraph(DAGraph, startnode, endnode)
        toSolve.SolveG()
        writef.write(str(int(toSolve.bestpathvalue)) + "\n")
        writef.write(toSolve.bestpathrev[-1])
        for i in reversed(toSolve.bestpathrev[:-1]):
            writef.write("->" + i)

    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()

