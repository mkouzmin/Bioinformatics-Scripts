# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys


class CommandLine():
    def __init__(self, inOpts=None):

        import argparse

        self.parser = argparse.ArgumentParser(
            description='Program finds hidden path that gives best probability of a given string',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-i', '--input', type=str, default=None, action='store',
                                 help='Input File')
        self.parser.add_argument('-o', '--output', type=str, default=None, action='store', help='Output File')
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
import numpy
import collections
class Markov_HPath_Matrix:# use to find probability of a given path
    def __init__(self, statelist):#creates matrix, fills with zeros
        self.statelist = statelist
        self.corrDict = dict() #dictionary of state to matrix position
        count = 0
        for i in statelist:
            self.corrDict[i] = count
            count+=1
        self.Mmatrix = numpy.zeros((count,count))

class Markov_Text_Matrix:
    def __init__(self, textstatelist, hiddenstatelist):  # creates matrix, fills with zeros
        self.textstatelist = textstatelist
        self.corrtextDict = dict()  # dictionary of text state to matrix position
        textcount = 0
        for i in textstatelist:
            self.corrtextDict[i] = textcount
            textcount += 1
        self.hiddenstatelist = hiddenstatelist
        self.corrhiddenDict = dict()  # dictionary of text state to matrix position
        hiddencount = 0
        for i in hiddenstatelist:
            self.corrhiddenDict[i] = hiddencount
            hiddencount += 1
        self.Mmatrix = numpy.zeros((textcount, hiddencount))



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
        self.graph[self.startnode].bestvalue = 1
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
                if self.graph[i.end].incount < 0:# if an edge is followed more than once, graph not acyclic
                    if self.graph[currnode].bestvalue != float("-inf"): # if loop in possible path location, return error
                       raise Usage("this is not an acyclic graph")
                    else:
                        continue # if  loop not in path section of graph, continue analysis with other edges from loop
                self.graph[i.end].incount-=1
                valuable = False # use to evaluate edges with positive values first
                endn = False #use to evaluate endnode immediately upon following all edges to it
                if self.graph[i.end].bestvalue < self.graph[currnode].bestvalue * i.weight:
                    if self.graph[i.end] == self.startnode:
                        continue #skip evaluation of edges going into startnode
                    else:
                        self.graph[i.end].bestvalue = self.graph[currnode].bestvalue * i.weight
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


class Markov_Matrices:
    def __init__(self, textstatelist, hiddenstatelist):  # creates matrix, fills with zeros
        self.textstatelist = textstatelist
        self.hiddenstatelist = hiddenstatelist
        self.M_P_Matrix = Markov_HPath_Matrix(hiddenstatelist)
        self.M_T_Matrix = Markov_Text_Matrix(textstatelist, hiddenstatelist)

    def Solve_Markov(self, text):
        startnode = 0
        DAGraph = dict()
        DAGraph[0] = Node()
        endnode = len(self.hiddenstatelist) * (len(text)) + 1
        for i in range(len(text) + 1):
            kcount = 1
            for k in self.hiddenstatelist:
                fromnode = len(self.hiddenstatelist) * (i - 1) + kcount
                if fromnode < 0:  # skip iteration over all but one repetition of start node
                    kcount += 1
                    continue
                lcount = 1
                for l in self.hiddenstatelist:
                    tonode = len(self.hiddenstatelist) * (i) + lcount
                    if tonode > len(self.hiddenstatelist) * (
                    len(text)) + 1:  # if at endnode, skip iteration over repetitions of edges
                        lcount += 1
                        continue
                    if fromnode <= 0:  # if at start node - should only happen once
                        fromnode = 0
                        edgeweight = 1 / len(self.hiddenstatelist) * \
                                     self.M_T_Matrix.Mmatrix.item(self.M_T_Matrix.corrtextDict[text[i]],
                                                                  self.M_T_Matrix.corrhiddenDict[l])
                    else:
                        if tonode >= len(self.hiddenstatelist) * (
                        len(text)) + 1:  # if at end node - should only happen once
                            edgeweight = 1
                        else:
                            edgeweight = self.M_P_Matrix.Mmatrix.item(self.M_P_Matrix.corrDict[l],
                                                                      self.M_P_Matrix.corrDict[k]) * \
                                         self.M_T_Matrix.Mmatrix.item(self.M_T_Matrix.corrtextDict[text[i]],
                                                                      self.M_T_Matrix.corrhiddenDict[l])
                    newedge = Edge(fromnode, tonode, edgeweight)
                    DAGraph[fromnode] = DAGraph.get(fromnode, Node())
                    DAGraph[fromnode].frontedges.append(newedge)
                    DAGraph[tonode] = DAGraph.get(tonode, Node())
                    DAGraph[tonode].incount += 1
                    lcount += 1
                kcount += 1
        toSolve = SolveGraph(DAGraph, startnode, endnode)
        toSolve.SolveG()
        b_path = str("")
        for i in reversed(toSolve.bestpathrev[1:-1]):
            b_path = b_path + self.hiddenstatelist[(i-1)%len(self.hiddenstatelist)]
        return b_path

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
        if readf is None:
            read = sys.stdin
        else:
            read = open(readf, "r")
        if myCommandLine.args.output is None:
            writef = sys.stdout
        else:
            writef = open(myCommandLine.args.output, "w")
        text = read.readline()[:-1]
        read.readline()
        text_state_string = read.readline()
        text_state_list = text_state_string.strip().split()
        read.readline()
        path_state_string = read.readline()
        path_state_list = path_state_string.strip().split()
        Mmatrices = Markov_Matrices(text_state_list,path_state_list)#generate empty matrices
        read.readline()
        read.readline()
        for vertpos in range(len(path_state_list)):#vertpos is vertical positionn in prevstate-state matrix
            line = read.readline()
            lineprob = line.strip().split()
            horizpos = 0
            for prob in lineprob[1:]:
                Mmatrices.M_P_Matrix.Mmatrix[horizpos, vertpos] = prob
                horizpos += 1

        line = read.readline()
        line = read.readline()
        for vertpos in range(len(path_state_list)):#read text-state Matrix into program
            line = read.readline()
            lineprob = line.strip().split()
            horizpos = 0
            for prob in lineprob[1:]:
                Mmatrices.M_T_Matrix.Mmatrix[horizpos, vertpos] = prob
                horizpos += 1

        b_path = Mmatrices.Solve_Markov(text)
        writef.write(b_path)

    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()

