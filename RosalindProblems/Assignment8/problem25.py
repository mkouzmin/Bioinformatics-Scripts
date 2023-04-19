# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys


class CommandLine():
    def __init__(self, inOpts=None):

        import argparse

        self.parser = argparse.ArgumentParser(
            description='Program implements forwards-backwards algorythm to solve the soft decoding problem',
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
import math
class Markov_Trans_Matrix:  # use to find probability of a given path
    def __init__(self, statelist):  # creates matrix, fills with zeros
        self.statelist = statelist
        self.corrDict = dict()  # dictionary of state to matrix position
        count = 0
        for i in statelist:
            self.corrDict[i] = count
            count += 1
        self.Mmatrix = numpy.zeros((count, count))

class Markov_Emit_Matrix:
    def __init__(self, emitstatelist, hiddenstatelist):  # creates matrix, fills with zeros
        self.emitstatelist = emitstatelist
        self.corremitDict = dict()  # dictionary of text state to matrix position
        emitcount = 0
        for i in emitstatelist:
            self.corremitDict[i] = emitcount
            emitcount += 1
        self.hiddenstatelist = hiddenstatelist
        self.corrhiddenDict = dict()  # dictionary of text state to matrix position
        hiddencount = 0
        for i in hiddenstatelist:
            self.corrhiddenDict[i] = hiddencount
            hiddencount += 1
        self.Mmatrix = numpy.zeros((hiddencount, emitcount))



class Edge:
    def __init__(self, start, end, weight):
        self.start = start
        self.end = end
        self.weight = weight

class Node:
    def __init__(self):
        self.frontedges = list()
        self.forwardtotalvalue = 0
        self.reversetotalvalue = 0
        self.backedges = list()
        self.incount = 0 # not reliable, used to keep track of edges unaccounted for
        self.outcount = 0

class SumGraph: # finds best path, returns list of nodes in reverse order
    def __init__(self, graph, startnode, endnode):
        self.graph = graph
        self.startnode = startnode
        self.endnode = endnode
        self.totalpathvalue = 0
    def SolveFG(self):
        self.graph[self.startnode].forwardtotalvalue = 1
        waiting = collections.deque()
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
                endn = False #use to evaluate endnode immediately upon following all edges to it
                if self.graph[i.end] == self.startnode:
                    continue #skip evaluation of edges going into startnode
                else:#add edge to next sum
                    self.graph[i.end].forwardtotalvalue = self.graph[i.end].forwardtotalvalue+\
                                                         (self.graph[currnode].forwardtotalvalue * i.weight)
                    if self.graph[i.end] == self.endnode:
                        endn = True
                if self.graph[i.end].incount == 0:
                    waiting.appendleft(i.end)
                    if endn:
                        break
        self.totalpathvalue = self.graph[self.endnode].forwardtotalvalue

    def SolveRG(self):
        self.graph[self.endnode].reversetotalvalue = 1
        waiting = collections.deque()
        waiting.appendleft(self.endnode)
        endwait = False
        while waiting and endwait == False: # iterate while nodes left without incoming edges unacounted for.
            currnode = waiting.popleft()
            if currnode == self.startnode:
                endwait = True
                break
            for i in self.graph[currnode].backedges:
                if self.graph[i.start].outcount < 0:# if an edge is followed more than once, graph not acyclic
                    if self.graph[currnode].bestvalue != float("-inf"): # if loop in possible path location, return error
                       raise Usage("this is not an acyclic graph")
                    else:
                        continue # if  loop not in path section of graph, continue analysis with other edges from loop
                self.graph[i.start].outcount-=1
                endn = False #use to evaluate endnode immediately upon following all edges to it
                if self.graph[i.start] == self.endnode:
                    continue #skip evaluation of edges going into startnode
                else:#add edge to next sum
                    self.graph[i.start].reversetotalvalue = self.graph[i.start].reversetotalvalue+\
                                                            (self.graph[currnode].reversetotalvalue * i.weight)
                    if self.graph[i.start] == self.startnode:
                        endn = True
                if self.graph[i.start].outcount == 0:
                    waiting.appendleft(i.start)
                    if endn:
                        break

    def Get_Norm_Value_Matrix(self, hiddenstatelist, text): # gets graph as matrix of normalized values
        pathlen = len(text)
        vertpos = 0
        horizpos = 0
        currnode = self.startnode
        hiddenlen = len(hiddenstatelist)
        V_Matrix = numpy.zeros((pathlen, hiddenlen))
        while self.graph[currnode].frontedges[0].end != self.endnode:
            for i in self.graph[currnode].frontedges:
                V_Matrix[vertpos][horizpos] = self.graph[i.end].forwardtotalvalue * self.graph[i.end].reversetotalvalue
                horizpos+=1
            V_Matrix[vertpos, :]/= self.totalpathvalue
            vertpos+=1
            horizpos = 0
            currnode = self.graph[currnode].frontedges[0].end
        return V_Matrix





class Markov_Matrices:
    def __init__(self, textstatelist, hiddenstatelist):  # creates matrix, fills with zeros
        self.textstatelist = textstatelist
        self.hiddenstatelist = hiddenstatelist
        self.M_T_Matrix = Markov_Trans_Matrix(hiddenstatelist)
        self.M_E_Matrix = Markov_Emit_Matrix(textstatelist, hiddenstatelist)

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
                                     self.M_E_Matrix.Mmatrix.item(self.M_E_Matrix.corrhiddenDict[l],
                                                                  self.M_E_Matrix.corremitDict[text[i]])
                    else:
                        if tonode >= len(self.hiddenstatelist) * (
                        len(text)) + 1:  # if at end node - should only happen once
                            edgeweight = 1
                        else:
                            edgeweight = self.M_T_Matrix.Mmatrix.item(self.M_T_Matrix.corrDict[k],
                                                                      self.M_T_Matrix.corrDict[l]) * \
                                         self.M_E_Matrix.Mmatrix.item(self.M_E_Matrix.corrhiddenDict[l],
                                                                      self.M_E_Matrix.corremitDict[text[i]])
                    newedge = Edge(fromnode, tonode, edgeweight)
                    DAGraph[fromnode] = DAGraph.get(fromnode, Node())
                    DAGraph[fromnode].frontedges.append(newedge)
                    DAGraph[fromnode].outcount +=1
                    DAGraph[tonode] = DAGraph.get(tonode, Node())
                    DAGraph[tonode].backedges.append(newedge)
                    DAGraph[tonode].incount += 1
                    lcount += 1
                kcount += 1
        toSolve = SumGraph(DAGraph, startnode, endnode)
        toSolve.SolveFG()
        toSolve.SolveRG()
        return toSolve


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
        hidden_state_string = read.readline()
        hidden_state_list = hidden_state_string.strip().split()
        MMatrices = Markov_Matrices(text_state_list,hidden_state_list)#generate empty matrices
        read.readline()
        read.readline()
        for vertpos in range(len(hidden_state_list)):#vertpos is vertical positionn in prevstate-state matrix
            line = read.readline()
            lineprob = line.strip().split()
            horizpos = 0
            for prob in lineprob[1:]:
                MMatrices.M_T_Matrix.Mmatrix[vertpos,horizpos] = prob
                horizpos += 1

        line = read.readline()
        line = read.readline()

        for vertpos in range(len(hidden_state_list)):#read text-state Matrix into program
            line = read.readline()
            lineprob = line.strip().split()
            horizpos = 0
            for prob in lineprob[1:]:
                MMatrices.M_E_Matrix.Mmatrix[vertpos, horizpos] = prob
                horizpos += 1


        S_Graph = MMatrices.Solve_Markov(text)
        V_Matrix = S_Graph.Get_Norm_Value_Matrix(hidden_state_list, text)


        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i)+ "\t")
        writef.write( "\n")
        for i in range(len(text)):
            row = V_Matrix[i]
            for j in numpy.nditer(row):
                writef.write('{0:<4.4g}'.format(float(j))+"\t")
            writef.write( "\n")



    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()

