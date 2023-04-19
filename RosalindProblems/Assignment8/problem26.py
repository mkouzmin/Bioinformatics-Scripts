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

class Resp_Matrices:
    def __init__(self, hiddenstatelist, text):
        self.hiddenstatelist = hiddenstatelist
        self.text = text
        self.pathlen = len(text)
        self.hiddenlen = len(hiddenstatelist)
        edgecount = 0
        self.edgedict = dict()
        self.edgelist = list()
        for i in hiddenstatelist:
            for j in hiddenstatelist:
                self.edgelist.append(str(i)+str(j))
                self.edgedict[str(i)+str(j)] = edgecount
                edgecount+=1
        self.R_Matrix = numpy.zeros((self.pathlen,self.hiddenlen))
        R_E_Matrix_len = len(self.text)-1
        self.R_E_Matrix = numpy.zeros((R_E_Matrix_len,edgecount))

    def Fill_Resp_Node_Matrix(self, Graph): # gets graph as matrix of normalized responsibility values
        vertpos = 0
        horizpos = 0
        currnode = Graph.startnode
        while Graph.graph[currnode].frontedges[0].end != Graph.endnode:
            for i in Graph.graph[currnode].frontedges:
                self.R_Matrix[vertpos][horizpos] = Graph.graph[i.end].forwardtotalvalue * Graph.graph[i.end].reversetotalvalue
                horizpos+=1
            self.R_Matrix[vertpos, :]/= Graph.totalpathvalue
            vertpos+=1
            horizpos = 0
            currnode = Graph.graph[currnode].frontedges[0].end


    def Fill_Resp_Edge_Matrix(self, Graph):
        currnode = Graph.startnode
        M_vert_count = 0
        while Graph.graph[Graph.graph[currnode].frontedges[0].end].frontedges[0].end != Graph.endnode:
            M_horiz_count = 0
            for i in Graph.graph[currnode].frontedges:
                for j in Graph.graph[i.end].frontedges:
                    self.R_E_Matrix[M_vert_count,M_horiz_count] = Graph.graph[j.start].forwardtotalvalue *\
                                                             Graph.graph[j.end].reversetotalvalue * j.weight
                    M_horiz_count+=1
            self.R_E_Matrix[M_vert_count, :]/= Graph.totalpathvalue
            M_vert_count+=1
            currnode = Graph.graph[currnode].frontedges[0].end





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

    def read_Resp_Matrix(self,R_Matrices,text):
        corrtextDict = dict()
        textcount = 0
        for i in self.textstatelist:
            corrtextDict[i] = textcount
            textcount += 1
        corrhiddenDict = dict()
        hiddencount = 0
        for i in self.hiddenstatelist:
            corrhiddenDict[i] = hiddencount
            hiddencount += 1
        textpos = 0
        for i in text:
            text_state_num = corrtextDict[i]
            state_pos = 0
            for j in R_Matrices.R_Matrix[textpos]:
                self.M_E_Matrix.Mmatrix[state_pos][text_state_num]+=j
                state_pos+=1
            textpos+=1
        for row in self.M_E_Matrix.Mmatrix: # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum!=0:
                row/=row_sum


        rownum = 0
        for row in R_Matrices.R_E_Matrix:
            colnum = 0
            for item in row:
                colid = R_Matrices.edgelist[colnum]
                i = corrhiddenDict[colid[0]]
                j = corrhiddenDict[colid[1]]
                self.M_T_Matrix.Mmatrix[i][j] += item
                colnum+=1
            rownum+=1

        for row in self.M_T_Matrix.Mmatrix:  # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum != 0:
                row /= row_sum






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
        reps = read.readline()[:-1]
        read.readline()
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

        for i in range(int(reps)):
            S_Graph = MMatrices.Solve_Markov(text)
            R_Matrices = Resp_Matrices(hidden_state_list,text)
            R_Matrices.Fill_Resp_Node_Matrix(S_Graph)
            R_Matrices.Fill_Resp_Edge_Matrix(S_Graph)
            MMatrices = Markov_Matrices(text_state_list, hidden_state_list)
            MMatrices.read_Resp_Matrix(R_Matrices,text)


        # print transition matrix
        writef.write('{0:<4s}'.format(""))
        for i in hidden_state_list:
            writef.write("\t" + '{0:<4s}'.format(i))
        writef.write( "\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = MMatrices.M_T_Matrix.Mmatrix[MMatrices.M_T_Matrix.corrDict[i]]
            for j in numpy.nditer(row):
                writef.write("\t"+'{0:<4.3g}'.format(math.trunc((float(j)*1000)+0.5)/1000))
            writef.write( "\n")

        writef.write("--------\n")

        # print emission matrix
        writef.write('{0:<4s}'.format(""))
        for i in text_state_list:
            writef.write("\t" + '{0:<4s}'.format(i))
        writef.write( "\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = MMatrices.M_E_Matrix.Mmatrix[MMatrices.M_E_Matrix.corrhiddenDict[i]]
            for j in numpy.nditer(row):
                writef.write("\t"+'{0:<4.3g}'.format(math.trunc((float(j)*1000)+0.5)/1000))
            writef.write( "\n")




    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()

