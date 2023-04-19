# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys


class CommandLine():
    def __init__(self, inOpts=None):

        import argparse

        self.parser = argparse.ArgumentParser(
            description='Program estimates emission and transmission matrices most likely to output a path or string',
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
class Markov_Matrices:
    def __init__(self, textstatelist, hiddenstatelist):  # creates matrix, fills with zeros
        self.textstatelist = textstatelist
        self.hiddenstatelist = hiddenstatelist
        self.M_T_Matrix = Markov_Trans_Matrix(hiddenstatelist)
        self.M_E_Matrix = Markov_Emit_Matrix(textstatelist, hiddenstatelist)

    def most_likely(self, text, path):
        l = len(text)
        for i in range(l):
            self.M_E_Matrix.Mmatrix[self.M_E_Matrix.corrhiddenDict[path[i]]][self.M_E_Matrix.corremitDict[text[i]]] += 1
            if i != 0:
                self.M_T_Matrix.Mmatrix[self.M_T_Matrix.corrDict[path[i-1]]][self.M_T_Matrix.corrDict[path[i]]] +=1

        for row in self.M_T_Matrix.Mmatrix: # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum!=0:
                row/=row_sum
            else:
                l = len(self.hiddenstatelist)
                for j in range(l):
                    row[j] = 1/l

        for row in self.M_E_Matrix.Mmatrix: # divide each number in emission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum!=0:
                row/=row_sum
            else:
                l = len(self.textstatelist)
                for j in range(l):
                    row[j] = 1 / l

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

########################################################################
# Main
# Here is the main program
#
#
########################################################################
import math

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
        text = str(read.readline().strip())
        read.readline()
        text_state_list = read.readline().strip().split()
        read.readline()
        path = str(read.readline().strip())
        read.readline()
        hidden_state_list = read.readline().strip().split()
        MMatrices = Markov_Matrices(text_state_list, hidden_state_list)
        MMatrices.most_likely(text, path)

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

