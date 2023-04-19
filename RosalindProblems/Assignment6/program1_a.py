# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys


class CommandLine():
    def __init__(self, inOpts=None):

        import argparse

        self.parser = argparse.ArgumentParser(
            description='Program finds probability of a hidden path',
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
class Markov_Matrix:
    def __init__(self, statelist):#creates matrix, fills with zeros
        self.statelist = statelist
        self.corrDict = dict() #dictionary of state to matrix position
        count = 0
        for i in statelist:
            self.corrDict[i] = count
            count+=1
        self.Mmatrix = numpy.zeros((count,count))
    def pathprob(self, path):
        totalprob = 1
        pos = 0
        for i in path:
            if pos == 0:
                totalprob = totalprob * 1/len(self.statelist)
            else:
                totalprob = totalprob * self.Mmatrix.item(self.corrDict[i],self.corrDict[path[pos-1]])
            pos+=1
        return totalprob



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
        path = read.readline()[:-1]
        read.readline()
        statestring = read.readline()
        statelist = statestring.strip().split()
        statelist = [x for x in statelist if x != ""]
        Mmatrix1 = Markov_Matrix(statelist)
        read.readline()
        read.readline()
        vertpos = 0
        for line in read:
            if line.strip() is not "":
                lineprob = line.strip().split()
                horizpos = 0
                for prob in lineprob[1:]:
                    Mmatrix1.Mmatrix[horizpos, vertpos] = prob
                    horizpos += 1
                vertpos += 1
            else:
                break
        writef.write(str(Mmatrix1.pathprob(path)))
    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()

