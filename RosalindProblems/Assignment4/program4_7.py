# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys
class CommandLine():
    def __init__(self, inOpts=None):

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program reads adjacency list of numbers and orders them into an Eulrian cicle using a DeBroijne graph',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-i', '--input', type=str, default=sys.stdin.readline().strip(), action='store',
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
class CalcPath:
    def __init__(self, dBGraph):
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
                raise Usage("A path is impossible, difference between number of exits and entrances in node >1 ")
            elif (osum == ival + 1):
                start = k
                startcount += 1
            elif (ival == osum + 1):
                end = k
                endcount += 1
        if startcount != endcount or startcount > 1:
            raise Usage("A path is impossible, more than one node each must be a start or an end")
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
        #  dBGraph stored in form dictionary{istring: [[ostring1, count][ostring2, count]], istringcount},
        #  where istringcount is # of times istring is referenced elsewhere
        #main program begins
        readf = myCommandLine.args.input
        read = open(readf, "r")
        if myCommandLine.args.output is None:
            writef = sys.stdout
        else:
            writef = open(myCommandLine.args.output, "w")
        dBGraph = dict()
        next(read)
        for line in read:
            firstn = line[0:-2]
            secondn = line[1:-1]
            #increase count of incoming edges for node secondn
            if secondn in dBGraph:
                dBGraph[secondn][1] += 1
            else:
                dBGraph[secondn] = [list(), 1]
            if firstn in dBGraph:#check if value of full edge in list, if so, inc. count,else add to dict
                inList = False
                for i in dBGraph[firstn][0]:
                    if secondn == i[0]:
                        inList = True
                        i[1] += 1
                if not inList:
                    dBGraph[firstn][0].append([secondn,1])
            else:
                dBGraph[firstn] = [[[secondn,1]], dBGraph.get(firstn, [list(),0])[1]]#form node in graph

        # end graph input

        cp = CalcPath(dBGraph)
        path = cp.findPath()
        writef.write(path[0])
        for i in path[1:]:
            writef.write(i[-1])
    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
    main()

