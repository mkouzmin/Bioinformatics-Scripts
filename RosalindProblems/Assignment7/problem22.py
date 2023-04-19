# scaffold taken from Canvas Website
########################################################################
# CommandLine
########################################################################
import sys


class CommandLine():
    def __init__(self, inOpts=None):

        import argparse

        self.parser = argparse.ArgumentParser(
            description='Program parses a string into its most likely hidden path',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-i', '--input', type=str, default="input.txt", action='store',
                                 help='Input File')
        self.parser.add_argument('-o', '--output', type=str, default="output.txt", action='store', help='Output File')
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

class Alignment:
    def __init__(self, alignmatrix, inprob):#takes matrix of allignments
        self.alignmatrix = alignmatrix
        self.threshold = float(inprob)
        alidim = self.alignmatrix.shape
        self.alilen = alidim[1]  # length of individual allignment
        self.alicount = alidim[0]  # number of given allignments
        self.insertlist = list()#create list of columns evaluated as insertions
        colnum = 0
        for col in self.alignmatrix.T:
            inscount = 0
            for i in numpy.nditer(col):
                if i == "-":
                    inscount += 1
            if (inscount / self.alicount) >= self.threshold:
                self.insertlist.append(colnum)
            colnum += 1
        self.vallen = self.alilen - len(self.insertlist)#get number of valued columns in alignment

    def fill_Matrix(self, MMatrices):
        for align in self.alignmatrix:
            colnum = 0#keep track of position in align
            colvalnum = 0#keep track of position in align without insert lines
            prevstate = "S"
            for col in numpy.nditer(align): # for each character in Alignment line
                if colnum in self.insertlist:
                    colnum+=1
                    if str(col) == "-":
                        continue
                    else:
                        currstate = "I" + str(colvalnum)
                else:
                    colnum += 1
                    colvalnum+=1
                    if col == "-":
                        currstate = "D" + str(colvalnum)
                    else:
                        currstate = "M" + str(colvalnum)
                MMatrices.M_T_Matrix.Mmatrix[MMatrices.M_T_Matrix.corrDict[prevstate],
                                             MMatrices.M_T_Matrix.corrDict[currstate]] += 1# add 1 to right space,
                                                                                            #  divide by rowsum later
                if str(col) != "-":
                    MMatrices.M_E_Matrix.Mmatrix[MMatrices.M_E_Matrix.corrhiddenDict[currstate],
                                                 MMatrices.M_E_Matrix.corremitDict[str(col)]] += 1

                prevstate = currstate
            MMatrices.M_T_Matrix.Mmatrix[MMatrices.M_T_Matrix.corrDict[prevstate],
                                         MMatrices.M_T_Matrix.corrDict["E"]] += 1

        for row in MMatrices.M_T_Matrix.Mmatrix: # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum!=0:
                row/=row_sum

        for row in MMatrices.M_E_Matrix.Mmatrix:  # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum!= 0:
                row/=row_sum

    def pseudo_norm_Matrix(self, MMatrices, pseudocount):

        for i,j in MMatrices.M_T_Matrix.corrDict.items():# pseudocount transmission matrix
            if i == "S" or i == "I0":

                MMatrices.M_T_Matrix.Mmatrix[int(j)] *= (1 - (3 * float(pseudocount)))
                for k in [MMatrices.M_T_Matrix.corrDict["I0"],
                          MMatrices.M_T_Matrix.corrDict["M1"],
                          MMatrices.M_T_Matrix.corrDict["D1"]]:
                    MMatrices.M_T_Matrix.Mmatrix[int(j)][int(k)] += float(pseudocount)

            elif i != "E":
                start_count = int(i[1:])
                if int(i[1:]) != self.vallen: # if 2nd value of key not equal to length of alignments
                    MMatrices.M_T_Matrix.Mmatrix[int(j)] *= 1-(3*float(pseudocount))
                    for k in [MMatrices.M_T_Matrix.corrDict[str("I"+ str(start_count))],
                              MMatrices.M_T_Matrix.corrDict[str("M"+ str(start_count+1))],
                              MMatrices.M_T_Matrix.corrDict[str("D"+str(start_count+1))]]:
                        MMatrices.M_T_Matrix.Mmatrix[int(j)][int(k)] += float(pseudocount)

                else:# if at last alignment character
                    MMatrices.M_T_Matrix.Mmatrix[int(j)] *= 1 - (2 * float(pseudocount))
                    for k in [MMatrices.M_T_Matrix.corrDict[str("I" + str(start_count))],
                              MMatrices.M_T_Matrix.corrDict["E"]]:
                        MMatrices.M_T_Matrix.Mmatrix[int(j)][int(k)] += float(pseudocount)

        for row in MMatrices.M_T_Matrix.Mmatrix: # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum!=0:
                row/=row_sum

        for i, j in MMatrices.M_E_Matrix.corrhiddenDict.items():
            if i[0] == "I" or i[0] == "M":
                MMatrices.M_E_Matrix.Mmatrix[int(j)] *= 1 - (float(pseudocount) * len(MMatrices.M_E_Matrix.emitstatelist))
                for k in MMatrices.M_E_Matrix.corremitDict.values():
                    MMatrices.M_E_Matrix.Mmatrix[int(j)][int(k)] += float(pseudocount)

        for row in MMatrices.M_E_Matrix.Mmatrix:  # divide each number in transmission matrix by sum of it's row
            row_sum = sum(row)
            if row_sum != 0:
                row /= row_sum

class Viterbi_G:
    def __init__(self, hidden_state_list, path):
        self.hidden_state_list = hidden_state_list
        self.path = str(path)
        self.path_len = len(path)
        self.corr_T_Dict = dict()  # dictionary of state to matrix position
        count = 0
        for i in hidden_state_list:
            self.corr_T_Dict[i] = count
            count += 1
        self.Vmatrix = numpy.empty((count, self.path_len + 2))#create -inf matrix for viterbi algorithm, path[0] corresponding to column 1
        self.Vmatrix.fill(0-numpy.inf)
        self.Pmatrix = numpy.zeros((count, self.path_len + 2))#use to keep track of last vertical position in matrix
        #self.Dmatrix = numpy.zeros((count, self.path_len + 2))#keeps track of best product to reach each location

    def get_T_Value(self, MMatrices, prevvalue, currvalue): #inputs are values in dictionary
        return MMatrices.M_T_Matrix.Mmatrix[prevvalue][currvalue]

    def get_E_Value(self, MMatrices, hiddenstatenum, letter):# hiddenstate is value in hidden state dictionary
        return MMatrices.M_E_Matrix.Mmatrix[hiddenstatenum][MMatrices.M_E_Matrix.corremitDict[letter]]


    def solve_Matrix(self,MMatrices):
        for i in range(self.path_len + 2):
            for k,v in self.corr_T_Dict.items():
                if k[0] == "S":
                    if i == 0:
                        self.Vmatrix[v][i] = numpy.float64(0)
                        self.Pmatrix[v][i] = -1
                elif k[0] == "I" and i > 0 and i < self.path_len+1:
                    if i == 1:#if in left column
                        self.Vmatrix[v][i] = self.Vmatrix[v - 1][i - 1] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, v - 1, v) \
                                             * self.get_E_Value(MMatrices, v, self.path[i - 1])))
                        self.Pmatrix[v][i] = v-1
                        #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, v-1, v)
                                                       #* self.get_E_Value(MMatrices, v, self.path[i-1]))
                    elif k[1:] == "0":#if in top row
                        self.Vmatrix[v][i] = self.Vmatrix[v][i-1] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, v, v)\
                                             * self.get_E_Value(MMatrices,v,self.path[i-1])))
                        self.Pmatrix[v][i] = v
                        #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, v, v)
                                                           #* self.get_E_Value(MMatrices, v, self.path[i-1]))
                    else:
                        for x in [v,v-1,v-2]:
                            num = self.Vmatrix[x][i-1] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, x, v)\
                                                 * self.get_E_Value(MMatrices,v,self.path[i-1])))
                            if self.Vmatrix[v][i] < num:
                                self.Vmatrix[v][i] = num
                                self.Pmatrix[v][i] = x
                                #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, x, v)\
                                                 #* self.get_E_Value(MMatrices,v,self.path[i-1]))
                elif k[0] == "M" and i>0 and i< self.path_len +1:
                    if i == 1:#if in left column
                        self.Vmatrix[v][i] = self.Vmatrix[v - 2][i - 1] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, v - 2, v) \
                                             * self.get_E_Value(MMatrices, v, self.path[i - 1])))
                        self.Pmatrix[v][i] = v-2
                        #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, v - 2, v) \
                                             #* self.get_E_Value(MMatrices, v, self.path[i - 1]))
                    elif k[1:] == "1":# if in top row
                        self.Vmatrix[v][i] = self.Vmatrix[v-1][i-1] + numpy.float64(numpy.log(self.get_T_Value(MMatrices, v-1, v)\
                                             * self.get_E_Value(MMatrices,v,self.path[i-1])))
                        self.Pmatrix[v][i] = v-1
                        #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, v-1, v)\
                                             #* self.get_E_Value(MMatrices,v,self.path[i-1]))
                    else:# for general case
                        for x in [v-1,v-2,v-3]:
                            num = self.Vmatrix[x][i-1] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, x, v)\
                                                 * self.get_E_Value(MMatrices,v,self.path[i-1])))
                            if self.Vmatrix[v][i] < num:
                                self.Vmatrix[v][i] = num
                                self.Pmatrix[v][i] = x
                                #self.Dmatrix[v][i] = 0-numpy.log(self.get_T_Value(MMatrices, x, v)\
                                                 #* self.get_E_Value(MMatrices,v,self.path[i-1]))
                elif k[0] =="D" and i< self.path_len +1:
                    if i == 0:# if in left column
                        self.Vmatrix[v][i] = self.Vmatrix[v - 3][i] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, v - 3, v)))
                        self.Pmatrix[v][i] = v-3
                        #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, v - 3, v))
                    elif k[1:] == "1":# if in top row
                        self.Vmatrix[v][i] = self.Vmatrix[v - 2][i] + numpy.float64(numpy.log( self.get_T_Value(MMatrices, v - 2, v)))
                        self.Pmatrix[v][i] = v-2
                        #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, v - 2, v))
                    else: # general case
                        for x in [v-2,v-3,v-4]:
                            num = self.Vmatrix[x][i] + numpy.float64(numpy.log(self.get_T_Value(MMatrices, x, v)))
                            if self.Vmatrix[v][i] < num:
                                self.Vmatrix[v][i] = num
                                self.Pmatrix[v][i] = x
                                #self.Dmatrix[v][i] = 0 - numpy.log(self.get_T_Value(MMatrices, x, v))
                elif k[0] == "E" and i == self.path_len + 1:
                    for x in [v - 1, v - 2, v - 3]:
                        num = self.Vmatrix[x][i-1]
                        if self.Vmatrix[v][i] < num:
                            self.Vmatrix[v][i] = num
                            self.Pmatrix[v][i] = x


        i = len(self.hidden_state_list)-1
        j = self.path_len + 1
        reverse_path = list()
        while self.Pmatrix[i][j] != -1:
            reverse_path.append(self.hidden_state_list[i])
            i2 = i
            i = int(self.Pmatrix[i][j])
            if self.hidden_state_list[i2][0] != "D":
                j-=1
        return reversed(reverse_path[1:])





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
        path = read.readline().strip().split()[0]
        read.readline()
        thresh_pseudo_line = read.readline().strip().split()
        inprob = thresh_pseudo_line[0]#get threshold
        pseudocount = thresh_pseudo_line[1]
        read.readline()
        dictstring = read.readline()
        dictlist = dictstring.strip().split()
        read.readline()
        alilist = list()
        for line in read:
            if line.strip() is not "":
                avalue = list()
                for l in line.strip():
                    avalue.append(l)
                alilist.append(avalue)
            else:
                break
        alignmatrix = numpy.matrix(alilist) # get allignment matrix, to put into allignment class
        alignment = Alignment(alignmatrix,inprob)
        hidden_state_list = list(["S", "I0"])  # create list of hidden states
        for i in range(alignment.vallen):
            for c in ["M", "D", "I"]:
                hidden_state_list.append(c + str(i + 1))
        hidden_state_list.append("E")
        MMatrices = Markov_Matrices(dictlist, hidden_state_list)
        alignment.fill_Matrix(MMatrices)
        alignment.pseudo_norm_Matrix(MMatrices, pseudocount)
        VGraph = Viterbi_G(hidden_state_list, path)
        path_iter = VGraph.solve_Matrix(MMatrices)
        for i in path_iter:
            writef.write(i + " ")
        writef.write("\n")
    except Usage as err:
        print(err.msg)
"""
        # print transition matrix
        
        writef.write('{0:<4s}'.format(""))
        for i in hidden_state_list:
            writef.write("\t" + '{0:<4s}'.format(i))
        writef.write( "\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = MMatrices.M_T_Matrix.Mmatrix[MMatrices.M_T_Matrix.corrDict[i]]
            for j in numpy.nditer(row):
                writef.write("\t"+'{0:<4.3g}'.format(float(j)))
            writef.write( "\n")

        writef.write("--------\n")

        # print emission matrix
        writef.write('{0:<4s}'.format(""))
        for i in dictlist:
            writef.write("\t" + '{0:<4s}'.format(i))
        writef.write( "\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = MMatrices.M_E_Matrix.Mmatrix[MMatrices.M_E_Matrix.corrhiddenDict[i]]
            for j in numpy.nditer(row):
                writef.write("\t"+'{0:<4.3g}'.format(float(j)))
            writef.write( "\n")



        writef.write('{0:<4s}'.format(""))
        for i in range(len(path)+2):
            writef.write("\t" + '{0:<4g}'.format(i))
        writef.write("\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = VGraph.Vmatrix[VGraph.corr_T_Dict[i]]
            for j in numpy.nditer(row):
                writef.write("\t" + '{0:<4.2g}'.format(float(j)))
            writef.write("\n")

        writef.write("--------------------\n")

        writef.write('{0:<4s}'.format(""))
        for i in range(len(path)+2):
            writef.write("\t" + '{0:<4g}'.format(i))
        writef.write("\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = VGraph.Pmatrix[VGraph.corr_T_Dict[i]]
            for j in numpy.nditer(row):
                writef.write("\t" + '{0:<4g}'.format(float(j)))
            writef.write("\n")

        writef.write("---------------------\n")

        writef.write('{0:<4s}'.format(""))
        for i in range(len(path)+2):
            writef.write("\t" + '{0:<4g}'.format(i))
        writef.write("\n")
        for i in hidden_state_list:
            writef.write('{0:<4s}'.format(i))
            row = VGraph.Dmatrix[VGraph.corr_T_Dict[i]]
            for j in numpy.nditer(row):
                writef.write("\t" + '{0:<4.2g}'.format(float(j)))
            writef.write("\n")
"""



if __name__ == "__main__":
    main()

