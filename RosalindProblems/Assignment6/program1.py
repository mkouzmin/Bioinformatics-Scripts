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
                totalprob = totalprob * 1/len(statelist)
            else:
                totalprob = totalprob * self.Mmatrix.item(self.corrDict[i],self.corrDict[path[pos-1]])
            pos+=1
        return totalprob


readf = "input.txt"
read = open(readf,"r")
writef = open("output.txt","w")
path = read.readline()[:-1]
read.readline()
statestring = read.readline()
statelist = statestring.strip().split(" ")
statelist = [x for x in statelist if x!= ""]
Mmatrix1 = Markov_Matrix(statelist)
read.readline()
read.readline()
vertpos = 0
for line in read:
    lineprob = line.strip().split()
    horizpos = 0
    for prob in lineprob[1:]:
        Mmatrix1.Mmatrix[horizpos,vertpos] = prob
        horizpos+=1
    vertpos+=1
writef.write (str(Mmatrix1.pathprob(path)))