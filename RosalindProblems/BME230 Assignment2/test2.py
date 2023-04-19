class SeedCluster:
    """ Represents a set of seeds between two strings.
    """

    def __init__(self, seeds):
        """ Seeds is a list of pairs [ (x_1, y_1), (x_2, y_2), ..., ], each is an instance of a seed
        (see static cluster seeds method below: static methods: https://realpython.com/blog/python/instance-class-and-static-methods-demystified/)
        """
        seeds = list(seeds)
        seeds.sort()
        self.seeds = seeds
        # Gather the minimum and maximum x and y coordinates
        self.minX = seeds[0][0]
        self.maxX = seeds[-1][0]
        ys = map(lambda (x, y): y, seeds)
        self.minY = min(ys)
        self.maxY = max(ys)

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.

        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of shared k-mer in both strings, termed a *seed*, such that the k-mer
        occurrence starts at position x in the first string and starts at position y_i in the second string.

        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.

        Consider a graph in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).

        The return value is a Python set of components, each component is a SeedCluster object.

        (QUESTION 1): The clustering of seeds is very simplistic. Can you suggest alternative strategies by
        which the seeds could be clustered, and what the potential benefits such alternative strategies could
        have? Consider the types of information you could use.
        """

        # Code to complete - you are free to define other functions as you like
        clusters = set()
        for xcoord,b in seeds:
            for ycoord in b:
                sub_clusters = set(a for a in clusters if a.can_append(xcoord,ycoord,l))
                if len(sub_clusters)>0:
                    elem = sub_clusters.pop()
                    clusters.difference_update(sub_clusters)
                    for a in sub_clusters:
                        elem.add_cluster(a)
                    elem.append_element(xcoord,ycoord)
                else:
                    clusters.add(SeedCluster([(xcoord,ycoord)]))

    def append_element(self,xcoord,ycoord):
        self.seeds.append((xcoord,ycoord))
        self.maxX = xcoord
        self.minY = min(self.minY,ycoord)
        self.maxY = max(self.maxY,ycoord)

    def add_cluster(self,other):
        i = 0
        j = 0
        selflen = len(self.seeds)
        otherlen = len(other.seeds)
        new_seeds = list()
        while i < selflen and j<otherlen:
            if self.seeds[i]<other.seeds[j]:
                new_seeds.append(self.seeds[i])
                i+=1
            else:
                new_seeds.append(other.seeds[j])
                j+=1
        new_seeds+=self.seeds[i:]
        new_seeds+=other.seeds[j:]
        self.seeds = new_seeds
        self.minX = min(self.minX,other.minX)
        self.maxX = max(self.maxX,other.maxX)
        self.minY = min(self.minY,other.minY)
        self.maxY = max(self.maxY,other.maxY)




    def can_append(self, xcoord,ycoord,l):
        for i in reversed(self.seeds):
            if i[0] < xcoord -l:
                return False
            else:
                if i[1]>=ycoord - l and i[1] <= ycoord + l:
                    return True
