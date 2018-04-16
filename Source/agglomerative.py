import itertools
from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

fasta = fasta_iter('dna.fa')
labels = []
for ff in fasta:
    headerStr, seq = ff
    labels.append(headerStr)
Z = [[0 for x in range(4)] for y in range(311)] 
import numpy
score = numpy.loadtxt("out.txt")
score = numpy.delete(score, (311), axis=0)
score = numpy.delete(score, (311), axis=1)
print score
class Cluster:
    def __init__(self):
        pass
    def __repr__(self):
        return '(%s,%s)' % (self.left, self.right)
    def add(self, clusters, grid, lefti, righti):
        self.left = clusters[lefti]
        self.right = clusters[righti]
        # merge columns grid[row][righti] and row grid[righti] into corresponding lefti
        for r in grid:
            r[lefti] = min(r[lefti], r.pop(righti))
        grid[lefti] = map(min, zip(grid[lefti], grid.pop(righti)))
        clusters.pop(righti)
        return (clusters, grid)

def agglomerate(labels, grid):
    clusters = labels
    allClusters = {}
    clusterNumbers = {}
    count = 0
    for c in clusters:
        allClusters[c] = count
        clusterNumbers[count] = 1
        count +=1
    it = 0
    while len(clusters) > 1:
        # find 2 closest clusters
        #print clusters
        distances = [(1, 0, grid[1][0])]
        for i,row in enumerate(grid[2:]):
            distances += [(i+2, j, c) for j,c in enumerate(row[:i+2])]
        j,i,_ = min(distances, key=lambda x:x[2])
        # merge i<-j
        c = Cluster()
        #print clusters[i]
        #print clusters[j]
        c1 = allClusters[clusters[i]]
        c2 = allClusters[clusters[j]]
        print c1 
        print c2
        Z[it][0] = c1
        Z[it][1] = c2
        Z[it][2] = _
        clusters, grid = c.add(clusters, grid, i, j)
        clusterNumbers[count] = clusterNumbers[c1]+clusterNumbers[c2]
        clusters[i] = c
        allClusters[clusters[i]] = count
        Z[it][3] = clusterNumbers[count]
        count+=1
        it+=1
    print Z
    print allClusters
    return clusters.pop()
import sys
sys.setrecursionlimit(100000)
import numpy
#random_matrix = numpy.random.randint(0,1000,(22,22))    
s = score.tolist()
agglomerate(labels,s)
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
plt.figure()

