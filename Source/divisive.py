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
Z = [[0 for x in range(4)] for y in range(309)] 
import numpy as np
score = np.loadtxt("out.txt")
score = np.delete(score, (311), axis=0)
score = np.delete(score, (311), axis=1)
print score
import pandas as pd
num_clusters = 0
dissimilarity_matrix = pd.DataFrame(score,index=labels, columns=labels)
print dissimilarity_matrix

def interCluster(ele, element_list):
    max_diameter = -np.inf
    sumd = 0
    for i in element_list:
        sumd += dissimilarity_matrix[ele][i]   
        if( dissimilarity_matrix[ele][i]  > max_diameter):
            max_diameter = dissimilarity_matrix[ele][i]
    if(len(element_list)>1):
        avg = sumd/(len(element_list)-1)
    else: 
        avg = 0
    return avg
def distance(cluster1,cluster2):
    sumd = 0
    for i in cluster1:
        for j in cluster2:
            sumd += dissimilarity_matrix[i][j]
    avg = sumd/(len(cluster1)*len(cluster2))
    return avg
            
def intraCluster(ele, main_list, splinter_list):
    if len(splinter_list) == 0:
        return 0
    sumd = 0
    for j in splinter_list:
        sumd = sumd + dissimilarity_matrix[ele][j]
    avg = sumd/(len(splinter_list))
    return avg
    
    
def splinter(main_list, splinters):
    max_val = -np.inf
    max_index = None
    for ele in main_list:
        x = interCluster(ele, main_list)
        y = intraCluster(ele, main_list, splinters)
        diff= x -y
        if diff > max_val:
            max_val = diff
            max_index = ele
    if(max_val>0):
        return  (max_index, 1)
    else:
        return (-1, -1)
    
def split(element_list):
    main_list = element_list
    splinters = []    
    (maxIndex,flag) = splinter(main_list, splinters)
    while(flag > 0):
        main_list.remove(maxIndex)
        splinters.append(maxIndex)
        (maxIndex,flag) = splinter(element_list, splinters)
    
    return (main_list, splinters)

def max_diameter(cluster_list):
    cluster_index = None
    cluster_value = -np.inf
    index = 0
    for element_list in cluster_list:
        for i in element_list:
            for j in element_list:
                if dissimilarity_matrix[i][j]  > cluster_value:
                    cluster_value = dissimilarity_matrix[i][j]
                    cluster_index = index
        
        index +=1
    
    if(cluster_value <= 0):
        return -1
    
    return cluster_index
    
def divisive(clusters):
    count = 0
    allClusters = {}
    clusterNumbers = {}
    invCluster = {}
    allClusters[count] = clusters
    clusterNumbers[count] = len(labels)
    count+=1
    
    index = 0
    while(index!=-1):  
        print(clusters)
        print
        (left, right) = split(clusters[index])
        del clusters[index]
        clusters.append(right)
        clusters.append(left)
        for every_list in clusters:
            temp = []
            flag = 1
            for every_element in every_list:
                temp.append(every_element)
            for key in allClusters:
                if allClusters[key]==temp:
                    flag = 0
            if flag==1:
                allClusters[count] = temp
                count+=1
        index = max_diameter(clusters)
    print(clusters)
    total = count
    for key in allClusters:
        invCluster[total-key-1] = allClusters[key]
        clusterNumbers[total-key-1] = len(allClusters[key])
    
    
    print 
    print invCluster
    print 
    print clusterNumbers
    it = 0
    key = 0
    while(key<=616):
        Z[it][0] = key
        Z[it][1] = key+1
        Z[it][2] = distance(invCluster[key],invCluster[key+1])
        Z[it][3] = clusterNumbers[key]+clusterNumbers[key+1]
        it+=1
        key+=2
    print Z
    


clusters = ([labels])
divisive(clusters)



from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
plt.figure()
import numpy 
n= numpy.array(Z)
n = n.astype(numpy.double)
print n
fig=plt.figure(figsize=(18, 12), dpi= 80, facecolor='w', edgecolor='k')

dn = hierarchy.dendrogram(n,p=25,truncate_mode='lastp')
plt.show()

fig=plt.figure(figsize=(18, 12), dpi= 80, facecolor='w', edgecolor='k')
dn = hierarchy.dendrogram(n)
plt.show()
