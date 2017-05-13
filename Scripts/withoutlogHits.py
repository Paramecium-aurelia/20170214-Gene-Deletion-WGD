from __future__ import division
import re, sys, math, operator, copy,collections
import numpy as np 

def ProcessCLI(args):
    #scaf_pb=readgff(args[1])
    bestHits= readpsl(args[1])
    for k in bestHits:
        print ' \t'.join([str(i) for i in list(k)])

#Read in and update the psl output and find hits
def readpsl(pslout):
    bhits=[]
    hits=collections.OrderedDict()
    phits=collections.OrderedDict()
    with open(pslout,'r') as file:
        for line in file:
            line=line.rstrip().split()
            ## Here are the filtering conditions: pident has to be than 50% and the minimum common subsequence overlap of  the querry and subject >=40%
            ## These criteria could be tuned and changed as desired
            if float(line[4])>40.0:# and max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])/min(float(line[5]),float(line[6]))>=0.20:
                #print ' '.join(line)
                if line[0]!=line[2]:
                    print ' '.join(line)
                    if line[1] not in hits.keys():
                        hits[line[1]]=[line[3]]
                    else:
                        hits[line[1]].append(line[3])
                    if line[0] not in phits.keys():
                        phits[line[0]]={line[2]:line[-2]}
                    else:
                        phits[line[0]][line[2]]=line[-2]
                        
#    print hits['scaffold_0001']
    return 

    for key in hits:
        #hits[key]=dict([a,x] for a, x in collections.Counter(hits[key]).iteritems())
        hits[key]=dict([a,-math.log(float(x)/sum(collections.Counter(hits[key]).values()))] for a, x in collections.Counter(hits[key]).iteritems\
())

#    print hits['scaffold_0269']
    for k1 in hits: 
        for k2 in hits[k1]:
            score1=0
            score2=0
            try:
                score1=hits[k1][k2]
                
            except:
                score1=10
            #print k1,k2,score1
            try:
                score2=hits[k2][k1]
            except:
                score2=10
            score=score1+score2
            
            bhits.append((k1,k2,score))
    for p in phits.keys():
        print p
        #print phits[p]
    return 
    return sorted(bhits, key=lambda tup: (tup[2]) )#bhits

#Read he gff3 files and create a dictionary: scaffolds: list of proteins                                   
def readgff(f):
    filename = open(f, "r")
    scaffolds={}
    pbi={}
    for line in filename:
        line=line.strip().split() 
        if line[0][:8]=='scaffold' and line[2]=='mRNA':
            pbi[re.split(r'[=;\s]',line[8])[1]]=line[0]
    return pbi 

if __name__ == '__main__':
     ProcessCLI(sys.argv)
