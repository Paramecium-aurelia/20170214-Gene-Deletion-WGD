import re, sys, math, operator, copy

def ProcessCLI(args):
    scaf_pb=readgff(args[1])
    bestHits,onewayhits= readpsl(args[2],scaf_pb)
    print 'Best Hits : reciprocal hits'
    print '------------'
    for hit in bestHits:
        print ' '.join(hit)
    print '\n***********************\n'
    print 'One way hits'
    print '------------'
    for hit1 in onewayhits:
        print ' '.join(hit1)

#Read in and update the psl output and find hits
def readpsl(pslout,pbidict):
    bhits=[]
    hits=[]
    with open(pslout,'r') as file:
        for line in file:
            line=line.rstrip().split()
            ## Here are the filtering conditions: pident has to be than 50% and the maximum common subsequence overlap of the querry and subject >=40%
            ## These criteria could be tuned and changed as desired
            if float(line[2])>50.0 and max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])/float(line[3])>=0.40:
                if line[0]!=line[1]:
                    maxover=max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])
                    overlap=maxover/float(line[3])
                    ### if the maximum coverage subsequence overlap % is 100%: this is the best hit
                    if overlap>=1.0:
                        bhits+=[(line[0],pbidict[line[0]],line[1],pbidict[line[1]])]
                    ### otherwise record the hit and then compare reciprocity later on
                    else:
                        hits+=[(line[0],pbidict[line[0]],line[1],pbidict[line[1]])]

    ###Check reciprocal hits
    onehits=copy.copy(hits)
    for h in range(len(hits)-1):
        for k in range(h,len(hits)):
            if hits[h][0]==hits[k][1]: ## if reciprocal ie (1,2)=(2,1) then best hit and remove from onehits list
                bhits.append(hits[h])
                onehits.remove(hits[h])

    #     print line[0],pbidict[line[0]],line[1],pbidict[line[1]],' '.join(line[2:-1]),maxover, maxover/float(line[3])
    return bhits,onehits
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

