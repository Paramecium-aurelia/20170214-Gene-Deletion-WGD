import re, sys, math, operator, copy, collections 

def ProcessCLI(args):
 
   #readFile(args[1],args[2])
    btree,utree=readOutput(args[1])
    best,ones=checkrbbh(readOutput(args[2])[0])
    print 'Ohno length: ', len(btree.keys())
    print 'Roy length: ' , len(best.keys())
    print 'Ohno ^ Roy ', len(set(btree.keys()).intersection(set(best.keys())))
    print 'Ohno - Roy ', len(set(btree.keys())-set(best.keys()))
    print 'Roy - Ohno ',len(set(best.keys())-set(btree.keys()))
    print 'Ohno 1. .1 ', len(utree)
    print 'Roy_owh ',len(ones.keys()) 
    print 'O1 ^ R1 ', len(set(utree).intersection(set(ones.keys())))

def checkrbbh(d):
    bestrbbh={}
    onewayhits={}
    seen=[]
    weirdhits={}
#    print len(d)
#    return 
    for k in d.keys():
        if d[k] in d.keys() and d[d[k]]==k and k not in seen and d[k] not in seen:
            bestrbbh[k]=d[k]
            seen.append(k)
            seen.append(d[k])
#            print 'best case: ', k, d[k]
        if d[k] not in d.keys():
            onewayhits[k]=d[k]
        if d[k] in d.keys() and d[d[k]]!=k:
            weirdhits[k]=d[k]
            onewayhits[k]=d[k]
            #print 'wierd', k, d[k]
            #print 'weird', d[k], d[d[k]]
            #if d[k] == d[d[d[k]]]:
             #   print "Roy is happy."
            #else:
            #    print 'super w', d[d[k]],d[d[d[k]]]
#    al=set(onewayhits.keys()).intersection(set(weirdhits.keys()))
#    print len(set(onewayhits)), len(set(weirdhits))
 #   print len(al)
#for p in weirdhits:
    #    print p 
        #if p in weirdhits:
      #      print 'sionfusing', bestrbbh[p]
    print len(bestrbbh)
    print len(onewayhits)
    return bestrbbh, onewayhits

#Read the gff3 files and create a dictionary: scaffolds: list of proteins                                   
def readFile(f1,f2):
    tree=readOutput(f1)
    output=readOutput(f2)
    same=0
    uniqtree=0
    diff=0
    print 'Input 1 length : ', len(tree)
    print 'Input 2 length : ', len(output)

    for k in tree.keys():
        if k in output:
            if tree[k]==output[k]:
                same+=1
            if tree[k]!=output[k]:
                diff+=1
                #print '1:',k, tree[k]
                #print '2:',k,output[k]
        elif k not in output.keys() and tree[k]!=".":
            uniqtree+=1
    #print 'Matched : ', same
    #print 'Differ  : ', diff 
    #print 'Unique to first input file : 'args[1], uniqtree

def readOutput(infile):
    d=collections.OrderedDict()
    d1=[]
    filename = open(infile,"r")
    for line in filename:
        if len(line)>0:
            line=line.strip().split()
            
            if GtoT(line[0]) not in d.keys() and "." not in GtoT(line[0]) and "." not in GtoT(line[1]):
                d[GtoT(line[0])]=GtoT(line[1])
            else:
                if "." not in GtoT(line[0]) and "." in GtoT(line[1]):
                    d1.append(GtoT(line[0]))
                if "." in GtoT(line[0]) and "." not in GtoT(line[1]):
                    d1.append(GtoT(line[1]))
    return d,d1 

def GtoT(pbid):
    pbid=list(pbid)
    if len(pbid)>5 and pbid[5]=='G':
        pbid[5]='T'
    return ''.join(pbid)

if __name__ == '__main__':
     ProcessCLI(sys.argv)
