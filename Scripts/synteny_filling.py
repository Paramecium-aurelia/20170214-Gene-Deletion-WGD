from __future__ import division
import re, sys, math, operator, copy,collections
import numpy as np 
from itertools import groupby

def ProcessCLI(args):
    ## comparison of the real McGrath output to the predicted and gap filled output
    pbi,scaf=readgff(args[1]) ## use gff file to find the protein and their respective scaffolds and verse versa
    htscd=createHTSC(args[2],pbi) ## This is to create the htscList in a Dict format
    htsc=[] ##This creates the input for the synteny-filling function (fillht), which is a List of Lists of Lists of Tuples (which are protein-protein hits, see below)
    for s1 in htscd.keys():
        LL=[]
        for s2 in htscd[s1].keys():
            LL+=[htscd[s1][s2]]
        htsc+=[LL]

    predicted_paralogs= fillht(htsc,args[3],pbi) 
    McGrath_paralogs,d=readOutput(args[4])  
    match=comp(McGrath_paralogs,predicted_paralogs)
    match2=comp(predicted_paralogs,McGrath_paralogs)
    McGrath_paralogsRev = [(t[1],t[0]) for t in McGrath_paralogs]#Consider each hit in both directions 

    print 'Pairs from McG that Roy recovered : ', match
    print 'Pairs from Roy that McG covered : ', match2
    print '-----Reversed pairs on McG-----------------'
    m=comp(CpairsRev,Rpairs)
    m2=comp(Rpairs,CpairsRev)
    print 'Pairs from McG that Roy recovered : ', m
    print 'Pairs from Roy that McG covered : ', m2



def comp(t1,t2):
    c=0
    for t in t1:
        if t in t2:
            c+=1
    return c 
    
### Read in the ohnolog input file
def readOutput(infile):
    d=collections.OrderedDict()
    pairs=[]
    filename = open(infile,"r")
    for line in filename:
        if len(line)>0:
            line=line.strip().split()

            if GtoT(line[0]) not in d.keys() and "." not in GtoT(line[0]):
                d[GtoT(line[0])]=GtoT(line[1])
                pairs.append((GtoT(line[0]),GtoT(line[1])))
    return pairs, d


def GtoT(pbid):
    pbid=list(pbid)
    if len(pbid)>5 and pbid[5]=='G':
        pbid[5]='T'
    return ''.join(pbid)



def getInt(pb):
    return [int(i) for i in re.split(r'(\d+)',pb) if len(i) >0 and i.isdigit()][0]
    
def fillht(htscL, en_pslout_file,pbi):#infer synteny if hits between pairs of scaffolds are uniformly decreasing
#or uniformly increasing; if there are gaps in the synteny find appropriate hits that fill the holes and mark the pairs as
#ohnologans
    holesList=[];l2hL=[];l2thhL=[];twohitsholesList=[]
    hole1way=[];alltg=[]#the structure of the 'htsc' is a [[[()]]]], list of lists of lists of tuples (also denoted LLLT);
                        #viz.,a List of scaffolds; for each scaffold a List of hit allo-scaffolds; for each scaffold-allo-scaffold pair, a List of tuples for which
                        #each Tuple represents a protein-protein hit
    newhtscList=[q for q in htscL]#make of copy of the incoming list
    for LT in newhtscList:
        for htsc in LT:
            b=len(htsc) ### count of htsc before gap filling
            print 'input htsc-->',len(htsc),sorted(htsc,key=lambda tup:tup[0])
            htsc=sorted(htsc,key=lambda tup:tup[0])#order each LT by gene number
            newhits=set()
            if len(htsc)>2:
                for i in range(1,len(htsc)-1):#here we identify any places in which any two adjacent good-hitting proteins from the first scaffold are not actually
#subsequent in number. This is explained in more detail in the L-534 paper Roy wrote. 
                    if getInt(htsc[i-1][0])< getInt(htsc[i][0]) and getInt(htsc[i][0])< getInt(htsc[i+1][0]):#case: both scaffolds are in ascending order
                        if getInt(htsc[i-1][1])< getInt(htsc[i][1]) and getInt(htsc[i][1])< getInt(htsc[i+1][1]):
                            if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])<(getInt(htsc[i][1])-1):
                                indivcol1=set(range(getInt(htsc[i-1][0])+1,getInt(htsc[i][0])))
                                indivcol2=set(range(getInt(htsc[i-1][1])+1,getInt(htsc[i][1])))
                                mycol=indivcol1.union(indivcol2)#we find every protein that should have had a good hit if there was perfect synteny
                                dummy=[]
                                with open(en_pslout_file,'r') as en_pslout:#for all the proteins that would appear in the case of perfect synteny, we try pairing them
#and choose the pairs with the best scores
                                    for line in en_pslout:
                                        line=line.rstrip().split()
                                        if getInt(line[0]) in indivcol1 and getInt(line[2]) in indivcol2 and pbi[line[0]]!=pbi[line[2]]:
                                            dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))
                                dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                                allowed=set()
                                kill_it=[]
                                for tup in dummy_s:
                                    if tup[0] not in allowed or tup[1] not in allowed:
                                        allowed.add(tup[0]);allowed.add(tup[1])
                                    else:
                                        kill_it.append(tup)
                                for tupk in kill_it:
                                    dummy_s.remove(tupk)
                                
                                length2=[(trip[0],trip[1]) for trip in dummy_s]
                                for trip in length2:
                                    newhits.add((trip[0],trip[1]))
                        
                                l2hL+=[length2]#we add new pairs to a list specifically for new pairs
                                holesList+=[dummy_s]
                        
                        elif getInt(htsc[i-1][1])> getInt(htsc[i][1]) and getInt(htsc[i][1])> getInt(htsc[i+1][1]):#case: scafs in reverse order
                            if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])>(getInt(htsc[i][1])+1):#if there is a hole on BOTH scafs
                                indivcol1=set(range(getInt(htsc[i-1][0])+1,getInt(htsc[i][0])))
                                indivcol2=set(range(getInt(htsc[i-1][1])-1,getInt(htsc[i][1]),-1))
                                mycol=indivcol1.union(indivcol2)
                                dummy=[]
                                with open(en_pslout_file,'r') as en_pslout:
                                    for line in en_pslout:
                                        line=line.rstrip().split()
                                        if getInt(line[0]) in list(indivcol1) and getInt(line[2]) in list(indivcol2) and pbi[line[0]]!=pbi[line[2]]:
                                            dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))#we now have a tuple with the ohnologan, & its bitscore&ev
                                dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                                allowed=set()
                                kill_it=[]
                                for tup in dummy_s:
                                    if tup[0] not in allowed or tup[1] not in allowed:
                                        allowed.add(tup[0]);allowed.add(tup[1])
                                    else:
                                        kill_it.append(tup)
                                for tupk in kill_it:
                                    dummy_s.remove(tupk)
                                length2=[(trip[0],trip[1]) for trip in dummy_s]
                                for trip in length2:
                                    newhits.add((trip[0],trip[1]))
                                l2hL+=[length2]
                                holesList+=[dummy_s]
                                
            elif len(htsc)==2:#special case in which there are two protein hits between the scaffolds
                for i in range(1,2):
                    if getInt(htsc[i-1][0])< getInt(htsc[i][0]):
                        if getInt(htsc[i-1][1])< getInt(htsc[i][1]):
                            if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])<(getInt(htsc[i][1])-1):
                                indivcol1=set(range(getInt(htsc[i-1][0])+1,getInt(htsc[i][0])))
                                indivcol2=set(range(getInt(htsc[i-1][1])+1,getInt(htsc[i][1])))
                                mycol=indivcol1.union(indivcol2)
                                dummy=[]
                                with open(en_pslout_file,'r') as en_pslout:
                                    for line in en_pslout:
                                        line=line.rstrip().split()
                                        if getInt(line[0]) in indivcol1 and getInt(line[2]) in indivcol2 and pbi[line[0]]!=pbi[line[2]]:
                                            dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))
                                dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                                allowed=set()
                                kill_it=[]
                                for tup in dummy_s:
                                    if tup[0] not in allowed or tup[1] not in allowed:
                                        allowed.add(tup[0]);allowed.add(tup[1])
                                    else:
                                        kill_it.append(tup)
                                for tupk in kill_it:
                                    dummy_s.remove(tupk)
                                
                                length2=[(trip[0],trip[1]) for trip in dummy_s]
                                for trip in length2:
                                    newhits.add((trip[0],trip[1]))
                                l2thhL+=[length2]
                                twohitsholesList+=[dummy_s]
                        
                        elif getInt(htsc[i-1][1])> getInt(htsc[i][1]):
                            if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])>(getInt(htsc[i][1])+1):
                                indivcol1=set(range(getInt(htsc[i-1][0])+1,getInt(htsc[i][0])))
                                indivcol2=set(range(getInt(htsc[i-1][1])-1,getInt(htsc[i][1]),-1))
                                mycol=indivcol1.union(indivcol2)
                                dummy=[]
                                with open(en_pslout_file,'r') as en_pslout:
                                    for line in en_pslout:
                                        line=line.rstrip().split()
                                        if getInt(line[0]) in indivcol1 and getInt(line[2]) in indivcol2 and pbi[line[0]]!=pbi[line[2]]:
                                            dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))
                                dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                                allowed=set()
                                kill_it=[]
                                for tup in dummy_s:
                                    if tup[0] not in allowed or tup[1] not in allowed:
                                        allowed.add(tup[0]);allowed.add(tup[1])
                                    else:
                                        kill_it.append(tup)
                                for tupk in kill_it:
                                    dummy_s.remove(tupk)
                                length2=[(trip[0],trip[1]) for trip in dummy_s]
                                for trip in length2:
                                    newhits.add((trip[0],trip[1]))
                                l2thhL+=[length2]
                                twohitsholesList+=[dummy_s] 
            for t in newhits:
                htsc+=[t]
                b+=1 ## update the number of new added htscs
            htsc=sorted(htsc,key=lambda tup: (tup[0]))
            alltg+=htsc
            print 'output htsc-->',len(htsc), htsc
        alltg=sorted(alltg,key=lambda tup:(tup[0]))
    return alltg, newhtscList #holesList, l2hL, twohitsholesList, l2thhL,newhtscList

## create htsc from best hits ouput
def createHTSC(bihits,pbidict):
    sca_sca_hits=collections.OrderedDict()
    with open(bihits,'r')as file:
        for line in file:
            line=line.rstrip().split()
            if pbidict[line[0]] not in sca_sca_hits.keys():
                sca_sca_hits[pbidict[line[0]]]=collections.OrderedDict()
                if pbidict[line[1]] not in sca_sca_hits[pbidict[line[0]]].keys(): 
                    sca_sca_hits[pbidict[line[0]]][pbidict[line[1]]]=[(line[0],line[1])]
                else: 
                    sca_sca_hits[pbidict[line[0]]][pbidict[line[1]]].append((line[0],line[1]))
            else:
                if pbidict[line[1]] not in sca_sca_hits[pbidict[line[0]]].keys():
                    sca_sca_hits[pbidict[line[0]]][pbidict[line[1]]]=[(line[0],line[1])]
                else:
                    sca_sca_hits[pbidict[line[0]]][pbidict[line[1]]].append((line[0],line[1]))
    return sca_sca_hits

#Read he gff3 files and create a dictionary: scaffolds: list of proteins                                   
def readgff(f):
    print 'Reading the gff3 file ...'
    filename = open(f, "r")
    scaffolds=collections.OrderedDict()
    pbi=collections.OrderedDict()
    for line in filename:
        line=line.strip().split()
        
        if line[0][:8]=='scaffold' and line[2]=='mRNA':
            pbi[re.split(r'[=;\s]',line[8])[1]]=line[0]
        if line[0] not in scaffolds.keys():
            scaffolds[line[0]]=[re.split(r'[=;\s]',line[8])[1]]
        else:
            scaffolds[line[0]].append(re.split(r'[=;\s]',line[8])[1])
    for k in scaffolds.keys():
        scaffolds[k]=sorted(scaffolds[k])
            
    return pbi,scaffolds  
 
if __name__ == '__main__':
    ProcessCLI(sys.argv)
