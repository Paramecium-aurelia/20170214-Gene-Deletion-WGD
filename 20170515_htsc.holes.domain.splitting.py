from __future__ import division
import re, sys, math, operator, copy,collections
import numpy as np 
from itertools import groupby

def ProcessCLI(args):
    scaf_pb=readgff(args[1])
    bestHits= readpsl(args[2],scaf_pb)
    #for k in bestHits:
    #    print ' \t'.join([str(i) for i in list(k)])

def getInt(pb):
    return [int(i) for i in re.split(r'(\d+)',pb) if len(i) >0 and i.isdigit()][0]

def getdup1(lst):
    d={}
    for key, value in lst:
        d[key] = d.get(key, []) + [value]
    
    for k in d:
        if len(d[k])==1:
            print k, ','.join(d[k])
        else:
            print 'dup:', k, '-->', ' '.join(d[k])
    #return pdup(d)
def domainsplit_tandemdup_ohno(htscList,en_psl,gff):#if a single protein is hit by multiple other proteins when we compare pairs
#of scaffolds, we infer that there is a tandem duplication, or a domain splitting of one large protein, and choose a single
#best mate to form the ohnologan
    bestmultL=[]
    noscorebestmultL=[]
    casesL=[]
    for htsc in htscList:
        doubles0=[]; doubles1=[]; mult={}
        for tup in htsc:
            if tup[0] not in doubles0:
                doubles0+=[tup[0]]
            else:
                mult.add(tup[0])
            if tup[1] not in doubles1:                
                doubles1+=[tup[1]]
            else:
                mult.add(tup[1])
        for tup in htsc:
            if tup[0] in mult:
                if tup[0] not in doublesD.keys():
                    doublesD[tup[0]]=set()
                else:
                    doublesD[tup[0]].add(tup[1])
            else:
                continue
            if tup[1] in mult:
                if tup[1] not in doublesD.keys():
                    doublesD[tup[1]]=set()
                else:
                    doublesD[tup[1]].add(tup[0])                
            else:
                continue

        for k in doublesD.keys():
            klen=0
            vlentot=0
            with open(gff,'r') as gfile:
                for line in gfile:
                    if k in re.split('\W+',line[8]):
                        klen+=float((math.fabs(line[4] - line[3])))
                    for v in doublesD[k]:
                        if v in re.split('\W+',line[8]):
                            vlentot+=float((math.fabs(line[4] - line[3])))
            if klen/vlentot > 1.5:
                casesL+=[(k,doublesD[k],"tandem duplication")]
            else:
                casesL+=[(k,doublesD[k]),"domain split"]
            #dummy=[]
            with open(en_psl,'r') as psl:
                for line in psl:
                    if ((line[0] == k and line[2] in doublesD[k]) or (line[0] in doublesD[k] and line[2]==k)) and line[1]!=line[3]:
                        bestmultL.append((line[0],line[2],(float(line[-3]),float(line[-2]))))#we now have a tuple with the ohnologan, & its bitscore&ev
                        noscorebestmultL.append((line[0],line[2]))
        #bestmultL+=bestdoubles; #noscorebestmultL+=
    return bestmultL, noscorebestmultL, casesL
    ##################################################################################################
    #what we need to do now:
    #1)  the items in bestmultL are (p1, p2, (bscore,escore)) triples
        #so, let's make 2-part tuples: for t in bmL: twomultL.append((t{0], t[1]))
        #then: for htsc in htscList:
            #build a list of mults
            #for t in htsc:
            #if t has a mult and t not in bmL and (t[1], t[0]) not in bmL:
                #delete the tuple from htsc---it was a tandem duplication or a split domain, so that isn't an ohnolog
                #if we need to know which was the case, we also have a list somewhere saying that (t[0], t[1], t[2]) was one of those cases
                                    #where t[0] is 1 protein,t[1] is the set of other proteins in the gff that it hits, and t[2] is a string indicated which case
def remove_ds_td_extras(htscList, nsbmL):#see comments above
    singly=htscList.copy()
    for htsc in singly:
        doubles0=[]; doubles1=[]; mult={}
        for tup in htsc:
            if tup[0] not in doubles0:
                doubles0+=[tup[0]]
            else:
                mult.add(tup[0])
            if tup[1] not in doubles1:
                doubles1+=[tup[1]]
            else:
                mult.add(tup[0])
        for tup in htsc:
            for p in tup:
                if p in mult:
                    if tup or (tup[1],tup[0]) in nsbmL:
                        continue
                    else:
                        #################################################
                        #DELETE THE TUPLE FROM THE htsc!!!!!!!!!!!!!!!!!!!!!!!!!
                        #But will this f**k up the indexing by which the for loop works? I hope not
                        #we will have to delete by value, not by index, that is for sure
                        #maybe use:
                        htsc.remove(tup)
    return singly
    
def fillht(htscList, en_pslout_file):#infer synteny if hits between pairs of scaffolds are uniformly increasing/incrasing
#or uniformly increasing; if there are gaps in the synteny find appropriate hits that fill the holes and mark the pairs as
#ohnologans
    holesList=[]
    l2hL=[]
    l2thhL=[]
    twohitsholesList=[]
    hole1way=[]
    for htsc in htscList:
        if len(htsc)>2:
            for i in range(1,len(htsc)-1):
                if getInt(htsc[i-1][0])< getInt(htsc[i][0]) and getInt(htsc[i][0])< getInt(htsc[i+1][0]):
                    if getInt(htsc[i-1][1])< getInt(htsc[i][1]) and getInt(htsc[i][1])< getInt(htsc[i+1][1]):
                        if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])<(getInt(htsc[i][1])-1):#if there is a hole on BOTH scafs
        #there was a gap on the scaffold. now we want to fill it in.
                            '''if getInt(htsc[i-1][0])-(getInt(htsc[i][0])-1) == getInt(htsc[i-1][1])-getInt(htsc[i][1]):
        #if the gaps are the same length on both scaffolds, we try to fill in the holes one at a time
                                for index,protein_int in enumerate(range((getInt(htsc[i-1][0])+1),getInt(htsc[i][0]))):
                                    for line in en_pslout:
                                        line=line.rstrip().split()
                                        if (getInt(line[0]==protein_int) and getInt(line[2])==(getInt(htsc[i-1][1])+(index+1))\#if the pair exists at all in the
                                        or getInt(line[2]==protein_int) and getInt(line[0])==(getInt(htsc[i-1][1])+(index+1))):#enriched psl output
                                            ##ADD THAT PAIR TO THE LIST OF OHNOLOGANS
                                            ##BUT IT DON'T KNOW WHAT THAT LIST IS CALLED---HELP ME, ETIENNE
                            elif getInt(htsc[i-1][0])-(getInt(htsc[i][0])-1) != getInt(htsc[i-1][1])-getInt(htsc[i][1]):'''
                            indivcol1=set(range(getInt(htsc[i-1][0])+1),getInt(htsc[i][0]))
                            indivcol2=set(range(getInt(htsc[i-1][1])+1),getInt(htsc[i][1]))
                            mycol=indivcol1.union(indivcol2)
                            
                            dummy=[]
                            with open(en_pslout_file,'r') as en_pslout:
                                for line in en_pslout:
                                    
                                    if getInt(line[0]) in mycol and getInt(line[2]) in mycol and line[1]!=line[3]:
                                        dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))#we now have a tuple with the ohnologan, & its bitscore&ev
                            dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                            allowed=set()#this code makes sure that we don't keep duplicates of a reciprocal hit
                            kill_it=[]
                            for tup in dummy_s:
                                if tup[0] not in allowed and tup[1] not in allowed:
                                    allowed.add(tup[0]);allowed.add(tup[1])
                                else:
                                    kill_it.append(tup)
                            for tupk in kill_it:
                                dummy_s.remove(tupk)
                            dummy_top=dummy_s[0:(min(len(indivcol1),len(indivcol2))+1)].
                            length2=[]
                            for trip in dummy_top:
                                length2.append((trip[0],trip[1]))
                            l2hL+=[length2]
                            holesList+=[dummy_top]
                    
                    elif getInt(htsc[i-1][1])> getInt(htsc[i][1]) and getInt(htsc[i][1])> getInt(htsc[i+1][1]):#case: scafs in reverse order
                        if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])>(getInt(htsc[i][1])+1):#if there is a hole on BOTH scafs
                            indivcol1=set(range(getInt(htsc[i-1][0])+1),getInt(htsc[i][0]))
                            indivcol2=set(range(getInt(htsc[i-1][1])-1),getInt(htsc[i][1]),-1)
                            mycol=indivcol1.union(indivcol2)
                            dummy=[]
                            with open(en_pslout_file,'r') as en_pslout:
                                for line in en_pslout:

                                    if getInt(line[0]) in mycol and getInt(line[2]) in mycol and line[1]!=line[3]:
                                        dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))#we now have a tuple with the ohnologan, & its bitscore&ev
                            dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                            allowed=set()
                            kill_it=[]
                            for tup in dummy_s:
                                if tup[0] not in allowed and tup[1] not in allowed:
                                    allowed.add(tup[0]);allowed.add(tup[1])
                                else:
                                    kill_it.append(tup)
                            for tupk in kill_it:
                                dummy_s.remove(tupk)
                            dummy_top=dummy_s[0:(min(len(indivcol1),len(indivcol2))+1)].
                            length2=[]
                            for trip in dummy_top:
                                length2.append((trip[0],trip[1]))
                            l2hL+=[length2]
                            holesList+=[dummy_top]
        elif len(htsc)==2:
            for i in range(1,2):
                if getInt(htsc[i-1][0])< getInt(htsc[i][0]):
                    if getInt(htsc[i-1][1])< getInt(htsc[i][1]):
                        if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])<(getInt(htsc[i][1])-1):#if there is a hole on BOTH scafs
        #there was a gap on the scaffold. now we want to fill it in.
                            indivcol1=set(range(getInt(htsc[i-1][0])+1),getInt(htsc[i][0]))
                            indivcol2=set(range(getInt(htsc[i-1][1])+1),getInt(htsc[i][1]))
                            mycol=indivcol1.union(indivcol2)
                            
                            dummy=[]
                            with open(en_pslout_file,'r') as en_pslout:
                                for line in en_pslout:
                                    if getInt(line[0]) in mycol and getInt(line[2]) in mycol and line[1]!=line[3]:
                                        dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))#we now have a tuple with the ohnologan, & its bitscore&ev
                            dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                            allowed=set()
                            kill_it=[]
                            for tup in dummy_s:
                                if tup[0] not in allowed and tup[1] not in allowed:
                                    allowed.add(tup[0]);allowed.add(tup[1])
                                else:
                                    kill_it.append(tup)
                            for tupk in kill_it:
                                dummy_s.remove(tupk)
                            dummy_top=dummy_s[0:(min(len(indivcol1),len(indivcol2))+1)].
                            length2=[]
                            for trip in dummy_top:
                                length2.append((trip[0],trip[1]))
                            l2thhL+=[length2]
                            twohitsholesList+=[dummy_top]
                    
                    elif getInt(htsc[i-1][1])> getInt(htsc[i][1]) and getInt(htsc[i][1])> getInt(htsc[i+1][1]):#case: scafs in reverse order
                        if getInt(htsc[i-1][0])<(getInt(htsc[i][0])-1) and getInt(htsc[i-1][1])>(getInt(htsc[i][1])+1):#if there is a hole on BOTH scafs
                            indivcol1=set(range(getInt(htsc[i-1][0])+1),getInt(htsc[i][0]))
                            indivcol2=set(range(getInt(htsc[i-1][1])-1),getInt(htsc[i][1]),-1)
                            mycol=indivcol1.union(indivcol2)
                            
                            dummy=[]
                            with open(en_pslout_file,'r') as en_pslout:
                                for line in en_pslout:
                                    if getInt(line[0]) in mycol and getInt(line[2]) in mycol and line[1]!=line[3]:
                                        dummy.append((line[0],line[2],(float(line[-3]),float(line[-2]))))#we now have a tuple with the ohnologan, & its bitscore&ev
                            dummy_s=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)
                            allowed=set()
                            kill_it=[]
                            for tup in dummy_s:
                                if tup[0] not in allowed and tup[1] not in allowed:
                                    allowed.add(tup[0]);allowed.add(tup[1])
                                else:
                                    kill_it.append(tup)
                            for tupk in kill_it:
                                dummy_s.remove(tupk)
                            dummy_top=dummy_s[0:(min(len(indivcol1),len(indivcol2))+1)].
                            length2=[]
                            for trip in dummy_top:
                                length2.append((trip[0],trip[1]))
                            l2thhL+=[length2]
                            twohitsholesList+=[dummy_top]  #i think this will make thhL a list of lists, one list for each sca-sca pair
                            #which may be good, or bad, but it is easy to change
                            #twohitsholesList+=sorted(dummy,key=lambda tup:(tup[2][0]),reverse=True)[0:(min(len(indivcol1),len(indivcol2))+1)]
    return holesList, l2hL, twohitsholesList, l2thhL

def printht(htsc):
    md=[]
    ml=[]
    if len(htsc)>2:                                                                        
        for i in range(1,len(htsc)-1):
            if getInt(htsc[i-1][0])< getInt(htsc[i][0]) and getInt(htsc[i][0])< getInt(htsc[i+1][0]):
                if getInt(htsc[i-1][1])< getInt(htsc[i][1]) and getInt(htsc[i][1])< getInt(htsc[i+1][1]):
                    if i==1:
                        ml+=[htsc[i-1],htsc[i],htsc[i+1]]
                    elif i>1:
                        ml+=[htsc[i+1]]
                elif getInt(htsc[i-1][1])> getInt(htsc[i][1]) and getInt(htsc[i][1])> getInt(htsc[i+1][1]):
                    if i==1:
                        ml+=[htsc[i-1],htsc[i],htsc[i+1]]
                    elif i>1:
                        ml+=[htsc[i+1]]
                else:
                    print 'failed:', htsc[i]
        print ml
    elif len(htsc)==2:                                                                             
        if getInt(htsc[0][0])< getInt(htsc[1][0]):
            if getInt(htsc[0][1])> getInt(htsc[1][1]):
                md+=[htsc[0]]
                md+=[htsc[1]]
            elif getInt(htsc[0][1])<getInt(htsc[1][1]):
                md+=[htsc[0]]
                md+=[htsc[1]]
        print md 
    else:
        print 'single : ', htsc 


def getdup2(lst):
    d={}
    dd={}
    for key, value in lst:
        d[value] = d.get(value, []) + [key]
    for k in d.keys():
        if len(d[k])==1:
            print ','.join(d[k]), k
        else:
            print 'dup: ',','.join( d[k]),"-->", k
def getdupboth(lst):
    d1={}
    d2={}
    
    for k,v in lst:
        if k in d1:
            d1[k]+=[v]
        else:
            d1[k]=[v]
        if v in d2 :
            d2[v]+=[k]
        else:
            d2[v]+=[k]
    for k1 in d1:
        print k1,':',d1[k1]
    for k2 in d2:
        #if len(d2[k2])==1:
        print d2[k2],':',k2
        
def pdup(d):
    for k in d:
        if len(d[k])==1:
            print (k,d[k])
        else:
            print 'dup:',(k,d[k])
def prTp(lst):
    for k in lst:
        print k[0],k[1]

def getCts(lst):
    recip=[]
    recipC=0
    oneway=[]
    onewayC=0
    clst=copy.copy(lst)
    recip_seen=[]
    for (k,v) in clst:
        if (v,k) in clst:
            if (v,k) not in recip_seen:
                recip_seen.append((k,v))
                recip_seen.append((v,k))
                recipC+=1
        elif (v,k) not in clst:
            oneway.append((k,v))
            onewayC+=1 
#    print 'Reciprocal hits : ', len(recip_seen),recipC 
#    print 'One way Hits : ', len(oneway), onewayC

def getDistrib(lst,diPL):
    d=[]
    seen=[]
    for k,v in lst:
        if (v,k) not in seen:
            
            t=float(diPL[k]) / float(diPL[v])
            if t > 1.0:
                t=1/t
            d.append(t)
            seen.append((k,v))
    return d            

#Read in and update the psl output and find hits
def readpsl(pslout,pbidict):
    bhits=[]
    hits=collections.OrderedDict()
    phits=collections.OrderedDict()
    scaffolds={}
    Lengthprot={}
    ps=collections.OrderedDict()
    with open(pslout,'r') as file:
        for line in file:
            line=line.rstrip().split()
            ov=max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])
            if float(line[4])>40.0 and ov/min(float(line[5]),float(line[6]))>=0.12:

                if line[0]!=line[2]:
                    if line[1] not in hits.keys():
                        hits[line[1]]=[line[3]]
                    else:
                        hits[line[1]].append(line[3])
                    if line[0] not in phits.keys():
                        phits[line[0]]={line[2]:line[-2]}
                    else:
                        phits[line[0]][line[2]]=line[-2]
                    if line[1] not in scaffolds.keys():
                        scaffolds[line[1]]=[line[0]]
                    else:
                        scaffolds[line[1]].append(line[0])
                    if line[3] not in scaffolds.keys():
                        scaffolds[line[3]]=[line[2]]
                    else:
                        scaffolds[line[3]].append(line[2])
                    if line[1]  in ps: 
                        if line[3] in ps[line[1]].keys():
                            ps[line[1]][line[3]]+=[(line[0],line[2],line[-2])]
                        else:
                            ps[line[1]][line[3]]=[(line[0],line[2],line[-2])]
                    else:
                        ps[line[1]]={}
                        ps[line[1]][line[3]]=[(line[0],line[2],line[-2])]
                    if line[0] not in Lengthprot:
                        Lengthprot[line[0]]=int(line[5])
                    if line[2] not in Lengthprot:
                        Lengthprot[line[2]]=int(line[6])
    for key in hits:
        kval=collections.Counter(hits[key])
        hits[key]=dict([a,-math.log(float(x)/sum(kval.values()))] for a, x in kval.iteritems())
    for k1 in hits: 
        for k2 in hits[k1]:
            score1=0
            score2=0
            try:
                score1=hits[k1][k2]
                
            except:
                score1=10
            try:
                score2=hits[k2][k1]
            except:
                score2=10
            score=score1+score2
            
            bhits.append((k1,k2,score))

    ghits=collections.OrderedDict()
    for p1 in phits.keys():
        minval=min([float(phits[p1][p2]) for p2 in phits[p1].keys()])
         
        uplimit=float(minval*10**10)
        check=0
        for p2 in phits[p1]:
            if float(phits[p1][p2]) <=uplimit:
                if p1 not in  ghits.keys():
                    ghits[p1]={}
                    ghits[p1][p2]=float(phits[p1][p2])
                else:
                    ghits[p1][p2]=float(phits[p1][p2])
            if float(phits[p1][p2])==minval and check==0:
                print p1,p2
                check=1

    betterHits=collections.OrderedDict()
    pairs=[]
#    scaffolds_score=sorted(bhits.items(),key=lambda kv: kv[1], reverse=True)
#    for k in bhits:
#        print k, bhits[k]
    return 
    for g in ghits.keys():
        score=15.0
        bh=None
        if pbidict[g] not in betterHits.keys():
            betterHits[pbidict[g]]={}
#        print 'P1: ', g
        for g2 in ghits[g].keys():
#            print 'P2',g2, ghits[g][g2],hits[pbidict[g]][pbidict[g2]]
            if hits[pbidict[g]][pbidict[g2]]<score:
                bh=g2
                score=hits[pbidict[g]][pbidict[g2]]
        #print 'Best Hit: ',g,bh, Scaffolds: ,pbidict[g],pbidict[bh]
#        print g,bh, 'Scaffolds: ',pbidict[g],pbidict[bh]
        pairs.append((g,bh))
#        print '          ',pbidict[g],pbidict[bh]
        if pbidict[bh] not in betterHits[pbidict[g]].keys():
            betterHits[pbidict[g]][pbidict[bh]]=[(g,bh)]
        else:
            betterHits[pbidict[g]][pbidict[bh]].append((g,bh))
 #       return
    getCts(pairs)
   # print getDistrib(pairs,Lengthprot) 
    for k in betterHits.keys():
        for g in betterHits[k].keys():
        #for g in ghits[k].keys():
            htsc=sorted(betterHits[k][g],key=lambda tup: (tup[0]))
            gh=[]
#            print k,g#,'\n'
            #print htsc
            '''
            s1=[i[0] for i in htsc]
            s2=[i[1] for i in htsc]
            if len(htsc)==len(set(s1)) and len(htsc)==len(set(s2)):
               # print len(htsc), len(set(s1)), len(set(s2))
                #print ' '.join(list(htsc))
                #prTp(htsc)
           #     print 'case 1'
                printht(htsc)
                #getscaf(htsc)
               # continue
               # if len(set(s1))==1 and len(set(s2))==1:
               #     print 'singleton: ', htsc 
               # if len(set(s1))==2 and len(set(s2))==2:
               #     print ' duo: ', htsc
            elif len(htsc)!=len(set(s1)) and len(htsc)!=len(set(s2)):
            #    print 'case 2'
                getdupboth(htsc)
                
            else:
                if len(htsc)==len(set(s1)):
             #       print 'case 3'
                    getdup2(htsc)
                elif len(htsc)==len(set(s2)):
              #      print 'case 4'
                    getdup1(htsc)
                #print len(htsc), len(set(s1)), len(set(s2))
                #print htsc
               '''
#                    print sorted(betterHits[k][g],key=lambda tup: (tup[0]))
#            print betterHits[g][k]
    #return
#    return sorted(bhits, key=lambda tup: (tup[2]) )#bhits

## get scaffolds dictionary
def readScaf(psl):
    scaffolds={}
    filename=open(psl,"r")
    for line in filename:
        line=line.rstrip().split()
        if line[1] in scaffolds and line[3] in scaffolds[line[1]]:
            scaffolds[line[1]][line[3]]+=[(line[0],line[2],line[-2])]
        else:
            scaffolds[line[1]]={line[3]:[(line[0],line[2],line[-2])]}

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
 
#if __name__ == '__main__':
     ProcessCLI(sys.argv)
