from __future__ import division 
import re, sys, math, itertools, string, subprocess, operator, copy,collections
from itertools import chain
import numpy as np
from itertools import groupby

def ProcessCLI(args):
    pslout=createpsl(args[1])
    pbidict=readgff(args[2])[1]
    enout=includeScaffTopsl(pslout,pbidict,args[3])
    sssd=getSSS(enout,pbidict)[1]
    anotssd=sssCheck(args[4],pbidict)
    printSSS(anotssd,'ohnologs5040')
    printSSS(sssd,'biscaffolds5040') 
    candidates=findCandidates(get_rbsssh(enout,pbidict,sssd))
    bestpairs,bpL=getbesthits(candidates)
    rbbh,oneways=getrbbsh(bestpairs)
    hits=rbbh+oneways
    printHits(hits,'bihits5040')
def printSSS(d,out):
    o_file=open(out,'w')
    for k in d.keys():
        if type(d[k])is dict:
            for k2 in d[k].keys():
                o_file.write(k +'\t'+k2+'\t'+str(d[k][k2])+'\n')
        else:
            o_file.write(k +'\t'+str(d[k])+'\n')
    o_file.close()
def printHits(hits,out):
    o_file=open(out,'w')
    for h in hits:
        o_file.write(h +'\n')
    o_file.close()
def getrbbsh(bestpairs):
    rbbsh=[]
    oneways=[]
    for t in bestpairs.keys(): 
        if (t[1],t[0]) in bestpairs.keys():
            rbbsh.append(' '.join([str(i) for i in list(list(t)+list(bestpairs[t]))]))
        else:
            oneways.append(' '.join([str(i) for i in list(list(t)+list(bestpairs[t]))]))
    return rbbsh,oneways
def getbesthits(candidates):
    bestpairs=[]
    bpairs={}
    for k in candidates:
        bestsss=200.0
        bestpair=None
        k2=None
        bestbit=0.0
        for v in candidates[k].keys():
            
            if float(candidates[k][v][2])<=bestsss and float(candidates[k][v][0])>bestbit:
                bestbit=float(candidates[k][v][0])
                bestsss=float(candidates[k][v][2])
                k2=v
        bestpairs.append((k,k2,candidates[k][k2]))
        bpairs[(k,k2)]=candidates[k][k2]
    return bpairs, bestpairs 
## Get candidates for each protein
def findCandidates(rbsssh):
    candidates={}
    for p1 in rbsssh.keys():
        minval=min([float(rbsssh[p1][p2][1]) for p2 in rbsssh[p1].keys()])
        uplimit=float(minval*10**15)
        check=0
        for p2 in rbsssh[p1].keys():
            if float(rbsssh[p1][p2][1])<=uplimit:
                if p1 not in  candidates.keys():
                    candidates[p1]={}
                    candidates[p1][p2]=rbsssh[p1][p2]
                else:
                    candidates[p1][p2]=rbsssh[p1][p2]
    return candidates
## Get Protein(key1):key2:(bscore,eval,sss)
def get_rbsssh(en_pslout,pbidict,sssd2):
    ppevsssd=collections.OrderedDict()
    for line in en_pslout:
        line=line.rstrip().split()
        jfident=(float(line[-4])/(min(float(line[5]),float(line[6]))))
        jfidentl=(float(line[-5])/(min(float(line[5]),float(line[6]))))
        ov=max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])
        if jfidentl>0.5 and float(line[4])>40.0:
            if line[0]!=line[2] and line[1]!=line[3]:
                sc=findSSS(pbidict[line[0]],pbidict[line[2]],sssd2)
                if line[0] not in ppevsssd.keys():
                    ppevsssd[line[0]]={} 
                    ppevsssd[line[0]][line[2]]=(float(line[-3]),float(line[-2]),sc)
                elif line[0] in ppevsssd.keys() and line[2] not in ppevsssd[line[0]].keys():
                    ppevsssd[line[0]][line[2]]=(float(line[-3]),float(line[-2]),sc)
                elif float(ppevsssd[line[0]][line[2]][0])> float(line[-3]):
                    continue 
                else:
                    ppevsssd[line[0]][line[2]]=(float(line[-3]),float(line[-2]),sc)
    return  ppevsssd
## Get ScaScaSc (SSS)::(sca_sca_list,Sca1(key):list of (sca1,sca2,sss),sssList)
def getSSS(en_pslout,pbidict):
    sssd=collections.OrderedDict()
    sca_List_sca=collections.OrderedDict()
    sssList=[]
    for line in en_pslout:
        line=line.rstrip().split()
        jfident=(float(line[-4])/(min(float(line[5]),float(line[6]))))
        jfidentl=(float(line[-5])/(min(float(line[5]),float(line[6]))))
        ov=max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])
        if jfidentl>0.5 and float(line[4])>40.0:
            if line[0]!=line[2] and line[1]!=line[3]: 
                if line[1] not in sca_List_sca.keys():
                    sca_List_sca[line[1]]=[line[3]]
                else:
                    sca_List_sca[line[1]].append(line[3])
#Find the scaffold_scaffold score pairs
    for key in sca_List_sca:
        kval=collections.Counter(sca_List_sca[key])
        sca_List_sca[key]=dict([a,-math.log(float(x)/sum(kval.values()))] for a, x in kval.iteritems())
    for k1 in sca_List_sca.keys():
        sssd[k1]=[]

        for k2 in sca_List_sca[k1].keys():
            score1=0
            score2=0
            try: 
                score1=sca_List_sca[k1][k2]
            except:
                score1=100
            try:
                score2=sca_List_sca[k2][k1]
            except:
                score2=100
            score=score1+score2
            sssd[k1].append((k1,k2,score))
            sssList.append((k1,k2,score))
    return sca_List_sca, makeDict(sortDictTu(sssd,2)),sortDictTu(sssList,2)
## Run command line 
def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout
def sortDictTu(indic,n):
    if type(indic)==dict:
        for k in indic:
            indic[k]=sorted(indic[k],key=lambda tup:(tup[n]))
    elif type(indic)==list:
        indic=sorted(indic,key=lambda tup:(tup[n]))
    return indic
def makeDict(sssd):
    newd={}
    for k in sssd.keys():
        newd[k]={}
        for k2 in sssd[k]:
           newd[k][k2[1]]=k2[2]
    return newd 
def findSSS(s1,s2,sssd2):
    try:
        x=float(sssd2[s1][s2])
        y=float(sssd2[s2][s1])
        return y
    except KeyError:
        return 200.0

def getInt(pb):
    return [int(i) for i in re.split(r'(\d+)',pb) if len(i) >0 and i.isdigit()][0]
## Get enriched output : scaffolds included
def includeScaffTopsl(pslout,pbidict,en_output):
    outpsl=open(en_output, 'w')
    output=[]
    pslout=pslout.split('\n')
    for line in pslout:
        line=line.rstrip().split()
        outpsl.write(' '.join([line[0],pbidict[line[0]],line[1],pbidict[line[1]],' \t'.join(line[2\
:])])+'\n')
        output.append(' '.join([line[0],pbidict[line[0]],line[1],pbidict[line[1]],' \t'.join(line[2:])]))
    outpsl.close()
    return output
### run blast and get blast output
def createpsl(infile):
    out=' -outfmt \"6 qseqid sseqid pident qlen slen qstart qend sstart send length nident bitscore evalue btop\"' 
    pslout=subprocess_cmd('blastp -query '+ infile + ' -subject '+ infile + out)   
    return pslout
### Read the gff3 files and create a dictionary -->scaffolds:list of proteins and pbi:scaffold, the scaff details info 
def readgff(gff):
    filename = open(gff, "r")
    scaffolds=collections.OrderedDict()
    pbi=collections.OrderedDict()
    scaff=collections.OrderedDict()
    for line in filename:
        line=line.strip().split() 
        if line[0][:8]=='scaffold' and line[2]=='mRNA':
            pbi[re.split(r'[=;\s]',line[8])[1]]=line[0]
            if line[0] not in scaffolds:
                scaffolds[line[0]]=[re.split(r'[=;\s]',line[8])[1]]
            else:
                scaffolds[line[0]].append(re.split(r'[=;\s]',line[8])[1])
        if line[0] not in scaff.keys():                                         
            scaff[line[0]]=[]
        else:
            scaff[line[0]]=[str(0),line[4]]
    return scaffolds, pbi, scaff
def sssCheck(infile,pbi):
    sca_List_sca=collections.OrderedDict()
    sssd={}
    sssList=[]
    filename = open(infile,"r")
    for line in filename:
        if len(line)>0 and 'hno' not in line and "." not in line:
            line=line.strip().split()
            if pbi[GtoT(line[0])] not in sca_List_sca.keys():
                sca_List_sca[pbi[GtoT(line[0])]]=[pbi[GtoT(line[1])]]
            else:
                sca_List_sca[pbi[GtoT(line[0])]].append(pbi[GtoT(line[1])])
    for key in sca_List_sca:
        kval=collections.Counter(sca_List_sca[key])
        sca_List_sca[key]=dict([a,-math.log(float(x)/sum(kval.values()))] for a, x in kval.iteritems())
    return sca_List_sca
def GtoT(pbid):
    pbid=list(pbid)
    if len(pbid)>5 and pbid[5]=='G':
        pbid[5]='T'
    return ''.join(pbid)
if __name__ == '__main__':
     ProcessCLI(sys.argv)
