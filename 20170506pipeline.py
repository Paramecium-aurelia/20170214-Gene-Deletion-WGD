import re, sys, math, itertools, string, subprocess, operator, copy,collections
from itertools import chain
from __future__ import division
import numpy as np
from itertools import groupby


def ProcessCLI(args):
    makesubfile(args[1])


def getSSS(pslout,pbidict):
    sssd=collections.OrderedDict()
    sca_List_sca=collections.OrderedDict()
    sssList=[]
    with open(pslout,'r') as file:
        for line in file:
            line=line.rstrip().split()
            ov=max([ int(i) for i in re.split(r'(\d+)',line[-1]) if len(i)>0 and i.isdigit()])
            if float(line[4])>40.0 and ov/min(float(line[5]),float(line[6]))>=0.12:

                if line[0]!=line[2]:
                    if line[1] not in hits.keys():
                        sca_List_sca[line[1]]=[line[3]]
                    else:
                        sca_List_sca[line[1]].append(line[3])

    for key in sca_List_sca:
        kval=collections.Counter(sca_List_sca[key])
        hits[key]=dict([a,-math.log(float(x)/sum(kval.values()))] for a, x in kval.iteritems())
    for k1 in sca_List_sca:
        for k2 in sca_List_sca[k1]:
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

    return sca_List_sca, sortDictTu(sssd,2),sortDictTu(sssList,2)

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


def getInt(pb):
    return [int(i) for i in re.split(r'(\d+)',pb) if len(i) >0 and i.isdigit()][0]

## Get enriched output : scaffolds included
def includeScaffTopsl(pslout,pbidict,en_output):
    outpsl=open(en_output, 'w')
    output=[]
    with open(pslout,'r') as file:
        for line in file:
            line=line.rstrip().split()
            outpsl.write(line[0],'\t',pbidict[line[0]],'\t',line[1],'\t',pbidict[line[1]],' \t'.join(line[2:]))
            output.append(' '.join([line[0],pbidict[line[0]],line[1],pbidict[line[1]],' \t'.join(line[2:])]))
    outpsl.close()
    return output

### run blast and get blast output
def createpsl(infile):
    out=' -outfmt \"6 qseqid sseqid pident qlen slen qstart qend sstart send bitscore evalue btop\"' 
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




if __name__ == '__main__':
     ProcessCLI(sys.argv)
