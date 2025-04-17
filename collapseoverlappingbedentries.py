import sys
import os
import pandas as pd
from skbio import DNA


def mergedegeneratebases(b1,b2):
    b1set=inversebasekey[b1]
    b2set=inversebasekey[b2]
    combineset=tuple(sorted(list(set(list(b1set) + list(b2set)))))
    bc=basekey[combineset]
    return(bc)

def mergewithdegseqalign(seq1,seq2):
    assert(len(seq1)==len(seq2))
    newseq=[]
    for pos in range(len(seq1)):
        if seq1[pos]==seq2[pos]: newseq.append(seq1[pos])
        else:newseq.append(mergedegeneratebases(seq1[pos],seq2[pos]))
    newseq="".join(newseq)
    return(newseq)

def mergetwoseqfrags(first,second,offset):
    stopold,offsetstart=len(first),offset
    startnew,stopnew,offsetstop=0,len(second),len(first)-offset#len(second)-(stopold-offset+1)
    res=first[0:offsetstart]+mergewithdegseqalign(first[offsetstart:stopold],second[startnew:offsetstop])+second[offsetstop:stopnew]
    return(res)



def writeentry(newdf,entry):
    return(newdf.append({0: entry[0], 1: entry[1], 2: entry[2],3: entry[3]}, ignore_index=True))

basekey = {tuple(sorted(values)):keys for keys,values in DNA.degenerate_map.items()}
basekey.update({('A',):'A',('G',):'G',('C',):'C',('T',):'T'})
inversebasekey= {basekey[values]:values for values in basekey}

file=sys.argv[1]
assert(os.path.isfile(file))

df = pd.read_csv(file, sep='\t', header=None)
df=df.sort_values(by=1)
assert(len(df.columns)==5)

currententry=[]
mergedentries = pd.DataFrame(columns=[0,1,2,3])

for pos in range(len(df)):
    if len(currententry)==0: currententry=list(df.iloc[pos,:])
    potfrag=df.iloc[pos,:]
    if (currententry[1]==potfrag[1]) and (currententry[2]==potfrag[2]):continue
#detect overlap between current interval and fragment
    if (currententry[2]>potfrag[1]) and (currententry[2]<potfrag[2]):
        currententry[2]=potfrag[2]
        mergedseq=mergetwoseqfrags(currententry[3],potfrag[3],potfrag[1]-currententry[1])
        currententry[3]=mergedseq
        #add seq to currententry
    else:
        mergedentries=writeentry(mergedentries,currententry)
        currententry=list(df.iloc[pos,:])

#write last interval
mergedentries=writeentry(mergedentries,currententry)

mergedentries.to_csv((file.split(".")[0])+"_merged.bed", sep='\t',header=False,index=False)
