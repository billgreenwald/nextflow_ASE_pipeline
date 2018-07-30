
# coding: utf-8

# In[21]:


import sys
import statsmodels.api as sm
from collections import Counter
import scipy.stats as ss
import numpy as np
import argparse
import os
import bgmodule as bg


# In[ ]:


parser=argparse.ArgumentParser()
parser._optionals.title = "Flag Arguments"
parser.add_argument('-i', help='Input mpileup file', required=True,default=5000)
parser.add_argument('-o', help='Output file name', required=True)
parser.add_argument('-vcf', help='Input vcf to check genotypes against.  Optional argument', required=False,default='%#$')
parser.add_argument('-vcfD', help='Folder of input vcfs to check genotypes against.  VCF must be named "chr#.vcf[.gz]". Optional argument', required=False,default='%#$')
parser.add_argument('-id', help='Name of individual in VCF file to check genotypes against', required=False,default="%#$")

args = vars(parser.parse_args())


# In[12]:


# #debugging
# args={'i':'',
#       'o':'',
#       'vcf':'%#$',
#       'vcfD':'',
#       'id':''}


# In[13]:


if args['vcf']!="%#$" and args['vcfD']!="%#$":
    print "Both vcf and vcfD were provided.  Using -vcf argument"


# In[14]:


def getZScore(x,mu,var):
    return float(x-mu)/np.sqrt(var)


# In[15]:


def cleanPileupString(pileup):
    new=""
    for x in pileup:
        if x== 'A' or x == 'C' or x == 'G' or x == 'T':
            new+=x
    return new


# In[16]:


def processPileup(pileup,refBase):
    pileup=cleanPileupString(pileup)
    if len(pileup)==0:
        return 'nan','nan','nan','nan'
    x=Counter(pileup)
    mostCommon=x.most_common()[0]
    mca='.' if mostCommon[0]==refBase else mostCommon[0]
    if len(x)==1:
        if mca=='.':
            return mca,mostCommon[1],0,'nan'
        else:
            return mca,0,mostCommon[1],'nan'
    else:
        mostCommon2=x.most_common()[1]
        if mostCommon[0]==refBase:
            mca=mostCommon2[0] 
            return mca,mostCommon[1],mostCommon2[1],float(mostCommon[1])/(mostCommon[1]+mostCommon2[1])
        else:
            mca=mostCommon[0]
            return mca,mostCommon2[1],mostCommon[1],float(mostCommon2[1])/(mostCommon[1]+mostCommon2[1])


# In[42]:


output=[]
ps=[]
for line in open(args['i']):
    line=line.split()
    chrom=line[0]
    pos=line[1]
    refBase=line[2].upper()
    pileup=line[4].replace(',',refBase).replace('.',refBase).upper()
    mca,rc,ac,mar=processPileup(pileup,refBase)
    tot=rc+ac
    if mar!='nan':
        z=getZScore(rc,0.5*tot,0.5*0.5*tot)
        p=ss.norm.sf(abs(z))*2
        ps.append(p)
        output.append([chrom,pos,refBase,mca,rc,ac,mar,z,p])
    else:
        output.append([chrom,pos,refBase,mca,rc,ac,mar,'nan','nan'])


# In[18]:


def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)


# In[19]:


genotypes=[]
if args['vcf']!="%#$":
    for line in open(args['vcf']):
        pass


# In[29]:


genotypes=[]
if args['vcfD']!="%#$" and args['vcf']=="%#$":
    whichSNP=0
    while whichSNP < len(output):
        whichChrom=output[whichSNP][0]
        fn=args['vcfD']+'/'+whichChrom+'.vcf'
        if not os.path.exists(fn):
            fn+='.gz'
            if not os.path.exists(fn):
                whichSNP+=1
                genotypes.append('nan')
                continue
        bg.eprint("Checking Genotypes in file {0}".format(fn))
        for line in openfile(fn):
            if whichSNP>=len(output):
                break
            curPos=int(output[whichSNP][1])
            if line[:2]=='##':
                continue
            elif line[0]=='#':
                try:
                    genCol=line.split().index(args['id'])
                except:
                    print "{0} not found in header of VCF. Exiting".format(args['id'])
                    asfsdf
                    exit()
            else:
                line=line.split()
                position=int(line[1])
                if position < curPos:
                    continue
                elif position==curPos:
                    geno=line[genCol][:3]
                    genotypes.append(geno)
                    whichSNP+=1
                    if whichSNP>=len(output):
                        break
                    if output[whichSNP][0] != whichChrom:
                        break
                else:
                    while position > curPos:
                        genotypes.append('nan')
                        whichSNP+=1
                        if whichSNP>=len(output):
                            break
                        if output[whichSNP][0] != whichChrom:
                            break
                        curPos=int(output[whichSNP][1])
                    if position==curPos:
                        geno=line[genCol][:3]
                        geno=line[genCol][:3]
                        genotypes.append(geno)
                        whichSNP+=1
                        if whichSNP>=len(output):
                            break
                        if output[whichSNP][0] != whichChrom:
                            break


# In[45]:


if args['vcf']!='%#$' or args['vcfD']!='%#$':
    ps=[]
    for i in range(len(output)):
        gen=genotypes[i]
        if gen[0]==gen[-1]:
            continue
        x=output[i]
        if type(x[-1])!=str:
            ps.append(x[-1])
    
qVals=sm.stats.multipletests(ps,method='fdr_bh')[1]


# In[53]:


i=0
cols=['#Chr','Pos','Ref_Genotype','Most_Common_Alt','Ref Counts','Alt Counts','Ref_Allele_Ratio','Z','p_value','q_value']
if args['vcf']!='%#$' or args['vcfD']!='%#$':
    cols.append('genotype')

with open(args['o'],'w+') as f:
    f.write("\t".join(cols)+"\n")
    if args['vcf']!='%#$' or args['vcfD']!='%#$':
        for x,g in zip(output,genotypes):
            if g[0]==g[-1]:
                x.append('nan')
            elif type(x[-1])!=str:
                x.append(qVals[i])
                i+=1
            x.append(g)
            f.write("\t".join([str(y) for y in x])+"\n")
    else:
        for x in output:
            if x[7]=='nan':
                x.append('nan')
            else:
                x.append(qVals[i])
                i+=1
            f.write("\t".join([str(y) for y in x])+"\n")


