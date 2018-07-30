
# coding: utf-8

# In[1]:


import argparse
import os


# In[ ]:


parser=argparse.ArgumentParser()
parser._optionals.title = "Flag Arguments"
parser.add_argument('-i',help="Input VCF or VCF folder.", required=True)
parser.add_argument('-o',help="Oupt SNP directory.", required=True)
parser.add_argument('-dir',help="Tells the program -i is specifying a VCF dir.", action='store_true')
parser.add_argument('-pref',help="VCF file name within folder with @ for chromsome number (ie some.name.chr@.vcf.", required=False)
parser.add_argument('-chr',help="Tells the program chromsomes are 'chr1' instead of '1'.", action='store_true')
args = vars(parser.parse_args())


# In[ ]:


if not os.path.exists(args['o']):
    os.mkdir(args['o'])


# In[ ]:


if args['i'][-3:]==".gz":
    import gzip
    openFile=gzip.open
else:
    openFile=open


# In[ ]:


if not args['dir']:
    for line in openFile(args['i']):
        if line[0]=="#":
            continue
        line=line.split()
        chrom=line[0] if args['chr'] else "chr"+line[0]
        with open("{0}/{1}.snps.txt".format(args['o'],chrom),'a+') as f:
            f.write("{0}\t{1}\t{2}\n".format(line[1],line[3],line[4]))
            
else:
    for c in range(1,23):
        chrom="chr"+str(c)
        with open("{0}/{1}.snps.txt".format(args['o'],chrom),'a+') as f:
            for line in openFile("{0}/{1}".format(args['i'],args['pref'].replace("@",c))):
                if line[0]=="#":
                    continue
                line=line.split()
                f.write("{0}\t{1}\t{2}\n".format(line[1],line[3],line[4]))
            


