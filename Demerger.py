#!/usr/bin/python3

# Use this simple parsing script to reformat the integrated Illumina callset VCF for ASE-SNP gene intersection with bedtools.
# Run on the command line as "$ python Demerger.py ALL_Illumina_Integrate_20170206.vcf"

import sys

infile = sys.argv[1]
svtype = infile.split('/')[-1].split('.')[1]
pref = '.'.join(infile.split('/')[-1].split('.')[0].split("_")[0:2])
outfy = pref+'.'+svtype+'_integrate.bed'

with open(infile) as SVs, open(outfy,'w') as outfile:
        for line in SVs.readlines():
                if line.startswith('#'):
                        continue
                else:
                        line = line.rstrip()
                        cols = line.split('\t')
                        start=int(cols[1])
                        end = cols[7].split(';')[1].split('=')[1]
                        cols.pop(2)
                        cols.insert(2,end)
                        trioinfo = []
                        info = cols[7].split(';')[7].split(',')
                        for genos in info:
                                geno = genos.split(':')
                                sampid = str(geno[4])
                                genetype = str(geno[3])
                                triotype = sampid + '\t' + genetype
                                SVdata=cols[0:7]
                                SVdata.append(triotype)
                                keep="\t".join(SVdata)
                                trioinfo.append(keep)
                        for tri in trioinfo:
                                outfile.write(tri+'\n')
