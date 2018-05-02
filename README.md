# HGSVC_ASE_Analysis

04/2018
Jia Wen, Conor Nodzak, Xinghua Mindy Shi 

If you have any question, please contact jwen6@uncc.edu, cnodzak@uncc.edu, X.Shi@uncc.edu

This repository includes results from SNP ASE analysis using strand-specific RNA-seq data and Whatshap strand-seq 10X phased SNPs nad SV ASE analysis as seen below:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160704_whatshap_strandseq_10X_phased_SNPs/

The index file of FASTQ files for strand-specific RNA-seq data of the 3 trios is located at:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/illumina_rna.sequence.index


We use a SNP ASE analysis pipeline based on mapping bias correction by WASP (Van de Geijn B et al. WASP: allele-specific software for robust discovery of molecular quantitative trait loci. bioRxiv (2014): 011221). The specific steps of this pipeline are as follows and results can be found in http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/201707_ASE_RES_Trios:

1. Map RNA-seq reads to the human reference GRCH38 using STAR v2.4.2a with default option. The reference we used here is from the link:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

```
   > command line: STAR --genomeDir ./genomedir --readFilesIn $indv_1.fastq $indv_2.fastq --runThreadN 4 --genomeLoad LoadAndKeep --   outFileNamePrefix $indv. --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0
```


2. The bam files output by STAR are corrected of mapping bias using WASP. The procedure involves remapping bam files that overlap with all SNPs, and discarding those reads that can’t be mapped to the same location with reference allele after flipped to alternative allele.

```
> command line: python $dir/WASP/mapping/find_intersecting_snps.py -p $indv.bam --snp_dir $dir/ --output_dir $dir/

> command line: STAR --genomeDir ./genomedir --readFilesIn $indv.remap.fq1 $indv.remap.fq2 --runThreadN 4 --genomeLoad LoadAndKeep —-outFileNamePrefix $indv_corrected. --outFilterScoreMin 0 --outFilterScoreMinOverLread 0 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0

> command line: python ./WASP/mapping/filter_remapped_reads.py $indv.to.remap.bam $indv_corrected.bam $indv_corrected.remap.keep.bam

> command line: samtools merge $indv.merged.bam $indv.keep.bam $indv_corrected.remap.keep.bam

```
3. Corrected bam files are sorted by Samtools and duplicate reads are then discarded using Picard.
```
> command line: java -jar picard.jar MarkDuplicates INPUT=$dir/$indv.merged.sorted.bam OUTPUT=$dir/$indv.Picard.bam METRICS_FILE=$indv.Picard.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
```

4. The reads with NM > 6 and MAPQ < 20 are further filtered to keep uniquely mapped reads using package Bamutils.

```
> command line: ./ngsutils-ngsutils-0.5.7/bin/bamutils filter $dir/$indv.Picard.bam $indv_filtered.bam -mapped -lte NM 6 -gte MAPQ 20
```

5. Use the Samtools mpileup to call variants from the filtered reads, and overlap them with heterozygous SNP (het-SNP) sites in each individual.

```
> command line: samtools mpileup -s -B -f ./GRCh38_full_analysis_set_plus_decoy_hla.fa $indv_filtered.bam -l $dir_snps/$indv/$indv_all_bisnp > $indv.pileup.txt
```

6. Filter het-SNP sites with total read count < 8 and only one allele seen at this het-SNP site. Only the reads with base quality > 10 are considered for counting reads.


7. Use binomial test with modified allelic ratio to test ASE for each filtered het-SNP site.
```
> command line: perl ./samase.pl -parse_pileup --sp $indv.pileup.txt --vcf $dir_snps/$indv/$indv_all_bisnp >  $indv.parsed_pileup.txt
> command line: perl ./samase.pl --calculate_ase_normal --pp $indv.parsed_pileup.txt -is $dir_snps/$indv/$indv_all_bisnp -r 0.5 -gf $indv.gene_snp_mapping_file.txt -mrs 8 -bq 10 > $indv_het.ase.txt
```
Perl script samase.pl used for generating reference and alternative read counts is from Kukurba et al. Allelic expression of deleterious protein-coding variants across human tissues. PLoS Genet 10.5 (2014): e1004304.

8. Finally, we conduct multi-test correction with Qvalue package in R.

9. Output significant ASE SNPs with qvalue < 0.05.

10. Bedtools is used to intersect genes with significant ASE SNPs (1bp overlap). The gene annotation file (Gencode v25) is downloaded from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz. 

```
> command line: bedtools intersect -a $indv_het.ase.txt -b gene_anno.bed -wa -wb > $indv.inter
```

Then intersect the heterozygous SV from PacBio callset for 3 children in trios with ASE genes. All the PacBio SV callset can be accessed in below link:
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20170109_UW_MSSM_Merged_PacBio/

We now begin a description of the SV-ASE analysis and results found in Table 5. SV-ASE.results.xlsx in the ftp: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/201707_ASE_RES_Trios.

The VCFs for Integrated Illumina SVs and Merged PacBio SVs were collected from the following: 
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/integration/20170206_Illumina_Integrate/
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20170109_UW_MSSM_Merged_PacBio/

The phased RNA-seq reads were gathered from: 
http:ftp://ftp-exchange.embl-heidelberg.de/pub/exchange/rausch/outgoing/haploRNA/

1. The Integrated Illumina VCF file was parsed to yield files of heterozygous SVs with a 'pass' value in the filter column, for each trio daughter.

```
> command line: grep 'PASS' ALL_Illumina_Integrate_20170206.vcf | grep -e '<DEL>' -e '<INS>' -  > PASS_Illumina_Integrate_20170206.DELINS.vcf
```
Demerger.py is a simple Python script to parse information from each sample's VCF INFO column into a "chrom start end" postion format output file. 
```
> command line: python Demerger.py PASS_Illumina_Integrate_20170206.DELINS.vcf > PASS.Illumina.DELINS_integrate.bed

> command line: sort -k1,1V -k2,2n  PASS.Illumina.DELINS_integrate.bed | uniq | awk '$9 == "0/1"' - | awk '{ $3 = ( $5 == "<INS>" ? $2+1 : $3) } 1' - | sed 's/\s/\t/g' - | grep -e "$sample" - > "$sample".PASS.Illumina.DELINS_heterozygous.bed
```

2. The heterozygous SVs from the Merged PacBio calls were extracted and formatted in a similar manner.
#
```
> command line: python PBSV.adjust.py 20170109_”$sample”.sv_calls.vcf  ## yields file called 20170109.”$sample".sv_calls_PBSV.bed
```
where PB.adjust.py is a Python script to parse information from each sample's VCF INFO column into a "chrom start end" postion format output file. 
```
> command line: grep -v '1|1' 20170109."$sample".sv_calls_PBSV.bed | awk '{ $3 = ( $5 == "<INS>" ? $2+1 : $3) } 1' - | sed 's/\s/\t/g' - > "$sample".hetsv_calls_PBSV.bed```		
```
3. Heterozygous SVs for each daughter were intersected with the SNP-ASE genes (genes identified by methods described above) to identify sets of SV-impacted genes to test for allele specific expression. The files “$sample".uniq.ASESNP.genes.txt contains the coordinates, gene IDs and gene names for ASE-SNP genes per sample identified by the SNP-ASE analysis above.
```
> command line: bedtools intersect -wa -wb -a “$sample".uniq.ASESNP.genes.txt  -b “$sample”.PASS.Illumina.DELINS_heterozygous.bed > $sample.intersect.hetILL.ASESNP.genes.txt

> command line: bedtools intersect -wa -wb -a "$sample".uniq.ASESNP.genes.txt -b "$sample".hetsv_calls_PBSV.bed > "$sample".intersect.hetPBSV.ASESNP.genes.txt
```
						 
4. Phased RNA-seq reads were then sorted with samtools and duplicates filtered with Picard. 
```
> command line: samtools sort "$sample".h1.bam > "$sample".sort.h1.bam  samtools sort "$sample".h2.bam > "$sample".sort.h2.bam

> command line: java -jar picard.jar MarkDuplicates INPUT="$sample".sort.h1.bam OUTPUT="$sample".filt.h1.bam M="$file".marked_dup_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
```

5. Phased RNA-seq reads were counted across coordinate regions of SNP-ASE genes that were intersected by heterozygous SVs using BEDtools multicov. This produces a count for each bam file, where “$sample”.filt.head.h2.bam and “$sample”.filt.head.h2.bam are the bam files for filtered RNA-seq reads phased to haplotype 1 or haplotype 2 for each daughter sample. 
```
> command line: bedtools multicov -bams “$sample”.filt.head.h1.bam “$sample”.filt.head.h2.bam  -bed “$sample”.intersect.hetILL.ASESNP.genes.txt > “$sample”.illumina.multicov.txt

> command line: bedtools multicov -bams “$sample”.filt.head.h1.bam “$sample”.filt.head.h2.bam  -bed $sample.intersect.hetPBSV.ASESNP.genes.txt > “$sample”.PBSV.multicov.txt
```

6. A binomial test was then applied using binom.test() in R. 
```
> R: binom.test(x, n, p = 0.5)
```
Where x = single RNA-seq haplotype counts, n = total RNA-seq counts for both haplotypes, p = hypothesized probability of success.

7. Multi-test correction was applied with p.adjust() in R.

8. Output significant SVs exhibit an ASE affect with FDR adjusted p-value <= 0.05.
