

#preapre bedgraph
samtools view pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.bam | sort -k10,10 > pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.sorted.sam

zcat pheno_apo.txt_k31_bimbam_all.assoc.txt.assoc.txt.gz.LRT0-01.txt.gz | sort -k2,2 > pheno_apo.txt_k31_bimbam_all.assoc.txt.assoc.txt.gz.LRT0-01.sorted.txt

#join and extract chr, start, stop, ltr_pval
join -1 10 -2 2 pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.sorted.sam pheno_apo.txt_k31_bimbam_all.assoc.txt.assoc.txt.gz.LRT0-01.sorted.txt | awk '{print $4 "\t" $5 "\t" $5+31 "\t" $33}' | sort -k1,1 -k2,2n > pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.sorted.join.bedgraph
#make pval into -log10(pval), also adjust seq name in this case
cat pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.sorted.join.bedgraph | sed 's/|/\t/g' |  awk '{printf "%s\t%.0f\t%.0f\t%.1f\n", $2,$4,$5,(-1*log($6))/log(10)}' > pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.sorted.join.log10.bedgraph



#preapre windows
samtools faidx Urochloa_basilisk_9B9D_20251016.fa
cut -f1,2 Urochloa_basilisk_9B9D_20251016.fa.fai > Urochloa_basilisk_9B9D_20251016.fa.fai.genome
bedtools makewindows -g Urochloa_basilisk_9B9D_20251016.fa.fai.genome -w 100000 | sort -k1,1 -k2,2n > Urochloa_basilisk_9B9D_20251016.100kbwindows.txt
cat Urochloa_basilisk_9B9D_20251016.100kbwindows.txt | sed 's/|/\t/g' | awk '{print $2 "\t" $4 "\t" $5}' > Urochloa_basilisk_9B9D_20251016.100kbwindows.renamed.txt

	##
#take the table that merged mapped kmer postions and kmer log10_pval and makes a bedgraph (chr, st, stop, pval)
 bedtools intersect -wa -wb -a pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.bowtie.sorted.join.log10.bedgraph -b Urochloa_basilisk_9B9D_20251016.100kbwindows.renamed.txt | bedtools groupby -g 5,6,7,4 -c 4 -o count > pheno_apo.k31.LRT0-01.Urochloa_basilisk_9B9D_20251016.100kbwindows.txt








