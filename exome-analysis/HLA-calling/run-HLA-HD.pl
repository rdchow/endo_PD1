#!/usr/bin/perl
# Filter bam files to MHC region
use strict;
use warnings;

my $fh;

#Starts with bam files already aligned to hg19 (g1kv37)
open ($fh, "bam-list.txt") || die "c";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @info = split ("/",$line);
    my $sample = $info[1];
    system("samtools view -h -b $line 6:28477797-33448354 > ./mhc-bam/$sample.mhc.bam"); #extract reads from MHC region
    system("samtools view -b -f 4 $line > ./mhc-bam/$sample.unmapped.bam"); #get unmapped reads
    system("samtools merge -@ 20 ./mhc-bam/$sample.merge.bam ./mhc-bam/$sample.unmapped.bam ./mhc-bam/$sample.mhc.bam"); #combine MHC-mapped and unaligned reads

    system("samtools sort -n -@ 20 -o ./mhc-bam/$sample.merge.qsort.bam ./mhc-bam/$sample.merge.bam"); #sort bam files
    system("bedtools bamtofastq -i ./mhc-bam/$sample.merge.qsort.bam -fq ./mhc-fq/$sample.hlatmp.1.fastq -fq2 ./mhc-fq/$sample.hlatmp.2.fastq"); #convert bam files to fastq
    
    #reformat fastq headers to show read 1 vs 2 labels
    system("cat ./mhc-fq/$sample.hlatmp.1.fastq | awk '{if(NR%4 == 1){O=\$0;gsub(\"/1\",\" 1\",O);print O}else{print \$0}}' > ./mhc-fq/$sample.hla.1.fastq");
    system("cat ./mhc-fq/$sample.hlatmp.2.fastq |awk '{if(NR%4 == 1){O=\$0;gsub(\"/2\",\" 2\",O);print O}else{print \$0}}' > ./mhc-fq/$sample.hla.2.fastq");

    #run HLA-HD
    system("hlahd.sh -t 20 -m 100 -f ~/gibbs/hlahd.1.4.0/freq_data ./mhc-fq/$sample.hla.1.fastq ./mhc-fq/$sample.hla.2.fastq ~/gibbs/hlahd.1.4.0/HLA_gene.split.txt ~/gibbs/hlahd.1.4.0/dictionary $sample ~/gibbs/pem/hlahd-out");

}
close $fh;