#!/usr/bin/perl
#grabs the VCF header from the OG vcf and adds onto the VEP-annotated files
use strict;
use warnings;

my $fh;
open ($fh, "vcf-files.txt") || die "c";
while (<$fh>){
	s/[\r\n]//g;
	my $line = $_;
	my @info = split ("_mut2",$line);
	my $sample = $info[0];

	system("vep --input_file ./vcf-annot/$sample.head.vcf --output_file ./vcf-annot/$sample.head.vep.vcf --fork 4 --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta /home/bioinfo/software/knightlab/genomes/hs37d5/human_g1k_v37_decoy.fasta --offline --cache --pick --plugin Frameshift --plugin Wildtype\n"); 

}


