#!/usr/bin/perl
# call allele-specific copy number with Sequenza

use strict;
use warnings;
my $fh;
open ($fh,"<samples.txt") || die "c"; # list of tumor-normal pairs
while (<$fh>){
	my $line = $_;
	$line =~ s/[\r\n]//g;
	my @lines = split ("\t",$line);
	my $tf = $lines[0];
	my $nf = $lines[1];
	#Run Sequenza on aligned bams
	system("sequenza-utils bam2seqz --normal /home/rdc55/gibbs/pem/exome-data/$nf/$nf.bam --tumor /home/rdc55/gibbs/pem/exome-data/$tf/$tf.bam --fasta /home/bioinfo/software/knightlab/genomes/hs37d5/human_g1k_v37_decoy.fasta -gc genome_gc50.wig.gz --parallel 23 -C 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X --output ./seqz-out/$tf.seqz.gz");	

	my $samples = "";

	#Merge chr-split seqz files
	for (my $i = 1; $i <= 22; $i++){
		$samples .= "./seqz-out/$tf"."_$i.seqz.gz ";
	}
	$samples .= "./seqz-out/$tf"."_X.seqz.gz";
	my $command = "zcat $samples | gawk ".'\'{if (NR!=1 && $1 != "chromosome") {print $0}}\''." | bgzip > ./seqz-out/$tf.merge.seqz.gz";
	system($command);

	#index w/ tabix
	system("tabix -f -s 1 -b 2 -e 2 -S 1 ./seqz-out/$tf.merge.seqz.gz");

	#50nt bins
	system("sequenza-utils seqz_binning --seqz ./seqz-out/$tf.merge.seqz.gz -w 50 -o ./seqz-out/$tf.bin50.seqz.gz");
}
close $fh;
