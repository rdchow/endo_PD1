#!/usr/bin/perl

use strict;
use warnings;

my $fh;
open ($fh, "<hla-calls-list.txt") || die "c";
<$fh>;
while (<$fh>){
	s/[\r\n]//g;
	my $line = $_;
	my @lines = split ("\t",$line);
	my $sample = $lines[0];
	my $genes = $lines[1];

	if ($sample eq "S16-18558"){
	system("singularity exec -B /gpfs/gibbs/pi/chen_sidi/rdc55/pem:/home /gpfs/gibbs/pi/chen_sidi/rdc55/pvactools_latest.sif pvacseq run --iedb-install-directory /opt/iedb --n-threads 20 /home/vcf-annot/$sample.head.vep.vcf TUMOR $genes MHCflurry /home/pvac-out-PEM19-flurry");
	}
}


