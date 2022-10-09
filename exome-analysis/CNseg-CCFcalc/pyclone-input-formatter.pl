#!/usr/bin/perl
# reformats Mutect2 variant calls into input files for pyclone-vi
use strict;
use warnings;

my %clinhash; # holds sample names
my $fh;
my $ofh;

# get sample names
open ($fh,"<sample-info-final.txt") || die "c";
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    $clinhash{$lines[0]} = $lines[4];
}
close $fh;

#for each sample
foreach my $k (sort keys %clinhash){

    my $ofh;
    open ($ofh,">pyclone-input/$k.pyclone.input.txt") || die "c";
    print $ofh "mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\ttumour_content\n";

    #save CN segments for this sample
    my %cnhash;
    my $fh2;
    my $infile = "all-segments/$k/$k"."_segments.txt";
    open ($fh2,"<$infile") || die "c";
    <$fh2>;
    my $counter = 0;
    while (<$fh2>){
        s/[\r\n]//g;
        my $l2 = $_;
        my @l2s = split("\t",$l2);
        my $c2 = $l2s[0];
        my $s2 = $l2s[1];
        my $e2 = $l2s[2];
        my $cn_a = $l2s[10];
        my $cn_b = $l2s[11];
        
        $cnhash{$counter} = [$c2,$s2,$e2,$cn_a,$cn_b];
        $counter++;
    }

    open ($fh,"<merged-all.vars-cts.tumors.txt") || die "a";
    <$fh>;
    while (<$fh>){
        s/[\r\n]//g;
        my $line = $_;
        my @lines = split("\t",$line);
        my $id = $lines[0];

        if ($id eq $k){ # if this variant belongs to the current sample
            my $sample = $lines[53];
            my $cellularity = $lines[54];
            my $refcts = $lines[51];
            my $altcts = $lines[52];
            
            my $chr = $lines[1];
            my $start = $lines[2];
            my $end = $lines[3];
            my $mid = ($start+$end)/2;
            my $aachg = $lines[10];

            # annotate with allele-specific CNs
            my $major_cn = "NA";
            my $minor_cn = "NA";

            foreach my $j (sort keys %cnhash){
                my @dd = @{$cnhash{$j}};
                if ($dd[0] eq $chr && ($dd[1] <= $mid && $dd[2] >= $mid)){ #variant is located within the segment
                    $major_cn = $dd[3];
                    $minor_cn = $dd[4];
                }
            }

            # if no CN calls at this region, drop the variant altogether
            if ($major_cn ne "NA" && $minor_cn ne "NA"){
                #mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\ttumour_content\n
                my $mutname = "$k:$sample:$chr:$start:$aachg";
                print $ofh "$mutname\t$k\t$refcts\t$altcts\t$major_cn\t$minor_cn\t2\t$cellularity\n";
            }

        }
    }
}


