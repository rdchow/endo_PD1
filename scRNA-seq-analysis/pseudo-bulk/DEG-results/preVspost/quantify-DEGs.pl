#!/usr/bin/perl

# NLR = epiR
# LR = mutR
# NR = NR

# quantify the # of DEGs (up / down) for each comparison, compile into table
use strict;
use warnings;
my $fh;

my @folders = ("LR","NLR","NR");

my %hash; #rows = cell types, columns = conditions.
# LR-up, LR-down, NLR-up, NLR-down, NR-up, NR-down

my $ofh;
open ($ofh,">preVspost.DEGcounts.txt") || die "c";

my $counter = 0;
foreach my $f (@folders){
    opendir my $dir, "./$f" || die "a";
    #my @files = readdir $dir;
    my @files = grep { /\.txt/ && -f "$f/$_" } readdir $dir;
    closedir $dir;

    foreach my $ff (@files){
        my @blah = split ("\Q.",$ff);
        my $cell = $blah[1];
        $cell =~ s/_deg//g;
        #print "$cell\n";

        open ($fh,"<./$f/$ff") || die "c";
        <$fh>;
        while (<$fh>){
            s/[\r\n]//g;
            my @lines = split ("\t",$_);
            my $padj =  $lines[5];
            my $lfc = $lines[2];

            if (!defined $hash{$cell}){
                $hash{$cell} = [0,0,0,0,0,0];
            }

            if ($padj < 0.05){
                if ($lfc > 0){
                    $hash{$cell}[$counter]++;
                }
                elsif ($lfc < 0){
                    $hash{$cell}[$counter+1]++;
                }
            }
        }
        close $fh;
    }
    $counter++;
    $counter++;
}

print $ofh "Celltype\tmutR_Post_up\tmutR_Post_down\tepiR_Post_up\tepiR_Post_down\tNR_Post_up\tNR_Post_down\n";

foreach my $k (sort keys %hash){
    my @dat = @{$hash{$k}};
    print $ofh "$k\t$dat[0]\t$dat[1]\t$dat[2]\t$dat[3]\t$dat[4]\t$dat[5]\n";
}

