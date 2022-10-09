#!/usr/bin/perl
use strict;
use warnings;
my $fh;

open ($fh,"<valero-ng-data.txt") || die "c";
my $head = <$fh>;
$head =~ s/[\r\n]//g;
my @headers = split ("\t",$head);
for (my $i = 0; $i < 17; $i++){
    print $headers[$i],"\t";
}
print "OS\tICB.binary\tMMRmut_canon\tMMRmut_any\tARID1Amut\n";

while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);

    
    my $mmr = "no"; # no MSI, canonical mut-MMRd, or other
    my $mutmmr = "no"; # no MSI, any mut-MMRd, or other
    my $aridmut = "no";
    if ($lines[14] =~ /ARID1A/){ $aridmut="yes"; }
    if ($lines[10] eq "Unstable"){
        if ($lines[14] =~ /MLH1/ || $lines[14] =~ /MSH2/ || $lines[14] =~ /MSH6/ || $lines[14] =~ /PMS2/){
            $mmr = "mut-MMRd";
        }
        else{
            $mmr = "other-MMRd";
        }

        if ($lines[14] ne ""){
            $mutmmr = "mut-MMRd";
        }
        else {
            $mutmmr = "other-MMRd";
        }
    } 

    for (my $i = 0; $i < 17; $i++){
            print $lines[$i],"\t";
    }

    my $os;
    if ($lines[19] == 0){ # if non-ICI treated, simply reprint OS
        print "$lines[17]\t$lines[19]\t$mmr\t$mutmmr\t$aridmut\n";
    }
    elsif ( ($lines[19] == 1)){ # if ICI-treated after diagnosis, subtract start->ICB and ICB->end
        my $diff = $lines[17] - $lines[18];
        print "$diff\t$lines[19]\t$mmr\t$mutmmr\t$aridmut\n"; #here, ICB treatment status is 0
    }
}
close $fh;
        