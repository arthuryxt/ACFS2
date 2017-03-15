#!/usr/bin/perl -w
use strict;
# concatenate the first N-bases to the end of the sequence, thus creating a pseudo-circle
die "Usage: $0  \"input_fasta\"    \"output\"  \"\(optional\) extend_N_bases double the sequence by default\"  \"\(optional\)only_junction==0 or 1\"" if (@ARGV < 2);
my $filein=$ARGV[0];    
my $fileout=$ARGV[1];
my $extend=0;
if (scalar(@ARGV) > 2) { $extend=$ARGV[2]; }
if ($extend < 0) { die "extend_N_bases must be larger than 0 and optimal to be the sequencing length, or set to 0 to double the sequence"; }
my $OJ=0;
if (scalar(@ARGV) > 3) { $OJ=$ARGV[3]; }
if (($OJ ne 0) and ($OJ ne 1)) { die "only_junction must be either 0 (output full sequence with first N_bases appended) or 1 (only junction of length <= 2*entend_N_bases)"; }

open(IN, $filein) or die "Cannot open input_fasta : $filein";
open(OUT, ">".$fileout) or die "Cannot open output file : $fileout";
my $count=0;
while(<IN>){
    chomp;
    my $id=$_;
    my $seq=<IN>;
    chomp $seq;
    print OUT $id,"\n";
    if ($OJ eq 0) {
        if ($extend > 0) { print OUT $seq.substr($seq,0,$extend),"\n";}
        else {print OUT $seq.$seq,"\n";}
    }
    elsif ($OJ eq 1) {
        if ((length($seq) > $extend) and ($extend > 0)) {
            print OUT substr($seq,(0-$extend)).substr($seq,0,$extend),"\n";
        }
        else {
            print OUT $seq.$seq,"\n";
        }
    }
    $count++;
}
if ($count eq 0) {
    print "Wangning!!! input fasta file is empty! Please check the file : $filein \n";
    print "Wangning!!! dummy PolyG sequence added to output. \n";
    print OUT ">nonsense\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n";
}
