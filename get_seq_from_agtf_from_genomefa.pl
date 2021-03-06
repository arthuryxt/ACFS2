#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output basename\"  \"agtf\"  \"genome.fa\"  \"\(optional\) extend N bases\"  \"\(optional\)remove_seq_with_N==0\" " if (@ARGV < 3);
# store postion of exons and orders.
# output mRNA sequences together with exon sequences
# put extended sequences in lower case
# treat transcript without strandness (mostly single-exonic) as from "+"
my $fileout=$ARGV[0];
my $anno=$ARGV[1];
my $genome=$ARGV[2];
my $Extend=0;
if (scalar(@ARGV) > 3) {$Extend=$ARGV[3];}
my $removeN=0;
if (scalar(@ARGV) > 4) {$removeN=$ARGV[4];}

my %uniq;   # store all exons according to chr
my %Exons;  # store the number of exons of each transcript
open IN,$anno;
open OUT,">".$fileout.".exon.fa";
open OUT1,">".$fileout.".gene.fa";
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    #if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
    #if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
    if ($a[2] eq "exon") {
    #if ((length($a[0]) <= 3) and ($a[2] eq "exon")) {
        my @b=split(/\_\_\_/,$a[8]);
        my $len=scalar(@b);
        my $gene_id=$b[0];
        my $transcript_id=$b[0];
        my $gene_name=$a[7];    # exon_number is ALWAYS counted from left to right
        my $exon_number=$b[1];
        if ($a[6] eq "-") { $exon_number=$b[2]-$b[1]+1; }
        my $strand="+";
        if ($a[6] eq "-") {$strand="-";}
        #                                           chr  strand   start end                           biotype
        #                   0        1              2     3       4     5     6          7            8
        my $info=join("___",$gene_id,$transcript_id,$a[0],$strand,$a[3],$a[4],$gene_name,$exon_number,$a[1]);
        $uniq{$a[0]}{$gene_id."___".$transcript_id}{$exon_number}=$info;
        $Exons{join("___",$a[0],$gene_id,$transcript_id)}++;
    }
}
close IN;

my %CHRSEQ;
open(INseq, $genome) or die "Cannot open genome fasta file : $genome";
my $id="";
my $seq="";
while (<INseq>) {
    chomp;
    my $line=$_;
    if ($line=~m/^>/) {
        if ($seq ne "") { $CHRSEQ{$id}=$seq; $seq="";}
        my @a=split(" ",$line);
        $a[0]=~s/^>//;
        $id=$a[0];
    }
    else {
        if ($seq eq "") { $seq=$line; }
        else { $seq=$seq.$line; }
    }
}
if (($id ne "") and ($seq ne "")) { $CHRSEQ{$id}=$seq;}
close INseq;

foreach my $chr (sort keys %CHRSEQ) {
    if (!exists $uniq{$chr}) { next; }
    print "Processing chr: $chr\n";
    my $seq=$CHRSEQ{$chr};
        foreach my $gene (sort keys %{$uniq{$chr}}) {
        my $mRNA="";
        my $left=9999999999;
        my $right=-1;
        my $biotype="";
        foreach my $exon (sort{$a <=> $b} keys %{$uniq{$chr}{$gene}}) {
            my @a=split("___",$uniq{$chr}{$gene}{$exon});
            if ($left > $a[4]) {$left=$a[4]}
            if ($right < $a[5]) {$right=$a[5]}
            $biotype=$a[8];
            my $exon_seq="";
            if ($a[3] eq "+") {
                if (($exon == 1) and ($Exons{$chr."___".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body;
                }
                elsif ($exon == $Exons{$chr."___".$gene}) {
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$down;
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$exon_seq,"\n";
                if ($mRNA eq "") {$mRNA=$exon_seq;}
                else {$mRNA=$mRNA.$exon_seq;}
            }
            else {
                if (($exon == 1) and ($Exons{$chr."___".$gene} == 1)) {
                    my $up=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $down=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$up.$body.$down;
                }
                elsif ($exon == 1) {
                    my $up=lc substr($seq,$a[5],$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$body.$up;
                }
                elsif ($exon == $Exons{$chr."___".$gene}) {
                    my $down=lc substr($seq,($a[4]-1-$Extend),$Extend);
                    my $body=substr($seq,($a[4]-1),($a[5]-$a[4]+1));
                    $exon_seq=$down.$body;
                }
                else { $exon_seq=substr($seq,($a[4]-1),($a[5]-$a[4]+1)); }
                $exon_seq=~tr/[atcgATCG]/[tagcTAGC]/;
                my $tmp=reverse scalar $exon_seq;
                print OUT ">".$uniq{$chr}{$gene}{$exon},"\n",$tmp,"\n";
                if ($mRNA eq "") {$mRNA=$tmp;}
                else {$mRNA=$mRNA.$tmp;}
            }
        }
        if ($removeN eq 1) {
            if ($mRNA!~m/N/i) { print OUT1 ">".$gene."___".$biotype,"\n",$mRNA,"\n"; }
        }
        else { print OUT1 ">".$gene."___".$biotype,"\n",$mRNA,"\n"; }
    }
}
close OUT;
