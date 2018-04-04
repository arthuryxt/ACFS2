#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"List\"  \"Database_to_screen_from\"  \"Output\"  \(select==1_or_-1\)" if (@ARGV < 3);
my $filein=$ARGV[0];    # List of header (single column file or select the first-column in a multi-column file)
my $DB=$ARGV[1];        # fasta or fastq file to choose from
my $fileout=$ARGV[2];   # de-/selected file
my $select=1;
if (scalar(@ARGV) > 3) {$select=$ARGV[3];}
if (($select ne 1) and ($select ne -1)) {
    die "select can be either 1 (select) or -1 (unselect) \n";
}

open IN,$filein;
my %uniq;
while(<IN>) {
    chomp $_;
    my @a=split(/\t/,$_);
    $uniq{$a[0]}=1;
}
close IN;
open OUT1,">".$fileout;

open IN2,$DB;
my $firstline=<IN2>;
chomp $firstline;
close IN2;

open IN1,$DB;
if ($firstline=~m/^>/) {
    while(<IN1>) {
        chomp $_;
        next unless m/^>/;
        s/^>//;
        my @a=split("\t",$_);
        my $seq=<IN1>;
        if ($select eq 1) {
            if (exists $uniq{$a[0]}) {
                print OUT1 ">".join("\t",@a),"\n",$seq;
            }
        }
        else {
            if (!exists $uniq{$a[0]}) {
                print OUT1 ">".join("\t",@a),"\n",$seq;
            }
        }
    }
}
else {
    while(<IN1>) {
        chomp $_;
        next unless m/^@/;
        s/^@//;
        my @a=split("\t",$_);
        my $seq=<IN1>;
        my $qua=<IN1>;
        $qua=<IN1>;
        if ($select eq 1) {
            if (exists $uniq{$a[0]}) {
                print OUT1 "@".join("\t",@a),"\n",$seq,"+\n",$qua;
            }
        }
        else {
            if (!exists $uniq{$a[0]}) {
                print OUT1 "@".join("\t",@a),"\n",$seq,"+\n",$qua;
            }
        }
    }
}

