#!/usr/bin/perl -w
use strict;
die "Usage: $0    \"output_bed\"   \"input_refFlat\"  \"expr\"   \"\(optional\)track_name\"    \"\(optional\)erase_CDS yes_or_no\"  \"\(optional\)add_\"chr\" yes_or_no\" " if (@ARGV < 2);
my $fileout=$ARGV[0];
my $filein1=$ARGV[1];
my $exprfile=$ARGV[2];
my $track_name=$filein1;
if (scalar(@ARGV) > 3) { $track_name=$ARGV[3]; }
my $erase="yes";
if (scalar(@ARGV) > 4) { $erase=$ARGV[4]; }
my $addchr="yes";
if (scalar(@ARGV) > 5) { $addchr=$ARGV[5]; }

my %expr;
open(INe, $exprfile);
while (<INe>) {
    chomp;
    my @a=split("\t",$_);
    if (scalar(@a) < 3) { next;}
    if (($a[0] eq "newid") or ($a[1] eq "Gname")) { next;}
    my $sum=$a[2];
    for(my $i=3; $i<scalar(@a); $i++) { $sum+=$a[$i]};
    $expr{$a[0]}=$sum;
}
close INe;

open(IN, $filein1) || die "Can't open $filein1 for reading!\n";
open(OUT, ">".$fileout) || die "Can't open $fileout for writing!\n";
# GTF/GFF and Ensembl : 1-based start, 1-based end;
# UCSC data           : 0-based start, 1-based end;
# UCSC display        : 1-based start, 1-based end;
# format examples: 
# I. refFlat   half-open zero-based. This means that the first 100 bases of a chromosome are represented as [0,100), i.e. 0-99.
# XLOC_005566 TCONS_00012033  chr6  + 171030441 171044890 171030441 171030441 3 171030441,171037018,171044838,  171030573,171037146,171044890,
# XLOC_005566 TCONS_00011428  chr6  + 171044838 171045633 171044838 171044838 2 171044838,171045435,  171044973,171045633,
# XLOC_006319 TCONS_00014242  chr7  + 158743206 158750045 158743206 158743206 3 158743206,158749273,158749907,  158743609,158749374,158750045,
#
# II. Bed:     half-open 1-based.
# chr1  11873 14409 uc001aaa.3  0 + 11873 11873 0 3 354,109,1189, 0,739,1347,
# chr1  11873 14409 uc010nxr.1  0 + 11873 11873 0 3 354,52,1189,  0,772,1347,
# chr1  11873 14409 uc010nxq.1  0 + 12189 13639 0 3 354,127,1007, 0,721,1529,

print OUT "track name=\"".$track_name."\" description=\"".$track_name."\" visibility=2 itemRgb=\"On\"\n";
#my $score=0;
my $col=0;
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    my $name=$a[1];
    $name = $a[0].'|'.$a[1] unless($a[0] eq $a[1]);
    my $cnt=1;
    if (exists $expr{$a[1]}) { $cnt=$expr{$a[1]};}
    my $SScore=50;
    my $col5="";
    if (scalar(@a) > 12) {
        $SScore = $a[11];
        $col5=$a[12]."_".$a[13];
        #for(my $i=12; $i<=16; $i++){ $col5=$col5."_".$a[$i]; }
    }
    if ($cnt <=0 ) { $cnt=0.1; }
    my $SSsign=1;
    if ($SScore > 1) { $SSsign=1; }
    else { $SSsign=-1; $SScore=2-$SScore; }
    my $score=sprintf("%.3f",log($cnt) + $SSsign * log($SScore));
    
    my $lengths="";
    my $starts="";
    my @Start=split(/\,/,$a[9]);
    my @End=split(/\,/,$a[10]);
    #$a[4]--;
    #$a[5]--;
    for (my $i=0; $i<$a[8]; $i++) {
        my $cur_s = $Start[$i] - $a[4];
        my $cur_l = $End[$i] - $Start[$i];
        $lengths = $lengths.$cur_l.',';
        $starts = $starts.$cur_s.',';
    }
    if ($erase eq "yes") {
        $a[6]=$a[4];
        $a[7]=$a[4];
    }
    else {
        $a[6]=$a[4];
        $a[7]=$a[5];
    }
    if (($addchr eq "yes") and ($a[2]!~m/chr/)) {
        $a[2]="chr".$a[2];
    }
    if ($a[2] eq "chrMT") { $a[2]="chrM"; }
    if (scalar(@a) eq 12) {
        my @b=split(/\,/,$a[11]);
        print OUT join("\t",$a[2],$a[4],$a[5],$name."|".scalar(@b),scalar(@b),$a[3],$a[6],$a[7],$col,$a[8],$lengths,$starts),"\n";
    }
    else {
        print OUT join("\t",$a[2],$a[4],$a[5],$a[1]."_".$col5."_".$a[0],$score,$a[3],$a[6],$a[7],$col,$a[8],$lengths,$starts),"\n";
    }
}
