#!/usr/bin/perl -w
use strict;
# convert bed_formatted pre_defined_circRNAs into proprietray sum format
die "Usage: $0   \"bed_input\"  \(optional\)output" if (@ARGV < 1);
my $filein=$ARGV[0];
my $fileout="temp.pre_defined_circRNA";
if (scalar(@ARGV) > 1) {$fileout=$ARGV[1];}
open IN,$filein;
open OUT,">".$fileout.".agtf";
open OUTs,">".$fileout.".stat";
my $count=0;
while(<IN>) {
	chomp;
	if (m/^>/){ next; }
	if (m/^@/){ next; }
	if (m/^#/){ next; }
	if (m/^track/){ next; }
	my @a=split("\t",$_);
	if (scalar(@a) >= 12) {
		my @b=split(/\_/,$a[3]);
		my $id=join("_",$b[0],$b[1],$b[2],$b[3]);
		print OUTs $id,"\n";
		my @L=split(/\,/,$a[10]);
		my @R=split(/\,/,$a[11]);
		for(my $i=0; $i<$a[9]; $i++) {
			if ($a[5] eq "+") {
				print OUT join("\t",$a[0],"split","exon",($a[1]+1+$R[$i]),($a[1]+$L[$i]+$R[$i]),"biotype",$a[5],$b[7],join("___",$id,(1+$i),$a[9])),"\n";
			}
			else {
				print OUT join("\t",$a[0],"split","exon",($a[1]+1+$R[$i]),($a[1]+$L[$i]+$R[$i]),"biotype",$a[5],$b[7],join("___",$id,(1+$a[9]-$i),$a[9])),"\n";
			}
		}
	}
}
