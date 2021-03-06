#!/usr/bin/perl -w
use strict;
# convert bed_formatted pre_defined_circRNAs into proprietray sum format
die "Usage: $0   \"bed_input\"  \(optional\)output" if (@ARGV < 1);
my $filein=$ARGV[0];
my $fileout="temp.pre_defined_circRNA";
if (scalar(@ARGV) > 1) {$fileout=$ARGV[1];}
open IN,$filein;
open OUT,">".$fileout;
my $count=0;
while(<IN>) {
	chomp;
	if (m/^>/){ next; }
	if (m/^@/){ next; }
	if (m/^#/){ next; }
	if (m/^track/){ next; }
	my @a=split("\t",$_);
	$count++;
	my $id="tpid-".$count."/1__1";
	my $chr=$a[0];
	my $LS=50;
	my $RS=50;
	my $SS=100;
	my $strand="+";
	if ($a[5] eq "+") {
		$strand="-";
		if (scalar(@a) >= 12) {
			my @b=split(/\_/,$a[3]);
			if (scalar(@b) >= 6) {
				$LS=$b[5];
				$RS=$b[6];
				$SS=$LS+$RS;
			}
			print OUT join("\t",$id,60,$chr,30,30,($a[2]-30+1),$a[2],$strand,0,30,($a[1]+1),($a[1]+30),$strand,0,($a[1]-$a[2]+1),$SS,$LS,$RS,$a[5],0,0,0,10,0),"\n";
		}
		elsif (scalar(@a) >= 6 ) {
			print OUT join("\t",$id,60,$chr,30,30,($a[2]-30+1),$a[2],$strand,0,30,($a[1]+1),($a[1]+30),$strand,0,($a[1]-$a[2]+1),$SS,$LS,$RS,$a[5],0,0,0,10,0),"\n";
		}
	}
	else {
		# $a[5] eq "-"
		if (scalar(@a) >= 12) {
			my @b=split(/\_/,$a[3]);
			if (scalar(@b) >= 6 ) {
				$LS=$b[5];
				$RS=$b[6];
				$SS=$LS+$RS;
			}
			print OUT join("\t",$id,60,$chr,0,30,($a[2]-30+1),($a[2]),$strand,30,30,($a[1]+1),($a[1]+30),$strand,0,($a[1]-$a[2]+1),$SS,$LS,$RS,$a[5],0,0,0,10,0),"\n";
		}
		elsif (scalar(@a) >= 6 ) {
			print OUT join("\t",$id,60,$chr,0,30,($a[2]-30+1),($a[2]),$strand,30,30,($a[1]+1),($a[1]+30),$strand,0,($a[1]-$a[2]+1),$SS,$LS,$RS,$a[5],0,0,0,10,0),"\n";
		}
	}
	
}
