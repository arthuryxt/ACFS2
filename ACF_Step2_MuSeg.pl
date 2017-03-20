#!/usr/bin/perl -w
use strict;
# rescue circ spanning very short exons
die "Usage: $0   \"multiple_segments\"   \"CB_splice_folder\"   \"genome_fasta\"    \"\(optional\) extend N bases\"   \"\(optional\)split_exon_gtf\"  \"\(optional\) fuzzy exon border recognition\"     \"\(DEBUG\)\"" if (@ARGV < 3);
my $fileout=$ARGV[0];
my $DIR=$ARGV[1];		# /home/arthur/CB_splice/
my $genome=$ARGV[2];    # /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/WholeGenomeFasta/genome.fa
my $Extend=15;           # 15nt by default. as 3' splice strength need 23nt, 20nt from intron and 3 from exon.
if (scalar(@ARGV) > 3) {$Extend=$ARGV[3];}
my $agtf="no";
if (scalar(@ARGV) > 4) {$agtf=$ARGV[4];}
my $Ext=10;
if (scalar(@ARGV) > 5) {$Ext=$ARGV[5];}
print "Fuzzy exon border recognition, extention = $Ext\n";
my $debug=-1;
if (scalar(@ARGV) > 6) {$debug=$ARGV[6];}
my $command="rm -f Step2_MuSeg_finished";
system($command);

my %uniqanno;
if ($agtf ne "no") {
    open IN1,$agtf;
    while(<IN1>) {
        chomp;
        my @a=split("\t",$_);
        #if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
        #if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
        my $id=$a[7]."___".$a[8];
        if ($a[6] eq "+") {
            # 5' of exon
            my $base=int($a[3]/1000);
            if (exists $uniqanno{$a[0]}{$base}{$a[3]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base}{$a[3]});
                if ($b[0] < 10) { $b[0]+=10; $uniqanno{$a[0]}{$base}{$a[3]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base}{$a[3]}=join("___",10,$id); }
            
            if (exists $uniqanno{$a[0]}{$base+1}{$a[3]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base+1}{$a[3]});
                if ($b[0] < 10) { $b[0]+=10; $uniqanno{$a[0]}{$base+1}{$a[3]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base+1}{$a[3]}=join("___",10,$id); }
            
            if (exists $uniqanno{$a[0]}{$base-1}{$a[3]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base-1}{$a[3]});
                if ($b[0] < 10) { $b[0]+=10; $uniqanno{$a[0]}{$base-1}{$a[3]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base-1}{$a[3]}=join("___",10,$id); }
            # 3' of exon
            $base=int($a[4]/1000);
            if (exists $uniqanno{$a[0]}{$base}{$a[4]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base}{$a[4]});
                if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniqanno{$a[0]}{$base}{$a[4]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base}{$a[4]}=join("___",1,$id); }
            
            if (exists $uniqanno{$a[0]}{$base+1}{$a[4]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base+1}{$a[4]});
                if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniqanno{$a[0]}{$base+1}{$a[4]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base+1}{$a[4]}=join("___",1,$id); }
            
            if (exists $uniqanno{$a[0]}{$base-1}{$a[4]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base-1}{$a[4]});
                if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniqanno{$a[0]}{$base-1}{$a[4]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base-1}{$a[4]}=join("___",1,$id); }
        }
        else {
            # 5' of exon
            my $base=int($a[4] / 1000);
            if (exists $uniqanno{$a[0]}{$base}{$a[4]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base}{$a[4]});
                if ($b[0] < 10) { $b[0]+=10; $uniqanno{$a[0]}{$base}{$a[4]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base}{$a[4]}=join("___",10,$id); }
            
            if (exists $uniqanno{$a[0]}{$base+1}{$a[4]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base+1}{$a[4]});
                if ($b[0] < 10) { $b[0]+=10; $uniqanno{$a[0]}{$base+1}{$a[4]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base+1}{$a[4]}=join("___",10,$id); }
            
            if (exists $uniqanno{$a[0]}{$base-1}{$a[4]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base-1}{$a[4]});
                if ($b[0] < 10) { $b[0]+=10; $uniqanno{$a[0]}{$base-1}{$a[4]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base-1}{$a[4]}=join("___",10,$id); }
            # 3' of exon
            $base=int($a[3] / 1000);
            if (exists $uniqanno{$a[0]}{$base}{$a[3]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base}{$a[3]});
                if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniqanno{$a[0]}{$base}{$a[3]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base}{$a[3]}=join("___",1,$id); }
            
            if (exists $uniqanno{$a[0]}{$base+1}{$a[3]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base+1}{$a[3]});
                if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniqanno{$a[0]}{$base+1}{$a[3]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base+1}{$a[3]}=join("___",1,$id); }
            
            if (exists $uniqanno{$a[0]}{$base-1}{$a[3]}) {
                my @b=split("___",$uniqanno{$a[0]}{$base-1}{$a[3]});
                if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniqanno{$a[0]}{$base-1}{$a[3]}=join("___",@b)}
            }
            else { $uniqanno{$a[0]}{$base-1}{$a[3]}=join("___",1,$id); }
        }
    }
    close IN1;
}

open IN2,$fileout;
open OUT1,">".$fileout.".S1.pp";
open OUT11,">".$fileout.".S1.ppparsed";
open OUT2,">".$fileout.".S1.pm";
open OUT22,">".$fileout.".S1.pmparsed";
while(<IN2>) {
    chomp;
    my @a=split("\t",$_);
    #if ($a[4]=~m/^chromosome/i) {$a[4]=~s/chromosome//i;}
    #if ($a[4]=~m/^chr/i) {$a[4]=~s/chr//i;}
    my $NR=int(scalar(@a)/5);
    my %strand;
	my $gap=-1;
    for(my $i=1; $i<$NR; $i++) {
        $a[5*$i+3]--;
		$strand{$a[5*$i + 4]}++;
		if ($debug eq 1) { print join("\t",$a[0],$a[4],$a[5*$i],$a[5*$i+1],$a[5*$i+2],$a[5*$i+3],$a[5*$i+4]),"\n"; }
		my $base=int(($a[5*$i + 2])/1000);
		if ($debug eq 1) { print $base,"\n"; }
		my $info="NA";
		if (exists $uniqanno{$a[4]}) {
		    if (exists $uniqanno{$a[4]}{$base}) {
				foreach my $pos (keys %{$uniqanno{$a[4]}{$base}}) {
				    if ($debug eq 1) {print join("\t",$a[0],$a[4],$a[5*$i],$a[5*$i+1],$a[5*$i+2],$a[5*$i+3],$a[5*$i+4],$pos,$uniqanno{$a[4]}{$base}{$pos}),"\n";}
				    if (abs($pos - $a[5*$i + 2]) <= $Ext) { $info=$uniqanno{$a[4]}{$base}{$pos}; last; }
				}
			}
		}
		$a[5*$i+2]=$a[5*$i + 2]."\t".$info;
	
		$base=int(($a[5*$i + 3])/1000);
		if ($debug eq 1) { print $base,"\n"; }
		$info="NA";
		if (exists $uniqanno{$a[4]}) {
		    if (exists $uniqanno{$a[4]}{$base}) {
				foreach my $pos (keys %{$uniqanno{$a[4]}{$base}}) {
				    if ($debug eq 1) {print join("\t",$a[0],$a[4],$a[5*$i],$a[5*$i+1],$a[5*$i+2],$a[5*$i+3],$a[5*$i+4],$pos,$uniqanno{$a[4]}{$base}{$pos}),"\n";}
				    if (abs($pos - $a[5*$i + 3]) <= $Ext) { $info=$uniqanno{$a[4]}{$base}{$pos}; last; }
				}
		    }
		}
		$a[5*$i+3]=$a[5*$i + 3]."\t".$info;
    }
    my $ss=0;
    foreach my $i (keys %strand) { $ss++; }
    if ($ss eq 1) {
		print OUT1 join("\t",@a),"\n";
		
		my @Outer_left=qw(999999999 999999999 NA);
		my @Outer_right=qw(-1 -1 NA);	# where the Back Splice took place
		my @Inner_left=qw(-1 -1 NA);
		my @Inner_right=qw(-1 -1 NA);	# should be the first and the last segments, but NO (left < right) relation assumed!!!
		my %POS;
        my $leftmostrank=-1;
        my $rightmostrank=-1;
		for(my $i=1; $i<$NR; $i++) {
		    my @tmpa=split("\t",$a[5*$i+2]);
		    my @tmpb=split("\t",$a[5*$i+3]);
		    $POS{$tmpa[0]}{$i}=join("____",$tmpa[1],$tmpb[0],$tmpb[1],$a[5*$i+4],$a[5*$i],$a[5*$i+1]);
		    if ($Inner_left[0] eq -1) { $Inner_left[0]=$tmpa[0]; $Inner_left[1]=$tmpb[0]; $Inner_left[2]=$tmpa[1];}
		    $Inner_right[0]=$tmpa[0]; $Inner_right[1]=$tmpb[0]; $Inner_right[2]=$tmpb[1];
		    if ($Outer_left[0] > $tmpa[0]) {$Outer_left[0] = $tmpa[0]; $Outer_left[1]=$tmpb[0]; $Outer_left[2]=$tmpa[1]; $leftmostrank=$i;}
		    if ($Outer_right[1] < $tmpb[0]) {$Outer_right[0] = $tmpa[0]; $Outer_right[1]=$tmpb[0]; $Outer_right[2]=$tmpb[1]; $rightmostrank=$i;}
            if ($debug eq 20) {
                print join("\t",@Outer_left),"\t\t",join("\t",@Outer_right),"\t\t",join("\t",@Inner_left),"\t\t",join("\t",@Inner_right),"\n";
            }
            
		}
	
		# check if Inner border is truely inner
		my $flag=0;	# 0:circle;	1:polymerase backslide or templateswitch;	10-or-100: backslide	110:linear
		#if ( (($Inner_right[0] <= $Inner_left[1]) and ($Inner_left[1] <= $Inner_right[1]))  or  (($Inner_left[0] <= $Inner_right[1]) and ($Inner_right[1] <= $Inner_left[1])) ) {$flag++;}
		#for(my $i=2; ($i<($NR-1))and($flag eq 0); $i++) {
		#    my @tmpa=split("\t",$a[5*$i+2]);
		#    my @tmpb=split("\t",$a[5*$i+3]);
		#    if ( $Inner_left[1]  <  $Inner_right[0]) {
		#		if (($Inner_left[1] <= $tmpa[0]) and ($tmpa[0] <= $Inner_right[0])) {$flag+=10;}
		#		if (($Inner_left[1] <= $tmpb[0]) and ($tmpb[0] <= $Inner_right[0])) {$flag+=100;}
		#    }
		#    else {
		#		if (($Inner_right[1] <= $tmpa[0]) and ($tmpa[0] <= $Inner_left[0])) {$flag+=10;}
		#		if (($Inner_right[1] <= $tmpb[0]) and ($tmpb[0] <= $Inner_left[0])) {$flag+=100;}
		#    }
		#}
		# check for circularity
		#if ($flag eq 0) {
		#	my $last=-1;
        #    my $breakpoint=0;
		#	foreach my $left (sort {$a <=> $b} keys %POS) {
		#		foreach my $rank (sort{$a <=> $b} keys %{$POS{$left}}) {
		#			if ($last eq -1) {}
		#			elsif (($last - $rank) > 1) {
		#				if ($a[-1] eq "+") {$breakpoint++;}
		#				else {$flag++;}
		#			}
		#			elsif (($rank - $last) > 1) {
		#				if ($a[-1] eq "-") {$breakpoint++;}
		#				else {$flag++;}
		#			}
		#			$last=$rank;
		#		}
		#	}
            #if($breakpoint ne 1){$flag=1000+$breakpoint;}
		#}
		# report the segmental order with sticky-end info
		my $info=join("\t",$a[0],$a[1],$a[2],$a[3],$a[4],$a[-1],$flag,$Outer_left[0],$Outer_right[1]);
		if ($debug eq 3) { $info=join("\t",$a[0],$a[1],$a[2],$a[3],$a[4],$a[-1],$flag,$Outer_left[0],$Outer_right[1],$Inner_left[0],$Inner_left[1],$Inner_right[0],$Inner_right[1]) }
		foreach my $left (sort {$a <=> $b} keys %POS) {
		    foreach my $rank (sort{$a <=> $b} keys %{$POS{$left}}) {
                if (($rank eq $rightmostrank) or ($rank eq $leftmostrank)) {
                    $info=$info."\t".$rank."____".$left."____".$POS{$left}{$rank};
                }
				#$info=$info."\t".$rank."____".$left."____".$POS{$left}{$rank};
		    }
		}
		# report
		print OUT11 $info,"\t",$a[-1],"\n";
    }
    else {
		print OUT2 join("\t",@a),"\n";
		
		my $Plus_left=999999999;
		my $Plus_right=-1;
		my $Minus_left=999999999;
		my $Minus_right=-1;
		my %POS;
		for(my $i=1; $i<$NR; $i++) {
		    my @tmpa=split("\t",$a[5*$i+2]);
		    my @tmpb=split("\t",$a[5*$i+3]);
		    $POS{$tmpa[0]}{$i}=join("____",$tmpa[1],$tmpb[0],$tmpb[1],$a[5*$i+4]);
		    if ($a[5*$i+4] eq "+") {
				if ($Plus_left >= $tmpa[0]) {$Plus_left = $tmpa[0];}
				if ($Plus_right <= $tmpb[0]) {$Plus_right = $tmpb[0];}
		    }
		    else {
				if ($Minus_left >= $tmpa[0]) {$Minus_left = $tmpa[0];}
				if ($Minus_right <= $tmpb[0]) {$Minus_right = $tmpb[0];}
		    }
		}
    	# Not necessary to check if Inner border is truely inner
		my $flag=1;	# 0:circle;	1:polymerase backslide or templateswitch;	10-or-100: backslide	110:linear
		# report the segmental order with sticky-end info
		my $info=join("\t",$a[0],$a[1],$a[2],$a[3],$a[4],$Plus_left,$Plus_right,$Minus_left,$Minus_right);
		foreach my $left (sort {$a <=> $b} keys %POS) {
		    foreach my $rank (sort{$a <=> $b} keys %{$POS{$left}}) {
				$info=$info."\t".$rank."____".$left."____".$POS{$left}{$rank};
		    }
		}
		# report
		print OUT22 $info,"\t",$a[-1],"\n";
    }
}
close IN2;
close OUT1;		# ".S1.pp"
close OUT11;	# ".S1.ppparsed"
close OUT2;		# ".S1.pm"
close OUT22;	# ".S1.pmparsed"

open INf, $fileout.".S1.ppparsed";
open(OUTf, ">".$fileout.".S1");

my %candy;
while(<INf>) {
	chomp;
	next if (m/^#/);
	my @a=split("\t",$_);
    next if (scalar(@a) eq 11);
    # newid-2161237__1	150	150	150	10	+	0	100012109	100019190	3____100012109____NA____100012131____NA____+____127____23	2____100019147____NA____100019190____NA____+____88____44
    my $strand=1;
    if ($a[5] eq "-") { $strand=-1;}
    my @tmpa=split(/\_\_\_\_/,$a[9]);
    my @tmpb=split(/\_\_\_\_/,$a[10]);
    my $ExonL=$tmpa[1].",".$tmpb[1];
    my $ExonR=$tmpa[3].",".$tmpb[3];
    my $min_gap=0;
    if ($a[5] eq "-") {
        # R-m+
        $min_gap=$tmpa[6]+$tmpa[7]-$tmpb[6];
    }
    else{
        # R+m-
        $min_gap=$tmpb[6]+$tmpb[7]-$tmpa[6];
    }
    #                                            0     1       2   3        4        5        6          7        8        9        10       11         12       13     14
    $candy{$a[4]}{$a[7]}{$a[8]}{$a[0]}=join("\t",$a[0],$a[3],$a[4],$tmpb[6],$tmpb[7],$tmpb[1],$tmpb[3]+1,$tmpb[5],$tmpa[6],$tmpa[7],$tmpa[1],$tmpa[3]+1,$tmpa[5],$a[11],$tmpa[1]-$tmpb[3]);
    print OUTf $candy{$a[4]}{$a[7]}{$a[8]}{$a[0]},"\n";
}
close INf;
close OUTf;

my %uniq;
my $cnt1=0;
open IN,$fileout.".S1";
open OUT,">".$fileout.".S2";
open OUT1,">".$fileout.".S2.0";	# splice-motif found				trust hit  CB prediction
open OUT2,">".$fileout.".S2.1";	# closure == 0,						trust best CB prediction
open OUTerr,">".$fileout.".S2_withN";
open OUTp,">".$fileout.".S2.pgap";
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[14] < 0) {
        #if ($a[2]=~m/^chromosome/i) {$a[2]=~s/chromosome//i;}
        #if ($a[2]=~m/^chr/i) {$a[2]=~s/chr//i;}
		$a[6]=$a[6]-1;
		$a[11]=$a[11]-1;
        $uniq{$a[2]}{$a[6]}{$a[10]}{$a[0]}=join("\t",@a);
        $cnt1++;
    }
    else {print OUTp join("\t",@a),"\n";}
}
close IN;

#foreach my $chr (sort keys %uniq) {
#    foreach my $start (sort keys %{$uniq{$chr}}) {
#        foreach my $end (sort keys %{$uniq{$chr}{$start}}) {
#            foreach my $seqid (keys %{$uniq{$chr}{$start}{$end}}) {
#                print $uniq{$chr}{$start}{$end}{$seqid},"\n";
#            }
#        }
#    }
#}

my $tmpfile1=$DIR."/me2x5";
my %me2x5 = &makescorematrix($tmpfile1);
my $tmpfile2=$DIR."/splicemodels/splice5sequences";
my %seq = &makesequencematrix($tmpfile2);
my %bgd;
$bgd{'A'} = 0.27;
$bgd{'C'} = 0.23;
$bgd{'G'} = 0.23;
$bgd{'T'} = 0.27; 
sub makesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub makescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub getrest5{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}

sub scoreconsensus5{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd; 
    $bgd{'A'} = 0.27; $bgd{'C'} = 0.23; $bgd{'G'} = 0.23; $bgd{'T'} = 0.27;  
    my %cons1;
    $cons1{'A'} = 0.004; $cons1{'C'} = 0.0032; $cons1{'G'} = 0.9896; $cons1{'T'} = 0.0032;
    my %cons2;
    $cons2{'A'} = 0.0034; $cons2{'C'} = 0.0039; $cons2{'G'} = 0.0042; $cons2{'T'} = 0.9884;
    my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
    return $addscore;
}

sub log2{
    my ($val) = @_;
    return log($val)/log(2);
}
sub hashseq{
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/ACGT/0123/;
    my @seqa = split(//,$seq);
    my $sum = 0;
    my $len = length($seq);
    my @four = (1,4,16,64,256,1024,4096,16384);
    my $i=0;
    while ($i<$len) {
        $sum+= $seqa[$i] * $four[$len - $i -1] ;
	$i++;
    }
    return $sum;
}

my @metables = &makemaxentscores;
sub makemaxentscores{
	my $dir = $DIR."/splicemodels/";
    my @list = ('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	$num++;
    }
    return @metables;
}
sub makewmmscores{
	my $dir = $DIR."/splicemodels/";
    my @list = ('me1s0acc1','me1s0acc2','me1s0acc3','me1s0acc4','me1s0acc5','me1s0acc6','me1s0acc7','me1s0acc8','me1s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	$num++;
    }
    return @metables;
}
sub makemmscores{
	my $dir = $DIR."/splicemodels/";
    my @list = ('me2s0acc1','me2s0acc2','me2s0acc3','me2s0acc4','me2s0acc5','me2s0acc6','me2s0acc7','me2s0acc8','me2s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	$num++;
    }
    return @metables;
}
sub maxentscore{
    my $seq = shift;
    my $table_ref = shift;
    my @metables = @$table_ref;
    my @sc;
    $sc[0] = $metables[0]{&hashseq(substr($seq,0,7))};
    $sc[1] = $metables[1]{&hashseq(substr($seq,7,7))};
    $sc[2] = $metables[2]{&hashseq(substr($seq,14,7))};
    $sc[3] = $metables[3]{&hashseq(substr($seq,4,7))};
    $sc[4] = $metables[4]{&hashseq(substr($seq,11,7))};
    $sc[5] = $metables[5]{&hashseq(substr($seq,4,3))};
    $sc[6] = $metables[6]{&hashseq(substr($seq,7,4))};
    $sc[7] = $metables[7]{&hashseq(substr($seq,11,3))};
    $sc[8] = $metables[8]{&hashseq(substr($seq,14,4))};
    my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] / ($sc[5] * $sc[6] * $sc[7] * $sc[8]);
    return $finalscore;
}    
    
sub getrest3{
    my $seq = shift;
    my $seq_noconsensus = substr($seq,0,18).substr($seq,20,3);
    return $seq_noconsensus;
}

sub scoreconsensus3{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd; 
    $bgd{'A'} = 0.27; $bgd{'C'} = 0.23; $bgd{'G'} = 0.23; $bgd{'T'} = 0.27;  
    my %cons1;
    $cons1{'A'} = 0.9903; $cons1{'C'} = 0.0032; $cons1{'G'} = 0.0034; $cons1{'T'} = 0.0030;
    my %cons2;
    $cons2{'A'} = 0.0027; $cons2{'C'} = 0.0037; $cons2{'G'} = 0.9905; $cons2{'T'} = 0.0030;
    my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]}); 
    return $addscore;
}

sub sanity {
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my $Nr=scalar(@seqa);
    my $flag=0;
    for(my $i=0; $i<$Nr; $i++) {
        if ($seqa[$i] !~m/[ATCG]/) {$flag++;}
    }
    return $flag;
}

my %SMotif;
$SMotif{"GTAG"}=3;
$SMotif{"GCAG"}=2;
$SMotif{"ATAC"}=1;
$SMotif{"CTAC"}=-3;
$SMotif{"CTGC"}=-2;
$SMotif{"GTAT"}=-1;

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
    my $SEQ=$CHRSEQ{$chr};
    foreach my $start (sort keys %{$uniq{$chr}}) {
        foreach my $end (sort keys %{$uniq{$chr}{$start}}) {
            foreach my $seqid (keys %{$uniq{$chr}{$start}{$end}}) {
                my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
				if  ($a[7] eq "+") {  EXIT_IF:{
					# R+m-
					my $overlap=$a[3]+$a[4]-$a[8];
					my $left_seq=substr($SEQ,($start-2*$Extend-1),(4*$Extend + 1));
					if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,"left_404"),"\n"; last EXIT_IF;}
					my $left_max=-99999;
					my $left_pos=0;
					my $left_flag=0;
					my $right_seq=substr($SEQ,($end-2*$Extend-1),(4*$Extend + 1));
					if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; last EXIT_IF; }
					if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,$end,$seqid,"right_404"),"\n";last EXIT_IF;}
					my $right_max=-99999;
					my $right_pos=0;
					my $right_flag=0;
					if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; last EXIT_IF;}
					if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; last EXIT_IF;}
					my $SMS=0;	
					my $PMS=0;
					my $SMscoren=-99999;
                    if ($debug eq 1){ print "\n"; }
                    if ($debug eq 1){print join("\t",@a),"\n",$left_seq,"\t",$right_seq,"\n",$overlap,"\n";}
					for(my $k=0; $k<=$overlap; $k++) {
						my $ml="";
						my $mr="";
						$ml=substr($left_seq,2*$Extend+1-$k,2);
						$mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
						if (($ml ne "") and ($mr ne "")) {
							my $motif=$ml.$mr;
                            if ($debug eq 1){print $k,"\t",$motif,"\n";}
							if ((exists $SMotif{$motif}) and ($SMotif{$motif} < 0)) {
								my $t=$left_seq;					
								$t=~tr/[atcgATCG]/[TAGCTAGC]/;
								my $rc_left_seq=scalar reverse $t;		
								$t=$right_seq;					
								$t=~tr/[atcgATCG]/[TAGCTAGC]/;
								my $rc_right_seq=scalar reverse $t;		
								my $str1="";
								$str1=uc substr($rc_left_seq,(2*$Extend+$k-20),23);
								my $str2="";
								$str2=uc substr($rc_right_seq,(2*$Extend+1-3-$overlap+$k),9);
								if ((length($str1) eq 23) and (length($str2) eq 9)) {
									my $left_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
									my $right_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
									my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
									if ($SMscoren < $sumt) { $SMscoren=$sumt; $SMS=$SMotif{$motif}; $PMS=$k; }
								}
							}
						}
					}
                    if ($debug eq 1){print "SMS = ",$SMS,"\tPSM = ",$PMS,"\tSMscoren = $SMscoren \n";}
					if ($SMS < 0) {
						my $t=$left_seq;					
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_left_seq=scalar reverse $t;		
						$t=$right_seq;					
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_right_seq=scalar reverse $t;		
						my $str1="";
						$str1=uc substr($rc_left_seq,(2*$Extend+$PMS-20),23);
						my $str2="";
						$str2=uc substr($rc_right_seq,(2*$Extend+1-3-$overlap+$PMS),9);
						if ((length($str1) eq 23) and (length($str2) eq 9)) {
							$left_max=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
							$right_max=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
							my $sum=sprintf("%.2f",($left_max + $right_max));
							$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."-5";
							$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."-3";
							$a[14]=$a[14]."\t".$sum."\t".$SMS."\t".$PMS;
							print OUT join("\t",@a),"\n";
							my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
							$a[4]=$a[4]-$PMS;
							$a[6]=$a[6]-$PMS;
							$a[8]=$a[8]+($overlap-$PMS);
							$a[9]=$a[9]-($overlap-$PMS);
							$a[10]=$a[10]+($overlap-$PMS);
							$a[14]=$a[14]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."-"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
							print OUT1 join("\t",@a),"\n";
							last EXIT_IF;
						}
					}
                    elsif (($overlap < 0) or ($a[13] > 0)) {
                        my $Overlap=$Ext;
                        $left_max=-99999;
                        $left_pos=0;
                        $right_max=-99999;
                        $right_pos=0;
                        my @SScoreL;
                        my @SScoreR;
                        for(my $k=0; $k<=$Overlap; $k++) {
                            my $t=$left_seq;					
                            $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                            my $rc_left_seq=scalar reverse $t;		
                            $t=$right_seq;					
                            $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                            my $rc_right_seq=scalar reverse $t;		
                            my $str1="";
                            $str1=uc substr($rc_left_seq,(2*$Extend+$k-20),23);
                            my $str2="";
                            $str2=uc substr($rc_right_seq,(2*$Extend+1-3-($Overlap-$k)),9);
                            if ((length($str1) eq 23) and (length($str2) eq 9)) {
                                $SScoreL[$k]=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
								$SScoreR[$Overlap-$k]=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
                                if ($debug > 0) { print join("\t",$k,$str1,$SScoreL[$k],$str2,$SScoreR[$Overlap-$k]),"\n"; }
                                if ($left_max < $SScoreL[$k]) {
                                    $left_max = $SScoreL[$k];
                                    $left_pos = $k;
                                }
                                if ($right_max < $SScoreR[$Overlap-$k]) {
                                    $right_max = $SScoreR[$Overlap-$k];
                                    $right_pos = $Overlap-$k;
                                }
                            }
                        }
                        $SMscoren=sprintf("%.2f",$left_max+$right_max);
                        if ($SMscoren > -99999) {
                            $left_flag=-5;
                            $right_flag=-3;
                        }
                        if ($debug eq 1){print $left_max,"\t",$left_pos,"\t",$right_max,"\t",$right_pos,"\n";}
						$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$left_pos."\t".$left_flag;
						$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".$right_pos."\t".$right_flag;
						$a[14]=$a[14]."\t".$SMscoren."\t".$SMS."\t".$PMS;
						print OUT join("\t",@a),"\n";
						my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
						my $strand="-";
						if ($left_flag > 0) {$strand="+";}
                        $a[4]=$a[4]-$left_pos;
						$a[6]=$a[6]-$left_pos;
						$a[8]=$a[8]+($right_pos);
						$a[9]=$a[9]-($right_pos);
						$a[10]=$a[10]+($right_pos);
                        $a[14]=$a[14]."\t".$SMscoren."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$left_pos."\t".($right_pos)."\t".$SMS."\t".$PMS;
                        print OUT2 join("\t",@a),"\n";
                    }
					else {
                        $SMscoren=-99999;
                        $left_max=0;
                        $left_pos=0;
                        $right_max=0;
                        $right_pos=0;
                        for(my $k=0; $k<=$overlap; $k++) {
                            my $t=$left_seq;					
                            $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                            my $rc_left_seq=scalar reverse $t;		
                            $t=$right_seq;					
                            $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                            my $rc_right_seq=scalar reverse $t;		
                            my $str1="";
                            $str1=uc substr($rc_left_seq,(2*$Extend+$k-20),23);
                            my $str2="";
                            $str2=uc substr($rc_right_seq,(2*$Extend+1-3-$overlap+$k),9);
                            if ((length($str1) eq 23) and (length($str2) eq 9)) {
                                my $left_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
								my $right_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
                                my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                                if ($SMscoren < $sumt) {
                                    $SMscoren=$sumt;
                                    $PMS=$k;
                                    $left_pos=$k;
                                    $left_max=$left_maxt;
                                    $right_pos=$overlap-$k;
                                    $right_max=$right_maxt;
                                    $left_flag=-5;
                                    $right_flag=-3;
                                }
                            }    
                        }
                        if ($debug eq 1){print $SMscoren,"\t",$PMS,"\t",($overlap-$PMS),"\t",$overlap,"\n";}
						$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$left_pos."\t".$left_flag;
						$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".$right_pos."\t".$right_flag;
						$a[14]=$a[14]."\t".$SMscoren."\t".$SMS."\t".$PMS;
						print OUT join("\t",@a),"\n";
						my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
						my $strand="-";
						if ($left_flag > 0) {$strand="+";}
                        $a[4]=$a[4]-$PMS;
						$a[6]=$a[6]-$PMS;
						$a[8]=$a[8]+($overlap-$PMS);
						$a[9]=$a[9]-($overlap-$PMS);
						$a[10]=$a[10]+($overlap-$PMS);
                        $a[14]=$a[14]."\t".$SMscoren."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
                        print OUT2 join("\t",@a),"\n";
					}
				}}
				if ($a[7] eq "-") {  EXIT_IF:{
					# R-m+
					my $overlap=$a[8]+$a[9]-$a[3];
					my $left_seq=substr($SEQ,($start-2*$Extend-1),(4*$Extend + 1));
					if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,"left_404"),"\n"; last EXIT_IF;}
					my $left_max=-99999;
					my $left_pos=0;
					my $left_flag=0;
					my $right_seq=substr($SEQ,($end-2*$Extend-1),(4*$Extend + 1));
					if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; last EXIT_IF; }
					if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,$end,$seqid,"right_404"),"\n";last EXIT_IF;}
					my $right_max=-99999;
					my $right_pos=0;
					my $right_flag=0;
					if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; last EXIT_IF;}
					if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; last EXIT_IF;}
					my $SMS=0;	
					my $PMS=0;
					my $SMscoren=-99999;
                    if ($debug eq 1){ print "\n"; }
                    if ($debug eq 1){print join("\t",@a),"\n",$left_seq,"\t",$right_seq,"\n",$overlap,"\n";}
					for(my $k=0; $k<=$overlap; $k++) {
						my $ml="";
						my $mr="";
						$ml=substr($left_seq,2*$Extend+1-$k,2);
						$mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
						if (($ml ne "") and ($mr ne "")) {
							my $motif=$ml.$mr;
                            if ($debug eq 1){print $k,"\t",$motif,"\n";}
							if ((exists $SMotif{$motif}) and ($SMotif{$motif} > 0)) {
								my $str1="";
								$str1=uc substr($left_seq,(2*$Extend+1-$k-3),9);
								my $str2="";
								$str2=uc substr($right_seq,(2*$Extend-2+$overlap-$k-20+2),23);
								if ((length($str1) eq 9) and (length($str2) eq 23)) {
									my $left_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
									my $right_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
									my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
									if ($SMscoren < $sumt) { $SMscoren=$sumt; $SMS=$SMotif{$motif}; $PMS=$k; }
								}
							}
						}
					}
                    if ($debug eq 1){print "SMS = ",$SMS,"\tPSM = ",$PMS,"\tSMscoren = $SMscoren \n";}
					if ($SMS > 0) {
						my $str1="";
						$str1=uc substr($left_seq,(2*$Extend+1-$PMS-3),9);
						my $str2="";
						$str2=uc substr($right_seq,(2*$Extend-2+$overlap-$PMS-20+2),23);
						if ((length($str1) eq 9) and (length($str2) eq 23)) {
							$left_max=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
							$right_max=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
							my $sum=sprintf("%.2f",($left_max + $right_max));
							$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."+5";
							$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."+3";
							$a[14]=$a[14]."\t".$sum."\t".$SMS."\t".$PMS;
							print OUT join("\t",@a),"\n";
							my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
							$a[3]=$a[3]+$PMS;
							$a[4]=$a[4]-$PMS;
							$a[6]=$a[6]-$PMS;
							$a[9]=$a[9]-($overlap-$PMS);
							$a[10]=$a[10]+($overlap-$PMS);
							$a[14]=$a[14]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."+"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
							print OUT1 join("\t",@a),"\n";
							last EXIT_IF;
						}
					}
                    elsif (($overlap < 0) or ($a[13] > 0)) {
                        my $Overlap=$Ext;
                        $left_max=-99999;
                        $left_pos=0;
                        $right_max=-99999;
                        $right_pos=0;
                        my @SScoreL;
                        my @SScoreR;
                        for(my $k=0; $k<=$Overlap; $k++) {
                            my $str1="";
                            $str1=uc substr($left_seq,(2*$Extend+1-$k-3),9);
                            my $str2="";
                            $str2=uc substr($right_seq,(2*$Extend-2+$Overlap-$k-20+2),23);
                            if ((length($str1) eq 9) and (length($str2) eq 23)) {
                                $SScoreL[$k]=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                                $SScoreR[$Overlap-$k]=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
                                if ($debug > 0) { print join("\t",$k,$str1,$SScoreL[$k],$str2,$SScoreR[$Overlap-$k]),"\n"; }
                                if ($left_max < $SScoreL[$k]) {
                                    $left_max = $SScoreL[$k];
                                    $left_pos = $k;
                                }
                                if ($right_max < $SScoreR[$Overlap-$k]) {
                                    $right_max = $SScoreR[$Overlap-$k];
                                    $right_pos = $Overlap-$k;
                                }
                            }
                        }
                        $SMscoren=sprintf("%.2f",$left_max+$right_max);
                        if ($SMscoren > -99999) {
                            $left_flag=+5;
                            $right_flag=+3;
                        }
                        if ($debug eq 1){print $left_max,"\t",$left_pos,"\t",$right_max,"\t",$right_pos,"\n";}
						$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$left_pos."\t".$left_flag;
						$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".$right_pos."\t".$right_flag;
						$a[14]=$a[14]."\t".$SMscoren."\t".$SMS."\t".$PMS;
						print OUT join("\t",@a),"\n";
						my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
						my $strand="+";
						if ($left_flag < 0) {$strand="-";}
                        $a[3]=$a[3]+$left_pos;
                        $a[4]=$a[4]-$left_pos;
                        $a[6]=$a[6]-$left_pos;
                        $a[9]=$a[9]-($right_pos);
                        $a[10]=$a[10]+($right_pos);
                        $a[14]=$a[14]."\t".$SMscoren."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$left_pos."\t".($right_pos)."\t".$SMS."\t".$PMS;
                        print OUT2 join("\t",@a),"\n";
                    }
					else {
                        $SMscoren=-99999;
                        $left_max=0;
                        $left_pos=0;
                        $right_max=0;
                        $right_pos=0;
                        for(my $k=0; $k<=$overlap; $k++) {
                            my $str1="";
                            $str1=uc substr($left_seq,(2*$Extend+1-$k-3),9);
                            my $str2="";
                            $str2=uc substr($right_seq,(2*$Extend-2+$overlap-$k-20+2),23);
                            if ((length($str1) eq 9) and (length($str2) eq 23)) {
                                my $left_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                                my $right_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
                                my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                                if ($SMscoren < $sumt) {
                                    $SMscoren=$sumt;
                                    $PMS=$k;
                                    $left_pos=$k;
                                    $left_max=$left_maxt;
                                    $right_pos=$overlap-$k;
                                    $right_max=$right_maxt;
                                    $left_flag=+5;
                                    $right_flag=+3;
                                }
                            }    
                        }
                        if ($debug eq 1){print $SMscoren,"\t",$PMS,"\t",($overlap-$PMS),"\t",$overlap,"\n";}
                        $a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$left_pos."\t".$left_flag;
                        $a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".$right_pos."\t".$right_flag;
                        $a[14]=$a[14]."\t".$SMscoren."\t".$SMS."\t".$PMS;
                        print OUT join("\t",@a),"\n";
                        my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
                        my $strand="+";
                        if ($left_flag < 0) {$strand="-";}
                        $a[3]=$a[3]+$PMS;
                        $a[4]=$a[4]-$PMS;
                        $a[6]=$a[6]-$PMS;
                        $a[9]=$a[9]-($overlap-$PMS);
                        $a[10]=$a[10]+($overlap-$PMS);
                        $a[14]=$a[14]."\t".$SMscoren."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
                        print OUT2 join("\t",@a),"\n";
					}
				}}
			}
	    }
    }
}
close OUT;
close OUT1;
close OUT2;

my $file0=$fileout.".S2.0";
my $file1=$fileout.".S2.1";
my $file2=$fileout.".S2.sum";
$command="cat $file0 $file1 | sort -k3,3 -k7,7n -k11,11n > $file2";
system($command);


