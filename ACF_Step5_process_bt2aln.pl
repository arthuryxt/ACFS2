#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output\"   \"input_bowtie2_alignment\"  \"UNMAP_expr\"  \"MEA_basename\"   \"CBR_basename\"   \"\(optional\)strandess = no\"    \"\(optional\)min_junction = 6\"    \"\(optional\)insertsize = 150\"    \"\(optional\)error = 0.05\"  \"\(optional\)debug\"   " if (@ARGV < 5);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];
my $UNMAPEXPR=$ARGV[2];
my $MEA=$ARGV[3];
my $CBR=$ARGV[4];
my $strandness="no";
if (scalar(@ARGV) > 5) {$strandness=$ARGV[5];}
if (($strandness ne "no") and ($strandness ne "+") and ($strandness ne "-")){
    die "strandess can only be one of the three: 1) no; 2) +; 3) -"
}
my $junc=6;
if (scalar(@ARGV) > 6) {$junc=$ARGV[6];}
my $insertSize=150;
if (scalar(@ARGV) > 7) {$insertSize=$ARGV[7];}
my $errorrate=0.1;
if (scalar(@ARGV) > 8) {$errorrate=$ARGV[8];}
if ($errorrate <= 0) { die "cutoff must be non-negative: (0,1] means error_rate; >1 means edit-distance." }
my $debug=0;
if (scalar(@ARGV) > 9) {$debug=$ARGV[9];}

my %usedreads;
my %scorereads;
my $tmpS3=$MEA;
if (-e $tmpS3) {
    open(IN1, $tmpS3);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        my @b=split(/\,/,$a[26]);
        for(@b){
            my $id=$_; if($id=~m/^rc\_R2\_/){$id=~s/^rc\_R2\_//;} if($id=~m/^rc\_/){$id=~s/^rc\_//;}  if($id=~m/^R2\_/){$id=~s/^R2\_//;}
            if(exists $usedreads{$id}){ if($scorereads{$id} < $a[17]){ $usedreads{$id}=$a[0]; $scorereads{$id}=$a[17];  } }
            else{ $usedreads{$id}=$a[0]; $scorereads{$id}=$a[17]; }
        }
    }
    close IN1;
}
$tmpS3=$CBR;
if (-e $tmpS3) {
    open(IN1, $tmpS3);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        my @b=split(/\,/,$a[26]);
        for(@b){
            my $id=$_; if($id=~m/^rc\_R2\_/){$id=~s/^rc\_R2\_//;} if($id=~m/^rc\_/){$id=~s/^rc\_//;}  if($id=~m/^R2\_/){$id=~s/^R2\_//;}
            if(exists $usedreads{$id}){ if($scorereads{$id} < $a[17]){ $usedreads{$id}=$a[0]; $scorereads{$id}=$a[17];  } }
            else{ $usedreads{$id}=$a[0]; $scorereads{$id}=$a[17]; }
        }
    }
    close IN1;
}
#foreach my $id (sort keys %usedreads){print $id,"\t",$usedreads{$id},"\n";}
my %refFlat;
my $fileintmp=$MEA.".refFlat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $refFlat{$a[1]}=join("\t",@a);
    }
    close IN1;
}
$fileintmp=$MEA.".ext.refFlat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $refFlat{$a[1]}=join("\t",@a);
    }
    close IN1;
}
$fileintmp=$CBR.".refFlat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $refFlat{$a[1]}=join("\t",@a);
    }
    close IN1;
}
$fileintmp=$CBR.".ext.refFlat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $refFlat{$a[1]}=join("\t",@a);
    }
    close IN1;
}

my %cidstat;
$fileintmp=$MEA.".stat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $cidstat{$a[0]}=1;
    }
    close IN1;
}
$fileintmp=$CBR.".stat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $cidstat{$a[0]}=1;
    }
    close IN1;
}
$fileintmp="temp.pre_defined_circRNA.stat";
if (-e $fileintmp) {
    open(IN1, $fileintmp);
    while (<IN1>) {
        chomp;
        my @a=split("\t",$_);
        $cidstat{$a[0]}=1;
    }
    close IN1;
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    #$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros return $str;
    my $decode = substr($str,-12);
    return $decode;
}

sub getSeqLenFromCIGAR ($) {
    my $MD=$_[0];
    my @CIGAR_op=($MD=~m/[MSHID]/g); 
    my @CIGAR_va=($MD=~m/\d+/g);
    my $Nr=scalar(@CIGAR_op);
    my $length=0;
    for(my $i=0; $i<$Nr; $i++){
        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; }
        elsif ($CIGAR_op[$i] eq "D") {}
        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; }
        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) { $length+=$CIGAR_va[$i]; }
    }
    return($length);
}

my %uniq;
my %uniq2;
my %reads;
my %reads2;
my %len;
open(IN, $filein);
while (<IN>) {
    chomp;
    if (m/^@/) {
        my @a=split("\t",$_);
        my @b=split(/\:/,$a[1]);
        if (($a[0] eq "\@SQ") and ($b[0] eq "SN")) {
            my @c=split(/\_\_\_/,$b[1]);
            my @d=split(/\:/,$a[2]);
            $len{$c[0]}=$d[1];
        }
        next;
    }
    next if (m/^#/);
    next if (m/^>/);
    my @a=split("\t",$_);
    $a[0]=~s/\/1$//;
    $a[0]=~s/\/2$//;
    my $cFlag=dec2bin($a[1]);
    my @aFlag=split('',$cFlag);
    if ($aFlag[-3] eq 1) { next; }
    if ($aFlag[-1] eq 1) {
        # PE reads
        my $seqid=$a[0];
        my $mate=1;
        if ($aFlag[-8] eq 1) { $mate=2; }
        if ((exists $usedreads{$a[0]}) and (!exists $reads2{$a[0]}{$mate}) ) { $reads2{$a[0]}{$mate}=$a[8];}
        
        if (($aFlag[-7] eq 1) and ($strandness eq "-") and ($aFlag[-5] eq 0)) { next; }
        if (($aFlag[-8] eq 1) and ($strandness eq "-") and ($aFlag[-5] eq 1)) { next; }
        if (($aFlag[-7] eq 1) and ($strandness eq "+") and ($aFlag[-5] eq 1)) { next; }
        if (($aFlag[-8] eq 1) and ($strandness eq "+") and ($aFlag[-5] eq 0)) { next; }
        my $NM=99999;
        my $MD="";
        for(my $i=11; $i<scalar(@a); $i++) {
            if ($a[$i]=~m/NM/) {
                my @b=split(/\:/,$a[$i]);
                $NM=$b[2];
            }
            if ($a[$i]=~m/MD/) { $MD=$a[$i]; }
        }
        # Number of Mismatches has to also count unmapped part
        my @CIGAR_op=($a[5]=~m/[MSHID]/g); 
        my @CIGAR_va=($a[5]=~m/\d+/g);
        my $Nr=scalar(@CIGAR_op);
        my $unmappedlen=0;
        for(my $i=0; $i<$Nr; $i++){
            if ($CIGAR_op[$i] eq "S") {$unmappedlen+=$CIGAR_va[$i]; }
            elsif ($CIGAR_op[$i] eq "H") {$unmappedlen+=$CIGAR_va[$i]; }
        }
        $NM+=$unmappedlen;
        my $errors=$errorrate;
        if ((0 < $errorrate) and ($errorrate < 1)) {  $errors=length($a[9])*$errorrate; }
        if ($NM > $errors) { next; }
        
        my @b=split(/\_\_\_/,$a[2]);
        my @c=split(/\_/,$b[0]);
        my @d=split(/\|/,$c[4]);
        my $length=0;
        if (exists $len{$b[0]}) { $length=$len{$b[0]}; }
        my $myinsert=$insertSize;
        if ($length < 2*$myinsert) { $myinsert=int($length/2); }
        if (exists $uniq2{$seqid}{$b[0]}{$mate}) {
            my @tmp=split("\t",$uniq2{$seqid}{$b[0]}{$mate});
            if ($NM < $tmp[1]) {
                $uniq2{$seqid}{$b[0]}{$mate}=join("\t",getSeqLenFromCIGAR($a[5]),scalar(@d),$NM,$a[3],$a[5],$MD,$length);
            }
        }
        else { $uniq2{$seqid}{$b[0]}{$mate}=join("\t",getSeqLenFromCIGAR($a[5]),scalar(@d),$NM,$a[3],$a[5],$MD,$length); }
        if ((exists $reads2{$seqid}) and (exists $reads2{$seqid}{$mate}) and (length($reads2{$seqid}{$mate}) < length($a[8]))){
            $reads2{$seqid}{$mate}=$a[8];
        }
        else {$reads2{$seqid}{$mate}=$a[8];}
        
    }
    elsif (($aFlag[-1] eq 0) or ($aFlag[-2] eq 0)) {
        # SE reads
        if ((exists $usedreads{$a[0]}) and (!exists $reads{$a[0]})) { $reads{$a[0]}=$a[8];}
        if (($strandness eq "+") and ($aFlag[-5] eq 1) and ($aFlag[-1] eq 0)) { next; }
        if (($strandness eq "-") and ($aFlag[-5] eq 0) and ($aFlag[-1] eq 0)) { next; }
        my $NM=99999;
        my $MD="";
        for(my $i=11; $i<scalar(@a); $i++) {
            if ($a[$i]=~m/NM/) {
                my @b=split(/\:/,$a[$i]);
                $NM=$b[2];
            }
            if ($a[$i]=~m/MD/) { $MD=$a[$i]; }
        }
        # Number of Mismatches has to also count unmapped part
        my @CIGAR_op=($a[5]=~m/[MSHID]/g); 
        my @CIGAR_va=($a[5]=~m/\d+/g);
        my $Nr=scalar(@CIGAR_op);
        my $unmappedlen=0;
        for(my $i=0; $i<$Nr; $i++){
            if ($CIGAR_op[$i] eq "S") {$unmappedlen+=$CIGAR_va[$i]; }
            elsif ($CIGAR_op[$i] eq "H") {$unmappedlen+=$CIGAR_va[$i]; }
        }
        $NM+=$unmappedlen;
        my $errors=$errorrate;
        if ((0 < $errorrate) and ($errorrate < 1)) {  $errors=length($a[9])*$errorrate; }
        if ($NM > $errors) { next; }
        my $seqid=$a[0];
        my @b=split(/\_\_\_/,$a[2]);
        my @c=split(/\_/,$b[0]);
        my @d=split(/\|/,$c[4]);
        my $length=0;
        if (exists $len{$b[0]}) { $length=$len{$b[0]}; }
        my $myinsert=$insertSize;
        if ($length < 2*$myinsert) { $myinsert=int($length/2); }
        my $seqlen=getSeqLenFromCIGAR($a[5]);
        if ((($a[3] + $junc) <= ($length - $myinsert)) and (($length - $myinsert) <= ($a[3] + $seqlen - $junc) )) {
            if (exists $uniq{$seqid}) {
                if ($debug > 0) { print $uniq{$seqid},"\n",join("\t",@a),"\n"; }
                my @tmpc=split("\t",$uniq{$seqid});
                if ($tmpc[3] > $NM) { # select the one with fewer mismatches
                    $uniq{$seqid}=join("\t",$b[0],$seqlen,scalar(@d),$NM,$a[3],$a[5],$MD,$length);
                }
                elsif ($tmpc[3] eq $NM) { 
                    if (exists $cidstat{$b[0]}) { 
                        my $Nr=int(scalar(@tmpc)/8);
                        my $flag=0;
                        for(my $i=0; $i<$Nr; $i++) { if(exists $cidstat{$tmpc[8*$i]}){ $flag++;} }
                        if ($debug > 0) { print "this is good\t",$flag,"\n"; }
                        if ($flag eq 0) { $uniq{$seqid}=join("\t",$b[0],$seqlen,scalar(@d),$NM,$a[3],$a[5],$MD,$length);}
                        else { $uniq{$seqid}=join("\t",$uniq{$seqid},$b[0],$seqlen,scalar(@d),$NM,$a[3],$a[5],$MD,$length); }
                    }
                    else {
                        my $Nr=int(scalar(@tmpc)/8);
                        my $flag=0;
                        for(my $i=0; $i<$Nr; $i++) { if(exists $cidstat{$tmpc[8*$i]}){ $flag++;} }
                        if ($debug > 0) { print "hash is good\t",$flag,"\n"; }
                        if ($flag eq 0) { $uniq{$seqid}=join("\t",$uniq{$seqid},$b[0],$seqlen,scalar(@d),$NM,$a[3],$a[5],$MD,$length); }
                    }
                }
                if ($debug > 0) { print $uniq{$seqid},"\n\n"; }
            }
            else { $uniq{$seqid}=join("\t",$b[0],$seqlen,scalar(@d),$NM,$a[3],$a[5],$MD,$length); }
            $reads{$seqid}=$a[8];
        }
    }
}
close IN;

my %circMap;
foreach my $id (keys %len) {
    my @a=split(/\_/,$id);
    my $newid=join("_",$a[0],$a[1],$a[2],$a[3]);
    $circMap{$newid}{$id}=1;
}

my %Anno;
my $header="";
my @Header;
if (($UNMAPEXPR ne "no") and (-e $UNMAPEXPR)) {
    open IN4, $UNMAPEXPR;
    $header=<IN4>;
    chomp $header;
    @Header=split("\t",$header);
    my $tmpid=$Header[0];
    $Header[0]=$Header[0]."\tGname";
    $header=join("\t",@Header);
    while(<IN4>) {
        chomp;
        my @a=split("\t",$_);
        $Anno{$a[0]}=$a[1];
        for(my $i=2; $i<scalar(@a); $i++) { $Anno{$a[0]}=$Anno{$a[0]}."\t".$a[$i]; }
    }
    close IN4;
}
else {
    $header=join("\t","newid","Sample");
    @Header=split("\t",$header);
    $Header[0]=$Header[0]."\tGname";
    $header=join("\t",@Header);
    foreach my $id(keys %usedreads) { $Anno{$id}=1; }
    foreach my $id(keys %uniq) { $Anno{$id}=1; }
    foreach my $id(keys %uniq2) { $Anno{$id}=1; }
}

my $Nr=scalar(@Header);
my $template=0;
for(my $i=2; $i<$Nr; $i++) {$template=$template."\t0";}



open(OUT1, ">".$fileout.".reads");
open(OUT2, ">".$fileout.".uncertain");
my %uncertain;
my %EXPR;
my %usedreads2;

foreach my $read (keys %uniq2){
    my %mismatches;
    my %exoncnt;
    my $good=0;
    foreach my $circ (keys %{$uniq2{$read}}) {
        my $p=0;
        foreach my $mate (keys %{$uniq2{$read}{$circ}}) { $p++; }
        if ($p eq 2) {
            my @mate1=split("\t",$uniq2{$read}{$circ}{1});
            my @mate2=split("\t",$uniq2{$read}{$circ}{2});
            my $left=$mate1[3] < $mate2[3] ? $mate1[3] : $mate2[3];
            if ($left eq $mate2[3]) {
                @mate1=split("\t",$uniq2{$read}{$circ}{2});
                @mate2=split("\t",$uniq2{$read}{$circ}{1});
            }
            my $right=$mate2[0]+$mate2[3];
            my $length=$mate1[6];
            my $myinsert=$insertSize;
            if ($length < 2*$myinsert) { $myinsert=int($length/2); }
            if ($debug > 0) { print $read,"\n",join("\t",$circ,@mate1,@mate2),"\n",join("\t",($left + $junc),($length - $myinsert),($right - $junc)),"\n";}
            if ((($left + $junc) < ($length - $myinsert)) and (($length - $myinsert) < ($right - $junc) )) {
                $mismatches{$circ}=$mate1[2]+$mate2[2];
                $exoncnt{$circ}=$mate1[1]+$mate2[1];
                if (exists $cidstat{$circ}) { $good++;}
            }
        }
        else {
            
        }
    }
    if ($debug > 0) { print $good,"\n";}
    my $Scirc="";
    if ($good > 0) {
        my $min=99999;
        foreach my $circ (sort{$mismatches{$a} <=> $mismatches{$b}} keys %mismatches) { if (exists $cidstat{$circ}){$min=$mismatches{$circ}; last;} }
        my $count=0;
        foreach my $circ (keys %mismatches) {
            if (($mismatches{$circ} eq $min) and (exists $cidstat{$circ})) { $count++; $Scirc=$circ; if ($debug > 0) { print $circ,"\n";}}
        }
        if ($count eq 1) {
            print OUT1 join("\t",$read,$Scirc,$uniq2{$read}{$Scirc}{1},$uniq2{$read}{$Scirc}{2}),"\n";
            $usedreads2{$read}=1;
            if (exists $EXPR{$Scirc}) {
                my @samples=split("\t",$Anno{$read});
                my @stored=split("\t",$EXPR{$Scirc});
                for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]; }
                $EXPR{$Scirc}=join("\t",@stored);
            }
            else{
                $EXPR{$Scirc}=$Anno{$read};
            }
        }
        elsif ($count > 1) {
            my $min_len=99999999;
            my $shortest_circ="";
            foreach my $circ (keys %mismatches) {
                if ($mismatches{$circ} eq $min) {
                    if (!exists $uncertain{$read}) {
                        $uncertain{$read}=join("\t",$circ,$uniq2{$read}{$circ}{1},$uniq2{$read}{$circ}{2});
                    }
                    else {
                        $uncertain{$read}=join("\t",$uncertain{$read},$circ,$uniq2{$read}{$circ}{1},$uniq2{$read}{$circ}{2});
                    }
                    print OUT2 join("\t",$read,$circ,$uniq2{$read}{$circ}{1},$uniq2{$read}{$circ}{2}),"\n";
                    my @tmp=split("\t",$uniq2{$read}{$circ}{1});
                    if ($tmp[6] < $min_len) {
                        $min_len=$tmp[6];
                        $shortest_circ=$circ;
                    }
                }
            }
            
            if ($shortest_circ ne "") {
                print OUT1 join("\t",$read,$Scirc,$uniq2{$read}{$shortest_circ}{1},$uniq2{$read}{$shortest_circ}{2}),"\n";
                $usedreads2{$read}=1;
                if (exists $EXPR{$shortest_circ}) {
                    my @samples=split("\t",$Anno{$read});
                    my @stored=split("\t",$EXPR{$shortest_circ});
                    for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]; }
                    $EXPR{$shortest_circ}=join("\t",@stored);
                }
                else{
                    $EXPR{$shortest_circ}=$Anno{$read};
                }
            }
        }
    }
    else {
        my $min=99999;
        foreach my $circ (sort{$mismatches{$a} <=> $mismatches{$b}} keys %mismatches) { $min=$mismatches{$circ}; last; }
        my $count=0;
        foreach my $circ (keys %mismatches) {
            if ($mismatches{$circ} eq $min) { $count++; $Scirc=$circ; }
        }
        if ($count eq 1) {
            print OUT1 join("\t",$read,$Scirc,$uniq2{$read}{$Scirc}{1},$uniq2{$read}{$Scirc}{2}),"\n";
            $usedreads2{$read}=1;
            if (exists $EXPR{$Scirc}) {
                my @samples=split("\t",$Anno{$read});
                my @stored=split("\t",$EXPR{$Scirc});
                for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]; }
                $EXPR{$Scirc}=join("\t",@stored);
            }
            else{
                $EXPR{$Scirc}=$Anno{$read};
            }
        }
        elsif ($count > 1) {
            my $min_len=99999999;
            my $shortest_circ="";
            foreach my $circ (keys %mismatches) {
                if ($mismatches{$circ} eq $min) {
                    if (!exists $uncertain{$read}) {
                        $uncertain{$read}=join("\t",$circ,$uniq2{$read}{$circ}{1},$uniq2{$read}{$circ}{2});
                    }
                    else {
                        $uncertain{$read}=join("\t",$uncertain{$read},$circ,$uniq2{$read}{$circ}{1},$uniq2{$read}{$circ}{2});
                    }
                    print OUT2 join("\t",$read,$circ,$uniq2{$read}{$circ}{1},$uniq2{$read}{$circ}{2}),"\n";
                    my @tmp=split("\t",$uniq2{$read}{$circ}{1});
                    if ($tmp[6] < $min_len) {
                        $min_len=$tmp[6];
                        $shortest_circ=$circ;
                    }
                }
            }
            
            if ($shortest_circ ne "") {
                print OUT1 join("\t",$read,$Scirc,$uniq2{$read}{$shortest_circ}{1},$uniq2{$read}{$shortest_circ}{2}),"\n";
                $usedreads2{$read}=1;
                if (exists $EXPR{$shortest_circ}) {
                    my @samples=split("\t",$Anno{$read});
                    my @stored=split("\t",$EXPR{$shortest_circ});
                    for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]; }
                    $EXPR{$shortest_circ}=join("\t",@stored);
                }
                else{
                    $EXPR{$shortest_circ}=$Anno{$read};
                }
            }
        }
    }
}

foreach my $read (keys %uniq){
    my @tmp=split("\t",$uniq{$read});
    my $Nr=int(scalar(@tmp)/8);
    if ($Nr > 1) {
        $uncertain{$read}=$uniq{$read};
        print OUT2 $read,"\t",$uniq{$read},"\n";
        # pick the shortest circRNA to report
        my $circ=0;
        my $min_len=99999999;
        for(my $i=0; $i<$Nr; $i++) {
            if ($tmp[8*$i + 7] < $min_len) {
                $min_len=$tmp[8*$i + 7];
                $circ=$i;
            }
        }
        
        print OUT1 $read,"\t",join("\t",$tmp[8*$circ],$tmp[8*$circ+1],$tmp[8*$circ+2],$tmp[8*$circ+3],$tmp[8*$circ+4],$tmp[8*$circ+5],$tmp[8*$circ+6],$tmp[8*$circ+7]),"\n";
        $usedreads2{$read}=1;
        if (exists $EXPR{$tmp[8*$circ]}) {
            my @samples=split("\t",$Anno{$read});
            my @stored=split("\t",$EXPR{$tmp[8*$circ]});
            for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]; }
            $EXPR{$tmp[0]}=join("\t",@stored);
        }
        else{
            $EXPR{$tmp[0]}=$Anno{$read};
        }
    }
    else {
        print OUT1 $read,"\t",$uniq{$read},"\n";
        $usedreads2{$read}=1;
        if (exists $EXPR{$tmp[0]}) {
            my @samples=split("\t",$Anno{$read});
            my @stored=split("\t",$EXPR{$tmp[0]});
            for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]; }
            $EXPR{$tmp[0]}=join("\t",@stored);
        }
        else{
            $EXPR{$tmp[0]}=$Anno{$read};
        }
    }
}

foreach my $read (keys %usedreads) {
    my $circ=$usedreads{$read};
    if ((exists $circMap{$circ}) and (!exists $usedreads2{$read})) {
        my $cnt=0;
        if (exists $reads{$read})  {
            print OUT1 join("\t",$read,$circ,0,0,0,0,0,0,0),"\n";
        }
        else {
            print OUT1 join("\t",$read,$circ,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"\n";
        }
        foreach my $mycirc (keys %{$circMap{$circ}}) {$cnt++;}
        foreach my $mycirc (keys %{$circMap{$circ}}) {
            if (exists $EXPR{$mycirc}) {
                my @samples=split("\t",$Anno{$read});
                my @stored=split("\t",$EXPR{$mycirc});
                for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]+=$samples[$i]/$cnt; }
                $EXPR{$mycirc}=join("\t",@stored);
            }
            else{
                my @samples=split("\t",$Anno{$read});
                my @stored;
                for (my $i=0; $i<scalar(@samples); $i++) { $stored[$i]=$samples[$i]/$cnt; }
                $EXPR{$mycirc}=join("\t",@stored);
            }
        }
    }
}

close OUT1;
close OUT2;

open(OUTexpr, ">".$fileout."_expr");
print OUTexpr join("\t",@Header),"\n";
open(OUTrefFlat, ">".$fileout."_refFlat");
foreach my $id (sort keys %EXPR) {
    my @ref=split("\t",$refFlat{$id});
    print OUTexpr join("\t",$id,$ref[0],$EXPR{$id}),"\n";
    print OUTrefFlat $refFlat{$id},"\n";
}

