#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output\"    \"BWA_sam\"   \"\(optional\)stranded=\"\-\"\"    \"\(optional\)min_AS default=30\"   \"\(optional\)coverage default=0.9\"" if (@ARGV < 2);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];
my $stranded="-";   # assume reads are anti-sense to mRNA as TruSeq-Stranded-RNASeq with default "-"; "+" otherwise; "no" assume no strandness and will report also the reverse-complemented
if (scalar(@ARGV) > 2) {$stranded=$ARGV[2];}
my $MAS=30;
if (scalar(@ARGV) > 3) {$MAS=$ARGV[3];}
my $cutoff=0.9;
if (scalar(@ARGV) > 4) {$cutoff=$ARGV[4];}
if (($cutoff <= 0) or ($cutoff >1)) { die "cutoff must be non-negative: (0,1] " }
my %uniq;   # store processed ID
open IN, $filein;
open OUT,">".$fileout.".tmp";
open OUT1,">".$fileout.".1";    # single hit
open OUT2,">".$fileout.".multi";
open OUT4,">".$fileout.".segs";
open OUT21,">".$fileout.".2pp"; # two hits on the same chr and same strand
open OUT22,">".$fileout.".2pm"; # two hits on the same chr BUT Different strand
open OUT3,">".$fileout.".unmap";
open OUTMT,">".$fileout.".MT";
open OUT2S,">".$fileout.".2pp.S1"; # selected candidates from two hits on the same chr and same strand
open OUTUID,">".$fileout.".UID";
my $command="rm -f Step1_finished";
system($command);

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    #$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros return $str;
    my $decode = substr($str,-12);
    return $decode;
}
sub reverseStrand($){
    my $s=shift;
    my $r;
    if ($s eq "+") { $r="-";}
    else{$r="+";}
    return $r;
}

while(<IN>) {
    chomp;
    if ((m/^#/) or (m/^@/)) {next;}
    my @a=split("\t",$_);
    #my $Len=length($a[9]);
    my $Len=0;
    my @len_va=($a[5]=~m/\d+/g);
    my @len_op=($a[5]=~m/[MHSID]/g); 
    for(my $i=0; $i<scalar(@len_va); $i++){ if($len_op[$i] ne "D" ){$Len+=$len_va[$i];} }
    
    my $cFlag=dec2bin($a[1]);
    my @aFlag=split('',$cFlag);

    if ($aFlag[-3] eq 1) {
        print OUT3 $a[0],"\t",$a[9],"\n";
        print OUTUID $a[0],"\n";
        next;
    }
    my $strand="+";
    if ($aFlag[-5] eq 1) {$strand="-";}
    if ($aFlag[-8] eq 1) {$strand=reverseStrand($strand); $a[0]="R2_".$a[0];}
    if (exists $uniq{$a[0]}) { next; }
    else {
        $uniq{$a[0]}=1;
        if (($stranded eq "-") or ($stranded eq "no")) {
            my @d=split(/\:/,$a[11]);
            my $info="";
            # $a[14] in bwa-0.7.3a : XP:Z:13,-91704253,49M52S,60,0;
            # $a[15] in bwa-0.7.15 : SA:Z:13,91704253,-,49M52S,60,0;
            for(my $i=11; $i<scalar(@a); $i++) {
                if ($a[$i]=~m/^XP:Z:/) {
                    my @b=split(/:/,$a[$i]);
                    my @c=split(/\;/,$b[2]);
                    if ($aFlag[-8] eq 0) {
                        $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
                    }
                    else {
                        for (my $j=0; $j<scalar(@c); $j++) {
                            my @tmpp=split(/\,/,$c[$j]);
                            my $tmpstrand=substr($tmpp[1],0,1);
                            my $tmppos=substr($tmpp[1],1);
                            $tmpp[1]=reverseStrand($tmpstrand).$tmppos;
                            $c[$j]=join(",",@tmpp);
                        }
                        $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
                    }
                    last;
                }
                elsif ($a[$i]=~m/^SA\:Z\:/) {
                    my @b=split(/:/,$a[$i]);
                    my @c=split(/\;/,$b[2]);
                    for(my $j=0; $j<scalar(@c); $j++) {
                        if (length($c[$j]) > 1) {
                            my @tmp=split(/\,/,$c[$j]);
                            if ($aFlag[-8] eq 0) { $c[$j]=join(",",$tmp[0],$tmp[2].$tmp[1],$tmp[3],$tmp[4],$tmp[5]); }
                            else { $c[$j]=join(",",$tmp[0],reverseStrand($tmp[2]).$tmp[1],$tmp[3],$tmp[4],$tmp[5]); }
                        }
                    }
                    $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
                    last;
                }
            }
            if ($info eq "") {
                $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]));
            }
            
            if (($a[2] eq "MT") or ($a[2] eq "M")) {
                print OUTMT join("\t",$a[0],$Len,$info),"\n";
                next;
            }
            my $indel=0;
            if(scalar(@a) <= 14){
                my @CIGAR_op=($a[5]=~m/[MHSID]/g); 
                my @CIGAR_va=($a[5]=~m/\d+/g);
                my $start=0;
                my $length=0;
                my $glength=0;
                my $Nr=scalar(@CIGAR_op);
                if ($strand eq "+") {
                    for(my $i=0; $i<$Nr; $i++){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$glength=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                else {
                    for(my $i=$Nr-1; $i>=0; $i--){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$glength+=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                print OUT1 join("\t",$a[0],$Len,$start,$length,$d[2],$a[2],$strand,$a[3],$a[4],$indel),"\n";
                print OUT join("\t",$a[0],$Len,$info),"\n";
                next;
            }
            print OUT join("\t",$a[0],$Len,$info),"\n";
            # process $info
            my %anno;
            my %Chr;
            my @p=split("\t",$info);
            if ($MAS > 0) {
                my $tmpmax=0;
                for(@p) { my @tmp=split(/\,/,$_); if ($tmp[3] > $tmpmax) {$tmpmax=$tmp[3]}; }
                if ($tmpmax < $MAS) { next; }
            }
            for(@p) {
                my @tmp=split(/\,/,$_);
                # ignore bad-quality alignments
                # if ($tmp[3] < $MAS) {next;}
                my $strandy="+";
                if ($tmp[1]=~m/^-/){
                    $strandy="-";
                    $tmp[1]=~s/^\-//;
                }
                else {
                    $tmp[1]=~s/^\+//;
                }
                my @CIGAR_op=($tmp[2]=~m/[MSID]/g); 
                my @CIGAR_va=($tmp[2]=~m/\d+/g);
                my $start=0;
                my $length=0;
                my $glength=0;
                my $Nr=scalar(@CIGAR_op);
                if ($strandy eq "+") {
                    for(my $i=0; $i<$Nr; $i++){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i]; }
                        elsif ($CIGAR_op[$i] eq "D") {$glength+=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                else {
                    for(my $i=$Nr-1; $i>=0; $i--){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$glength+=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                if (exists $anno{$tmp[0]}{$start}) {}
                else {
                    $Chr{$tmp[0]}++;
                    $anno{$tmp[0]}{$start}=join("\t",$length,$tmp[1],($tmp[1]+$glength),$strandy);
                }
            }
            # pick the most-hit chromosome
            my $chromo="";
            foreach my $id(sort{$Chr{$b} <=> $Chr{$a}} keys %Chr) {
                if ($Chr{$id} > 1) {$chromo=$id; last;}
            }
            # each part hit once
            if ($chromo eq "") {
                #my $result=$Len;
                #my $cov1=0; 
                #my $cov2=0; 
                #my @Cov;
                #for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                #foreach my $id(sort keys %anno) {
                #    foreach my $pos (sort{$a <=> $b} keys %{$anno{$id}}) {
                #        $result=$result."\t".$id."\t".$pos."\t".$anno{$id}{$pos};
                #        my @t=split("\t",$anno{$id}{$pos});
                #        for(my $i=0; $i<$t[0]; $i++) {$Cov[$i+$pos]=1;}
                #    }
                #}
                #my $first=-1;
                #my $last=-1;
                #for(my $i=0; $i<$Len; $i++) {
                #    if ($Cov[$i] eq 1) {
                #        $cov1++;
                #        if ($first eq -1) {$first = $i;}
                #        $last=$i;
                #    }
                #}
                #$cov2=$last - $first + 1;
                #if ($cov1 > 0) { print OUT2 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\n"; }
            }
            # exactly two-hitter, interesting !!!
            elsif ($Chr{$chromo} eq 2) {
                my $flag=1;
                my $gap=0;
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    if ($tmp[3] eq "-") {$flag+=2;}
                }
                if ($flag eq 1) {
                    # both segments on the plus strand
                    my $result=$Len."\t".$chromo;
                    my $cov1=0; 
                    my $cov2=0; 
                    my @Cov;
                    for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                    foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                        $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                        my @tmp=split("\t",$anno{$chromo}{$pos});
                        for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                        if ($gap eq 0) {$gap = $tmp[2];}
                        else {$gap=$tmp[1] - $gap;}
                    }
                    my $first=-1;
                    my $last=-1;
                    for(my $i=0; $i<$Len; $i++) {
                        if ($Cov[$i] eq 1) {
                            $cov1++;
                            if ($first eq -1) {$first = $i;}
                            $last=$i;
                        }
                    }
                    $cov2=$last - $first + 1;
                    print OUT21 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    print OUTUID $a[0],"\n";
                    if (($cov2 eq  $cov1) and ($cov2 >= $Len * $cutoff)) {
                        print OUT2S $a[0],"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    }
                }
                elsif ($flag eq 5) {
                    # both segments on the minus strand
                    my $result=$Len."\t".$chromo;
                    my $cov1=0; 
                    my $cov2=0; 
                    my @Cov;
                    for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                    foreach my $pos (sort{$b <=> $a} keys %{$anno{$chromo}}) {
                        $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                        my @tmp=split("\t",$anno{$chromo}{$pos});
                        for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                        if ($gap eq 0) {$gap = $tmp[2];}
                        else {$gap=$tmp[1] - $gap;}
                    }
                    my $first=-1;
                    my $last=-1;
                    for(my $i=0; $i<$Len; $i++) {
                        if ($Cov[$i] eq 1) {
                            $cov1++;
                            if ($first eq -1) {$first = $i;}
                            $last=$i;
                        }
                    }
                    $cov2=$last - $first + 1;
                    print OUT21 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    print OUTUID $a[0],"\n";
                    if (($cov2 eq  $cov1) and ($cov2 >= $Len * $cutoff)) {
                        print OUT2S $a[0],"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    }
                }
                else {
                    # two segments on different strands
                    my $result=$Len."\t".$chromo;
                    my $cov1=0; 
                    my $cov2=0; 
                    my @Cov;
                    for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                    foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                        $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                        my @tmp=split("\t",$anno{$chromo}{$pos});
                        for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                        if ($gap eq 0) {$gap = $tmp[1];}
                        else {$gap=$tmp[1] - $gap;}
                    }
                    my $first=-1;
                    my $last=-1;
                    for(my $i=0; $i<$Len; $i++) {
                        if ($Cov[$i] eq 1) {
                            $cov1++;
                            if ($first eq -1) {$first = $i;}
                            $last=$i;
                        }
                    }
                    $cov2=$last - $first + 1;
                    print OUT22 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\t",$gap,"\n";
                }
            }
            # more than two-hitter, whassap ??
            elsif ($Chr{$chromo} > 2) {
                my $result=$Len."\t".$chromo;
                my $cov1=0;
                my $cov2=0; 
                my @Cov;
                for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                    $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                }
                my $first=-1;
                my $last=-1;
                for(my $i=0; $i<$Len; $i++) {
                    if ($Cov[$i] eq 1) {
                        $cov1++;
                        if ($first eq -1) {$first = $i;}
                        $last=$i;
                    }
                }
                $cov2=$last - $first + 1;
                print OUT4 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\n";
                print OUTUID $a[0],"\n";
            }
        }
        if (($stranded eq "+") or ($stranded eq "no")) {
            $strand=reverseStrand($strand);
            my @d=split(/\:/,$a[11]);
            my $info="";
            # $a[14] in bwa-0.7.3a : XP:Z:13,-91704253,49M52S,60,0;
            # $a[15] in bwa-0.7.15 : SA:Z:13,91704253,-,49M52S,60,0;
            for(my $i=11; $i<scalar(@a); $i++) {
                if ($a[$i]=~m/^XP:Z:/) {
                    my @b=split(/:/,$a[$i]);
                    my @c=split(/\;/,$b[2]);
                    if ($aFlag[-8] eq 1) {
                        $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
                    }
                    else {
                        for (my $j=0; $j<scalar(@c); $j++) {
                            my @tmpp=split(/\,/,$c[$j]);
                            my $tmpstrand=substr($tmpp[1],0,1);
                            my $tmppos=substr($tmpp[1],1);
                            $tmpp[1]=reverseStrand($tmpstrand).$tmppos;
                            $c[$j]=join(",",@tmpp);
                        }
                        $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
                    }
                    last;
                }
                elsif ($a[$i]=~m/^SA\:Z\:/) {
                    my @b=split(/:/,$a[$i]);
                    my @c=split(/\;/,$b[2]);
                    for(my $j=0; $j<scalar(@c); $j++) {
                        if (length($c[$j]) > 1) {
                            my @tmp=split(/\,/,$c[$j]);
                            if ($aFlag[-8] eq 1) { $c[$j]=join(",",$tmp[0],$tmp[2].$tmp[1],$tmp[3],$tmp[4],$tmp[5]); }
                            else { $c[$j]=join(",",$tmp[0],reverseStrand($tmp[2]).$tmp[1],$tmp[3],$tmp[4],$tmp[5]); }
                        }
                    }
                    $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
                    last;
                }
            }
            if ($info eq "") {
                $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]));
            }
            
            if (($a[2] eq "MT") or ($a[2] eq "M")) {
                print OUTMT join("\t","rc_".$a[0],$Len,$info),"\n";
                next;
            }
            my $indel=0;
            if(scalar(@a) <= 14){
                my @CIGAR_op=($a[5]=~m/[MHSID]/g); 
                my @CIGAR_va=($a[5]=~m/\d+/g);
                my $start=0;
                my $length=0;
                my $glength=0;
                my $Nr=scalar(@CIGAR_op);
                if ($strand eq "+") {
                    for(my $i=0; $i<$Nr; $i++){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$glength=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                else {
                    for(my $i=$Nr-1; $i>=0; $i--){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$glength+=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                print OUT1 join("\t","rc_".$a[0],$Len,$start,$length,$d[2],$a[2],$strand,$a[3],$a[4],$indel),"\n";
                print OUT join("\t","rc_".$a[0],$Len,$info),"\n";
                next;
            }
            print OUT join("\t","rc_".$a[0],$Len,$info),"\n";
            # process $info
            my %anno;
            my %Chr;
            my @p=split("\t",$info);
            if ($MAS > 0) {
                my $tmpmax=0;
                for(@p) { my @tmp=split(/\,/,$_); if ($tmp[3] > $tmpmax) {$tmpmax=$tmp[3]}; }
                if ($tmpmax < $MAS) { next; }
            }
            for(@p) {
                my @tmp=split(/\,/,$_);
                # ignore bad-quality alignments
                # if ($tmp[3] < $MAS) {next;}
                my $strandy="+";
                if ($tmp[1]=~m/^-/){
                    $strandy="-";
                    $tmp[1]=~s/^\-//;
                }
                else {
                    $tmp[1]=~s/^\+//;
                }
                my @CIGAR_op=($tmp[2]=~m/[MSID]/g); 
                my @CIGAR_va=($tmp[2]=~m/\d+/g);
                my $start=0;
                my $length=0;
                my $glength=0;
                my $Nr=scalar(@CIGAR_op);
                if ($strandy eq "+") {
                    for(my $i=0; $i<$Nr; $i++){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i]; }
                        elsif ($CIGAR_op[$i] eq "D") {$glength+=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                else {
                    for(my $i=$Nr-1; $i>=0; $i--){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; $glength+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$glength+=$CIGAR_va[$i]; $indel++;}
                        elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; $indel++;}
                        elsif (($CIGAR_op[$i] eq "S") or ($CIGAR_op[$i] eq "H")) {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                if (exists $anno{$tmp[0]}{$start}) {}
                else {
                    $Chr{$tmp[0]}++;
                    $anno{$tmp[0]}{$start}=join("\t",$length,$tmp[1],($tmp[1]+$glength),$strandy);
                }
            }
            # pick the most-hit chromosome
            my $chromo="";
            foreach my $id(sort{$Chr{$b} <=> $Chr{$a}} keys %Chr) {
                if ($Chr{$id} > 1) {$chromo=$id; last;}
            }
            # each part hit once
            if ($chromo eq "") {
                #my $result=$Len;
                #my $cov1=0; 
                #my $cov2=0; 
                #my @Cov;
                #for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                #foreach my $id(sort keys %anno) {
                #    foreach my $pos (sort{$a <=> $b} keys %{$anno{$id}}) {
                #        $result=$result."\t".$id."\t".$pos."\t".$anno{$id}{$pos};
                #        my @t=split("\t",$anno{$id}{$pos});
                #        for(my $i=0; $i<$t[0]; $i++) {$Cov[$i+$pos]=1;}
                #    }
                #}
                #my $first=-1;
                #my $last=-1;
                #for(my $i=0; $i<$Len; $i++) {
                #    if ($Cov[$i] eq 1) {
                #        $cov1++;
                #        if ($first eq -1) {$first = $i;}
                #        $last=$i;
                #    }
                #}
                #$cov2=$last - $first + 1;
                #if ($cov1 > 0) { print OUT2 "rc_".$a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\n"; }
            }
            # exactly two-hitter, interesting !!!
            elsif ($Chr{$chromo} eq 2) {
                my $flag=1;
                my $gap=0;
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    if ($tmp[3] eq "-") {$flag+=2;}
                }
                if ($flag eq 1) {
                    # both segments on the plus strand
                    my $result=$Len."\t".$chromo;
                    my $cov1=0; 
                    my $cov2=0; 
                    my @Cov;
                    for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                    foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                        $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                        my @tmp=split("\t",$anno{$chromo}{$pos});
                        for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                        if ($gap eq 0) {$gap = $tmp[2];}
                        else {$gap=$tmp[1] - $gap;}
                    }
                    my $first=-1;
                    my $last=-1;
                    for(my $i=0; $i<$Len; $i++) {
                        if ($Cov[$i] eq 1) {
                            $cov1++;
                            if ($first eq -1) {$first = $i;}
                            $last=$i;
                        }
                    }
                    $cov2=$last - $first + 1;
                    print OUT21 "rc_".$a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    print OUTUID "rc_".$a[0],"\n";
                    if (($cov2 eq  $cov1) and ($cov2 >= $Len * $cutoff)) {
                        print OUT2S "rc_".$a[0],"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    }
                }
                elsif ($flag eq 5) {
                    # both segments on the minus strand
                    my $result=$Len."\t".$chromo;
                    my $cov1=0; 
                    my $cov2=0; 
                    my @Cov;
                    for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                    foreach my $pos (sort{$b <=> $a} keys %{$anno{$chromo}}) {
                        $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                        my @tmp=split("\t",$anno{$chromo}{$pos});
                        for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                        if ($gap eq 0) {$gap = $tmp[2];}
                        else {$gap=$tmp[1] - $gap;}
                    }
                    my $first=-1;
                    my $last=-1;
                    for(my $i=0; $i<$Len; $i++) {
                        if ($Cov[$i] eq 1) {
                            $cov1++;
                            if ($first eq -1) {$first = $i;}
                            $last=$i;
                        }
                    }
                    $cov2=$last - $first + 1;
                    print OUT21 "rc_".$a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    print OUTUID "rc_".$a[0],"\n";
                    if (($cov2 eq  $cov1) and ($cov2 >= $Len * $cutoff)) {
                        print OUT2S "rc_".$a[0],"\t",$result,"\t",$indel,"\t",$gap,"\n";
                    }
                }
                else {
                    # two segments on different strands
                    my $result=$Len."\t".$chromo;
                    my $cov1=0; 
                    my $cov2=0; 
                    my @Cov;
                    for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                    foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                        $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                        my @tmp=split("\t",$anno{$chromo}{$pos});
                        for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                        if ($gap eq 0) {$gap = $tmp[1];}
                        else {$gap=$tmp[1] - $gap;}
                    }
                    my $first=-1;
                    my $last=-1;
                    for(my $i=0; $i<$Len; $i++) {
                        if ($Cov[$i] eq 1) {
                            $cov1++;
                            if ($first eq -1) {$first = $i;}
                            $last=$i;
                        }
                    }
                    $cov2=$last - $first + 1;
                    print OUT22 "rc_".$a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\t",$gap,"\n";
                }
            }
            # more than two-hitter, whassap ??
            elsif ($Chr{$chromo} > 2) {
                my $result=$Len."\t".$chromo;
                my $cov1=0;
                my $cov2=0; 
                my @Cov;
                for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                    $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                }
                my $first=-1;
                my $last=-1;
                for(my $i=0; $i<$Len; $i++) {
                    if ($Cov[$i] eq 1) {
                        $cov1++;
                        if ($first eq -1) {$first = $i;}
                        $last=$i;
                    }
                }
                $cov2=$last - $first + 1;
                print OUT4 "rc_".$a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$indel,"\n";
                print OUTUID $a[0],"\n";
            }
        }

    }
}
close IN;


open(OUTFLAG,">Step1_finished");
print OUTFLAG "Step1_finished\n";
close OUTFLAG;
