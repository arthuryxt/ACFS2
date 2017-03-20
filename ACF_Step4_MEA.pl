#!/usr/bin/perl -w
use strict;
# convert circRNA structure file into raw_gtf format 
die "Usage: $0  \"output_basename\"   \"circRNA\"  \"split_exon_gtf\"   \"\(optional\) make combination default=0_or_1\"   \"\(optional\) insert_size default=150\"    \"\(optional\) extend N bases\"  \"\(optional\)do_NOT_fix_border_using_annotation==1_or_0\"    \"\(optional\)debug\"  " if (@ARGV < 3);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];    
my $gtf=$ARGV[2];       
my $make_combi=0;
if (scalar(@ARGV) > 3) {$make_combi=$ARGV[3];}
my $libInsertSize=150;
if (scalar(@ARGV) > 4) {$libInsertSize=$ARGV[4];}
my $Extend=50;         # 100nt by default.
if (scalar(@ARGV) > 5) {$Extend=$ARGV[5];}
my $do_not_fix=1;           # fix the junctional border using annotation
if (scalar(@ARGV) > 6) {$do_not_fix=$ARGV[6];}
my $debug=0;
if (scalar(@ARGV) > 7) {$debug=$ARGV[7];}
my $command="rm -f Step4_MEA_finished";
system($command);

my %anno;
my %ExonNr;
if ($gtf ne "no") {
    open IN,$gtf;
    while(<IN>) {
        chomp;
        my @a=split("\t",$_);
        if ($a[2] eq "exon") {
            #if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
            #if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
            my @b=split(/\_\_\_/,$a[8]);
            #     ENSG   exon_number
            $anno{$b[0]}{$b[1]}=join("\t",@a);
            $ExonNr{$b[0]}=$b[2];
        }
    }
    close IN;
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    my $cnt= shift;
    #$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros return $str;
    my $decode = substr($str,0-$cnt);
    return $decode;
}
sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

open IN,$filein;
open OUT,">".$fileout.".gtf";               # close zero-based. This means that the first 100 bases of a chromosome are represented as [0,99]
open OUTrefFlat,">".$fileout.".refFlat";    # half-open zero-based. This means that the first 100 bases of a chromosome are represented as [0,100), i.e. 0-99. The second 100 bases are represented as [100,200), i.e. 100-199.
open OUT1gtf,">".$fileout.".ext.gtf";
open OUT1refFlat,">".$fileout.".ext.refFlat";
open OUT2,">".$fileout.".gtf_2G";
open OUT3,">".$fileout.".err";
open OUTstat,">".$fileout.".stat";
my %usedcid;
my %cidstat;    # 1 for circRNAs identified by reads
while(<IN>) {
    chomp;
    if (m/^#/) { next; }
    my @a=split("\t",$_);
    #if ($a[1]=~m/^chromosome/i) {$a[1]=~s/chromosome//i;}
    #if ($a[1]=~m/^chr/i) {$a[1]=~s/chr//i;}
    my @left=split(/\_\_\_/,$a[5]);
    my @right=split(/\_\_\_/,$a[11]);
    if ($left[0] ne $right[0]) { print OUT2 join("\t",@a),"\n"; next;}
    #if ($a[6] ne $a[12]) { print OUT2 join("\t",@a),"\n"; next;}
    my $start=$left[1] < $right[1] ? $left[1] : $right[1];
    my $end=$left[1] > $right[1] ? $left[1] : $right[1];
    my $Nr=$end - $start +1;
    my @ExonL;
    my @ExonR;
    my @ExonLen;
    my $strand=$a[7];
    if ($a[7] ne $a[13]) { print OUT3 "discrepancy in two-part-annotation\t",join("\t",@a),"\n"; next; }
    if ($a[7] ne $a[20]) { print OUT3 "discrepancy in annotation and read strandness\t",join("\t",@a),"\n"; next; }
    if ($a[2] ne $a[10]) { print OUT3 "discrepancy in left-border\t",join("\t",@a),"\n"; }
    if ($a[3] ne $a[15]) { print OUT3 "discrepancy in right-border\t",join("\t",@a),"\n"; }
    my $cid=$a[0]."_".$start;
    for(my $i=$start+1; $i<=$end; $i++){ $cid=$cid."|".$i; }
    if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
    $cidstat{$cid}=1;
    if ($debug > 0) { print $cid,"\n"; }
    my $biotype=$a[8];
    my $gene_name=$a[6];
    my $gene_name2=$left[0];
    for(my $i=1; $i<=$Nr; $i++) {
        my @b=split("\t",$anno{$left[0]}{$i+$start-1});
        if ($debug > 0) {print join("\t",@b),"\n";}
        my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$i\"");
        my $tmp_s=$b[3];
        my $tmp_e=$b[4];
        if (($tmp_s <= $a[3]) and ($a[3] <= $tmp_e)) { $tmp_s=$a[3]; }
        elsif (($tmp_s <= $a[2]) and ($a[2] <= $tmp_e)) { $tmp_e=$a[2]; }
        if ($do_not_fix eq 1) {
            if ($i eq 1) { $tmp_s=$a[3];}
            elsif($i eq $Nr) { $tmp_e=$a[2]; }
        }
        if ($debug > 0) {print join("\t",$b[0],$b[1],$b[2],$tmp_s,$tmp_e,$b[5],$b[6],$b[7]),"\n";}
        print OUT join("\t",$b[0],$b[1],$b[2],$tmp_s,$tmp_e,$b[5],$b[6],$b[7],$info),"\n";
        if ($strand eq "+") {
            $ExonL[$i-1]=$tmp_s;
            $ExonR[$i-1]=$tmp_e;
            $ExonLen[$i-1]=$tmp_e+1-$tmp_s;
        }
        else {
            $ExonL[$Nr-$i]=$tmp_s;
            $ExonR[$Nr-$i]=$tmp_e;
            $ExonLen[$Nr-$i]=$tmp_e+1-$tmp_s;
        }
    }
    if (($do_not_fix eq 1) and (($ExonL[0] ne $a[3]) or ($ExonR[$Nr-1] ne ($a[2])))) {
        print OUT3 "fix-border\t",join("\t",@a),"\n";
        $ExonL[0]=$a[3];
        $ExonR[$Nr-1]=$a[2];
        $ExonLen[0]=$ExonR[0]+1-$ExonL[0];
        $ExonLen[$Nr-1]=$ExonR[$Nr-1]+1-$ExonL[$Nr-1];
    }
    print OUTrefFlat join("\t",$a[6],$cid,$a[1],$strand,$a[3],$a[2]+1,$a[3],$a[3],$Nr,join(",",@ExonL),join(",",@ExonR),$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
    if (($make_combi > 0) and ($Nr > 2)) {
        # generate all possible alternative splicing events that can be determined given $libInsertSize
        my @Reachable;
        for(my $i=0; $i<$Nr; $i++){ $Reachable[$i]=0; if ($debug > 0) {print $i,"\t",$ExonLen[$i],"\n";} }
        if ($debug > 0) { my $stdout=$Reachable[0]; for(my $i=1; $i<$Nr; $i++){ $stdout=$stdout.$Reachable[$i]; }  print $stdout,"\n";}
        my $posLeft=0;
        my $reached=0;
        for(my $i=0; $i<$Nr; $i++) { $reached=$reached+$ExonLen[$i]; if($reached > $libInsertSize){ $posLeft=$i; last; } }
        for(my $i=0; $i<=$posLeft; $i++) { $Reachable[$i]=1; }
        if ($debug > 0) { my $stdout=$Reachable[0]; for(my $i=1; $i<$Nr; $i++){ $stdout=$stdout.$Reachable[$i]; }  print $stdout,"\n";}
        my $posRight=0;
        $reached=0;
        for(my $i=0; $i<$Nr; $i++) { $reached=$reached+$ExonLen[$Nr-1-$i]; if($reached > $libInsertSize){ $posRight=$Nr-1-$i; last; } }
        for(my $i=$Nr-1; $i>=$posRight; $i--) { $Reachable[$i]=1; }
        if ($debug > 0) { my $stdout=$Reachable[0]; for(my $i=1; $i<$Nr; $i++){ $stdout=$stdout.$Reachable[$i]; }  print $stdout,"\n";}
        my $apase=1;
        for(my $k=0; $k<($Nr-3); $k++){ $apase=10*$apase + 1; }
        my $patterns=bin2dec($apase);
        for(my $k=0; $k<$patterns; $k++) {
            my @usage=split("",dec2bin($k,$Nr-2));
            my $cid=$a[0]."_".$start;
            for(my $j=0; $j<scalar(@usage); $j++) { if(($usage[$j] eq 1) and ($Reachable[1+$j] eq 1)){ $cid=$cid."|".($start+1+$j); }   }
            $cid=$cid."|".$end;
            if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
            if ($debug > 0) { print $cid,"\t",join("",@usage),"\n"; }
            my $exonL=$ExonL[0]-1;
            my $exonR=$ExonR[0];
            if ($strand eq "-") { $exonL=$ExonL[$Nr-1]-1; $exonR=$ExonR[$Nr-1]; }
            my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"1\"");
            print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[0],$ExonR[0],$biotype,$a[20],$gene_name2,$info),"\n";
            my $exoncnt=1;
            for(my $j=0; $j<scalar(@usage); $j++) {
                if(($usage[$j] eq 1) and ($Reachable[1+$j] eq 1)){
                    if ($strand eq "+") { $exonL=$exonL.",".($ExonL[1+$j]-1); $exonR=$exonR.",".$ExonR[1+$j]; }
                    else { $exonL=($ExonL[1+$j]-1).",".$exonL; $exonR=$ExonR[1+$j].",".$exonR; } 
                    $exoncnt++;
                    $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$exoncnt\"");
                    print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[1+$j],$ExonR[1+$j],$biotype,$a[20],$gene_name2,$info),"\n";
                }
            }
            if ($strand eq "+") { $exonL=$exonL.",".($ExonL[$Nr-1]-1); $exonR=$exonR.",".$ExonR[$Nr-1]; }
            else { $exonL=($ExonL[0]-1).",".$exonL; $exonR=$ExonR[0].",".$exonR; } 
            $exoncnt++;
            $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$exoncnt\"");
            print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[$Nr-1],$ExonR[$Nr-1],$biotype,$a[20],$gene_name2,$info),"\n";
            print OUT1refFlat join("\t",$a[6],$cid,$a[1],$strand,$a[3]-1,$a[2],$a[3]-1,$a[3]-1,$exoncnt,$exonL,$exonR,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
        }
    }
    
    my $cumlen=0;
    # check for extension on left side
    my $start_move=0;
    for (my $i=$start-1; $i>0; $i--) {
        my @tmp=split("\t",$anno{$left[0]}{$i});
        my $tmplen=$tmp[4] - $tmp[3];
        $cumlen+=$tmplen;
        if ($debug > 0) { print join("\t","start_pos:",$i,$cumlen,@tmp),"\n";}
        if ($cumlen < $Extend) {$start_move++;}
        else {last;}
    }
    $cumlen=0;
    # check for extension on right side
    my $end_move=0;
    for (my $i=$end+1; $i<=$ExonNr{$left[0]}; $i++) {
        my @tmp=split("\t",$anno{$left[0]}{$i});
        my $tmplen=$tmp[4] - $tmp[3];
        $cumlen+=$tmplen;
        if ($debug > 0) { print join("\t","end_pos:",$i,$cumlen,@tmp),"\n";}
        if ($cumlen < $Extend) {$end_move++;}
        else {last;}
    }
    if ($debug > 0) { print join("\t",$start,$end,$start_move,$end_move,$Nr),"\n";}
    for (my $i=0; $i<=$start_move; $i++) {
        for (my $j=0; $j<=$end_move; $j++) {
            # output all possible combinations, except the original ones
            if (($i eq 0) and ($j eq 0)) {next;}
            my $nstart=$start-$i;
            my $nend=$end+$j;
            my @Anstart=split("\t",$anno{$left[0]}{$nstart});
            my @Anend=split("\t",$anno{$left[0]}{$nend});
            my $nleftmost=$Anstart[3] < $Anend[3] ? $Anstart[3] : $Anend[3];
            my $nrightmost=$Anstart[4] > $Anend[4] ? $Anstart[4] : $Anend[4];
            
            # output the extended version
            my $exonL="";
            my $exonR="";
            my $cid0=join("_",$a[1],$nrightmost,$nleftmost,$Anstart[6].abs($nrightmost-$nleftmost));
            my $cid=$cid0."_".$nstart;
            for(my $tmp=$nstart+1; $tmp<=$nend; $tmp++) { $cid=$cid."|".$tmp; }
            if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
            $cidstat{$cid}=1;
            if ($debug > 0) { print $cid,"\n",$i,"\t",$j,"\n";}
            for(my $k=1; $k<=($nend - $nstart +1); $k++) {
                my @b=split("\t",$anno{$left[0]}{$k+$nstart-1});
                my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$k\"");
                print OUT join("\t",$b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$info),"\n";
                if ($strand eq "+") {
                    if ($exonL eq "") { $exonL=$b[3]; } else { $exonL=$exonL.",".$b[3]; }
                    if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$exonR.",".$b[4]; }
                }
                else {
                    if ($exonL eq "") { $exonL=$b[3]; } else { $exonL=$b[3].",".$exonL; }
                    if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$b[4].",".$exonR; }
                }
            }
            print OUTrefFlat join("\t",$a[6],$cid,$a[1],$strand,$nleftmost,$nrightmost,$nleftmost,$nleftmost,($nend - $nstart +1),$exonL,$exonR,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";

            # generate all possible alternative splicing events that can be determined given $libInsertSize
            if (($make_combi > 0) and ($Nr > 2)) {
                my @Reachable;
                for(my $i=0; $i<$Nr; $i++){ $Reachable[$i]=0; if ($debug > 0) {print $i,"\t",$ExonLen[$i],"\n";} }
                if ($debug > 0) { my $stdout=$Reachable[0]; for(my $i=1; $i<$Nr; $i++){ $stdout=$stdout.$Reachable[$i]; }  print $stdout,"\n";}
                my $posLeft=0;
                my $reached=0;
                for(my $i=0; $i<$Nr; $i++) { $reached=$reached+$ExonLen[$i]; if($reached > $libInsertSize){ $posLeft=$i; last; } }
                for(my $i=0; $i<=$posLeft; $i++) { $Reachable[$i]=1; }
                if ($debug > 0) { my $stdout=$Reachable[0]; for(my $i=1; $i<$Nr; $i++){ $stdout=$stdout.$Reachable[$i]; }  print $stdout,"\n";}
                my $posRight=0;
                $reached=0;
                for(my $i=0; $i<$Nr; $i++) { $reached=$reached+$ExonLen[$Nr-1-$i]; if($reached > $libInsertSize){ $posRight=$Nr-1-$i; last; } }
                for(my $i=$Nr-1; $i>=$posRight; $i--) { $Reachable[$i]=1; }
                if ($debug > 0) { my $stdout=$Reachable[0]; for(my $i=1; $i<$Nr; $i++){ $stdout=$stdout.$Reachable[$i]; }  print $stdout,"\n";}
            
                my $apase=1;
                for(my $k=0; $k<($Nr-3); $k++){ $apase=10*$apase + 1; }
                my $patterns=bin2dec($apase);
                for(my $k=0; $k<$patterns; $k++) {
                    my @usage=split("",dec2bin($k,$Nr-2));
                    my $cid=$cid0."_".($start-$i);
                    for(my $tmp=$i-1; $tmp>=0; $tmp--) { $cid=$cid."|".($start-$tmp); }
                    for(my $tmp=0; $tmp<scalar(@usage); $tmp++) { if(($usage[$tmp] eq 1) and ($Reachable[1+$tmp] eq 1)){ $cid=$cid."|".($start+1+$tmp); }   }
                    $cid=$cid."|".$end;
                    for(my $tmp=1; $tmp<=$j; $tmp++) { $cid=$cid."|".($end+$tmp); }
                    if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
                    if ($debug > 0) { print $cid,"\n"; }
                    my $exonL="";
                    my $exonR="";
                    my $exoncnt=0;
                    for(my $tmp=$i; $tmp>=0; $tmp--) {
                        my @b=split("\t",$anno{$left[0]}{$start-$tmp});
                        my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$exoncnt\"");
                        print OUT1gtf join("\t",$b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$info),"\n";
                        if ($strand eq "+") {
                            if ($exonL eq "") { $exonL=$b[3]-1; } else { $exonL=$exonL.",".($b[3]-1); }
                            if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$exonR.",".$b[4]; }
                        }
                        else {
                            if ($exonL eq "") { $exonL=$b[3]-1; } else { $exonL=($b[3]-1).",".$exonL; }
                            if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$b[4].",".$exonR; }
                        }
                        $exoncnt++;
                    }
                    for(my $tmp=0; $tmp<scalar(@usage); $tmp++) {
                        if(($usage[$tmp] eq 1) and ($Reachable[1+$tmp] eq 1)){
                            my @b=split("\t",$anno{$left[0]}{$start+$tmp});
                            my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$exoncnt\"");
                            print OUT1gtf join("\t",$b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$info),"\n";
                            if ($strand eq "+") {
                                if ($exonL eq "") { $exonL=$b[3]-1; } else { $exonL=$exonL.",".($b[3]-1); }
                                if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$exonR.",".$b[4]; }
                            }
                            else {
                                if ($exonL eq "") { $exonL=$b[3]-1; } else { $exonL=($b[3]-1).",".$exonL; }
                                if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$b[4].",".$exonR; }
                            }
                            $exoncnt++;
                        }
                    }
                    for(my $tmp=0; $tmp<=$j; $tmp++) {
                        my @b=split("\t",$anno{$left[0]}{$end+$tmp});
                        my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$exoncnt\"");
                        print OUT1gtf join("\t",$b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$info),"\n";
                        if ($strand eq "+") {
                            if ($exonL eq "") { $exonL=$b[3]-1; } else { $exonL=$exonL.",".($b[3]-1); }
                            if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$exonR.",".$b[4]; }
                        }
                        else {
                            if ($exonL eq "") { $exonL=$b[3]-1; } else { $exonL=($b[3]-1).",".$exonL; }
                            if ($exonR eq "") { $exonR=$b[4]; } else { $exonR=$b[4].",".$exonR; }
                        }
                        $exoncnt++;
                    }
                    print OUT1refFlat join("\t",$a[6],$cid,$a[1],$strand,$nleftmost-1,$nrightmost,$nleftmost-1,$nleftmost-1,$exoncnt,$exonL,$exonR,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
                }
            }
        }
    }
}

foreach my $id (keys %cidstat) { print OUTstat $id,"\n"; }

open(OUTFLAG,">Step4_MEA_finished");
print OUTFLAG "Step4_MEA_finished\n";
close OUTFLAG;

