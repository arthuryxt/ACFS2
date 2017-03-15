#!/usr/bin/perl -w
use strict;
# convert circRNA structure file into raw_gtf format 
die "Usage: $0  \"output_basename\"   \"circRNA\"  \"split_exon_gtf\"   \"\(optional\) make combination default=0_or_1\"   \"\(optional\) insert_size default=150\"   \"\(optional\) extend N bases\"  \"\(optional\)do_NOT_fix_border_using_annotation==1_or_0\"  " if (@ARGV < 3);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];    
my $gtf=$ARGV[2];       
my $make_combi=0;
if (scalar(@ARGV) > 3) {$make_combi=$ARGV[3];}
my $libInsertSize=150;
if (scalar(@ARGV) > 4) {$libInsertSize=$ARGV[4];}
my $Extend=50;         
if (scalar(@ARGV) > 5) {$Extend=$ARGV[5];}
my $do_not_fix=1;           # fix the junctional border using annotation
if (scalar(@ARGV) > 6) {$do_not_fix=$ARGV[6];}
my $debug=0;
if (scalar(@ARGV) > 7) {$debug=$ARGV[7];}
my $command="rm -f Step4_CBR_finished";
system($command);

my %anno;
my %ExonNr;
if ($gtf ne "no") {
    open IN,$gtf;
    while(<IN>) {
        chomp;
        my @a=split("\t",$_);
        if ($a[2] eq "exon") {
            my @b=split(/\_\_\_/,$a[8]);
            $anno{$a[0]}{$a[6]}{$a[3]."\t".$a[4]}=join("\t",@a);
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
open OUTstat,">".$fileout.".stat";
my %usedcid;
my %cidstat;
while(<IN>) {
    chomp;
    if (m/^#/) { next; }
    my @a=split("\t",$_);
    if ($debug > 0) { print join("\t",@a),"\n";}
    my $left=$a[2] < $a[3] ? $a[2] : $a[3];
    my $right=$a[2] > $a[3] ? $a[2] : $a[3];
    my $rank++;
    if ((exists $anno{$a[1]}) and (exists $anno{$a[1]}{$a[20]})) {
        my %record;
        my %Gname;
        foreach my $info (keys %{$anno{$a[1]}{$a[20]}}) {
            my @b=split("\t",$anno{$a[1]}{$a[20]}{$info});
            if (($left <= $b[3]) and ($b[4] <= $right)) {
                $Gname{$b[7]}++;
                $record{$b[3]}=$b[4];
            }
        }
        foreach my $pos (sort{$a <=> $b} keys %record) {
            my $tmp=$record{$pos};
            delete $record{$pos};
            $record{$left}=$tmp;
            last;
        }
        foreach my $pos (sort{$b <=> $a} keys %record) {
            $record{$pos}=$right;
            last;
        }
        my $gene_name="na";
        foreach my $id (sort{$Gname{$b} <=> $Gname{$a}} keys %Gname) { $gene_name=$id; last; }
        
        my @ExonL;
        my @ExonR;
        my @ExonLen;
        my $Nr=0;
        my $cid=$a[0]."_1";
        foreach my $pos (sort{$a <=> $b} keys %record) { $Nr++; if ($debug > 0) { print join("\t",$Nr, $pos, $record{$pos}),"\n";}  }
        for(my $i=2; $i<=$Nr; $i++){ $cid=$cid."|".$i; }
        if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
        $cidstat{$cid}=1;
        if ($Nr > 0) {            
            if ($a[20] eq "+") {
                my $count=1;
                foreach my $pos (sort {$a <=> $b} keys %record) {
                    my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$count\"");
                    print OUT join("\t",$a[1],"split","exon",$pos,$record{$pos},"na",$a[20],$gene_name,$info),"\n";
                    $ExonL[$count-1]=$pos;
                    $ExonR[$count-1]=$record{$pos};
                    $ExonLen[$count-1]=$record{$pos}-$pos;
                    $count++;
                }
            }
            else {
                my $count=$Nr;
                foreach my $pos (sort {$a <=> $b} keys %record) {
                    my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$count\"");
                    print OUT join("\t",$a[1],"split","exon",$pos,$record{$pos},"na",$a[20],$gene_name,$info),"\n";
                    $ExonL[$Nr - $count]=$pos;
                    $ExonR[$Nr - $count]=$record{$pos};
                    $ExonLen[$Nr - $count]=$record{$pos}-$pos;
                    $count--;
                }
            }
            print OUTrefFlat join("\t",$gene_name,$cid,$a[1],$a[20],$left,$right+1,$left,$left,($Nr),join(",",@ExonL),join(",",@ExonR),$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
            
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
        
                # generate all possible alternative splicing events
                my $apase=1;
                for(my $k=0; $k<($Nr-3); $k++){ $apase=10*$apase + 1; }
                my $patterns=bin2dec($apase);
                for(my $k=0; $k<$patterns; $k++) {
                    my @usage=split("",dec2bin($k,$Nr-2));
                    if ($debug > 0) { print join("\t",$apase,@usage),"\n";} 
                    my $cid=$a[0]."_1";
                    for(my $j=0; $j<scalar(@usage); $j++) { if(($usage[$j] eq 1) and ($Reachable[$j+1] eq 1)){ $cid=$cid."|".(2+$j); }   }
                    $cid=$cid."|".$Nr;
                    if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
                    if ($debug > 0) { print $cid,"\n"; }
                    my $exonL=$ExonL[0];
                    my $exonR=$ExonR[0];
                    my $exoncnt=1;
                    my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$exoncnt\"");
                    print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[0],$ExonR[0],"na",$a[20],$gene_name,$info),"\n";
                    for(my $j=0; $j<scalar(@usage); $j++) {
                        if(($usage[$j] eq 1) and ($Reachable[$j+1] eq 1)){
                            $exonL=$exonL.",".$ExonL[1+$j];
                            $exonR=$exonR.",".$ExonR[1+$j];
                            $exoncnt++;
                            $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$exoncnt\"");
                            print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[1+$j],$ExonR[1+$j],"na",$a[20],$gene_name,$info),"\n";
                        }
                    }
                    $exonL=$exonL.",".$ExonL[$Nr-1];
                    $exonR=$exonR.",".$ExonR[$Nr-1];
                    $exoncnt++;
                    $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"$exoncnt\"");
                    print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[$Nr-1],$ExonR[$Nr-1],"na",$a[20],$gene_name,$info),"\n";
                    print OUT1refFlat join("\t",$a[6],$cid,$a[1],$a[20],$a[3],$a[2]+1,$a[3],$a[3],$exoncnt,$exonL,$exonR,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
                }
            }
        }
        else {
            my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"1\"");
            print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],$gene_name,$info),"\n";
            print OUTrefFlat join("\t","na",$cid,$a[1],$a[20],$left,$right+1,$left,$left,1,$left,$right+1,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
        }
    }
    else {
        my $cid=$a[0]."_1";
        if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
        $cidstat{$cid}=1;
        my $gene_name="na";
        my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"1\"");
        print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],$gene_name,$info),"\n";
        print OUTrefFlat join("\t","na",$cid,$a[1],$a[20],$left,$right+1,$left,$left,1,$left,$right+1,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
    }
}

foreach my $id (keys %cidstat) { print OUTstat $id,"\n"; }

open(OUTFLAG,">Step4_CBR_finished");
print OUTFLAG "Step4_CBR_finished\n";
close OUTFLAG;

