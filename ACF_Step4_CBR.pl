#!/usr/bin/perl -w
use strict;
use Math::BigInt;
# convert circRNA structure file into raw_gtf format 
die "Usage: $0  \"output_basename\"   \"circRNA\"  \"split_exon_gtf\"   \"\(optional\) make combination default=0_or_1\"  \"\(optional\) max_AS default=10\"   \"\(optional\) insert_size default=150\"   \"\(optional\) extend N bases\"  \"\(optional\)do_fix_border_using_annotation==1_or_0\"  \"\(optional\)debug\" " if (@ARGV < 3);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];    
my $gtf=$ARGV[2];       
my $make_combi=0;
if (scalar(@ARGV) > 3) {$make_combi=$ARGV[3];}
my $max_AS=10;
if (scalar(@ARGV) > 4) {$max_AS=$ARGV[4];}
my $libInsertSize=150;
if (scalar(@ARGV) > 5) {$libInsertSize=$ARGV[5];}
my $Extend=50;          # not used anyway
if (scalar(@ARGV) > 6) {$Extend=$ARGV[6];}
my $do_fix=0;           # will NOT fix the junctional border using annotation
if (scalar(@ARGV) > 7) {$do_fix=$ARGV[7];}
$do_fix=0;
my $debug=0;
if (scalar(@ARGV) > 8) {$debug=$ARGV[8];}
my $command="rm -f Step4_CBR_finished";
system($command);

my %usedcid;
my %cidstat;
my %anno;
my %anno3;
my %ExonNr;
if ($gtf ne "no") {
    open IN,$gtf;
    while(<IN>) {
        chomp;
        my @a=split("\t",$_);
        if ($a[2] eq "exon") {
            my @b=split(/\_\_\_/,$a[8]);
            # 1	split	exon	1330895	1330917	protein_coding	-	CCNL2	ENSG00000221978___45___58___01
            $anno3{$a[0]}{$a[6]}{$a[3]."\t".$a[4]}=join("\t",@a);
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

sub bigbin2dec {
  my $bin = shift;
  return Math::BigInt->new("0b$bin");
}

sub bigdec2bin {
  my $dec = shift;
  my $cnt = shift;
  my $mybin = Math::BigInt->new($dec);
  my $mybin2 = substr($mybin->as_bin(), 2);
  my $mybin3 = "";
  if (length($mybin2) >= $cnt) {
    $mybin3=substr($mybin2,0-$cnt);
  }
  else {
    my $heading=0;
    for(my $i=1; $i<($cnt-length($mybin2)); $i++) {$heading=$heading."0";}
    $mybin3=$heading.$mybin2;
  }
  return $mybin3;
}


sub joinexon($$$$$$$) {
    my $Gene=$_[0];
    my $mleft=$_[1];
    my $mright=$_[2];
    my @exoncnt=split("",join("","1",$_[3],"1"));
    my $mycid=$_[4];
    my $myi=$_[5];
    my $myj=$_[6];
    if ($debug > 0) {print join("\t",@_),"\n";}
    my @exonlen;
    for(my $i=$mleft; $i<=$mright; $i++) {
        if ($exoncnt[$i-$mleft] eq 1) {
            my @a=split("\t",$anno{$Gene}{$i});
            $exonlen[$i-$mleft]=$a[4]-$a[3]+1;
        }
        else { $exonlen[$i-$mleft]=0; }
    }
    my @FLeft;
    for(my $i=$mleft; $i<=$mright; $i++) { $FLeft[$i-$mleft]=0; }
    my $length=0;
    for(my $i=$mleft; $i<=$mright; $i++) {
        if($length < $libInsertSize) {
            if ($exoncnt[$i-$mleft] eq 1) {$FLeft[$i-$mleft]=1; $length+=$exonlen[$i-$mleft]; if($debug > 0){ print "lef\t",$i-$mleft,"\t",$exonlen[$i-$mleft],"\t",$length,"\n"; } }
        }else{last;}
    }
    my @FRight;
    for(my $i=$mleft; $i<=$mright; $i++) { $FRight[$i-$mleft]=0; }
    $length=0;
    for(my $i=$mright; $i>=$mleft; $i--) {
        if($length < $libInsertSize) {
            if ($exoncnt[$i-$mleft] eq 1) {$FRight[$i-$mleft]=1; $length+=$exonlen[$i-$mleft]; if($debug > 0){ print "rig\t",$i-$mleft,"\t",$exonlen[$i-$mleft],"\t",$length,"\n"; } }
        }else{last;}
    }
    if ($debug > 0) { print "lef result \t= ",join("",@FLeft),"\n"; }
    if ($debug > 0) { print "rig result \t= ",join("",@FRight),"\n"; }
    my $tmpid=$mleft;
    for(my $i=$mleft+1; $i<=$mright; $i++) {
        if (($FLeft[$i-$mleft] eq 1) or ($FRight[$i-$mleft] eq 1)) { $tmpid=$tmpid."|".$i; }
    }
    if ($myi > 0) {
        for(my $i=1; $i<=$myi; $i++) { $tmpid=($mleft-$i)."|".$tmpid; }
    }
    if ($myj > 0) {
        for(my $i=1; $i<=$myj; $i++) { $tmpid=$tmpid."|".($mright+$i) }
    }
    my $cid=$mycid."_".$tmpid;
    if (exists $usedcid{$cid}) { return(-1); } 
    else{
    my $result="";
    my $lastl="";
    my $lastr="";
    my $lastexon="";
    for(my $i=$mleft; $i<=$mright; $i++) {
        if ($result eq -1) { last; }
        if ($exoncnt[$i-$mleft] eq 0) { $result=$result."0"; next; }
        if (($FLeft[$i-$mleft] eq 0) and ($FRight[$i-$mleft] eq 0)) { $result=$result."0"; $lastl=""; next; }
        my @a=split("\t",$anno{$Gene}{$i});
        my @b=split(/\_\_\_/,$a[8]);
        my $lb=substr($b[3],0,1);
        my $rb=substr($b[3],1,1);
        if ($lastl eq "") { # the first exon is guaranteed to be included
            if ($result eq ""){$result=1;}else{$result=$result."1";} $lastl=1; $lastr=$rb; $lastexon=$i;
        }
        else {
            my @posA=split("\t",$anno{$Gene}{$lastexon});
            if (($lastl eq 1) and ($lastr eq 0)) {
                if (($a[3] - $posA[4]) <= 1) { $result=$result."1"; $lastr=$rb; $lastexon=$i; }
                else { $result=-1; }
            }
            elsif (($lastl eq 1) and ($lastr eq 1)) {
                if ($lb eq 0) {
                    if (($a[3] - $posA[4]) <= 1) { $result=$result."1"; $lastr=$rb; $lastexon=$i; }
                    else { $result=-1; }
                }
                elsif ($lb eq 1)  {
                    $result=$result."1"; $lastl=$lb; $lastr=$rb; $lastexon=$i;
                }
                else  {
                    $result=-1;
                }
            }
            elsif (($lastl eq 1) and ($lastr eq 2)) {
                if ($lb eq 0) {
                    if (($a[3] - $posA[4]) <= 1) { $result=$result."1"; $lastr=$rb; $lastexon=$i; }
                    else { $result=-1; }
                }
                #elsif ($lb eq 1)  {
                #    $result=-1;
                #}
                else  {
                    $result=-1; 
                }
            }
        }
    }
    if ($debug > 0) { print "raw result \t= ",$result,"\n"; }
    if ($result ne -1) {
        my @toUse=split("",$result);
        
        for(my $i=$mleft; $i<=$mright; $i++) {
            if (($FLeft[$i-$mleft] eq 0) and ($FRight[$i-$mleft] eq 0)) { $toUse[$i-$mleft]=0; }
        }
        $result=join("",@toUse);
        if ($debug > 0) { print "final result \t= ",$result,"\n"; }
        return(substr($result,1,(length($result)-2)));
    }
    else { return($result); }
    }
}

open IN,$filein;
open OUT,">".$fileout.".gtf";               # close zero-based. This means that the first 100 bases of a chromosome are represented as [0,99]
open OUTrefFlat,">".$fileout.".refFlat";    # half-open zero-based. This means that the first 100 bases of a chromosome are represented as [0,100), i.e. 0-99. The second 100 bases are represented as [100,200), i.e. 100-199.
open OUT1gtf,">".$fileout.".ext.gtf";
open OUT1refFlat,">".$fileout.".ext.refFlat";
open OUTstat,">".$fileout.".stat";
while(<IN>) {
    chomp;
    if (m/^#/) { next; }
    my @a=split("\t",$_);
    my @tmp=split(/\_/,$a[0]);
    my $strand="+";
    if (substr($tmp[3],0,1) eq "-") { $strand="-";}
    if ($debug > 0) { print join("\t",@a),"\n";}
    my $left=$a[2] < $a[3] ? $a[2] : $a[3];
    my $right=$a[2] > $a[3] ? $a[2] : $a[3];
    if ((exists $anno3{$a[1]}) and (exists $anno3{$a[1]}{$a[20]})) {
        my %record;
        my %Gid;
        my %Gname;
        foreach my $info (keys %{$anno3{$a[1]}{$a[20]}}) {
            my @b=split("\t",$anno3{$a[1]}{$a[20]}{$info});
            if (($left <= $b[3]) and ($b[4] <= $right)) {
                my @c=split(/\_\_\_/,$b[8]);
                $Gid{$c[0]}++;
                $Gname{$c[0]}=$b[7];
                $record{$b[3]}=$b[4];
                if ($debug > 0) { print join("\t",@b),"\n"; }
            }
        }
        my $gene_name="na";
        my $gene_name2="na";
        foreach my $id (sort{$Gid{$b} <=> $Gid{$a}} keys %Gid) { $gene_name=$id; $gene_name2=$Gname{$gene_name}; last; }
        if ($debug > 0) { print join("\t",$gene_name,$gene_name2),"\n"; }
        if ($gene_name eq "na") {
            my $cid=$a[0]."_1";
            if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
            $cidstat{$cid}=1;
            my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"1\"");
            print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],$gene_name,$info),"\n";
            print OUTrefFlat join("\t","na",$cid,$a[1],$a[20],$left-1,$right,$left-1,$left-1,1,$left-1,$right,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
        }
        else {
            my $start=-1;
            my $end=-1;
            foreach my $i (sort{$a <=> $b} keys %{$anno{$gene_name}}) {
                my @b=split("\t",$anno{$gene_name}{$i});
                if ($b[4] < $left) { $start=$i; }
                else { $start=$i; last; }
            }
            foreach my $i (sort{$a <=> $b} keys %{$anno{$gene_name}}) {
                my @b=split("\t",$anno{$gene_name}{$i});
                if ($b[3] < $right) { $end=$i; }
                else { last; }
            }
            if ($debug > 0) { print join("\t",$start,$end),"\n"; }
            if (($end > 1) and ($start > 0) and ($end - $start) > 0) {
                my $cid=$a[0]."_".$start;
                for(my $i=$start+1; $i<=$end; $i++){ $cid=$cid."|".$i; }
                if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
                $cidstat{$cid}=1;
                if ($debug > 0) { print $cid,"\n"; }
                my @ExonL;
                my @ExonR;
                my @ExonLen;
                my $exoncnt=0;
                my $exonL=-1;
                my $exonR=-1;
                for(my $i=$start; $i<=$end; $i++) {
                    my @b=split("\t",$anno{$gene_name}{$i});
                    if ($debug > 0) {print join("\t",@b),"\n";}
                    $exoncnt++;
                    my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$exoncnt\"");
                    my $tmp_s=$b[3];
                    my $tmp_e=$b[4];
                    if ($i eq $start) { $tmp_s=$a[3]; }
                    if ($i eq $end) { $tmp_e=$a[2]; }
                    if ($debug > 0) {print join("\t",$b[0],$b[1],$b[2],$tmp_s,$tmp_e,$b[5],$b[6],$b[7]),"\n";}
                    print OUT join("\t",$b[0],$b[1],$b[2],$tmp_s,$tmp_e,$b[5],$b[6],$b[7],$info),"\n";
                    $ExonL[$i-$start]=$tmp_s;
                    $ExonR[$i-$start]=$tmp_e;
                    $ExonLen[$i-$start]=$tmp_e+1-$tmp_s;
                    if ($exonL eq -1) { $exonL=$tmp_s-1; $exonR=$tmp_e }
                    else { $exonL=$exonL.",".($tmp_s-1); $exonR=$exonR.",".$tmp_e; }
                }
                my $Nr=$end-$start+1;
                print OUTrefFlat join("\t",$a[6],$cid,$a[1],$strand,$a[3]-1,$a[2],$a[3]-1,$a[3]-1,$Nr,$exonL,$exonR,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
                
                if (($make_combi > 0) and ($Nr > 2) and ($max_AS >= $Nr)) {
                    my $apase="1";
                    for(my $k=0; $k<($Nr-3); $k++){ $apase=$apase."1"; }
                    my $patterns=bigbin2dec($apase);
                    for(my $k=0; $k<$patterns; $k++) {
                        my $bigbin=bigdec2bin($k,$Nr-2);
                        my $myreachable=joinexon($gene_name,$start,$end,$bigbin,$a[0],0,0);
                        if ($myreachable eq -1) { next;}
                        my @usage=split("",$myreachable);
                        my $cid=$a[0]."_".$start;
                        for(my $j=0; $j<scalar(@usage); $j++) { if($usage[$j] eq 1) { $cid=$cid."|".($start+1+$j); }   }
                        $cid=$cid."|".$end;
                        if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
                        if ($debug > 0) { print $cid,"\t",join("",@usage),"\n"; }
                        my $exonL=$ExonL[0]-1;
                        my $exonR=$ExonR[0];
                        my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"1\"");
                        print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[0],$ExonR[0],".",$a[20],$gene_name,$info),"\n";
                        my $exoncnt=1;
                        for(my $j=0; $j<scalar(@usage); $j++) {
                            if($usage[$j] eq 1){
                                $exonL=$exonL.",".($ExonL[1+$j]-1);
                                $exonR=$exonR.",".$ExonR[1+$j];
                                $exoncnt++;
                                $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$exoncnt\"");
                                print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[1+$j],$ExonR[1+$j],".",$a[20],$gene_name,$info),"\n";
                            }
                        }
                        $exonL=$exonL.",".($ExonL[$Nr-1]-1);
                        $exonR=$exonR.",".$ExonR[$Nr-1];
                        $exoncnt++;
                        $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"$exoncnt\"");
                        print OUT1gtf join("\t",$a[1],"split","exon",$ExonL[$Nr-1],$ExonR[$Nr-1],".",$a[20],$gene_name,$info),"\n";
                        print OUT1refFlat join("\t",$a[6],$cid,$a[1],$strand,$a[3]-1,$a[2],$a[3]-1,$a[3]-1,$exoncnt,$exonL,$exonR,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
                    }
                }
            }
            else {
                my $cid=$a[0]."_1";
                if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
                $cidstat{$cid}=1;
                my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name2\"\;","exon_number","\"1\"");
                print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],$gene_name,$info),"\n";
                print OUTrefFlat join("\t","na",$cid,$a[1],$a[20],$left-1,$right,$left-1,$left-1,1,$left-1,$right,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
            }
        }
    }
    else {
        my $cid=$a[0]."_1";
        if (exists $usedcid{$cid}) { next;} else { $usedcid{$cid}=1; }
        $cidstat{$cid}=1;
        my $gene_name="na";
        my $info=join(" ","gene_id","\"$cid\"\;","transcript_id","\"$cid\"\;","gene_name","\"$gene_name\"\;","gene_name2","\"$gene_name\"\;","exon_number","\"1\"");
        print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],$gene_name,$info),"\n";
        print OUTrefFlat join("\t","na",$cid,$a[1],$a[20],$left-1,$right,$left-1,$left-1,1,$left-1,$right,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
    }
}

foreach my $id (keys %cidstat) { print OUTstat $id,"\n"; }

open(OUTFLAG,">Step4_CBR_finished");
print OUTFLAG "Step4_CBR_finished\n";
close OUTFLAG;

