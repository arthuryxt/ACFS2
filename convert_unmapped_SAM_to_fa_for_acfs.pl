#!/usr/bin/perl -w
use strict;
# samtools sort -m 20G --threads 4 -n -o unmapped.sorted.bam unmapped.bam
# samtools view unmapped.sorted.bam > unmapped.sam
# convert paired-end sam into fasta complying with acfs header format, and report also reverse-complementary sequences for the second read, still assuming the stranded RNA-Seq protocol : read-1 antisense to mRNA; read-2 sense to mRNA
# 2015-08-20 by Arthur
die "Usage: $0   \"output\"   \"input_sam\"   \"newName\"       \"enforce_acfs_header==1_or_0\"    \"min_length\(default==50nt\)\"    \"stranded_mRNA==1_or_0\"   \"verbose=1\" " if (@ARGV < 3);
my $fileout=$ARGV[0];   # Ctrl21.fa
my $filein1=$ARGV[1];   # SRA1234679.unmapped.sam
my $newName=$ARGV[2];   # Ctrl21
my $enforce=1;
if (scalar(@ARGV) > 3) { $enforce=$ARGV[3]; }
if (($enforce ne 0) and ($enforce ne 1)) { die "enable acfs header by setting to 1, disable by setting to 0.\n";}
my $min_len=50;
if (scalar(@ARGV) > 4) { $min_len=$ARGV[4]; }
my $stranded=1;     # 1 means TruSeq stranded mRNA-Seq; 0 means unstraned mRNA-Seq
if (($stranded ne 1) and ($stranded ne 0)) { die "stranded_mRNA could only be 1 or 0: 1 means TruSeq stranded mRNA-Seq; 0 means unstraned mRNA-Seq\;" }
if (scalar(@ARGV) > 5) { $stranded=$ARGV[5]; }
my $verbose=0;
if (scalar(@ARGV) > 6) { $verbose=$ARGV[6]; }
my $Error_rate=0.15;

my %uniq1;
my %uniq2;
open(IN, $filein1) or die "Cannot open input_sam : $filein1\n";
my $max_len=0;
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[2] ne "*") {next;}
    my $len=length($a[9]);
    if ($len < $min_len) { next; }
    if ($len > $max_len) { $max_len=$len; }
    if ($a[1] eq 4) {
        # single-end data
        $uniq1{$a[0]}=$a[9];
    }
    elsif ($a[1] eq 69){
        # pair-end data, read-1
        $uniq1{$a[0]}=$a[9];
    }
    elsif ($a[1] eq 133){
        # pair-end data, read-2
        $uniq2{$a[0]}=$a[9];
    }
    #if (exists $uniq1{$a[0]}) { $uniq2{$a[0]}=$a[9];}
    #else { $uniq1{$a[0]}=$a[9]; }
}
open(OUT, ">".$fileout) or die "Cannot write to output : $fileout\n";

sub rev_cpm($){
    my $seq=$_[0];
    $seq=~tr/[atcgATCG]/[TAGCTAGC]/;
    my $rc=scalar reverse $seq;
    return $rc;
}

sub min($$$){
    if ($_[0] < $_[2]){ pop @_; } else { shift @_; }
    return $_[0] < $_[1]? $_[0]:$_[1];
}

sub levenshtein($$){
    my @A=split //, uc shift; # upper case
    my @B=split //, uc shift;
    my @W=(0..@B);
    my ($i, $j, $cur, $next);
    for $i (0..$#A){
        $cur=$i+1;
        for $j (0..$#B){
            $next=min( $W[$j+1]+1, $cur+1, ($A[$i] ne $B[$j])+$W[$j] );
            $W[$j]=$cur;
            $cur=$next;
        }
        $W[@B]=$next;
    }
    return $next;
}

sub overlap_pairs($$){
    my $seq1=uc shift;
    my $seq2=uc shift;
    my $seq3=rev_cpm($seq2);
    my %uniq;
    my $kmer=6;
    for(my $i=0; $i<(length($seq1)-$kmer); $i=$i+1) {$uniq{$i}=substr($seq1,$i,$kmer)}
    my $f=0;
    my $hit='';
    my %hits;
    foreach my $anchor(sort {$a <=> $b} keys %uniq) {
        #print $anchor,"\t",$uniq{$anchor},"\n";
        my $offset=-1;
        my $found=index($seq3, $uniq{$anchor}, $offset);
        #if ($found < $anchor) {next;}
        if ($found != -1) {if($hit=~m/$uniq{$anchor}/){}else{$hit=$hit.$uniq{$anchor}}}
        while($found != -1) {
            my $flag=0;
            for (my $p=1; ($p <= $f)and($flag == 0); $p++){
                my @info=split("\t",$hits{$p});
                if (abs($anchor - $info[0] - $found + $info[1]) < $kmer) {$info[2]++; $hits{$p}=join("\t",@info); $flag=1;}
            }
            if ($flag eq 0) {$f++; $hits{$f}=join("\t",$anchor,$found,1);}
            $offset=$found+1;
            $found=index($seq3, $uniq{$anchor}, $offset);
        }
    }
    if ($verbose) { foreach my $id (keys %hits) { print $hits{$id},"\n"; } }
    my $result='';
    if ($f ne 0) {
        my %result;
        for(my $i=1; $i<=$f; $i++ ) {
            #print $hits{$i},"\n";
            my @tmp=split("\t",$hits{$i});
            $result{$tmp[2]}{$tmp[0]}{$tmp[1]}=join("\t",@tmp);
        }
        $f=0;
        foreach my $kmer (sort{$b <=> $a} keys %result) {
            if ($f eq 0) {
                foreach my $start (sort{$a <=> $b} keys %{$result{$kmer}}) {
                    if ($f eq 0) {
                        foreach my $q_start (sort{$a <=> $b} keys %{$result{$kmer}{$start}}) {
                            if ($f eq 0) {
                                my @tmp=split("\t",$result{$kmer}{$start}{$q_start});
                                my $Kmer=$kmer;
                                my $Klen1=length($seq1) - $start;
                                my $Klen3=length($seq3) - $q_start;
                                if (($Kmer > $Klen1) | ($Kmer > $Klen3)) {
                                    $Kmer=$Klen1 < $Klen3 ? $Klen1 : $Klen3;
                                }
                                if ($verbose) {  print "inner loop : common_length = $Kmer  start_on_seq1 = $start  start_on_seq3 = $q_start \n"; }
                                if ($start > $q_start) {
                                    my $tmpseq1=substr($seq1,$start);
                                    my $tmplen=length($tmpseq1);
                                    my $tmpseq3=substr($seq3,0,$tmplen);
                                    my $dist=levenshtein($tmpseq1, $tmpseq3);
                                    if ($verbose) { 
                                        print "inner loop : tmpseq1 = $tmpseq1\n";
                                        print "inner loop : tmpseq3 = $tmpseq3\n";
                                        print "inner loop : dist = $dist \n";
                                    }
                                    if ($dist/$tmplen <= $Error_rate ) {
                                        $f=1;
                                        if (($tmplen > length($seq3)) and ($verbose)) {
                                            print "inner loop : tmplen = $tmplen  start_on_seq1 = $start  start_on_seq3 = $q_start \n";
                                            print "seq1 = $seq1\n";
                                            print "seq3 = $seq3\n";
                                        }
                                        if ($tmplen >= length($seq3)) {
                                            $result=$f."\t".$seq1;
                                        }
                                        else {
                                            $result=$f."\t".$seq1.substr($seq3,$tmplen);
                                        } 
                                    }
                                }
                                elsif ($start < $q_start) {
                                    my $tmpseq3=substr($seq3,$q_start);
                                    my $tmplen=length($tmpseq3);
                                    my $tmpseq1=substr($seq1,0,$tmplen);
                                    my $dist=levenshtein($tmpseq1, $tmpseq3);
                                    if ($verbose) { 
                                        print "inner loop : tmpseq1 = $tmpseq1\n";
                                        print "inner loop : tmpseq3 = $tmpseq3\n";
                                        print "inner loop : dist = $dist \n";
                                    }
                                    if ($dist/$tmplen <= $Error_rate ) {
                                        $f=1;
                                        if (($tmplen > length($seq3)) and ($verbose)) {
                                            print "inner loop : tmplen = $tmplen   start_on_seq1 = $start  start_on_seq3 = $q_start \n";
                                            print "seq1 = $seq1\n";
                                            print "seq3 = $seq3\n";
                                        }
                                        if ($tmplen >= length($seq1)) {
                                            $result=$f."\t".$seq3;
                                        }
                                        else {
                                            $result=$f."\t".$seq3.substr($seq1,$tmplen);
                                        } 
                                    }
                                }
                                else {
                                    my $tmplen1=length($seq1);
                                    my $tmplen3=length($seq3);
                                    my $tmplen=$tmplen1 < $tmplen3 ? $tmplen1 : $tmplen3;
                                    my $tmpseq1=substr($seq1,0,$tmplen);
                                    my $tmpseq3=substr($seq3,0,$tmplen);
                                    my $dist=levenshtein($tmpseq1, $tmpseq3);
                                    if ($verbose) { 
                                        print "inner loop : tmpseq1 = $tmpseq1\n";
                                        print "inner loop : tmpseq3 = $tmpseq3\n";
                                        print "inner loop : dist = $dist \n";
                                    }
                                    if ($dist/$tmplen <= $Error_rate ) {
                                        $f=1;
                                        $result=$f."\t".$tmpseq1;
                                    }
                                }
                            }
                        } 
                    }
                }
            }
        }
    }
    return $result;
}

if ($enforce eq 1) {
    my %reported;
    foreach my $id (keys %uniq1) {
        if (exists $uniq2{$id}) {
            # there is a counterpart, check if the two overlap
            if (length($uniq1{$id}) < ($max_len-2)) {
                # the two reads must overlap, and read-1 is reported
                print OUT ">Truseq_".$newName."_".$id.".1\n",$uniq1{$id},"\n";
                if($stranded eq 0 ) { print OUT ">Truseq_".$newName."_".$id.".2\n",$uniq2{$id},"\n"; }
            }
            else {
                my $overlap=overlap_pairs($uniq1{$id},$uniq2{$id});
                if ($verbose) { print $overlap,"\n"; }
                my @Overlap=split("\t",$overlap);
                if (scalar(@Overlap) > 1) {
                    # the two reads overlaps
                    print OUT ">Truseq_".$newName."_".$id.".1\n",$Overlap[1],"\n";
                    if($stranded eq 0 ) { print OUT ">Truseq_".$newName."_".$id.".2\n",rev_cpm($Overlap[1]),"\n"; }
                }
                else{
                    # the two reads are far away from each other
                    print OUT ">Truseq_".$newName."_".$id.".11\n",$uniq1{$id},"\n";
                    if($stranded eq 0 ) { print OUT ">Truseq_".$newName."_".$id.".12\n",rev_cpm($uniq1{$id}),"\n"; }
                    if($stranded eq 0 ) { print OUT ">Truseq_".$newName."_".$id.".21\n",$uniq2{$id},"\n"; }
                    print OUT ">Truseq_".$newName."_".$id.".22\n",rev_cpm($uniq2{$id}),"\n";
                }
            }
        }
        else {
            print OUT ">Truseq_".$newName."_".$id.".11\n",$uniq1{$id},"\n";
            if($stranded eq 0 ) { print OUT ">Truseq_".$newName."_".$id.".12\n",rev_cpm($uniq1{$id}),"\n"; }
        }
        $reported{$id}=1;
    }
    foreach my $id (keys %uniq2){
        if (!exists $reported{$id}) {
            if($stranded eq 0 ) { print OUT ">Truseq_".$newName."_".$id.".21\n",$uniq2{$id},"\n"; }
            print OUT ">Truseq_".$newName."_".$id.".22\n",rev_cpm($uniq2{$id}),"\n";
        }
        
    }
}
else {
    my %reported;
    foreach my $id (keys %uniq1) {
        if (exists $uniq2{$id}) {
            # there is a counterpart, check if the two overlap
            if (length($uniq1{$id}) < $max_len) {
                # the two reads must overlap, and one of the two is NOT junction-spanning due to the wrong strandedness
                print OUT ">".$id.".1\n",$uniq1{$id},"\n";
                if($stranded eq 0 ) { print OUT ">".$id.".2\n",$uniq2{$id},"\n"; }
            }
            else {
                my $overlap=overlap_pairs($uniq1{$id},$uniq2{$id});
                if ($verbose) { print $overlap,"\n"; }
                my @Overlap=split("\t",$overlap);
                if (scalar(@Overlap) > 1) {
                    # the two reads overlaps
                    print OUT ">".$id.".1\n",$Overlap[1],"\n";
                    if($stranded eq 0 ) { print OUT ">".$id.".2\n",rev_cpm($Overlap[1]),"\n"; }
                }
                else{
                    # the two reads are far away from each other
                    print OUT ">".$id.".11\n",$uniq1{$id},"\n";
                    if($stranded eq 0 ) { print OUT ">".$id.".12\n",rev_cpm($uniq1{$id}),"\n"; }
                    if($stranded eq 0 ) { print OUT ">".$id.".21\n",$uniq2{$id},"\n"; }
                    print OUT ">".$id.".22\n",rev_cpm($uniq2{$id}),"\n";
                }
            }
        }
        else {
            print OUT ">".$id.".11\n",$uniq1{$id},"\n";
            if($stranded eq 0 ) { print OUT ">".$id.".12\n",rev_cpm($uniq1{$id}),"\n"; }
        }
        $reported{$id}=1;
    }
    foreach my $id (keys %uniq2){
        if (!exists $reported{$id}) {
            if($stranded eq 0 ) { print OUT ">".$id.".21\n",$uniq2{$id},"\n"; }
            print OUT ">".$id.".22\n",rev_cpm($uniq2{$id}),"\n";
        }
        
    }
}


#>Truseq_ERR861805_ERR861805.1269249.1
#GAGTACGGGAGGATCAATGTGTAGATGCCATTCTTGTAGCTCTTGCCTGCGAGTTTAGCTTGTGTTAAGGCATTTGGGGGGGGGACCATTTAACTTACTG
#                     TAGATGCCATTCTTGTAGCTCTTGCCTGCGAGTTTAGCTTGTGTTAAGGCATTTGGGGGGGGGACCATTTAACTTACTGCTAGTTTGGGATGTCTTTGTG
#>Truseq_ERR861805_ERR861805.1269249.2
#CACAAAGACATCCCAAACTAGCAGTAAGTTAAATGGTCCCCCCCCCAAATGCCTTAACACAAGCTAAACTCGCAGGCAAGAGCTACAAGAATGGCATCTA
#                     CAGTAAGTTAAATGGTCCCCCCCCCAAATGCCTTAACACAAGCTAAACTCGCAGGCAAGAGCTACAAGAATGGCATCTACACATTGATCCTCCCGTACTC
