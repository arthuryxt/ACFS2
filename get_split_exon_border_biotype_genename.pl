#!/usr/bin/perl -w
use strict;
# make stupid split exons regardless of junctions
die "Usage: $0  \"output\"   \"input\"   \"\(optional\)input_2 ... \" " if (@ARGV < 2);
my $fileout=$ARGV[0];    
my $fileins=scalar(@ARGV);

my $debug=0;
#if (scalar(@ARGV) > 2) {$debug=$ARGV[2];}
open OUT,">".$fileout;
open OUT2,">".$fileout."_gene";
my %uniq;
my %biotype;
my %Gname;
my %Terminals;
my %Internals;

for(my $fileNr=1; $fileNr<$fileins; $fileNr++) {
    my $filein=$ARGV[$fileNr];
    print "processing $filein...\n";
    if (-e $filein) {
        open(IN, $ARGV[$fileNr]);
        while(<IN>) {
            chomp;
            if (m/^#/) {next;}
            my @a=split("\t",$_);
            # MT      protein_coding  exon    3307    4262    .       +       .        gene_id "ENSG00000198888"; transcript_id "ENST00000361390"; exon_number "1"; gene_name "MT-ND1"; gene_biotype "protein_coding"; transcript_name "MT-ND1-201";
            if ($a[2] eq "exon") {
                my @b=split(/\"/,$a[8]);
                my $Nr=scalar(@b);
                my $gene_id="";
                my $transcript_id="";
                #my $exon_number="";
                my $gene_name="";
                my $gene_biotype=$a[5];
                for(my $i=0; $i<$Nr; $i=$i+2){
                    $b[$i]=~s/\s//g; $b[$i]=~s/\;//g;
                    if ($b[$i] eq "gene_id") {$gene_id=$b[$i+1];}
                    elsif ($b[$i] eq "transcript_id") {$transcript_id=$b[$i+1];}
                    elsif ($b[$i] eq "gene_name") {$gene_name=$b[$i+1];}
                    elsif ($b[$i] eq "gene_biotype") {$gene_biotype=$b[$i+1];}  # ensembl
                    elsif ($b[$i] eq "gene_type") {$gene_biotype=$b[$i+1];}     # gencode
                }
                if ($gene_name eq "") {$gene_name=$gene_id;}
                my $id=$gene_id."\t".$a[0]."\t".$a[6];
                if (!exists $uniq{$id}) {
                    $biotype{$id}=$gene_biotype;
                    $Gname{$id}=$gene_name;
                }
                $uniq{$id}{$a[3]}{$a[4]}=join("\t",@a);
                if ((exists $Terminals{$id}) and (exists $Terminals{$id}{$transcript_id})) {
                    my @T=split("\t",$Terminals{$id}{$transcript_id});
                    if ($a[3] < $T[0]) { $T[0]=$a[3]; }
                    if ($a[4] > $T[1]) { $T[1]=$a[4]; }
                    $Terminals{$id}{$transcript_id}=join("\t",@T);
                    $Internals{$id}{$transcript_id}=$Internals{$id}{$transcript_id}."\t".join("\t",$a[3],$a[4]);
                }
                else {
                    $Terminals{$id}{$transcript_id}=join("\t",$a[3],$a[4]);
                    $Internals{$id}{$transcript_id}=join("\t",$a[3],$a[4]);
                }
            }
        }
        close IN;
    }
}

foreach my $id (keys %uniq) {
    my %L;
    my %R;
    my %pair;
    my %internalborder;
    foreach my $tid (keys %{$Internals{$id}}) {
        my @Ttmp=split("\t",$Internals{$id}{$tid});
        for(my $i=1; $i<(scalar(@Ttmp)-1); $i++) { $internalborder{$Ttmp[$i]}=1; }
        my $Nr=scalar(@Ttmp)/2;
        my %tmp;
        for(my $i=0; $i<$Nr; $i++) { $tmp{$Ttmp[2*$i]}=$Ttmp[2*$i+1]; }
        my @Stmp;
        $Nr=0;
        foreach my $pos (sort{$a <=> $b} keys %tmp){ $Stmp[$Nr]=$pos; $Nr++; $Stmp[$Nr]=$tmp{$pos}; $Nr++;}
        undef %tmp;
        $Nr=scalar(@Stmp)/2;
        for(my $i=0; $i<$Nr; $i++) {
            if (($i>0) and (($Stmp[2*$i] - $Stmp[2*$i-1]) < 1)) {
                print "Warning1: Two exons overlaps. Gene = ".$id."   Transcript = ".$tid,"\n";
                print join("\t",$Stmp[2*$i-2],$Stmp[2*$i-1],$Stmp[2*$i],$Stmp[2*$i+1]),"\n";
                delete $pair{$Stmp[2*$i-2]}{$Stmp[2*$i-1]};
                $pair{$Ttmp[2*$i-2]}{$Ttmp[2*$i+1]}=1;
                delete $R{$Stmp[2*$i-1]};
                $R{$Stmp[2*$i+1]}=2;
            }
            else{
                $pair{$Ttmp[2*$i]}{$Ttmp[2*$i+1]}=1;
                $L{$Ttmp[2*$i]}=1;
                $R{$Ttmp[2*$i+1]}=2;
            }
        }
    }
    my %border;
    foreach my $pos (sort{$a <=> $b} keys %L) {
        if (exists $border{$pos}) { $border{$pos}+=1; }
        else {  $border{$pos}=1; }
    }
    foreach my $pos (sort{$a <=> $b} keys %R) {
        if (exists $border{$pos}) { $border{$pos}+=2; }
        else {  $border{$pos}=2; }
    }
    my $cond="L";     # L for sit on a left-border and is connecting to a right-border; R for sit on a right border and will break
    my @Start;
    my @End;
    my $exoncnt=0;
    my %splitpos;
    foreach my $pos (sort{$a <=> $b} keys %border) {
        if ($debug > 0) {
            print $pos,"\t",$border{$pos},"\t",$exoncnt,"\n";
            my $tmp="Left";
            for(my $i=1 > ($exoncnt-5) ? 1 : ($exoncnt-5); $i<$exoncnt; $i++){ $tmp=$tmp."\t".$Start[$i]; }
            print $tmp,"\n";
            $tmp="Rght";
            for(my $i=1 > ($exoncnt-5) ? 1 : ($exoncnt-5); $i<$exoncnt; $i++){ $tmp=$tmp."\t".$End[$i]; }
            print $tmp,"\n";
        }
        
        if ($exoncnt eq 0) {
            $exoncnt++;
            $cond="L";
            $Start[$exoncnt]=$pos;
        }
        else {
            if ($cond eq "L") {
                if ($border{$pos} eq 1) {
                    # tandem exon
                    $End[$exoncnt]=$pos-1;
                    $splitpos{$pos-1}=1;
                    $exoncnt++;
                    $Start[$exoncnt]=$pos;
                    $cond="L";
                }
                elsif ($border{$pos} eq 2) {
                    $End[$exoncnt]=$pos;
                    $cond="R";
                }
                elsif ($border{$pos} eq 3) {
                    print "Warning2: Two exons overlaps. Gene = ".$id."   exon border = ".$pos,"\n";
                    $End[$exoncnt]=$pos-1;
                    $exoncnt++;
                    $Start[$exoncnt]=$pos;
                    $End[$exoncnt]=$pos;
                    $exoncnt++;
                    $Start[$exoncnt]=$pos+1;
                    $splitpos{$pos-1}=1;
                    $splitpos{$pos}=1;
                    $splitpos{$pos+1}=1;
                    $cond="L";
                }
            }
            else {
                if ($border{$pos} eq 1) {
                    $exoncnt++;
                    $Start[$exoncnt]=$pos;
                    $cond="L";
                }
                elsif ($border{$pos} eq 2) {
                    $exoncnt++;
                    $Start[$exoncnt]=$End[$exoncnt-1]+1;
                    $splitpos{$End[$exoncnt-1]+1}=1;
                    $End[$exoncnt]=$pos;
                    $cond="R";
                }
                else {
                    print "Warning3: Two exons overlaps. Gene = ".$id."   exon border = ".$pos,"\n";
                    $exoncnt++;
                    $Start[$exoncnt]=$pos;
                    $End[$exoncnt]=$pos;
                    $splitpos{$pos}=1;
                    $cond="R";
                }
            }
        }
    }
    # check all exon borders to handel retained intron
    my %newpairs;
    for(my $i=1; $i<=$exoncnt; $i++) { $newpairs{$Start[$i]}=$End[$i]; }
    foreach my $lpos (sort{$a <=> $b} keys %pair) {
        foreach my $rpos (sort{$a <=> $b} keys %{$pair{$lpos}}) {
            my $newcnt=0;
            my @newL;
            my @newR;
            foreach my $pos (sort{$a <=> $b} keys %newpairs) {
                if (($lpos <= $pos) and ($newpairs{$pos} <= $rpos)) {
                    if ($pos eq $lpos) {
                        $newcnt=1;
                        $newL[$newcnt]=$pos;
                        $newR[$newcnt]=$newpairs{$pos};
                    }
                    else {
                        $newcnt++;
                        $newL[$newcnt]=$pos;
                        $newR[$newcnt]=$newpairs{$pos};
                    }
                }
            }
            for(my $i=1; $i<$newcnt; $i++) {
                if (($newL[$i+1] - $newR[$i]) > 1) {
                    $newpairs{$newR[$i]+1}=$newL[$i+1]-1;
                    $splitpos{$newR[$i]+1}=1;
                    $splitpos{$newL[$i+1]-1}=1;
                }
            }
        }
    }
    my %Term;
    foreach my $tid (keys %{$Terminals{$id}}) {
        my @Ttmp=split("\t",$Terminals{$id}{$tid});
        $Term{$Ttmp[0]}=1;
        $Term{$Ttmp[1]}=1;
        my $Nr=scalar(@Ttmp)/2;
    }
    my @chrinfo=split("\t",$id);
    my @FStart;
    my @FEnd;
    my $Fr=0;
    foreach my $pos (sort{$a <=> $b} keys %newpairs) {
        $Fr++;
        $FStart[$Fr]=$pos;
        $FEnd[$Fr]=$newpairs{$pos};
    }
    for(my $i=1; $i<=$Fr; $i++) {
        my $splicing="";
        if ((exists $Term{$FStart[$i]}) and (!exists $internalborder{$FStart[$i]})) { $splicing="2"; }
        elsif (exists $splitpos{$FStart[$i]}) { $splicing="0"; }
        else { $splicing="1"; }
        if ((exists $Term{$FEnd[$i]}) and (!exists $internalborder{$FEnd[$i]})) { $splicing=$splicing."2"; }
        elsif (exists $splitpos{$FEnd[$i]}) { $splicing=$splicing."0"; }
        else { $splicing=$splicing."1"; }
        print OUT join("\t",$chrinfo[1],"split","exon",$FStart[$i],$FEnd[$i],$biotype{$id},$chrinfo[2],$Gname{$id},$chrinfo[0]."___".$i."___".$Fr."___".$splicing),"\n";
    }
    #print OUT2 join("\t",$chrinfo[1],"split","exon",$FStart[1],$FEnd[$exoncnt],$biotype{$id},$chrinfo[2],$Gname{$id},$chrinfo[0]),"\n";
    print OUT2 join("\t",$chrinfo[1],"split","exon",$FStart[1],$FEnd[$Fr],$biotype{$id},$chrinfo[2],$Gname{$id},$chrinfo[0]),"\n";
}
