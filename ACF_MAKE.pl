#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"config_file\"   \"output_sh\"  " if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my %SPEC;
open(IN,$filein) or die "Cannot open config_file $filein";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ((scalar(@a) < 2) or ($a[0] eq "")) { next; }
    $SPEC{$a[0]}=$a[1];
}
close IN;
my $command="rm -rf $fileout";
system($command);

my $bwa_seed_len=16;
my $bwa_min_score=20;
my $thread=1;
my $minJump=100;
my $maxJump=2500000;
my $minSSSum=10;
my $minSamplecnt=1;
my $minReadcnt=1;
my $MAS=20;
my $coverage=0.9;
my $Junc=6;
my $ER=0.05;
my $stranded="no";
my $pre_defined_circRNA="no";
my $search_trans_splicing="no";
my $blat_path="";
my $do_blat_search="no";
my $ts_coverage=$coverage;
my $ts_MAS=0;
my $ts_minSSSum=$minSSSum;
my $ts_maxSpan=$maxJump;
my $remove_temp="yes";
my $make_AS="1";
my $max_AS=10;
my $UNMAP2="no";
# check if all parameters are set
if (!exists $SPEC{"BWA_folder"}) { die "ERROR: BWA_folder must by specified in the config_file $filein";}
if (!exists $SPEC{"BWA_genome_Index"}) { die "ERROR: BWA_genome_Index must by specified in the config_file $filein";}
if (!exists $SPEC{"genome_fasta"}) { $SPEC{"genome_fasta"}=$SPEC{"BWA_genome_Index"}}
if (!exists $SPEC{"ACF_folder"}) { die "ERROR: ACF_folder must by specified in the config_file $filein";}
$SPEC{"CBR_folder"}=$SPEC{"ACF_folder"}."/CB_splice/";
if (!exists $SPEC{"Agtf"}) { $SPEC{"Agtf"}="no"; print "No gene annotation file is provided. Providing annotation would enhance the performance.\n";}
if (!exists $SPEC{"UNMAP"}) { die "UNMAP must by specified in the config_file $filein";}
#if (!exists $SPEC{"UNMAP_expr"}) { die "UNMAP_expr must by specified in the config_file $filein";}
my $UNMAP_expr_filename=$SPEC{"UNMAP"}."_expr";
$SPEC{"UNMAP_expr"}="no";
if (-e $UNMAP_expr_filename) { $SPEC{"UNMAP_expr"}=$UNMAP_expr_filename; }
if (!exists $SPEC{"Seq_len"}) { die "ERROR: length of SE reads or insert_size(mean+SD) of PE reads must by specified in the config_file $filein";}
if (exists $SPEC{"bowtie2_folder"}) { print "Bowtie2 will be used, insteand of BWA, for abundance quantification.\n";}
if (exists $SPEC{"remove_temp"}) { $remove_temp=$SPEC{"remove_temp"}; }
if (exists $SPEC{"make_AS"}) { $make_AS=$SPEC{"make_AS"}; }
if (exists $SPEC{"max_AS"}) { $max_AS=$SPEC{"max_AS"}; }
if (($make_AS ne 0) and ($make_AS ne 1)) { die "make_AS can be either 1 or 0\n"; }
if ($make_AS eq 1 ) { print "Potential Alternative splicing events will be enumerated for BSJs containing up to $max_AS split-exons. This can be very slow.\n"; }
if (exists $SPEC{"Thread"}) { $thread=$SPEC{"Thread"}; }
if (exists $SPEC{"minJump"}) { $minJump=$SPEC{"minJump"}; }
if (exists $SPEC{"maxJump"}) { $maxJump=$SPEC{"maxJump"}; }
if (exists $SPEC{"minSplicingScore"}) { $minSSSum=$SPEC{"minSplicingScore"}; }
if (exists $SPEC{"minSampleCnt"}) { $minSamplecnt=$SPEC{"minSampleCnt"}; }
if (exists $SPEC{"minReadCnt"}) { $minReadcnt=$SPEC{"minReadCnt"}; }
if (exists $SPEC{"minMappingQuality"}) { $MAS=$SPEC{"minMappingQuality"}; }
if (exists $SPEC{"Coverage"}) { $coverage=$SPEC{"Coverage"}; }
if (exists $SPEC{"minSpanJunc"}) { $Junc=$SPEC{"minSpanJunc"}; }
if (exists $SPEC{"ErrorRate"}) { $ER=$SPEC{"ErrorRate"}; }
if (exists $SPEC{"Strandness"}) { $stranded=$SPEC{"Strandness"}; }
if (exists $SPEC{"pre_defined_circle_bed"}) { $pre_defined_circRNA=$SPEC{"pre_defined_circle_bed"}; }
if (exists $SPEC{"BWA_seed_length"}) { $bwa_seed_len=$SPEC{"BWA_seed_length"}; }
if (exists $SPEC{"BWA_min_score"}) { $bwa_min_score=$SPEC{"BWA_min_score"}; }
if (exists $SPEC{"Search_trans_splicing"}) { $search_trans_splicing=$SPEC{"Search_trans_splicing"}; }
if (exists $SPEC{"blat_path"}) { $blat_path=$SPEC{"blat_path"}; }
if (exists $SPEC{"blat_search"}) { $do_blat_search=$SPEC{"blat_search"}; }
if ($blat_path eq "") { $do_blat_search="no"; }
if (exists $SPEC{"trans_splicing_coverage"}) { $ts_coverage=$SPEC{"trans_splicing_coverage"}; }
if (exists $SPEC{"trans_splicing_minMappingQuality"}) { $ts_MAS=$SPEC{"trans_splicing_minMappingQuality"}; }
$ts_minSSSum=$minSSSum;
if (exists $SPEC{"trans_splicing_minSplicingScore"}) { $ts_minSSSum=$SPEC{"trans_splicing_minSplicingScore"}; }
if (exists $SPEC{"trans_splicing_maxSpan"}) { $ts_maxSpan=$SPEC{"trans_splicing_maxSpan"}; }

if (exists $SPEC{"UNMAP2"}) {
    if ((exists $SPEC{"bowtie2_folder"}) and ($SPEC{"UNMAP_expr"} eq "no")) {
        $UNMAP2=$SPEC{"UNMAP2"};
        print "PE information will be used for abundance quantification.\n";
        #$stranded="no";
        #print "No strandness information is used.\n";
    }
    elsif (!exists $SPEC{"bowtie2_folder"})  {
        die "ERROR: Bowtie2 is needed for quantification using PE reads. Please provide the path to bowtie2.\n";
    }
    elsif ($SPEC{"UNMAP_expr"} ne "no") {
        $SPEC{"UNMAP_expr"}="no";
        warn "WARNing: No UNMAP_expr is needed when processing PE reads.  Changing UNMAP_expr to \"no\"\n";
    }
}

my $is_fasta=1;
open(IN1, $SPEC{"UNMAP"});
my $line=<IN1>;
chomp $line;
if ($line=~/^@/) { $is_fasta=0;}
close IN1;

open(OUT, ">".$fileout) or die "ERROR: Cannot open output_sh file $fileout";
print OUT "#!/bin/bash\n\n";
print OUT "date\n";
print OUT "#Step1\n";
print OUT "echo \"Step1 maping_the_unmapped_reads_to_genome Started\" \n";
if (($SPEC{"UNMAP_expr"} eq "no") and (exists $SPEC{"UNMAP2"})) {
    $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." ".$SPEC{"BWA_genome_Index"}." ".$SPEC{"UNMAP"}." ".$SPEC{"UNMAP2"}." \> temp.unmap.sam";
}
else {
    $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." ".$SPEC{"BWA_genome_Index"}." ".$SPEC{"UNMAP"}." \> temp.unmap.sam";
}
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step1.pl temp.unmap.parsed temp.unmap.sam $stranded $MAS $coverage";
print OUT $command,"\n";
#$command="rm -rf temp.unmap.parsed.UID.fa";
#print OUT $command,"\n";
#$command="ln -s ".$SPEC{"UNMAP"}." temp.unmap.parsed.UID.fa";
#print OUT $command,"\n";
print OUT "echo \"Step1 maping_the_unmapped_reads_to_genome Finished\" \n\n\n";


print OUT "#Step2\n";
print OUT "date\n";
print OUT "echo \"Step2 find_circle_supporting_sequences Started\" \n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step2.pl temp.unmap.parsed.2pp.S2 ".$SPEC{"CBR_folder"}."/ temp.unmap.parsed.2pp.S1 ".$SPEC{"genome_fasta"};
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step2_MuSeg.pl temp.unmap.parsed.segs ".$SPEC{"CBR_folder"}."/ ".$SPEC{"genome_fasta"}." 15 ".$SPEC{"Agtf"};
print OUT $command,"\n";
print OUT "echo \"Step2 find_circle_supporting_sequences Finished\" \n\n\n";


print OUT "#Step3\n";
print OUT "date\n";
print OUT "echo \"Step3 define_circle Started\" \n";
if ($pre_defined_circRNA eq "no"){
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step3.pl temp.unmap.parsed.2pp.S3 temp.unmap.parsed.2pp.S2.sum temp.unmap.parsed.segs.S2.sum";
    print OUT $command,"\n";
}
else {
    $command="perl ".$SPEC{"ACF_folder"}."/get_circRNA_BSJ_from_bed.pl ".$pre_defined_circRNA;
    print OUT $command,"\n";
    $command="perl ".$SPEC{"ACF_folder"}."/get_circRNA_agtf_from_bed.pl ".$pre_defined_circRNA;
    print OUT $command,"\n";
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step3.pl temp.unmap.parsed.2pp.S3 temp.unmap.parsed.2pp.S2.sum temp.unmap.parsed.segs.S2.sum temp.pre_defined_circRNA";
    print OUT $command,"\n";
}
print OUT "echo \"Step3 define_circle Finished\" \n\n\n";


print OUT "#Step4\n";
print OUT "date\n";
print OUT "echo \"Step4 annotate_select_and_make_pseudo_sequences_for_circles Started\" \n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step4.pl temp_circ temp.unmap.parsed.2pp.S3 ".$SPEC{"Agtf"}." 0 $minJump $maxJump $minSSSum";
print OUT $command,"\n";

$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step4_MEA.pl temp_circ_MEA temp_circ_MEA ".$SPEC{"Agtf"}." ".$make_AS." ".$max_AS." ".$SPEC{"Seq_len"};
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step4_CBR.pl temp_circ_CBR temp_circ_CBR ".$SPEC{"Agtf"}." ".$make_AS." ".$max_AS." ".$SPEC{"Seq_len"};
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_split_exon_border_biotype_genename.pl temp_circ.agtf temp_circ_MEA.gtf temp_circ_MEA.ext.gtf temp_circ_CBR.gtf temp_circ_CBR.ext.gtf temp.pre_defined_circRNA.agtf";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_seq_from_agtf_from_genomefa.pl temp_circ.pseudo temp_circ.agtf ".$SPEC{"genome_fasta"};
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_pseudo_circle.pl temp_circ.pseudo.gene.fa temp_circ.CL ".$SPEC{"Seq_len"};
print OUT $command,"\n";
print OUT "echo \"Step4 define_circle Finished\" \n\n\n";



print OUT "#Step5\n";
print OUT "date\n";
print OUT "echo \"Step5 caliberate_the_expression_of_circles Started\" \n";
if (($SPEC{"UNMAP_expr"} eq "no") and (exists $SPEC{"UNMAP2"})) {
    if (exists $SPEC{"bowtie2_folder"}) {
        $command=$SPEC{"bowtie2_folder"}."/bowtie2-build -q -o 1 temp_circ.CL temp_circ.CL.bt2index";
        print OUT $command,"\n";
        if ($is_fasta eq 1) {
            $command=$SPEC{"bowtie2_folder"}."/bowtie2 -p ".$thread." -f -a -x temp_circ.CL.bt2index -1 ".$SPEC{"UNMAP"}." -2 ".$SPEC{"UNMAP2"}." > temp_circ.bt2map2";
        }
        else {
            $command=$SPEC{"bowtie2_folder"}."/bowtie2 -p ".$thread." -a -x temp_circ.CL.bt2index -1 ".$SPEC{"UNMAP"}." -2 ".$SPEC{"UNMAP2"}." > temp_circ.bt2map2";
        }
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5_process_bt2aln.pl circRNA_candidates temp_circ.bt2map2 ".$SPEC{"UNMAP_expr"}." temp_circ_MEA temp_circ_CBR $stranded $Junc ".$SPEC{"Seq_len"}." $ER";
        print OUT $command,"\n";
    }
    else{
        $command=$SPEC{"BWA_folder"}."/bwa index temp_circ.CL";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." temp_circ.CL ".$SPEC{"UNMAP"}." ".$SPEC{"UNMAP2"}." > temp_circ.bwamap2";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5_process_bt2aln.pl circRNA_candidates temp_circ.bwamap2 ".$SPEC{"UNMAP_expr"}." temp_circ_MEA temp_circ_CBR $stranded $Junc ".$SPEC{"Seq_len"}." $ER";
        print OUT $command,"\n";
    }
} else {
    if (exists $SPEC{"bowtie2_folder"}) {
        $command=$SPEC{"bowtie2_folder"}."/bowtie2-build -q -o 1 temp_circ.CL temp_circ.CL.bt2index";
        print OUT $command,"\n";
        if ($is_fasta eq 1) {
            $command=$SPEC{"bowtie2_folder"}."/bowtie2 -p ".$thread." -f -a -x temp_circ.CL.bt2index -U ".$SPEC{"UNMAP"}." > temp_circ.bt2map";
        }
        else {
            $command=$SPEC{"bowtie2_folder"}."/bowtie2 -p ".$thread." -a -x temp_circ.CL.bt2index -U ".$SPEC{"UNMAP"}." > temp_circ.bt2map";
        }
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5_process_bt2aln.pl circRNA_candidates temp_circ.bt2map ".$SPEC{"UNMAP_expr"}." temp_circ_MEA temp_circ_CBR $stranded $Junc ".$SPEC{"Seq_len"}." $ER";
        print OUT $command,"\n";
    }
    else{
        $command=$SPEC{"BWA_folder"}."/bwa index temp_circ.CL";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." temp_circ.CL ".$SPEC{"UNMAP"}." > temp_circ.bwamap";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5_process_bt2aln.pl circRNA_candidates temp_circ.bwamap ".$SPEC{"UNMAP_expr"}." temp_circ_MEA temp_circ_CBR $stranded $Junc ".$SPEC{"Seq_len"}." $ER";
        print OUT $command,"\n";
    }
}


$command="perl ".$SPEC{"ACF_folder"}."/get_bed12_from_refFlat.pl circRNA_candidates.bed12 circRNA_candidates_refFlat circRNA_candidates_expr circRNA_candidates.bed12";
print OUT $command,"\n";


print OUT "date\n";


if (($search_trans_splicing eq "yes") and (!exists $SPEC{"UNMAP2"})) {
    print OUT "\n#Extra Step\: finding trans_splicing events\n";
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step1.pl ".$SPEC{"CBR_folder"}."/ temp.unmap.parsed.tmp ".$SPEC{"genome_fasta"}." temp.unmap.trans.splicing $ts_coverage 15 $ts_MAS $maxJump";
    print OUT $command,"\n";
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step2.pl temp.unmap.trans.splicing ".$SPEC{"genome_fasta"}." temp.unmap.trans.splicing ".$SPEC{"Seq_len"}." $ts_minSSSum ".$SPEC{"Agtf"};
    print OUT $command,"\n";
    if (($do_blat_search eq "yes") and (-e $SPEC{"BWA_genome_Index"})) {
        #increasing the -minScore parameter value beyond one-half of the query size has no further effect
        #my $newlen=int(1.3*($SPEC{"Seq_len"}));
        my $newlen=$SPEC{"Seq_len"};
        $command=$blat_path." -minScore=".$newlen." -noHead ".$SPEC{"BWA_genome_Index"}." temp.unmap.trans.splicing.tsloci.fa temp.unmap.trans.splicing.tsloci.psl";
        print OUT $command,"\n";
        #$command="perl -ne \'chomp; my \@a=split(\"\\t\",\$_); if(abs(\$a[11]-\$a[12]) > $newlen){print \$a[9],\"\\n\"\;}\' unmap.trans.splicing.tsloci.psl \> unmap.trans.splicing.tsloci.psl.badid";
        #print OUT $command,"\n";
        #$command="sort \-u unmap.trans.splicing.tsloci.psl.badid \> unmap.trans.splicing.tsloci.psl.badids";
        #print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_linear_id_from_psl.pl temp.unmap.trans.splicing.tsloci.psl temp.unmap.trans.splicing.tsloci.psl.badids";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_selected_from_pool_singleline.pl temp.unmap.trans.splicing.tsloci.psl.badids temp.unmap.trans.splicing.tsloci temp.unmap.trans.splicing.tsloci.good 0 -1";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_selected_from_pool_singleline.pl temp.unmap.trans.splicing.tsloci.psl.badids temp.unmap.trans.splicing.tsloci.anno temp.unmap.trans.splicing.tsloci.anno.good 0 -1";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_selected_fa_from_pool.pl temp.unmap.trans.splicing.tsloci.psl.badids temp.unmap.trans.splicing.tsloci.fa temp.unmap.trans.splicing.tsloci.fa.good -1";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa index temp.unmap.trans.splicing.tsloci.fa.good";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." temp.unmap.trans.splicing.tsloci.fa.good ".$SPEC{"UNMAP"}." \> temp.unmap.trans.splicing.sam";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5.pl temp.unmap.trans.splicing.sam temp.unmap.trans.splicing.tsloci.fa temp.unmap.trans.splicing.p1 ".$SPEC{"Seq_len"}." $Junc $ER $stranded";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step3.pl temp.unmap.parsed.tmp temp.unmap.trans.splicing ".$SPEC{"UNMAP_expr"};
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_fusion_circRNAs.pl fusion_circRNAs temp.unmap.trans.splicing ".$ts_maxSpan." temp.unmap.trans.splicing.expr";
        print OUT $command,"\n";
        $command="mv temp.unmap.trans.splicing.p1.2 trans.splicing.reads";
        print OUT $command,"\n";
        $command="mv temp.unmap.trans.splicing.expr trans.splicing.expr";
        print OUT $command,"\n";
    }
    else {
        $command=$SPEC{"BWA_folder"}."/bwa index temp.unmap.trans.splicing.tsloci.fa";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." temp.unmap.trans.splicing.tsloci.fa ".$SPEC{"UNMAP"}." \> temp.unmap.trans.splicing.sam";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5.pl temp.unmap.trans.splicing.sam temp.unmap.trans.splicing.tsloci.fa temp.unmap.trans.splicing.p1 ".$SPEC{"Seq_len"}." $Junc $ER $stranded";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step3.pl temp.unmap.parsed.tmp temp.unmap.trans.splicing ".$SPEC{"UNMAP_expr"};
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_fusion_circRNAs.pl fusion_circRNAs temp.unmap.trans.splicing ".$ts_maxSpan." temp.unmap.trans.splicing.expr";
        print OUT $command,"\n";
        $command="mv temp.unmap.trans.splicing.p1.2 trans.splicing.reads";
        print OUT $command,"\n";
        $command="mv temp.unmap.trans.splicing.expr trans.splicing.expr";
        print OUT $command,"\n";
    }
}

print OUT "date\n";

if ($remove_temp eq "yes") {
    $command="rm -rf temp\*";
    print OUT $command,"\n";
    $command="rm -rf Step*finished";
    print OUT $command,"\n";
}

close OUT;

