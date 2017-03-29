#!/bin/bash

date
#Step1
echo "Step1 maping_the_unmapped_reads_to_genome Started" 
/Users/arthur/Desktop/Dev/bwa0715//bwa mem -t 1 -k 16 -T 20 myGenome.fa simu_circ_SE_150.fa > temp.unmap.sam
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step1.pl temp.unmap.parsed temp.unmap.sam no 20 0.9
rm -rf temp.unmap.parsed.UID.fa
ln -s simu_circ_SE_150.fa temp.unmap.parsed.UID.fa
echo "Step1 maping_the_unmapped_reads_to_genome Finished" 


#Step2
date
echo "Step2 find_circle_supporting_sequences Started" 
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step2.pl temp.unmap.parsed.2pp.S2 /Users/arthur/Desktop/Dev/acfs2//CB_splice// temp.unmap.parsed.2pp.S1 myGenome.fa
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step2_MuSeg.pl temp.unmap.parsed.segs /Users/arthur/Desktop/Dev/acfs2//CB_splice// myGenome.fa 15 myGenes.agtf
echo "Step2 find_circle_supporting_sequences Finished" 


#Step3
date
echo "Step3 define_circle Started" 
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step3.pl temp.unmap.parsed.2pp.S3 temp.unmap.parsed.2pp.S2.sum temp.unmap.parsed.segs.S2.sum
echo "Step3 define_circle Finished" 


#Step4
date
echo "Step4 annotate_select_and_make_pseudo_sequences_for_circles Started" 
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step4.pl temp_circ temp.unmap.parsed.2pp.S3 myGenes.agtf 0 100 2500000 0
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step4_MEA.pl temp_circ_MEA temp_circ_MEA myGenes.agtf 0 150
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step4_CBR.pl temp_circ_CBR temp_circ_CBR myGenes.agtf 0 150
perl /Users/arthur/Desktop/Dev/acfs2//get_split_exon_border_biotype_genename.pl temp_circ.agtf temp_circ_MEA.gtf temp_circ_MEA.ext.gtf temp_circ_CBR.gtf temp_circ_CBR.ext.gtf temp.pre_defined_circRNA.agtf
perl /Users/arthur/Desktop/Dev/acfs2//get_seq_from_agtf_from_genomefa.pl temp_circ.pseudo temp_circ.agtf myGenome.fa
perl /Users/arthur/Desktop/Dev/acfs2//get_pseudo_circle.pl temp_circ.pseudo.gene.fa temp_circ.CL 150
echo "Step4 define_circle Finished" 


#Step5
date
echo "Step5 caliberate_the_expression_of_circles Started" 
/Users/arthur/Desktop/Dev/bwa0715//bwa index temp_circ.CL
/Users/arthur/Desktop/Dev/bwa0715//bwa mem -t 1 -k 16 -T 20 temp_circ.CL simu_circ_SE_150.fa > temp_circ.bwamap
perl /Users/arthur/Desktop/Dev/acfs2//ACF_Step5_process_bt2aln.pl circRNA_candidates temp_circ.bwamap no temp_circ_MEA temp_circ_CBR no 6 150 0.05
perl /Users/arthur/Desktop/Dev/acfs2//get_bed12_from_refFlat.pl circRNA_candidates.bed12 circRNA_candidates_refFlat circRNA_candidates_expr no no 
date
date
rm -rf temp*
rm -rf Step*finished
