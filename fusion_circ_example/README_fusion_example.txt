We have:
a toy genome file : myGenome.fa
a teo gene gtf file : myGenome.gtf1, which file fusion circRNAs are generated
a toy gene agtf file : myGenes.agtf

#_Step-1_# pre-processing
# generate fusion circRNAs and sampling reads (SE 100nt) from them
perl ../simulate_gtf_for_fusion_circRNA.pl myGenome.gtf1 simu_fusion_circ 2 5 10
perl ../get_split_exon_border_biotype_genename.pl simu_fusion_circ.agtf simu_fusion_circ
perl ../get_seq_from_agtf_from_genomefa.pl simu_fusion_circ_seq simu_fusion_circ.agtf myGenome.fa
perl ../get_id_from_simulate_fusion_circRNA_gtf.pl simu_fusion_circ_id simu_fusion_circ
perl ../simulate_reads_for_fusion_circRNA.pl simu_fusion_circ_seq.gene.fa simu_fusion_circ_SE.fa 100

# make genome index for BWA
/Users/arthur/Desktop/Dev/bwa0715/bwa index myGenome.fa


#_Step-2_# configure parameters according to your settings. An example is shown in below
rm -rf SPEC.txt
perl -e 'print join("\t","BWA_folder","/Users/arthur/Desktop/Dev/bwa0715/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","BWA_genome_Index","myGenome.fa"),"\n";' >> SPEC.txt
perl -e 'print join("\t","ACF_folder","/Users/arthur/Desktop/Dev/acfs2/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Agtf","myGenes.agtf"),"\n";' >> SPEC.txt
perl -e 'print join("\t","UNMAP","simu_fusion_circ_SE.fa"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Seq_len","100"),"\n";' >> SPEC.txt
perl -e 'print join("\t","make_AS","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","minSplicingScore","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","remove_temp","yes"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Search_trans_splicing","yes"),"\n";' >> SPEC.txt


#_Step-3_# configure generate pipeline
perl ../ACF_MAKE.pl SPEC.txt run_acfs2.sh


#_Step-4_# run ACFS
bash run_acfs2.sh
The result is : fusion_circRNAs
