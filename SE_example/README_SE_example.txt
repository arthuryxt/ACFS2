We have:
a toy genome file : myGenome.fa
a toy gene agtf file : myGenes.agtf
Single-end reads: simu_circ_SE_150.fa 

#_Step-1_# pre-processing
# change fasta header
perl ../change_fastq_header.pl simu_circ_SE_150.fa simu_circ_SE_150.fasta Truseq_SE
# collapse reads
perl ../Truseq_merge_unique_fa.pl UNMAP simu_circ_SE_150.fasta
# make genome index for BWA
/Users/arthur/Desktop/Dev/bwa0715/bwa index myGenome.fa


#_Step-2_# configure parameters according to your settings. An example is shown in below
rm -rf SPEC.txt
perl -e 'print join("\t","BWA_folder","/Users/arthur/Desktop/Dev/bwa0715/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","BWA_genome_Index","myGenome.fa"),"\n";' >> SPEC.txt
perl -e 'print join("\t","ACF_folder","/Users/arthur/Desktop/Dev/acfs2/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Agtf","myGenes.agtf"),"\n";' >> SPEC.txt
perl -e 'print join("\t","UNMAP","UNMAP"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Seq_len","150"),"\n";' >> SPEC.txt
perl -e 'print join("\t","make_AS","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","minSplicingScore","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","remove_temp","yes"),"\n";' >> SPEC.txt


#_Step-3_# configure generate pipeline
perl ../ACF_MAKE.pl SPEC.txt run_acfs2.sh


#_Step-4_# run ACFS
bash run_acfs2.sh






############## alternatively, since we have only one sample, doing read-collapsing would probably not speed up the whole process. We can therefore start with the raw reads.

#_Step-1_# pre-processing, nothing to be done


#_Step-2_# configure parameters according to your settings. An example is shown in below. Note changes to the parameters : "UNMAP" and "UNMAP_expr"
rm -rf SPEC.txt
perl -e 'print join("\t","BWA_folder","/Users/arthur/Desktop/Dev/bwa0715/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","BWA_genome_Index","myGenome.fa"),"\n";' >> SPEC.txt
perl -e 'print join("\t","ACF_folder","/Users/arthur/Desktop/Dev/acfs2/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Agtf","myGenes.agtf"),"\n";' >> SPEC.txt
perl -e 'print join("\t","UNMAP","simu_circ_SE_150.fa"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Seq_len","150"),"\n";' >> SPEC.txt
perl -e 'print join("\t","make_AS","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","minSplicingScore","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","remove_temp","yes"),"\n";' >> SPEC.txt
perl -e 'print join("\t","UNMAP_expr","no"),"\n";' >> SPEC.txt


#_Step-3_# configure generate pipeline
perl ../ACF_MAKE.pl SPEC.txt run_acfs2.sh

#_Step-4_# run ACFS
bash run_acfs2.sh