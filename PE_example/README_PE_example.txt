We have:
a toy genome file : myGenome.fa
a toy gene agtf file : myGenes.agtf
Paired-end reads: simu_circ_PE_150.1.fa , simu_circ_PE_150.2.fa, which were simulated with insert size mean=250, sd=50. So we will set the Seq_len as 250+50.

#_Step-1_# pre-processing
# reverse complement READ-2
perl ../reverse_complement.pl simu_circ_PE_150.2.fa simu_circ_PE_150.2.farc
# change fasta header
perl ../change_fastq_header.pl simu_circ_PE_150.1.fa simu_circ_PE_150.1.fasta Truseq_PE1
perl ../change_fastq_header.pl simu_circ_PE_150.2.farc simu_circ_PE_150.2.fasta Truseq_PE2
# collapse reads
perl ../Truseq_merge_unique_fa.pl UNMAP simu_circ_PE_150.1.fasta simu_circ_PE_150.2.fasta
# make genome index for BWA
/Users/arthur/Desktop/Dev/bwa0715/bwa index myGenome.fa


#_Step-2_# configure parameters according to your settings. An example is shown in below
rm -rf SPEC.txt
perl -e 'print join("\t","BWA_folder","/Users/arthur/Desktop/Dev/bwa0715/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","BWA_genome_Index","myGenome.fa"),"\n";' >> SPEC.txt
perl -e 'print join("\t","ACF_folder","/Users/arthur/Desktop/Dev/acfs2/"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Agtf","myGenes.agtf"),"\n";' >> SPEC.txt
perl -e 'print join("\t","UNMAP","UNMAP"),"\n";' >> SPEC.txt
perl -e 'print join("\t","Seq_len","300"),"\n";' >> SPEC.txt
perl -e 'print join("\t","make_AS","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","minSplicingScore","0"),"\n";' >> SPEC.txt
perl -e 'print join("\t","remove_temp","yes"),"\n";' >> SPEC.txt


#_Step-3_# configure generate pipeline
perl ../ACF_MAKE.pl SPEC.txt run_acfs2.sh


#_Step-4_# run ACFS
bash run_acfs2.sh



As we treated PE reads as two individual SE reasd, the final read count could be slightly inflated if both mates of a read span the BSJ.