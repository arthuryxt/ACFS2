# ACFS2
Accurate CircRNA Finder Suite. Discovering circRNAs from RNA-Seq data, with alternative splicing.


# Overview
CircRNAs are generated through splicing, and alternative splicing within circRNAs have been observed. While ACFS aims at detecting the Back-Splice Junction (BSJ), ACFS2 fully utilize existing gene annotation to infer the internal exon structure of circRNAs.


# Change Log
- First release on 2017-03-14


# Installation
Simply unpack the ACFS2 package.


# Requirement
- bwa-0.7.15 (included in the package, but you **_need_** to "make")
- perl
- blat (not necessary, sometime it helps to rule out false positive fusion-circRNAs when there is a gene-duplication or gene-pseudogene dilemma)


# Pipeline scheme
1. Prepare genome- and transcriptome- unmapped reads or read-pairs.
2. Preprocess reads.
3. Set up parameters file <SPEC_example.txt>
4. Run ACFS2 using the preprocessed reads.
- make make Pipeline Bash file
    ```
    perl ACF_MAKE.pl <SPEC_example.txt> <BASH_example.sh>
    ```
- find circles
    ```
    bash <BASH_example.sh>
    ```


# Parameter file <SPEC_example.txt>
There are eight mandatory parameters to run ACFS2 in a basic mode. Searching for fusion-circRNAs is disabled by default. Please modify the config file "SPEC_example.txt" according to your specific organism, experimental design and sequence specs. The config file "SPEC_example.txt" is a two-column tab-delimited file : \<name\>\t\<value\>

Mandatory paramters:

| Parameter | value | Note | 
| --------- | ----- | ---- | 
| BWA_folder | /home/bin/bwa037a/ | path of the folder of bwa |
| BWA_genome_Index | /data/.../BWAIndex/genome.fa | full path to the index files | 
| BWA_genome_folder | /data/.../genome.fa | can be ignored if is the same as BWA_genome_Index | 
| ACF_folder | /home/bin/ACFS2/ | path of the folder of ACFS2 | 
| Agtf | /data/.../Homo_sapiens.GRCh37.71_split_exon.gtf | processed annotation, see the next section point-4 | 
| UNMAP | UNMAP | the collapsed fasta file | 
| Seq_len | 150 | length of sequencing reads | 

Optional parameters, the values in below are set as default:

| Parameter | value | Note |
| --------- | ----- | ---- |
| bowtie2_folder | /home/bin/bowtie2/  | path of the folder of bowtie2. Bowtie2 can be more sensitive in dealing with multiple mapping. Bowtie2 will be used for abundance quantification insteand of BWA if provided |
| remove_temp | yes | set to "no" to keep all temp files |
| make_AS | 1 | set to 0 to assume all internal exons are presented in the circRNA |
| Thread | 16 | number of threads used in bwa |
| BWA_seed_length | 16 | bwa seed length  |
| BWA_min_score | 20 | bwa min score to trigger report. For shorter reads, e.g. 50nt, set to 10 or lower could report more circRNAs at risk of higher FDR |
| minJump | 100 | the minimum distance of a back-splice. The smaller, the more likely you can find circles from short exons |
| maxJump | 2500000 | the maximum distance of a back-splice. The larger, the more likely you can find circles from long genes. The longest human gene is CNTNAP2 which spans 2.3M bp |
| minSplicingScore | 10 | the minimum score for the sum of splicing strength at both splice site, 10 corresponds to 95% of all human/mouse splice site pairs. One could also set it to a lower value, e.g. zero, and do a post-filtering after running acfs |
| minSampleCnt | 1 | the minimum number of samples that detect any given circle |
| minReadCnt | 1 | the minimum number of reads (from all samples) that detect any given circle |
| minMappingQuality | 20 | the minimum mapping quality of any given sequence |
| minSpanJunc | 6 | the minimum number of bases reach beyond the back-splice-site. The larger the less likely of false-positive |
| Coverage | 0.9 | the minimum percentage of any given read is aligned. The larger the more conserved the results are |
| ErrorRate | 0.05 | the maximum error rate for re-alignment. The smaller the better |
| Strandness | - | the strand information of sequencing, must be one of the {+, -, no}  |
| pre_defined_circle_bed | no | pre-defined circle annotation in bed6 or bed12 format (to increase sensitivity for lowly expressed circRNAs please include bed12 files of annotated one, e.g. merge the bed files from GSE61991) |
| Search_trans_splicing | no | set to "yes" to seach for trans-splicing reads |
| blat_search | yes | use blat to discard false positives results from gene duplication, turn off by "no" |
| blat_path | blat | full path to the executable, such as "/usr/bin/blat/blat", it is ignored if the blat_search option is set to no |
| trans_splicing_coverage | 0.9 | see Coverage |
| trans_splicing_minMappingQuality | 0 | see minMappingQuality |
| trans_splicing_minSplicingScore | 10 | see minSplicingScore |
| trans_splicing_maxSpan | 2500000 | the maximum distance between the junctions on the same gene for fusion circRNAs |




# Before running ACFS2, a few pre-process (Same as ACFS)
0. Map the RNA-Seq reads to genome and transcriptom, and extract the unmapped reads. This is **_recommended_** as those mapped reads will NOT span the back-splice sites, and therefore do NOT contribute to circRNA discovery. 

1. Change fasta/fastq header format to allow processing multiple samples in one run.
    This is **_IMPORTANT_** ! ACFS expects a special header format so that multiple samples can be processed in one run. Do change the default header such as ">HWUSI-EAS100R:6:73:941:1973" into ">Truseq_sample1_HWUSI-EAS100R:6:73:941:1973", where the "sample1" is the name of your choice describing the sample. Do the conversion as:
    ```
    perl change_fastq_header.pl SRR650317_1.fasta SRR650317_1.fa Truseq_SRR650317left
    perl change_fastq_header.pl SRR650317_2.fasta SRR650317_2.fa Truseq_SRR650317right
    ```
    Make sure there is **_No underline_** within the sample name. e.g. ">Truseq_ctrl.1_Default_header" and ">Truseq_AGO2KO.1_Default_header" are OK; ">Truseq_ctrl_1_Default_header" is BAD because ACF will register "ctrl" as the sample name instead of "ctrl_1".

2. Merge sequences from multiple fasta/fastq files into one fasta file, which saves time for mapping.
    ```
    perl Truseq_merge_unique_fa.pl UNMAP SRR650317_1.fa SRR650317_2.fa
    ```
    Alternatively, if there are many files to merge, generate a file (named filelist for example) contains the full-path of each file
    ```
    perl Truseq_merge_unique_filelist.pl <UNMAP> <filelist>
    ```
    ```UNMAP``` is the collasped fasta file which will be processed by ACFS, and ```UNMAP_expr``` contains the readcount per sequence in all the samples. If you change the name of ```UNMAP``` to ```SomethingElse```, then the readcount file will be automatically named as ```SomethingElse_expr```
    
    However, one **can** bypass the previous and this step to run ACFS **sample by sample**. This way, no fasta header reformatting and reads collapsing is needed. For each sample, set the value of ```UNMAP``` to the name of fasta/fastq in the SPEC_example.txt file, and set the value of ```UNMAP_expr``` to "no".
    
3. Build BWA index, using verion 0.73a (currently not support for other versions as the output format changes between versions)
    ```
    /bin/bwa073a/bwa index /data/iGenome/human/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
    ```
    
4. Prepare for annotation (recommended)
    Download the gtf file from iGenome package    or    download ensembl gtf here : ftp://ftp.ensembl.org/pub/current_gtf/
    Then run:
    ```
    perl get_split_exon_border_biotype_genename.pl </data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/genes.gtf> </data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf> ...
    ```
    The first argument is the output agtf file, the second argument (and possibly other arguments) is the input gtf file(s).
    Note, there are **differences** between this script in ACFS2 with the one in Acfs. Make sure you don't mix up the annotation files.

5. Strandedness is assumed as from the Truseq Stranded RNA-Seq, so the reads are reverse-complementary to mRNAs.
    - If the reads are actually sense (the same 5'->3' direction as mRNA), please reverse-complement all reads.
    - If the reads are actually stransless or pair-ended, please run in parallel original reads and reverse-complemented reads.
    - No pair-end information is used, as the exact junction-site must be supported by a single read. Seeing is believing.
    
6. For Paired-end data, it is highly recommended to align the reads to genome+transcriptome first (e.g. using Tophat2), and extract the unmapped read (read pairs) using the following script (the fasta header line is also modified)
    ```
    perl convert_unmapped_SAM_to_fa_for_acfs.pl <output_file_name> <unmapped.sam> sample_id
    ```



# Results
1. circRNAs are stored in bed12 format, which can be visualzed using UCSC Genome Browser:
    - circRNA_candidates.bed12
    The 1-st column of the bed12 file might need to be prefix by "chr" to visualize on Genome Browser.
    The higher the value in 5th column, the more "likely" that circle is true.  
    The name in 4th column can be seperated by "_" into six segments: chromosome, SD site, SA site, internal exon combination, strength of left Splice site and strength of right Splice site.  

2. The circRNA expression is stored here:
    - circRNA_candidates_expr
    The second column denote the name of the gene from which this circRNA is derived, whenever possible.

3. Fusion-circRNAs, if enabled:  
    - fusion_circRNAs  
    The Junctional sequences are stored in "unmap.trans.splicing.tsloci.fa".





# Tutorial
We provide two toy examples in the ```SE_example``` folder and ```PE_example```. Follow the README file inside each folder and run ACFS2.




# A few useful scripts for simulation  
To see the usage, simply run the perl scripts with no arguments.  

1. simulating SE reads from linear transcripts
    ```
    simulate_SE_reads_from_linear.pl
    ```

2. simulating PE reads from linear transcripts
    ```
    simulate_PE_reads_from_linear.pl
    ```

3. simulating circRNAs
    ```
    simulate_gtf_for_circRNA.pl
    get_split_exon_border_biotype_genename.pl
    get_seq_from_agtf.pl
    ```

4. simulating SE reads from circRNAs
    ```
    simulate_SE_reads_from_circRNA.pl
    ```

5. simulating PE reads from circRNAs
    ```
    simulate_PE_reads_from_circRNA.pl
    ```

6. simulating fusion-circRNAs
    ```
    simulate_gtf_for_fusion_circRNA.pl
    get_split_exon_border_biotype_genename.pl
    get_seq_from_agtf.pl
    get_id_from_simulate_fusion_circRNA_gtf.pl
    ```

7. simulating SE reads from fusion-circRNAs (or PE reads where the insert-size is the read length)
    ```
    simulate_reads_for_fusion_circRNA.pl
    ```



    
# Contact
This pipeline is developed and is maintained by Arthur Xintian You: arthur.yxt@gmail.com. I will do my best to respond in a timely manner.

# Cite
Scientific Reports 6, Article number: 38820 (2016) doi:10.1038/srep38820




