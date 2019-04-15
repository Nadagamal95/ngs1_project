
#############################################
##   - NGS1_course                         ##
##   - NGS1_Assignment                     ##
##   - Bash script                         ##
##   - Assignment_Steps:-                  ##
##   - Date: 9 April 2019                  ##
##   - copyright Nada Gamal                ##
##   - Nile University                     ##
#############################################

#############################################

source activate ngs1

#download whole sra data

mkdir ~/ngs1_project/ && cd ~/ngs1_project/
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra 

#Download the SRAtoolkit

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz

#unzip the SRAtoolkit

tar -xzf sratoolkit.2.9.6-ubuntu64.tar.gz

#add SRAtoolkit to the PATH:

PATH=$PATH: ~/ngs1_project/sratoolkit.2.9.6-ubuntu64/bin

##########################################################

#Extract the First 5M Reads From basic file (PAIRED ends)and convert them from SRA file to FASTQ file

source activate ngs1

cd ~/ngs1_project/
fastq-dump --split-files -X 5000000 SRR8797509.sra 

#counting number of sequences in the fastq file

grep '@' SRR8797509.fastq | wc -l

####################################################

#install seqkit

conda install -c bioconda seqkit
source activate ngs1

#Splitting the 5M reads inTo 5 samples each one has 1M read for the two paired ends

seqkit split2 -1 SRR8797509_1.fastq -2 SRR8797509_2.fastq -p 5 -O split_main_reads -f

###################################################

#shuffleing the main fastq file (paired-end reads r1 & r2) each containing 5000000 reads

mkdir ~/ngs1_project/shuffled_reads 

seqkit shuffle SRR8797509_1.fastq > cd ~/ngs1_project/shuffled_reads/shuffled_SRR8797509_1.part_001.fastq 
seqkit shuffle SRR8797509_2.fastq > cd ~/ngs1_project/shuffled_reads/shuffled_SRR8797509_2.part_001.fastq

#Splitting 5M shuffled reads into 1M (two paired)  after shuffling

seqkit split2 -1 shuffled_SRR8797509_1.part_001.fastq  -2 shuffled_SRR8797509_2.part_001.fastq  -p 5 -O split_shuffled_reads -f

##################################################

#fastaqc:

#installation

source activate ngs1
conda install -c bioconda fastqc
conda install -c bioconda multiqc

#fastqc for sample one
 
#fastqc assessment for the first sample before shuffling

mkdir ~/ngs1_project/fastqc/fastqc_unshuffled && cd ~/ngs1_project/fastqc/fastqc_unshuffled

for f in ~/ngs1_project/split_main_reads/SRR8797509_*.part_001.fastq;do fastqc -t 1 -f fastq -noextract -o . $f;done

#merging the fastQC reports of the two reads of sample 1 before shuffling

multiqc -z -o . .

####################################################

#fastqc assessment for the first sample after shuffling

mkdir ~/ngs1_project/fastqc/fastqc_shuffled && cd ~/ngs1_project/fastqc/fastqc_shuffled 

for f in ~/ngs1_project/shuffled_reads/split_shuffled_reads/shuffled_SRR8797509_*.part_001.part_001.fastq;do fastqc -t 1 -f fastq -noextract -o . $f;done

#merging the fastQC reports of the two reads of sample 1 after shuffling

multiqc -z -o . .

################################

#trimming

#mild trimming:(unshuffled)

mkdir ~/ngs1_project/trimming/Unshuffled && cd ~/ngs1_project/trimming/Unshuffled

for ((i=1;i<=5;i++));do

f1="Unshuffled/SRR8797509_1.part_00$i.fastq"
f2="Unshuffled/SRR8797509_2.part_00$i.fastq"

newf1="Unshuffled/AfterTrimming/SRR8797509_1.part_00$i.pe.trim.fq"
newf2="Unshuffled/AfterTrimming/SRR8797509_2.part_00$i.pe.trim.fq"

newf1U="Unshuffled/AfterTrimming/SRR8797509_1.part_00$i.se.trim.fq"
newf2U="Unshuffled/AfterTrimming/SRR8797509_2.part_00$i.se.trim.fq"

adap="~/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile \
$f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36

done

########################################

#Aggressive Trimmimg

mkdir ~/ngs1_project/trimming/Shuffled && cd ~/ngs1_project/trimming/Shuffled

for ((i=1;i<=5;i++));do

f1="Shuffled/shuffled_SRR8797509_1.part_001.part_00$i.fastq"
f2="Shuffled/shuffled_SRR8797509_2.part_001.part_00$i.fastq"

newf1="Shuffled/AfterTrimming/shuffled_SRR8797509_1.part_00$i.pe.trim.fq"
newf2="Shuffled/AfterTrimming/shuffled_SRR8797509_2.part_00$i.pe.trim.fq"

newf1U="Shuffled/AfterTrimming/shuffled_SRR8797509_1.part_00$i.se.trim.fq"
newf2U="Shuffled/AfterTrimming/shuffled_SRR8797509_2.part_00$i.se.trim.fq"

adap="~/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile \
$f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:25 MINLEN:36

done

########################################

#Alignment 

#Download human reference

mkdir -p ~/ngs1_project/sample_data && cd ~/ngs1_project/sample_data

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
gunzip gencode.v29.pc_transcripts.fa.gz

# Download the Transcriptome Annotation File

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

# Select the transcripts of Chr22

source activate ngs1
cd ~/ngs1_project/sample_data/
READS=$(grep "^chr22" gencode.v29.annotation.gtf | awk -F'\t' '{print $9}' | awk -F';' '{print $1}' | awk -F' ' '{print $2}' | awk -F'"' '{print $2}' | sort | uniq)

for value in $READS
    do  
        echo "Processing: $value"
        seqkit grep -r -p ${value} gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' >> gencode.v29.pc_transcripts.chr22.simplified.fa
    done

################################################

#Alignment for the 5 unshuffled samples by BWA
 
#install BWA

source activate ngs1
conda install -c bioconda bwa

#################################################

#index the genome

mkdir -p ~/ngs1_project/bwa_align/bwaIndex && cd ~/ngs1_project/bwa_align/bwaIndex
ln -s ~/ngs1_project/sample_data/gencode.v29.pc_transcripts.chr22.simplified.fa .
bwa index -a bwtsw gencode.v29.pc_transcripts.chr22.simplified.fa

###################################################

#sequence alignment

cd ~/ngs1_project/bwa_align
for ((i=1;i<=5;i++));
    do	
    R1="~/ngs1_project/trimming/Unshuffled/AfterTrimming/SRR8797509_1.part_00${i}.pe.trim.fq"
    R2="~/ngs1_project/trimming/Unshuffled/AfterTrimming/SRR8797509_2.part_00${i}.pe.trim.fq"
    /usr/bin/time -v bwa mem bwaIndex/gencode.v29.pc_transcripts.chr22.simplified.fa $R1 $R2 > bwaAlign_unshuffled_sample_part_00${i}.sam
    done

cd ~/ngs1_project/bwa_align/

for ((i=1;i<=5;i++));
do
    samtools flagstat unshuffled_sample_part_00${i}.sam > bwaAlign_unshuffled_sample_part_00${i}_stats.out
done

###############################################
	
#Alignment for the 5 shuffled samples by HISAT

#Download Data

cd ~/ngs1_project/sample_data
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

# gunzip chr22_with_ERCC92.fa.gz

wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

###################################################### 

#install Hisat

source activate ngs1
conda install -c bioconda hisat2

#####################################################

#Indexing

mkdir -p ~/ngs1_project/hisat_align/hisatIndex && cd ~/ngs1_project/hisat_align/hisatIndex
ln -s ~/ngs1_project/sample_data/chr22_with_ERCC92.fa .
hisat2_extract_splice_sites.py ~/ngs1_project/sample_data/chr22_with_ERCC92.gtf > splicesites.tsv
hisat2_extract_exons.py ~/ngs1_project/sample_data/chr22_with_ERCC92.gtf > exons.tsv
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92 

#########################################################
#sequence alignment

cd ~/ngs1_project/hisat_align

for ((i=1;i<=5;i++));
do

R1="~/ngs1_project/trimming/Shuffled/AfterTrimming/shuffled_SRR8797509_1.part_00${i}.pe.trim.fq"
R2="~/ngs1_project/trimming/Shuffled/AfterTrimming/shuffled_SRR8797509_2.part_00${i}.pe.trim.fq"
hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S shuffeled_hisatAlign${i}.sam

done

cd ~/ngs1_project/hisat_align/
for ((i=1;i<=5;i++));
do
    samtools flagstat shuffeled_hisatAlign${i}.sam > hisatlign_shuffled_sample_part_00${i}_stats.out
done

################################################

#prepare the SAM files for assembly

#install Samtools

source activate ngs1
conda install samtools

#################################################

#convert the SAM files to BAM files for the 5 unshuffled samples (BWA alignment)

cd ~/ngs1_project/bwa_align/

for ((i=1;i<=5;i++));
do

samtools view -bS bwaAlign_unshuffled_sample_part_00${i}.sam > bwaAlign_unshuffled_sample_part_00${i}.bam

done

###################################################
 
#convert the BAM file to a sorted BAM file for the 5 unshuffled samples. 

cd ~/ngs1_project/bwa_align/

for ((i=1;i<=5;i++));
do

samtools sort bwaAlign_unshuffled_sample_part_00${i}.bam -o bwaAlign_unshuffled_sample_part_00${i}.sorted.bam

done
#####################################################

#Export some useful statistics report for each sample indvidually

cd ~/ngs1_project/bwa_align/

for f in 1 2 3 4 5; do samtools flagstat bwaAlign_unshuffled_sample_part_00$f.sorted.bam  > useful_stat_bam$f.txt; done

#####################################################

#convert the SAM files to BAM files for the 5 shuffled samples (Hisat alignment)

cd ~/ngs-01/ngs1_project/hisat_align/
for ((i=1;i<=5;i++));
do

samtools view -bS shuffeled_hisatAlign${i}.sam > shuffeled_hisatAlign${i}.bam

done

#####################################################
 
#convert the BAM file to a sorted BAM file for the 5 shuffled samples (Hisat align). 

cd ~/ngs1_project/hisat_align/
for ((i=1;i<=5;i++));
do

samtools sort shuffeled_hisatAlign${i}.bam -o shuffeled_hisatAlign${i}.sorted.bam

done

#####################################################

#Export some useful statistics report for each sample indvidually

cd ~/ngs-01/ngs1_project/hisat_align/
for f in 1 2 3 4 5; do samtools flagstat shuffeled_hisatAlign$f.sorted.bam > shuff_useful_stat_$f.txt; done

#############################################

#Assembly  using stringTie

#install stringtie

source activate ngs1
conda install stringtie

###############################################

#Assembly with known previous annotations

#assembly for bwa alignment 
  
cd ~/ngs1_project/bwa_align/

for SAMPLE in 1 2 3 4 5;
    do
        stringtie bwaAlign_unshuffled_sample_part_00${SAMPLE}.sorted.bam --rf -l ref_sup_${SAMPLE} -G ~/workdir/sample_data/chr22_with_ERCC92.gtf -o ref_sup_${SAMPLE}.gtf 
done

## how many transcript do you have?

cd ~/ngs1_project/bwa_align
for SAMPLE in 1 2 3 4 5;
do
cat ref_sup_${SAMPLE}.gtf  | grep -v "^@" | awk '$3=="transcript"' | wc -l

done

##############################################

#assembly for haisat alignment

cd ~/ngs1_project/hisat_align/
for ((i=1;i<=5;i++));

do

stringtie shuffeled_hisatAlign${i}.sorted.bam --rf -l ref_sup_${i} -G ~/ngs1_project/sample_data/chr22_with_ERCC92.gtf -o ref_sup_${i}.gtf 

done

## how many transcript do you have?

cd ~/ngs1_project/hisat_align 
for ((i=1;i<=5;i++));

do

cat ref_sup_${i}.gtf  | grep -v "^@" | awk '$3=="transcript"' | wc -l

done

############################################## 

#Assembly without known annotations

#assembly for bwa alignment 

cd ~/ngs1_project/bwa_align/ 
for ((i=1;i<=5;i++));

do

stringtie bwaAlign_unshuffled_sample_part_00${i}.sorted.bam --rf -l ref_free_${i} -o ref_free_${i}.gtf 

done

## how many transcript do you have?

cd ~/ngs1_project/bwa_align
for ((i=1;i<=5;i++));
do

cat ref_free_${i}.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l

done

##########################################

#assembly for haisat alignment
 
 cd ~/ngs1_project/hisat_align/
for ((i=1;i<=5;i++));

do

stringtie shuffeled_hisatAlign${i}.sorted.bam --rf -l ref_free_${i} -o ref_free_${i}.gtf 

done

## how many transcript do you have?

cd ~/ngs1_project/hisat_align 

for ((i=1;i<=5;i++));

do

cat ref_free_${i}.gtf  | grep -v "^@" | awk '$3=="transcript"' | wc -l

done

############################################## 

#GTF-Compare

#Using GTF-Compare to Compare the Generated Annotation Files to a Reference Annotation.

#Create virtual evironment with conda

conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py

###############################################

#Install prerequisites

pypy3 -m pip install gffutils numpy tqdm 'intervaltree<3.0'

#################################################

#GtF compare bwa aligment

mkdir -p ~/ngs1_project/gtf-compare/gtfs && cd ~/ngs1_project/gtf-compare/gtfs
    
ln -s ~/ngs1_project/bwa_align/ref_sup_*.gtf .
ln -s ~/ngs1_project/bwa_align/ref_free_*.gtf .

# ln -s ~/workdir/sample_data/chr22_with_ERCC92.gtf .

mkdir -p ~/ngs1_project/gtf-compare/method_one && cd ~/ngs1_project/gtf-compare/method_one
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py

###################################################  

#Run

source activate ngs-gtf

cd ~/ngs1_project/bwa_gtf-compare/method_one
for f in 1 2 3 4 5;
        do
                pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/ref_free_*.gtf 
                pypy3 stat.py
done

# pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/chr22_with_ERCC92.gtf

##################################################################

#GtF compare Hisat aligment

mkdir -p ~/ngs1_project/gtf_compare/gtfs && cd ~/ngs1_project/gtf_compare/gtfs

ln -s ~/ngs1_project/hisat_Align/ref_sup_*.gtf .
ln -s ~/ngs1_project/hisat_Align/ref_free_*.gtf .

# ln -s ~/workdir/sample_data/chr22_with_ERCC92.gtf .

mkdir -p ~/ngs1_project/hisat_gtf-compare/method_one && cd ~/ngs1_project/hisat_gtf-compare/method_one
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py

######################################################################
  
#Run
source activate ngs-gtf

cd ~/ngs1_project/gtf-compare/method_one

for f in 1 2 3 4 5;
        do
                pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/ref_free_*.gtf 
                pypy3 stat.py
done

# pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/chr22_with_ERCC92.gtf

####################################################################################################

#Differential Expression

#Download data 

mkdir -p ~/ngs1_project/diff_exp && cd ~/ngs1_project/diff_exp/
wget -c https://0x0.st/zK57.gz -O ref.tar.gz
tar xvzf ref.tar.gz
wget -c https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/Day-6/deseq1.r
wget -c https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/Day-6/draw-heatmap.r
# Extract sample_data

#############################################################################

#Setup enviornemnt

conda activate ngs1
conda install kallisto
conda install samtools

# Install subread, we will use featureCount : a software program developed for counting reads to genomic features such as genes, exons, promoters and genomic bins.
conda install subread

# install r and dependicies
conda install r
conda install -y bioconductor-deseq r-gplots

##########################################################################

#Step 2 (Quantification)

mkdir -p ~/ngs1_project/diff_exp && cd ~/ngs1_project/diff_exp/
GTF=~/ngs1_project/sample_data/chr22_with_ERCC92.gtf

# Generate the counts.

featureCounts -a $GTF -g gene_name -o counts.txt ~/ngs1_project/bwa_align/bwaAlign_unshuffled_sample_part_*.sorted.bam ~/ngs1_project/hisat_align/shuffeled_hisatAlign*.sorted.bam

# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt

# Simplify the file to keep only the count columns.
#cat counts.txt | cut -f 1,7-16 > simple_counts.txt

###################################################################################  

# Analyze the counts with DESeq1.

cat simple_counts.txt | Rscript deseq1.r 5x5 > results_deseq1.tsv

#cat simple_counts.txt | Rscript deseq1.r 3x3 > results2_deseq1.tsv
#####################################################################################

#View only rows with pval < 0.05

cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > final_output.pdf

#cat results2_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
#cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > final2_output.pdf

######################################################################
