#!/bin/bash
#EXAMPLE
# rnaseq_human -s3_in rstr-raw-fastq -s3_subdir Uganda-HIVneg-orig \   #Get fastq from AWS S3
#              -out ~/project/rstr.orig -name rstr.HIVneg.orig \       #Set output directory and filename
#              -ref_download false -s3_ref kadm-ref -release 102 \     #Get reference genome
#              -threads 90 \ 
#              -s3_out kadm-rstr-hiv-neg/discovery                     #Set output directory fo AWS S3
# 
# s3_in="rstr-raw-fastq"
# s3_subdir="Uganda-HIVneg-orig"
# out=~/project/rstr.orig
# name="rstr.HIVneg.orig"
# ref_download="false"
# s3_ref="kadm-ref"
# release=102
# threads=90
# s3_out="kadm-rstr-hiv-neg/discovery"
# adapter_length=7

rnaseq_human () {

##### Set parameters #####
#Defaults
subdir_default="fastq"
subdir_default="fastq"
name_default="rnaseq"
out_default=~/rnaseq
ref_default="false"
threads_default=1

#From user input
## If input not given, use default

while getopts s3_in:s3_subdir:out:name:ref_download:s3_ref:release:threads:s3_out: flag
do
    case "${flag}" in
        s3_in) s3_in=${OPTARG};;
        s3_subdir) s3_subdir=${OPTARG:-$subdir_default};;
        out) out=${OPTARG:-$out_default};;
        name) name=${OPTARG:-$name_default};;
        ref_download) ref_download=${OPTARG:-$ref_default};;
        s3_ref) s3_ref=${OPTARG};;
        release) release=${OPTARG};;
        threads) threads=${OPTARG:-$threads_default};;
        s3_out) s3_out=${OPTARG};;
    esac
done

##### Make directories #####
#Raw fastq data
mkdir -p "$out"_data
#All results
mkdir -p "$out"_results
#Result subdirectories
mkdir -p "$out"_results/fastqc_raw/
mkdir -p "$out"_results/fastq_trim/
mkdir -p "$out"_results/fastqc_trim/
mkdir -p "$out"_results/bam/
mkdir -p "$out"_results/bam_metrics/
mkdir -p "$out"_results/bam_filter/
mkdir -p "$out"_results/counts/
#reference genome
mkdir -p "$out"_ref

##### Std out and error settings #####
#Save stdout and stderr
command > "$out"_results/"$name".stdout.txt 2>&1

#Stop if any errors
set -e

##### Check software #####
fastqc --version
AdapterRemoval --version
samtools --version
bedtools --version
echo "Picard"
java -jar ~/apps/anaconda3/share/picard*/picard.jar CollectRnaSeqMetrics --version
echo "STAR"
STAR --version
featureCounts -v

##### Fuse data bucket #####
s3fs $s3_in "$out"_data -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007

##### Quality assessment 1 ##### 
#Include all fastq in all subdirectories
for filename in $(find "$out"_data/"$s3_subdir"/ -type f -name "*.fastq.gz" -print) ;
do
echo "Starting FastQC analysis of" $filename

fastqc $filename \
       -o "$out"_results/fastqc_raw/ \
       -t $threads
done

#Save to S3
aws s3 sync "$out"_results/ s3://$s3_out

##### Adapter removal ##### 
## Remove adapters 
## Remove reads with > 1 ambiguous base
## Trim ends until reach base witrstr_datah quality 30+
## Remove reads < 15 bp

#Wait for user input of adapter length determined from FastQC results

echo "Input adapter length for 5' trimming"

read -n 1 -p "Length (bp):" adapter_length

#Run trimming
paste <(find "$out"_data/"$s3_subdir"/ -type f -name "*R1*.fastq.gz" -print) \
      <(find "$out"_data/"$s3_subdir"/ -type f -name "*R2*.fastq.gz" -print) |

while read file1 file2;
do
  name=$(paste -d '\0' \
            <(awk -F'[_]S' '{print $1}' <(basename $file1)))
  
  AdapterRemoval --file1 $file1 --file2 $file2 \
    --basename "$out"_results/fastq_trim/$name --gzip \
    --trim5p $adapter_length --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --threads $threads
done

##### Quality assessment 2 ##### 
## Assess trimmed read quality using FastQC.
for filename in "$out"_results/fastq_trim/*pair[12].truncated.gz ;
do
echo "Starting FastQC analysis of" $filename
fastqc $filename \
       -o "$out"_results/fastqc_trim/ \
       -t $threads
done
    
#Save to S3
aws s3 sync "$out"_results/ s3://$s3_out

##### Alignment ##### 
# Download human ref genome if indicated
if [$ref_download = "true"]; then

  mkdir "$out"_ref/STARref
  cd "$out"_ref/STARref
  sudo curl -O ftp://ftp.ensembl.org/pub/release-$release/gtf/homo_sapiens/Homo_sapiens.GRCh38.$release.gtf.gz
  sudo curl -O ftp://ftp.ensembl.org/pub/release-$release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

  gunzip *

# Make genome index
  mkdir "$out"_ref/STARindex
  STAR --runMode genomeGenerate \
     --genomeDir "$out"_ref/STARindex \
     --genomeFastaFiles "$out"_ref/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile "$out"_ref/STARref/Homo_sapiens.GRCh38.$release.gtf \
     --sjdbOverhang 99 \
     --runThreadN $threads

  cp ./Log.out "$out"_ref/STARindex

#Save to S3
  aws s3 sync "$out"_ref/ s3://kadm-ref/release$release/
  cd

else
#Fuse index files from S3
  s3fs $s3_ref "$out"_ref -o passwd_file=~/.passwd-s3fs \
      -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007

fi

#Increase number of RAM files limit
if [[ $(ulimit -n) = 1048 ]]; then

  sudo printf "* soft nofile 100000\n* hard nofile 100000\nroot hard nofile 100000\nroot soft nofile 100000\n" >> /etc/security/limits.d/nofile.conf
  sudo printf "fs.file-max = 100000" >> /etc/sysctl.conf
  sudo printf "session required pam_limits.so" >> etc/pam.d/sshd
  sudo printf "session required pam_limits.so" >> /etc/pam.d/login
  sudo printf "UsePAM yes" >> /etc/ssh/sshd_config

fi

## Align with STAR
paste <(ls "$out"_results/fastq_trim/*pair1.truncated.gz) \
      <(ls "$out"_results/fastq_trim/*pair2.truncated.gz) |
      
while read file1 file2;
do
    echo "Aligning" $(basename  -- "$file1");
    
    name=$(paste -d '\0' \
            <(awk -F'[.]pair' '{print $1}' <(basename $file1)) \
            <(echo '_'))
    
    STAR --genomeDir "$out"_ref/release$release/STARindex \
         --readFilesIn $file1 $file2 \
         --readFilesCommand zcat \
         --outFileNamePrefix "$out"_results/bam/$name \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN $threads \
         --runRNGseed 8756
done

#Save to S3
aws s3 sync "$out"_results/ s3://$s3_out

##### PICARD assess alignments ##### 
#Get Picard ref if haven't already fused ref S3 bucket
if[$ref_download = "true"]; then

  mkdir "$out"_ref/PICARDref
  cd "$out"_ref/PICARDref
  sudo curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
  gunzip refFlat.txt.gz
  ## Remove chr in chromosome name to match ensembl alignment
  sed 's/chr//' refFlat.txt > refFlat.ensembl.txt
  cd ../..
  
fi

# median CV of gene model coverage
for bam_file in "$out"_results/bam/*sortedByCoord.out.bam ;
do
    java -XX:ParallelGCThreads=$threads \
        -jar ~/project/apps/anaconda3/share/picard*/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT="$out"_ref/PICARDref/refFlat.ensembl.txt \
        I=$bam_file  \
        O="$out"_results/bam_metrics/temp.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500

    #Append results
    echo $bam_file >> "$out"_results/bam_metrics/bam.metrics.tsv
    cat "$out"_results/bam_metrics/temp.tsv >> "$out"_results/bam_metrics/bam.metrics.tsv
    #Remove this iteration
    rm "$out"_results/bam_metrics/temp.tsv
done

##### SAMTOOLS assess alignments #####
## mapped_reads_w_dups aka alignment percentage

for bam_file in "$out"_results/bam/*sortedByCoord.out.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> "$out"_results/bam_metrics/bam.summary.tsv
    
    samtools flagstat -@ $threads $bam_file \
    >> "$out"_results/bam_metrics/bam.summary.tsv
done

#Save to S3
aws s3 sync "$out"_results/ s3://$s3_out

##### Quality filter BAM ##### 
## -h: keep header
## -f 3: keeps paired reads where both mapped
## -F 1284: removes unmapped reads, non-primary alignments, and PCR duplicates
## -q 30: min MAPQ of 30

## Paired reads
for bam_file in "$out"_results/bam/*sortedByCoord.out.bam ;
do
  filename=$(paste -d '\0' \
            <(awk -F'[_]Aligned' '{print $1}' <(basename $bam_file)) \
            <(echo '_filter.bam'))
  
  echo "Filtering" $bam_file          
  samtools view $bam_file \
      -h -f 3 -F 1284 -q 30 \
      -@ $threads \
      > "$out"_results/bam_filter/$filename
done

#Save to S3
aws s3 sync "$out"_results/ s3://$s3_out

##### SAMTOOLS assess filtered alignments #####

for bam_file in "$out"_results/bam_filter/*_filter.bam ;
do
    echo "Processing" $bam_file
    echo $bam_file >> "$out"_results/bam_metrics/bam.filter.summary.tsv
    
    samtools flagstat -@ $threads $bam_file \
    >> "$out"_results/bam_metrics/bam.filter.summary.tsv
done

##### Count reads in genes ##### 
#Function max is 64 threads
if[$threads < 65]; then
featureCounts -T $$threads -g gene_id -t exon -p \
  -a "$out"_ref/release$release/STARref/Homo_sapiens.GRCh38.$release.gtf \
  -o "$out"_results/counts/"$name".featurecounts.paired.tsv \
  "$out"_results/bam_filter/*_filter.bam
  else
featureCounts -T 64 -g gene_id -t exon -p \
  -a "$out"_ref/release$release/STARref/Homo_sapiens.GRCh38.$release.gtf \
  -o "$out"_results/counts/"$name".featurecounts.paired.tsv \
  "$out"_results/bam_filter/*_filter.bam
fi

##### Save to S3 #####
aws s3 sync "$out"_results/ s3://$s3_out

#unmount data
fusermount -u "$out"_data
if [$ref_download = "false"]; then fusermount -u "$out"_ref fi

}

##### END ##### 