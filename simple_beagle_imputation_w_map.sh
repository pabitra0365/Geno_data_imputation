#!/bin/bash
set -e

## Simple Beagle missing data imputation with map
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
## 
## This script takes a (possibly gzipped) .vcf file containing missing genotype
## calls, and runs BEAGLE imputation on this file chromosome-by-chromosome.
## As the name suggests, it uses a genetic map file (in PLINK format) with genetic
## distance in cM. This map is currently hardcoded as the variable "map_file"
## 
## The script returns the input file, but with all calls imputed and phased. This script
## DOES NOT perform imputation in the sense of adding additional SNPs by comparing
## against a high-density reference panel.
##
## Positional inputs are read in from the command line. These are (in order)
##   1) input VCF (absolute path)
##   2) Genetic map file (absolute path)
##   3) number of threads to use
##
## The output imputed VCF file will be bgzipped, and will have the same name as 
## the input file, except with "_imp" added before the .vcf.gz extension.
##
## NOTES ON BEAGLE PARAMETERS
## 
## The default size of the sliding window is 50,000 SNPs, which should in
## general encompass whole chromosomes in wheat datasets. According to a
## personal communication with Beagle's author (Brian Browning), Beagle should
## perform perfectly well in this way. However, note that Beagle's options can be
## customized below, with the line that starts "java -jar $beajar"
##
##
## DEPENDENCIES:
##
## 1. java 8: Open a terminal and type "java -showversion". If you see
##    something like: openjdk version "1.8.X_XXX" at the top of the output,
##    you're all set. If not, you can type:
##
##      sudo apt-get update
##      sudo apt-get install default-jre
##
## 2. Beagle .jar file, with bash alias "beajar" set in ~/.profile
##    (note this will likely instead be ~/.bashrc in non-Ubuntu distributions)
##    
##    To set this up, cd to the desired directory for Beagle, then run:
##
## 	  wget http://faculty.washington.edu/browning/beagle/beagle.08Jun17.d8b.jar
##    echo export alias beajar="${PWD}/beagle.08Jun17.d8b.jar" >> ~/.profile
##
##    Note that the version number of Beagle may change:
##    check http://faculty.washington.edu/browning/beagle
##
## 3. bgzip, included with samtools, bcftools, or htslib: 
##      http://www.htslib.org/download/
##
## IMPORTANT: After everything is installed, either restart or logout/login to 
## update the user's $PATH
###############################################################################


#### Read input parameters ####
vcf_in=$1
map_file=$2
n_threads=$3


#### Executable ####

## Generate names for output directory and output VCF (no extension)
out_dir=$(dirname "${vcf_in}")
vcf_out=$(echo "${vcf_in}" | sed 's/.gz$//')
vcf_out=$(echo "$vcf_out" | sed 's/.vcf$//')_imp

## Sanity check on user-supplied number of cores
if (( $n_threads > $(nproc) )); then
	echo "User specified more cores (real or hyperthreaded) than supported on machine"
	echo "Setting n_cores to number on machine: $(nproc)"
	n_threads=$(nproc)
fi

## Pull contig info out of input VCF if present
gunzip -c $vcf_in | head -n 200 | grep "^##contig" > ${out_dir}/temp_contigs.txt

## Perform the imputation
## Beagle 5.0 works differently from 4.1. In 4.1, the user specified a window
## size in number of markers. In 5.0 the window size is specified in cM. If
## no genetic map is supplied, Beagle assumes 1cM = 1Mb. For our purposes,
## we probably just want to supply each chromosome as a window, so setting
## window to 1,000 will accomodate the size of the largest wheat chrom. (~ 1Gb).
java -jar $beajar gt="$vcf_in" out="$vcf_out" map="$map_file" nthreads=$n_threads window=205

## Add contig info to output VCF
zcat "${vcf_out}.vcf.gz" | head -n 200 | grep "^##" > "${vcf_out}_tempheadupdate.vcf"
cat "${out_dir}/temp_contigs.txt" >> "${vcf_out}_tempheadupdate.vcf"
zgrep -v "^##" "${vcf_out}.vcf.gz" >> "${vcf_out}_tempheadupdate.vcf"
bgzip "${vcf_out}_tempheadupdate.vcf"
mv "${vcf_out}_tempheadupdate.vcf.gz" "${vcf_out}.vcf.gz"

bcftools index -c "${vcf_out}.vcf.gz"

rm "${out_dir}/temp_contigs.txt"
exit 0;
