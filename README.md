# Host-Parasite
# SNP Calling Pipeline

## Requirements
The following modules/software are required to run the pipeline:
- **FastQC**
- **cutadapt**
- **Stacks (version 2.67)**
- **VCFtools (version 0.1.15)**

## Pipeline Steps
1. **Quality Control**: Remove adapter contamination and assess quality using FastQC.
2. **Demultiplexing**: Separate samples using barcodes.
3. **SNP Calling**: Generate SNPs using `denovo_map.pl` (de novo assembly).
4. **Filtering**: Apply quality filters to the resulting VCF files.


# Quality Control Script

module load FastQC
fastqc *fq

module load cutadapt
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 75 --length 75 -q 25 \
    -o trimmed_S1_R1_001.fastq -p trimmed_S1_R2_001.fastq \
    Vahid_GBS_S1_R1_001.fastq.gz Vahid_GBS_S1_R2_001.fastq.gz

# Check the quality of trimmed files
fastqc trimmed_S1_R1_001.fastq trimmed_S1_R2_001.fastq

#!/bin/bash
# Demultiplexing Script

# Create required directories
mkdir -p raw samples

# Copy trimmed files into the raw folder
cp trimmed_S1_R1_001.fastq raw/
cp trimmed_S1_R2_001.fastq raw/

# Run process_radtags
module load Stacks
process_radtags -P -p raw/ -o samples/ -b Barcodes.txt -e pstI -r -c -q --inline-inline

# Calculate read counts
output_file="read_counts.txt"
for file in samples/*.fq.gz; do
    read_count=$(zcat "$file" | wc -l)
    read_count=$((read_count / 4))
    echo "File: $file, Read Count: $read_count" >> "$output_file"
done

# SNP Calling Script

module load Stacks
denovo_map.pl --samples samples/ --paired --popmap Fish.txt -o Fish.with.het.filter \
    -M 2 -n 3 -p 0.8 -X "populations: --max-obs-het 0.6 --write-random-snp --vcf"


# Filter based on missing data and minor allele frequency
module load VCFtools
vcftools --vcf populations.snps.vcf --max-missing 0.5 --maf 0.04 --minDP 3 \
    --recode --recode-INFO-all --out miss50

# Remove individuals with more than 80% missing data
vcftools --vcf miss50.recode.vcf --missing-indv
awk '$5 > 0.80' out.imiss | cut -f1 > lowDP-80.indv
vcftools --vcf miss50.recode.vcf --remove lowDP-80.indv --recode --recode-INFO-all --out miss50.INDV

# Apply additional filtering
vcftools --vcf miss50.INDV.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out final_snps

# Hetrosigosity measurment
# Load required libraries
library(vcfR)
library(ggplot2)

# Function to calculate heterozygosity
calculate_heterozygosity <- function(vcf_file_path) {
  vcf <- read.vcfR(vcf_file_path)
  geno <- extract.gt(vcf)
  het <- apply(geno, 2, function(x) sum(x == "0/1", na.rm = TRUE) / length(x))
  return(het)
}

# File paths
fish_vcf_path <- "data/Fish/final.50.recode.vcf"
copepod_vcf_path <- "data/Copepod/final.50.recode.vcf"

# Calculate heterozygosity for fish and copepods
fish_het <- calculate_heterozygosity(fish_vcf_path)
copepod_het <- calculate_heterozygosity(copepod_vcf_path)

# Shapiro-Wilk Normality Test
cat("Shapiro-Wilk Normality Test Results:\n")
fish_shapiro <- shapiro.test(fish_het)
copepod_shapiro <- shapiro.test(copepod_het)
cat("Fish heterozygosity p-value:", fish_shapiro$p.value, "\n")
cat("Copepod heterozygosity p-value:", copepod_shapiro$p.value, "\n")

# Non-parametric Wilcoxon Rank-Sum Test
cat("\nWilcoxon Rank-Sum Test Results:\n")
wilcox_test <- wilcox.test(fish_het, copepod_het, exact = FALSE)
print(wilcox_test)

# Combine the data for plotting
group <- c(rep("Fish", length(fish_het)), rep("Copepod", length(copepod_het)))
heterozygosity <- c(fish_het, copepod_het)
data <- data.frame(Group = group, Heterozygosity = heterozygosity)

# Violin plot with boxplot overlay
ggplot(data, aes(x = Group, y = Heterozygosity, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Heterozygosity Distribution: Fish vs. Copepods",
       x = "Group",
       y = "Heterozygosity") +
  scale_fill_manual(values = c("Fish" = "skyblue", "Copepod" = "salmon"))

# Boxplot
ggplot(data, aes(x = Group, y = Heterozygosity, fill = Group)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, notch = TRUE) +
  theme_minimal() +
  labs(title = "Boxplot: Heterozygosity of Fish vs. Copepods",
       x = "Group",
       y = "Heterozygosity") +
  scale_fill_manual(values = c("Fish" = "skyblue", "Copepod" = "salmon"))






