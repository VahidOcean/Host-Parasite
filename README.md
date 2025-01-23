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

### Hetrozygosity 
setwd("/Users/sepva552/Desktop/GBS/Fish and Copepod/Fish/Fish.with.het.filter/Whole.50")

library(vcfR)
vcf <- read.vcfR("final.50.recode.vcf")
geno <- extract.gt(vcf)
het <- apply(geno, 2, function(x) sum(x == "0/1", na.rm = TRUE) / length(x))
print(het)


setwd("/Users/sepva552/Desktop/GBS/Fish and Copepod/copepod/LEA/with.het.filter/50whole")
library(vcfR)
vcf <- read.vcfR("final.50.recode.vcf")
geno <- extract.gt(vcf)
het <- apply(geno, 2, function(x) sum(x == "0/1", na.rm = TRUE) / length(x))
print(het)

fish_het <- c(0.16393443, 0.19262295, 0.12522769, 0.14845173, 0.17349727, 0.13615665, 0.16621129, 
              0.08060109, 0.13661202, 0.19535519, 0.11475410, 0.16621129, 0.10428051, 0.18897996, 
              0.20582878, 0.09836066, 0.29872495, 0.15619308, 0.15437158, 0.10519126, 0.13752277, 
              0.14708561, 0.17941712, 0.06785064, 0.07969035, 0.10473588, 0.12841530, 0.13387978, 
              0.12295082, 0.14207650, 0.13752277, 0.10063752, 0.09790528, 0.13570128, 0.10382514, 
              0.06147541, 0.07422587, 0.12613843)

copepod_het <- c(0.11523478, 0.17838835, 0.16567562, 0.16239491, 0.18412959, 0.14230059, 0.15501333, 
                 0.18823047, 0.12876769, 0.07176543, 0.20094320, 0.21160550, 0.21242567, 0.21878204, 
                 0.19315153, 0.19335657, 0.21673160, 0.21837195)
het_data <- data.frame(
  Heterozygosity = c(fish_het, copepod_het),
  Group = c(rep("Fish", length(fish_het)), rep("Copepod", length(copepod_het)))
)
boxplot(Heterozygosity ~ Group, data = het_data, col = c("skyblue", "orange"),
        main = "Heterozygosity Comparison", ylab = "Heterozygosity"
        
        ######## Test for Normality
shapiro.test(fish_het)
shapiro.test(copepod_het)

The results of the Shapiro-Wilk normality test show:
  
  Fish heterozygosity (fish_het): 
  𝑝
=
  0.02662
p=0.02662, indicating the data does not follow a normal distribution (at a significance level of 0.05).
Copepod heterozygosity (copepod_het): 
  𝑝
=
  0.06006
p=0.06006, suggesting the data might follow a normal distribution, though it’s borderline.

####non-parametric statistical test, since data wasnt normal

wilcox.test(fish_het, copepod_het, exact = FALSE)

The results of the Wilcoxon rank-sum test indicate:
  
  W = 155: This is the test statistic for the Wilcoxon rank-sum test.
p-value = 0.001068: This is highly significant (
  𝑝
  <
    0.05
  p<0.05), indicating that there is a statistically significant difference in heterozygosity between the fish and copepods.
Alternative hypothesis: The true location shift (median difference) between the two groups is not equal to 0

######plot
library(ggplot2)

# Combine the data(Violin plot)
group <- c(rep("Fish", length(fish_het)), rep("Copepod", length(copepod_het)))
heterozygosity <- c(fish_het, copepod_het)

# Create a data frame
data <- data.frame(Group = group, Heterozygosity = heterozygosity)
ggplot(data, aes(x = Group, y = Heterozygosity, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Heterozygosity Distribution: Fish vs. Copepods",
       x = "Group",
       y = "Heterozygosity") +
  scale_fill_manual(values = c("Fish" = "skyblue", "Copepod" = "salmon"))


# Combine the data(Box plot)
group <- c(rep("Fish", length(fish_het)), rep("Copepod", length(copepod_het)))
heterozygosity <- c(fish_het, copepod_het)

# Create a data frame
data <- data.frame(Group = group, Heterozygosity = heterozygosity)

# Boxplot
library(ggplot2)
ggplot(data, aes(x = Group, y = Heterozygosity, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Comparison of Heterozygosity Between Fish and Copepods",
       x = "Group",
       y = "Heterozygosity") +
  scale_fill_manual(values = c("Fish" = "skyblue", "Copepod" = "salmon"))

