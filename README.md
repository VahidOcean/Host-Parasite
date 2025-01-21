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

## Usage
Each step is run using shell scripts found in the `scripts/` directory.


