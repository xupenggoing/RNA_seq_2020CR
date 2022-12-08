# Cancer_Res_2020
This is the learning of 2020 Cancer Res paper (doi: 10.1158/0008-5472)

## 1. Data preparation
According to the description of paper: RNA-sequencing (RNA-seq) data were deposited in the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus database under accession number GSE121105. Searching GSE121105, we find there are 21 files in this dataset. From the description here and the original paper, we need to figure out the experimental design.

## 1.1 Download Datasets
In GEO page, we can acquire the **SRP164949**, also we can click the SRA Run Selector, from which we can get more detailed information about the experiment design. Also, there are two important information: 1) the cells come from human; 2) this is the single-ended sequencing; 3) sequencing platform is Illumina HiSeq 2000; 4) the AvgSpotLen is 100.

To download these original data, I prefer to use SRA Explorer. When I search SRP164949, I can get 21 results. Click the box on the left side of Title, `add 21 to collection`. Click upright `21 saved datasets` and click `Raw FastQ Download URLs`.

```
mkdir 2020_cancer_res
cd 2020_cancer_res
mkdir 1_raw_data
cd 1_raw_data

vi download_data.sh

#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 10
#SBATCH --mem=16G
#SBATCH --mail-type=ALL

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/005/SRR7997165/SRR7997165.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/002/SRR7997162/SRR7997162.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/006/SRR7997166/SRR7997166.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/004/SRR7997164/SRR7997164.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/001/SRR7997161/SRR7997161.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/003/SRR7997163/SRR7997163.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/000/SRR7997170/SRR7997170.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/008/SRR7997168/SRR7997168.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/002/SRR7997172/SRR7997172.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/001/SRR7997171/SRR7997171.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/007/SRR7997167/SRR7997167.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/009/SRR7997169/SRR7997169.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/003/SRR7997173/SRR7997173.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/004/SRR7997174/SRR7997174.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/005/SRR7997175/SRR7997175.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/007/SRR7997177/SRR7997177.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/008/SRR7997178/SRR7997178.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/006/SRR7997176/SRR7997176.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/009/SRR7997179/SRR7997179.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/001/SRR7997181/SRR7997181.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR799/000/SRR7997180/SRR7997180.fastq.gz

sbatch download_data.sh
```

## 1.2 Unzip and Rename Datasets
I prefer to unzip and rename these datasets according to the experiment design. Below it the experiment design
<img width="1207" alt="image" src="https://user-images.githubusercontent.com/24839999/206330149-cf3df222-d7d4-4ab3-a3e5-81b5aa93e5da.png">

We can perform unzip and rename in the same step. TRA is short for trastuzumab, PER is short for pertuzumab.
gunzip -c SRR7997161.fastq.gz > BT474_Rep1
gunzip -c SRR7997162.fastq.gz > BT474_Rep2
gunzip -c SRR7997163.fastq.gz > BT474_Rep3
gunzip -c SRR7997164.fastq.gz > BT474_TRA_Rep1
gunzip -c SRR7997165.fastq.gz > BT474_TRA_Rep2
gunzip -c SRR7997166.fastq.gz > BT474_TRA_Rep3
