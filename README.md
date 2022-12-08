# Cancer_Res_2020
This is the learning of 2020 Cancer Res paper (doi: 10.1158/0008-5472)

## 1. Data preparation
According to the description of paper: RNA-sequencing (RNA-seq) data were deposited in the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus database under accession number GSE121105. Searching GSE121105, we find there are 21 files in this dataset. From the description here and the original paper, we need to figure out the experimental design.

### 1.1 Download Datasets
In GEO page, we can acquire the **SRP164949**, also we can click the **SRA Run Selector**, from which we can get more detailed information about the experiment design. Also, there are two important information: 1) the cells come from human; 2) this is the single-ended sequencing; 3) sequencing platform is Illumina HiSeq 2000; 4) the AvgSpotLen is 100.

To download these original data, I prefer to use SRA Explorer. When I search SRP164949, I can get 21 results. Click the box on the left side of Title, `add 21 to collection`. Click upright `21 saved datasets` and click `Raw FastQ Download URLs`.

```
mkdir 2020_cancer_res
cd 2020_cancer_res
mkdir 1_raw_data
cd 1_raw_data

##creat a data downloading script - start line##

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

##creat a data downloading script - end line##

sbatch download_data.sh
```

### 1.2 Unzip and Rename Datasets
I prefer to unzip and rename these datasets according to the experiment design. Below it’s the experiment design:

<img width="1207" alt="image" src="https://user-images.githubusercontent.com/24839999/206330149-cf3df222-d7d4-4ab3-a3e5-81b5aa93e5da.png">

We can perform unzip and rename in a single step. Click the Metadata in **SRA Run Selector**, by which you can download it. Open it with Microsoft Excel, choose Data->Text to Columns->Delimited->Next->Comma->Next->Finish. Then we can get the metadata and can easily sort out the file format we need.

<img width="1208" alt="image" src="https://user-images.githubusercontent.com/24839999/206332508-36b927e6-cc88-4fcb-9859-a814031229fb.png">

```
##Creat a sample information file under 2020_cancer_res
vi sample_info.txt

SRR7997161	BT474_Rep1
SRR7997162	BT474_Rep2
SRR7997163	BT474_Rep3
SRR7997164	BT474_TRA_Rep1
SRR7997165	BT474_TRA_Rep2
SRR7997166	BT474_TRA_Rep3
SRR7997167	BT474_TRA_PER_Rep1
SRR7997168	BT474_TRA_PER_Rep2
SRR7997169	BT474_TRA_PER_Rep3
SRR7997170	BT_TR1_Rep1
SRR7997171	BT_TR1_Rep2
SRR7997172	BT_TR1_Rep3
SRR7997173	BT_TR2_Rep1
SRR7997174	BT_TR2_Rep2
SRR7997175	BT_TR2_Rep3
SRR7997176	BT_TPR1_Rep1
SRR7997177	BT_TPR1_Rep2
SRR7997178	BT_TPR1_Rep3
SRR7997179	BT_TPR2_Rep1
SRR7997180	BT_TPR2_Rep2
SRR7997181	BT_TPR2_Rep3
```
We generate a sample_info.txt file here and copy content from metadata to this file. **TRA** is short for trastuzumab, **PER** is short for pertuzumab, **TR** stands for trastuzumab-resistant, **TPR** stands for trastuzumab + pertuzumab-resistant.

Then we can use `sample_info.txt` to simply the downstream analysis. The final command line of unzip and rename is `gunzip -c raw_data_name.fastq.gz > novel_name.fastq`
```
awk '{print "gunzip -c "$1".fastq.gz > "$2".fastq"}' ../sample_info.txt > unzip_rename.sh

##creat a data run_unzip_rename script - start line##

vi run_unzip_rename.sh

#!/bin/bash
#SBATCH --job-name=run_unzip_rename
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

source unzip_rename.sh

##creat a data run_unzip_rename script - end line##

sbatch run_unzip_rename.sh
```
- **Trouble shooting**: 1) in `run_unzip_rename.sh` script, I can't call `unzip_rename.sh` directly, it turned out `unzip_rename.sh: command not found`; 2) when I add work directory for `unzip_rename.sh` and put `./unzip_rename.sh`, it turned out `./unzip_rename.sh: Permission denied`;3) I used `source` and put `source unzip_rename.sh` to get this issue fixed.
## 2. QC and trimming
Previously, we can perform QC and Trimming seperately. However, the novel softwares are emerging, some are powerful enough to perform QC and Trimming in a single step, such as Fastp.

I installed Fastp in my environment before. You can create a new environment called `RNA_Seq` (specified with the `-n` option) and install `fastp` into it with the following command lines.
```
module load miniconda/4.9.2 #load miniconda in the module
conda create -n RNA_Seq fastp

##to use the applications in your environment, run the following##
module load miniconda
conda activate RNA_Seq
```
The sample usage of Fastp can refer to https://github.com/OpenGene/fastp
- for single end data (not compressed)
```
fastp -i in.fq -o out.fq
```
- for paired end data (gzip compressed)
```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```
To begin the QC and trimming with Fastp, I make a novel directory under `2020_cancer_res` named `2_QC_trimming_via_fastp` and move all the unzipped and renamed files into it.
```
mkdir 2_QC_trimming_via_fastp
mv 1_raw_data/BT* 2_QC_trimming_via_fastp/
```
For the QC and trimming, taking `BT474_Rep1.fastq` as an example, the command line is `fastp -i BT474_Rep1.fastq -o BT474_Rep1_trimmed.fastq -h BT474_Rep1.html -j BT474_Rep1.json`.

```
awk '{print "fastp -i "$2".fastq -o "$2"_trimmed.fastq -h "$2".html -j "$2".json"}' ../sample_info.txt > fastp.sh

##creat a run_fastp script - start line##

vi run_fastp.sh

#!/bin/bash
#SBATCH --job-name=run_fastp
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

source fastp.sh

##creat a run_fastp script - end line##

module load miniconda
conda activate RNA_Seq
sbatch run_fastp.sh
```
After the QC and trimming by Fastp, I will get three types of files: one ended up with `trimmed.fastq`, one with `.html`, and one with `.json`. Make novel directories for `.html` and `.json` files.
```
mkdir html_files
mkdir json_files
mv *.html ./html_files
mv *.json ./json_files
```
## 3. Mapping
Mapping is the most important step in the RNA-seq analysis. Before mapping, make a directory `3_mapping_via_hisat2` under `2020_cancer_res` and move all the *trimmed.fastq` files in `2020_cancer_res/2_QC_trimming_via_fastp` into it. 
```
mkdir 3_mapping_via_hisat2
cd 3_mapping_via_hisat2
mv ../2_QC_trimming_via_fastp/*trimmed.fastq ./
```
### 3.1 Build human reference genome
You can download the human reference genome from the following link: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40. Note that you need to download both reference genome data (.fna) and annotation data (.gft). You can download them to you laptop and then tranfer then to the web server. Also, you can download them to the web server directly via `wget link_address`.

Before building the human reference genome, make a directory `human_ref_genome` for it, then perform the following steps.
```
wget https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_6390e4bdef461250ef424de0&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true
wget https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_6390e4bdef461250ef424de0&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_GTF&Flat=true

nohup gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz &
nohup gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz &
```
Then you need to build the reference genome because they can't be used directly.
```
##creat a ref_genome_build - start line##
vi ref_genome_build.sh

#!/bin/bash
#SBATCH --job-name=ref_genome_build
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

hisat2-build GCF_000001405.40_GRCh38.p14_genomic.fna hisat2_built_Genome 1>hisat2-build.log 2>&1

##creat a ref_genome_build - end line##

module load HISAT2/2.2.1-gompi-2020b
sbatch ref_genome_build.sh

```
After the building, you will get 8 files. Previously, I have perform the same steps in other projects. So I don't need to re-run the code. I can use what I generated before. Find the location of the previous `Human_Ref_Genome` folder, then copy and paste it to the `human_ref_genome` folder directly.`cp` command can't copy directories, you need to add `-r` or `-R`.
```
cp -r Human_Ref_Genome/ ~/scratch60/2020_cancer_res/3_mapping_via_hisat2/
mv Human_Ref_Genome/ human_ref_genome
```
### 3.2 Mapping to the reference genome
We will use HISAT2 for RNA-seq reads mapping.

```
hisat2 --new-summary -p 10 -x ./grch38_snp_tran/genome_snp_tran -U G1_Rep1_R1_val_1.fq -S G1_Rep1_R1.sam --rna-strandness RF --summary-file G1_Rep1_R1_alignment_summary.txt

awk '{print "fastp -i "$2".fq -o "$2"_trimmed.fq"}' sample_info.txt > unpaired_hisat2.sh
```
Here `-U` means unpaired, `--ran-strandness`, RF. Then we generate `run_unpaired_hisat2_mapping.sh`
```
vi run_unpaired_hisat2_mapping.sh

##creat a run_unpaired_hisat2_mapping script - start line##

#!/bin/bash
#SBATCH --job-name=run_unpaired_hisat2_mapping
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

unpaired_hisat2.sh

##creat a run_unpaired_hisat2_mapping script - start line##

module load HISAT2/2.2.1-gompi-2020b
sbatch run_unpaired_hisat2_mapping.sh

### 3.3 SAM files to BAM files

We use Samtools to convert SAM files to BAM files.
```
samtools sort G1_Rep1.sam -o G1_Rep1.bam
awk '{print "fastp -i "$2".fq -o "$2"_trimmed.fq"}' sample_info.txt > sam2bam.sh

```
vi run_sam2bam.sh

##creat a run_sam2bam script - start line##

#!/bin/bash
#SBATCH --job-name=run_sam2bam
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

sam2bam.sh
##creat a run_sam2bam script - end line##

module load SAMtools/1.16-GCCcore-10.2.0
sbatch run_sam2bam.sh
```

### 3.4 Indexing BAM files

We use Samtools to generate .bam.bai files.
```
samtools index G1_Rep1.bam
awk '{print "fastp -i "$2".fq -o "$2"_trimmed.fq"}' sample_info.txt > indexing_bam.sh
```

```
vi run_indexing_bam.sh

##creat a run_indexing_bam script - start line##

#!/bin/bash
#SBATCH --job-name=run_sam2bam
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

indexing_bam.sh
##creat a run_indexing_bam script - end line##

module load SAMtools/1.16-GCCcore-10.2.0
sbatch run_indexing_bam.sh
```

## 4. Counting reads by Subread
We use `featureCounts` in `Subread` to perform the reads counting.

```
#!/bin/bash
#SBATCH --job-name=featureCounts_via_Subread_selfIndexRefGenome
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

featureCounts -t exon -g gene_id -a ../Human_Ref_Genome/GCF_000001405.40_GRCh38.p14_genomic.gtf -o Yao_RNA_counts.txt G1_Rep1_selfIndex.bam G1_Rep2_selfIndex.bam G1_Rep3_selfIndex.bam G2_Rep1_selfIndex.bam G2_Rep2_selfIndex.bam G2_Rep3_selfIndex.bam

module load Subread/2.0.3-GCC-10.2.0
sbatch featureCounts_via_Subread_selfIndexRefGenome.sh
```

For paired-ended reads,
```
-p, paired-end reads
-s, strand-specific read counting
-T, threads

featureCounts -p -T 24 -t exon -g gene_id -a ../Human_Ref_Genome/GCF_000001405.40_GRCh38.p14_genomic.gtf -o Yao_RNA_counts.txt G1_Rep1_selfIndex.bam G1_Rep2_selfIndex.bam G1_Rep3_selfIndex.bam G2_Rep1_selfIndex.bam G2_Rep2_selfIndex.bam G2_Rep3_selfIndex.bam
```
Then we use other tutorials to perform the downstream analysis.
