--------------------------------------------------------------------------------ViralScout_tools
---------DVSA (Differential Viral Segment Abundance) analysis tool:
# Using log2(FC) and log10(TPM+1) values to calculate the Euclidean distance of
contigs to benchmark viral contigs to identify virus-like contigs.

---------VSFV (Viral Small-RNAs Feature Vector) analysis tool:
Using 13-, 17-, 26-, 52-, and 104-dimensional feature vectors of all or uniquely
mapped sRNAs of contigs to identify virus-like contigs based on Spearman
correlation to benchmark viral contigs or self-defined reference sRNA feature vectors.


--------------------------------------------------------------------------------Running_dependencies
python=3.9
bowtie2=2.5.1
samtools=1.17
pysam=0.21.0
pandas=1.5.3
numpy=1.24.3
scipy=1.10.1
matplotlib=3.7.1
seaborn=0.12.2
scikit-learn=1.3.0
joblib=1.3.2
biopython=1.79

Environment management system: conda 25.5.1
Tested platform: WSL (Windows Subsystem for Linux) with Ubuntu 22.04.3 LTS


--------------------------------------------------------------------------------Installation
---------Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc

# This will create a Conda environment to support ViralScout tools
conda create -n ViralScout_env -c conda-forge -c bioconda -c defaults \
    python=3.9 \
    bowtie2=2.5.1 \
    samtools=1.17 \
    pysam=0.21.0 \
    pandas=1.5.3 \
    numpy=1.24.3 \
    scipy=1.10.1 \
    matplotlib=3.7.1 \
    seaborn=0.12.2 \
    scikit-learn=1.3.0 \
    joblib=1.3.2 \
    biopython=1.79

# Activate the Conda environment before use
conda activate ViralScout_env 

---------Install ViralScout
sudo apt update

sudo apt install git

git clone https://github.com/qq371260/ViralScout.git

cd ViralScout

--------------------------------------------------------------------------------Running
---------Run DVSA tool
# Please place all required files in this folder and set it as the working directory
cd DVSA

# Simply run for test files: single reads
bash DVSA.sh \
  --contigs contig_test.fa \
  --positive sRNA-P_test.fastq \
  --negative sRNA-N_test.fastq \
  --benchmark benchmark_test.txt \
  --output ./result_1 \
  --threads 32

# Simply run for test files: paired reads
bash DVSA.sh \
  --contigs contig_test.fa \
  --positive RNA-P_test_R1.fastq,RNA-P_test_R2.fastq \
  --negative RNA-N_test_R1.fastq,RNA-N_test_R2.fastq \
  --benchmark benchmark_test.txt \
  --output ./result_2 \
  --threads 32

# For more usage details, please see the usage file

---------Run VSFV tool
# Please place all required files in this folder and set it as the working directory
cd VSFV

# Simply run for test files: using benchmark viral contigs
bash unisRNA_vfv.sh \
  -f contig_test.fa \
  -b benchmark_test.txt \
  -r sRNA_test.fastq \
  -o ./result_1 \
  -d 104

# Simply run for test files: using self-defined vsiRNA_simulant file
bash unisRNA_vfv.sh \
  -f contig_test.fa \
  -r sRNA_test.fastq \
  -o ./result_2 \
  -d 104 \

# Simply run for test files: using all sRNAs
bash unisRNA_vfv.sh \
  -f contig_test.fa \
  -b benchmark_test.txt \
  -r sRNA_test.fastq \
  -o ./result_3 \
  -d 104 \
  --keep-non-unique

# For more usage details, please see the usage file


