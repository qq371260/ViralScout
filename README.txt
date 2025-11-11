--------------------------------------------------------------------------------ViralScout_tools
---------DVSA (Differential Viral Segment Abundance) analysis tool:
log2(FC) and log10(TPM+1)

---------VSFV (Viral Small-RNAs Feature Vector) analysis tool:
13, 17, 26, 52, 104 dimensional feature vector of 

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

conda activate ViralScout

---------Install ViralScout
sudo apt update
sudo apt install git
git clone https://github.com/qq371260/ViralScout.git
cd ViralScout

--------------------------------------------------------------------------------Running
cd DVSA # Please provide all needed files in this folder


cd VSFV # Please provide all needed files in this folder, 


