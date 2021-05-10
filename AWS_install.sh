#!/bin/bash

#### Basic ####

sudo yum update -y
sudo yum upgrade -y

#link volume
sudo mkfs -t ext4 /dev/nvme2n1
sudo mkdir ~/apps
sudo mount /dev/nvme2n1 ~/apps/

sudo chmod 777 -R ~/apps/

#### AWS CLI ####

sudo printf "[default]\naws_access_key_id = ##########\naws_secret_access_key = ############" \
> ~/.aws/credentials

sudo printf "[default]\noutput = txt\nregion = us-west-2" \
> ~/.aws/credentials

#### Fuse ####

# sudo amazon-linux-extras install -y epel
# sudo yum install -y s3fs-fuse
# 
# ### Setup key
# echo ###########:############## > ~/.passwd-s3fs
# chmod 600  ~/.passwd-s3fs

#### Conda ####

sudo mkdir ~/apps
cd ~/apps
sudo curl -O https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
sudo bash Anaconda3-2020.11-Linux-x86_64.sh
### Save to /home/ec2-user/apps/anaconda3
sudo chmod 777 -R ~/apps/anaconda3
eval "$(/home/ec2-user/apps/anaconda3/bin/conda shell.bash hook)"
conda init

## Configure channel priority
conda config --add channels bioconda 
conda config --add channels conda-forge
conda config --add channels biobakery
conda config --set channel_priority false
conda config --set allow_conda_downgrades true

#### python3 ####
# conda install -c anaconda python=3.7 -y
# 
# alias python='/usr/bin/python3.7'
# . ~/.bashrc

#### RNA-seq tools ####
conda install -c conda-forge fastqc samtools -y
conda install adapterremoval bedtools -y
conda install -c bioconda/label/cf201901 picard star subread -y

#### Humann3 ####
#Install humann3
conda install -c biobakery humann -y
#Downgrade to get bowtie to work
conda install tbb=2020.2 -y
#Check that diamond is at least 0.9.36
#If not upgrade
# conda config --set allow_conda_downgrades false
# conda install diamond=2.0.8 -y

#### R ####
cd ~/apps/  
  wget https://cran.r-project.org/src/base/R-4/R-4.0.5.tar.gz
tar xf R-4.0.5.tar.gz
cd R-4.0.5/
  
  ### Dependencies
  sudo yum install -y gcc gcc-c++ gcc-gfortran readline-devel \
zlib-devel bzip2 bzip2-devel xz xz-devel \
libcurl libcurl.i686 libcurl-devel.x86_64 \
openssl-devel findutils libffi-devel \
libxml2-devel

### Compile and install
./configure --prefix=$HOME/R-4.0.5/ --with-x=no
make

## Set default path to R
export PATH=~/apps/R-4.0.5/bin:$PATH
## Close and reopen to take effect

#### R packages ####
R

options(Ncpus = 50)
install.packages(c("tidyverse","data.table","readxl","BiocManager","devtools",
                   "lme4","coxme","car","broom","foreach","doParallel"),
                 repos='http://cran.us.r-project.org')
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene","annotatr","GenomicRanges",
                       "org.Hs.eg.db", "limma"), ask=FALSE)
devtools::install_github("mixOmicsTeam/mixOmics")