# SRAminer
**Authors: Era Sharma and Amanda M. V. Brown** 

**Brown Lab, Department of Biological Sciences, Texas Tech University, Texas, US**

SRAminer, developed in Python, is specifically tailored to extract data from the NCBI SRA (National Center for Biotechnology Information Sequence Read Archive), the primary repository housing raw sequences from diverse samples, studies, and projects. Its objective is to identify targets within this vast data set efficiently. The targets can be genes or genomes. 
## Purpose:
SRAminer has been crafted to facilitate seamless retrieval, extraction, and identification of targets from the extensive petabytes of sequencing data available within the NCBI SRA. Written in Python, SRAminer leverages a combination of Python libraries and pre-existing bioinformatics tools to achieve its functionality. It offers versatility in operation, supporting command-line usage and integration into bash script for High-performance clustering (HPC) environments. 
## Features:
- **Parallel processing**- Enabling efficient utilization of multiple CPU cores. Distributing SRA runs across individual CPU cores ensures optimal performance, preventing any CPU from remaining idle and allowing them to operate independently.
- **Overcoming Access restrictions**- SRAminer can bypass the conventional limitation of accessing a maximum of 10,000 SRA Runs per search in one go. In addition to accessing the list of SRA Runs, SRAMiner provides access to associated metadata for the specific SRA Runs. 
- **Streamline search process**- With SRAminer, users can streamline the process of target identification as it eliminates the need for manual retrieval of data.
- **Efficient storage management**- SRAminer automatically removes unnecessary files and folders and frees up valuable storage space ensuring that the storage resources are utilized optimally.

## To install and run SRAminer, follow these steps:
1.	Clone the SRAminer GitHub repository.
2.	Run the workflow using the pre-built Apptainer container image: sraminer_pipeline.sif

The container includes all required software dependencies, so no manual tool installation is necessary.

### The sraminer_pipeline.sif container bundles the following tools:
- SRA Toolkit (prefetch v3.1.1): /opt/sratoolkit/bin/prefetch
- BLAST+ (blastn v2.14.0+): /opt/ncbi-blast/bin/blastn
- MEGAHIT: /opt/megahit/bin/megahit
- BWA: /usr/bin/bwa
- SAMtools (v1.13): /usr/bin/samtools
- NCBI Datasets CLI: /usr/local/bin/datasets
  
### The only required user downloads are
1.	BLAST core_nt database
2.	NCBI taxonomy database (taxdb)

#### If you already have the databases
- Add the path of the core_nt database in the config.ini file.
- Place the files taxdb.bti and taxdb.btd in the directory where you will run the SRAminer scripts.
No additional downloads are required.

#### If you donot have the databases
1. Download the BLAST core_nt Database
   
First, check how many core_nt archive files are available:
###
     wget -q -O - https://ftp.ncbi.nlm.nih.gov/blast/db/ | grep core_nt | grep .tar.gz
Note the file number range shown in the output (for example, 00..01).
Then download the corresponding files:
###
     wget https://ftp.ncbi.nlm.nih.gov/blast/db/core_nt.{00..01}.tar.gz
Adjust the range (00..01) according to the number of files listed.
Extract the downloaded files:
###
     tar -xvzf core_nt.*.tar.gz

3. Download the NCBI Taxonomy Database (taxdb)
Download the taxonomy database:
###
     wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
Extract the archive:
###
     tar -xvzf taxdb.tar.gz

After completing these steps, you will have all the required databases. Remember to add the path to the core_nt database in config.ini and place the taxdb.bti and taxdb.btd in the directory where you will run the SRAminer scripts.
