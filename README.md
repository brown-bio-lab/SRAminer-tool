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

### SRAminer execution 
SRAminer is designed to run the three stages (Retrieval, Filtering and Extraction, Validation) in two separate scripts: sraminer_screen.py and sraminer_validate.py.

The first script is designed to identify potential hits (i.e., SRA runs) containing the specified target (i.e., top blastn matches to the user-input multi-fasta file), while the second script verifies accuracy of the target hits (matches) through de novo assembling the metagenomes from each initially positive SRA run, then screening these for blastn matches again.

There are two ways to execute the script: run it in a Unix/Linux shell in Terminal or submit it as a job on a high-performance computing cluster. When running it via the terminal, be aware that if the terminal session is closed, the job will terminate. To prevent this, you can run the script in the background using the screen command.

The user input command line argument required to run the script 
| Command-line arguments | Required/Optional | Default | Description |
|-----------------------|------------------|---------|-------------|
| `-np` | Required | — | Number of MPI processes (CPUs) to use for the run. This should match the number of available CPU cores on a local machine, or the number of CPU cores allocated by the scheduler on an HPC/HPCC system. |
| `-m` | Required | — | `mpi4py`: In Python, `-m` tells the interpreter to execute a module. |
| `--keyword` | Required | — | Specifies a mining keyword used to select a subset of NCBI SRA runs and must be enclosed in quotes on the command line. |
| `--threads` | Optional | 1 | Specifies the number of threads to be used by BLAST during sequence alignment. Increasing the number of threads can improve performance on multicore systems but should not exceed available CPU resources. |
| `--pi` | Required | — | Sets the minimum percent identity threshold used to post-filter BLAST tabular results produced by the pipeline. |
| `--len` | Required | — | Minimum base pair length used to post-filter BLAST tabular results produced by the pipeline. |
| `--grep_keyword` | Required | — | Post-filters BLAST results by matching a keyword pattern, allowing restriction to specific organisms or terms. |
| `--db` | Required | sra | Specifies the database by default (Sequence Read Archive). |
| `--mail` | Required | — | Email address required for NCBI server. |
| `--taxon` | Optional | — | Downloads genomes for a taxon using NCBI datasets if `user_db` is not provided. |
| `--user_db` | Required | — | The target FASTA file defining the biological signal of interest. Used to construct a BLAST database for screening SRA reads. |
| `--wd` | Optional | current directory | Defines a user-managed working directory created and populated by the pipeline. Users should provide a unique directory per analysis to avoid overwriting. |

### Running SRAminer
A. Running sraminer screen.py in a terminal with a user database:
sraminer_screen.py starts by taking the input as a command-line argument. Here is a generalized example of how command-line arguments are used:

###
     apptainer exec \
       --bind "<path_of_your_working_location>:/work" \
       "<path_to_sraminer_pipeline.sif>" \
       mpirun -np <number_of_CPUs> python3 /work/sraminer_screen.py \
         --keyword "<dataset_to_search>" \
         --threads 1 \
         --wd /work \
         --user_db <user_database.fasta> \
         --db sra \
         --mail your@email.com \
         --pi <minimum_percentage_identity> \
         --len <minimum_base_pair_length> \
         --grep_keyword "<Target>"

Upon execution, the potential positive hits are stored in a folder Potential_positive_SRA_Runs, and later you can assign the CPU number equal to the number of Potential positive hits for running the sraminer_validate.py.

B. Running sraminer_validate.py in a terminal with a user database:
###
     apptainer exec \
       --bind "<path_of_your_working_location>:/work" \
       "<path_to_sraminer_pipeline.sif>" \
       mpirun -np <number_of_CPUs> python3 /work/sraminer_validate.py \
         --keyword "<dataset_to_search>" \
         --threads 1 \
         --wd /work \
         --user_db <user_database.fasta> \
         --db sra \
         --mail your@email.com \
         --pi <minimum_percentage_identity> \
         --len <minimum_base_pair_length> \
         --grep_keyword "<Target>"

The positive target sequence will be stored in a directory named Positive_SRA, while metadata information will be stored in a file named metadata_tables.txt.

C. Running sraminer_screen.py in a terminal without the user database:
###
     apptainer exec \
       --bind "<path_of_your_working_location>:/work" \
       "<path_to_sraminer_pipeline.sif>" \
       mpirun -np <number_of_CPUs> python3 /work/sraminer_screen.py \
         --keyword "<dataset_to_search>" \
         --threads 1 \
         --wd /work \
         --taxon <taxon> \
         --db sra \
         --mail your@email.com \
         --pi <minimum_percentage_identity> \
         --len <minimum_base_pair_length> \
         --grep_keyword "<Target>"

D. Running sraminer_validate.py in a terminal without a user database:
###
     apptainer exec \
       --bind "<path_of_your_working_location>:/work" \
       "<path_to_sraminer_pipeline.sif>" \
       mpirun -np <number_of_CPUs> python3 /work/sraminer_validate.py \
         --keyword "<dataset_to_search>" \
         --threads 1 \
         --wd /work \
         --taxon <taxon> \
         --db sra \
         --mail your@email.com \
         --pi <minimum_percentage_identity> \
         --len <minimum_base_pair_length> \
         --grep_keyword "<Target>"

