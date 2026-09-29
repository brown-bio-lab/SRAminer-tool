####################################################################################
####################################################################################
##                                                                                ##
##                                                                                ##
##      _ _ _ _                      _      _                                     ##
##     |  _ _ _|                    | \    / |  _                                 ##
##     | |        _    __  _ _ _ _  |  \  /  | |_|  _ _ _ _   _ _ _ _   _    __   ##
##     | |_ _ _  | |  / / |_ _ _  | |   \/   |  _  |  _ _  | |  _ _  | | |  / /   ##
##     |_ _ _  | | | / /   _ _ _| | | |\__/| | | | | |   | | | |_ _| | | | / /    ##
##           | | | |/_/   |  _ _| | | |    | | | | | |   | | | |_ _ _| | |/_/     ##
##      _ _ _| | | |      | |_ _| | | |    | | | | | |   | | | |_ _ _  | |        ##
##     |_ _ _ _| |_|      |_ _ _ _| |_|    |_| |_| |_|   |_| |_ _ _ _| |_|        ##
##                                                                                ##
##                                                                                ##
##                                                                                ##
##                                                                                ##
####################################################################################
####################################################################################
"""
   SraMiner - version 0.4(09/27/2026)-Part1

   Author - Era Sharma, Amanda M.V. Brown
   Organization - Brown Lab, Department of Biology, Texas Tech University
   Contact - ersharma@ttu.edu, amanda.mv.brown@ttu.edu

   Purpose- The purpose of this tool is to provide user a command line interface that will
            search over the NCBI Databse, with user defined keyword and find hidden targets
"""
import psutil
import configparser
from math import ceil
import re
import fileinput
import sys
import os
from Bio import Entrez
from urllib.request import urlopen
from urllib.error import HTTPError
import json
import csv
from collections import defaultdict
import subprocess
import argparse
import importlib
import xml.etree.ElementTree as ET
import glob
import fileinput
from Bio import SeqIO
import multiprocessing
import io
from multiprocessing import Pool
import time
import shutil
import shlex
import math
import re
import threading
import platform
import concurrent.futures
import zipfile
from prettytable import PrettyTable
from mpi4py import MPI
from prettytable import PrettyTable

parser = argparse.ArgumentParser()
parser.add_argument('--wd', type=str, default=os.getcwd())
parser.add_argument('--mail', type=str)
parser.add_argument('--keyword', type=str)
parser.add_argument('--db', type=str)
parser.add_argument('--user_db', type = str, help = "Custom databse required to do Blast (database format - fasta)")
parser.add_argument('--taxon', type = str, help = "Taxon to make a custom database")
parser.add_argument('--threads', help = "No of threads to use")
parser.add_argument('--grep_keyword', help = "target to search for as 'Wolbachia|endosymbiont'")
parser.add_argument('--pi', help = "Percentage identity")
parser.add_argument('--len', help = "Length in basepair")
parser.add_argument('--field', help = "field")
parser.add_argument('--maximum_returned_items', help = "mri")
parser.add_argument('--downloader', choices=["sratools", "sracha"], default="sratools",
                    help="SRA download tool: 'sratools' (default, all platforms) or "
                         "'sracha' (faster; not for 454/SOLiD/Ion Torrent)")
args = parser.parse_args()


''' Function that checks if a path exists off of the working directory. If there is no path, one is created. '''

def pathExist(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"The new directory '{path}' is created!")
    else:
        print(f"The directory '{path}' already exists!")

# Define the list of directories to create
directories = [
    args.wd + "/SRR_FILES",
    args.wd + "/FASTQ_FILES",
    args.wd + "/BLAST",
    args.wd + "/BLAST_NT_RESULTS",
    args.wd + "/BLAST_RESULT",
    args.wd + "/Potential_positive_SRA_Runs",
    args.wd + "/blast_table",
    args.wd + "/BLAST_SPLIT_NT"
]

#Function-1
# Function: get_ncbi_ids
# Description: uses the Entrez.efetch to return ids based on search term
# Parameters:
#   - email
#   - search_term
#   - database
# Returns:
# id_list

def get_ncbi_ids(email, database, search_term, field, maximum_returned_items, threads):

    start_time_op1 = time.time()

    Entrez.email = email
    search_results_xml_1 = Entrez.esearch(db=database, term=search_term, field=field, retmax=maximum_returned_items, retstart=0)
    ncbi_dict = Entrez.read(search_results_xml_1)
    count = ncbi_dict["Count"]

    iterations = int(count)/10000

    total_iterations= math.ceil(iterations)

    retmax = 10000

    retstarts = [i * retmax for i in range(int(total_iterations))]
    id_list = []

    # Iterate through the retstarts and retrieve the IDs
    for offset in retstarts:
        search_results_xml = Entrez.esearch(db=database, term=search_term, field=field, retmax=retmax, retstart=offset)
        ncbi_dict = Entrez.read(search_results_xml)
        id_list.extend(ncbi_dict["IdList"])

    end_time_op1 = time.time()
    time_taken_op1 = end_time_op1 - start_time_op1
    time_taken_op1_minutes = int(time_taken_op1 // 60)
    time_taken_op1_seconds = time_taken_op1 % 60
    print("Time taken for fetching SRA Id list:", time_taken_op1_minutes, "minutes and {:.2f} seconds".format(time_taken_op1_seconds))
    return(id_list)

#Function-2
# Function: get_sra_runs
# Description: uses the Entrez.efetch to take id_list as an input to retrieve the SRA Run associated with the id
# Parameters:
#   - id_list
# Returns:
# SRA Run Ids

def get_sra_runs(id_list, database):

    start_time_op2 = time.time()
    size = int(len(id_list))

    retmax =10000
    iterations = math.ceil(size/retmax)

    sra_list = []
    sra_final_list = []

    for i in range(0, size, retmax):
        batch = id_list[i:i + retmax]
        batch_ids = ",".join(batch)

        try:
           handle= Entrez.efetch(db="sra", id=batch_ids, rettype="runinfo", retmode="txt")
           result= str(handle.read()).split("\n")
           sra_list.append(result)

           for item in result :
               field1 = item.split(",")
               sra_final_list.append(field1[0])

        except Exception as e:
            print(f"Error fetching batch {i}-{i + retmax}: {str(e)}")

    end_time_op2 = time.time()
    time_taken_op2 = end_time_op2 - start_time_op2
    time_taken_op2_minutes = int(time_taken_op2 // 60)
    time_taken_op2_seconds = time_taken_op2 % 60
    print("Time taken to retrieve the SRA Runs:", time_taken_op2_minutes, "minutes and {:.2f} seconds".format(time_taken_op2_seconds))
    return(sra_list, sra_final_list)

#Function-3
# Function: download_fasta_by_taxon
# Description: dataset download is used to make database incase the user does not provide the user_database
# Parameters:
#   - args.taxon.fasta

def download_fasta_by_taxon(taxon, datasets_path):
    # Execute datasets command to get fasta sequences based on taxon
    command = f"{datasets_path} download genome taxon {shlex.quote(taxon)}"
    subprocess.call(command, shell = True)
    # Unzip the downloaded file
    with zipfile.ZipFile('ncbi_dataset.zip', 'r') as zip_ref:
        zip_ref.extractall('ncbi_dataset')

    # Construct the filename based on taxon name, replacing spaces with underscores
    filename = f'{args.taxon.replace(" ", "_")}.fasta' if ' ' in args.taxon else f'{args.taxon}.fasta'

    # Concatenate .fna files into taxon.fasta
    with open(filename, 'w') as outfile:
        for root, _, files in os.walk('ncbi_dataset/ncbi_dataset/data/'):
            for filename in files:
                if filename.endswith('.fna'):
                    filepath = os.path.join(root, filename)
                    with open(filepath, 'r') as infile:
                        outfile.write(infile.read())
#Function-4
# Function: create_blast_db
# Description: This function creates a BLAST database using either a user-provided FASTA file or a FASTA file formed after function-3.
# Parameters:
#   - args.user_db.fasta
#   - args.taxon.fasta

def create_blast_db(user_db=None):

  # Load configuration
    config = configparser.ConfigParser()
    config.read("config.ini")

  # Extract path to makeblastdb
    make_blast_db_path = config["make_blast_db"]["path"]

    if user_db:
       # Construct the makeblastdb command
       make_blast_command = make_blast_db_path + ' ' + "-in" + ' ' + args.user_db + ' ' + "-dbtype nucl"
       subprocess.call(make_blast_command, shell=True)
    else:
        # Check if the taxon is a single word or multiple words and construct the filename accordingly
        filename = f"{args.taxon.replace(' ', '_')}.fasta" if ' ' in args.taxon else f"{args.taxon}.fasta"

        # Check if the fasta file exists
        if not os.path.exists(filename):
            print(f"Taxon fasta file {filename} not found.")
            return

        # Construct the makeblastdb command using the taxon FASTA file
        make_blast_command = f"{make_blast_db_path} -in {filename} -dbtype nucl"
        subprocess.call(make_blast_command, shell=True)

#Function-5
# Function: analysis_fun
# Description: This function Runs in parallel for each sra_id to retrieve the raw reads in form of fastq format, later use the raw reads to find the target sequence using alignemnt tool BLAST
# Parameters:
#   - sra_id

def analysis_fun(sra_id):

    start_time = time.time() # Record the start time
    config = configparser.ConfigParser() # Load configuration file
    config.read("config.ini")
    use_sratools = True

# Download with sracha (fast), if selected
    if args.downloader == "sracha":
        sracha_path = config["sracha"]["path"]
        sracha_command = sracha_path + ' ' + "get" + ' ' + sra_id + ' ' + "-O" + ' ' + args.wd + "/FASTQ_FILES" + ' ' + "--split split-3 --no-gzip --no-progress"
        status = subprocess.call(sracha_command, shell=True)

        if status == 0 and (os.path.exists(args.wd + "/FASTQ_FILES/" + sra_id + "_1.fastq") or os.path.exists(args.wd + "/FASTQ_FILES/" + sra_id + ".fastq")):
            use_sratools = False
        else:
            print(sra_id + ": sracha failed, falling back to sra-tools")
# Prefetch the SRA data (sra-tools)
    if use_sratools:
        prefetch_path = config["prefetch"]["path"]
        prefetch_command = prefetch_path +  ' ' + "--max-size u" + ' ' + sra_id +' ' + "-O"+ ' '+ args.wd + "/SRR_FILES"
        subprocess.call(prefetch_command, shell=True)

# Convert the .sra file to a fastq file
        fastq_dump_path = config["fastq-dump"]["path"]
        fastqdump_command = fastq_dump_path + ' '+ "--split-3 --outdir" + ' ' + args.wd + "/FASTQ_FILES" + ' ' + args.wd + "/SRR_FILES/" + sra_id
        subprocess.call(fastqdump_command, shell=True)

# Concatenate the forward and reverse reads(Paired end reads)
    if os.path.exists(args.wd + "/FASTQ_FILES/" + sra_id + "_1.fastq"):
        concat = "cat" +  ' ' + args.wd+ "/FASTQ_FILES/" + sra_id + "_1.fastq" + ' ' + args.wd + "/FASTQ_FILES/" + sra_id + "_2.fastq"+ ' '  +">"+ ' ' + args.wd + "/BLAST/" + sra_id + "_readscat.fastq"
        subprocess.call(concat, shell = True)

# Copy fastq file to the same location as the paired-end file(Single end reads)
    else:
        fastqfasta = "cp" + ' ' + args.wd + "/FASTQ_FILES/" + sra_id + ".fastq"+ ' ' + args.wd + "/BLAST/" + sra_id + "_readscat.fastq"
        subprocess.call(fastqfasta, shell = True)
# Convert the fastq file to a fasta file
    fastq_fasta_cat= "sed -n '1~4s/^@/>/p;2~4p'" + ' ' + args.wd + "/BLAST/" + sra_id +"_readscat.fastq"+ ' '+ ">" + ' ' + args.wd + "/BLAST/" + sra_id + "_readscat.fasta"
    subprocess.call(fastq_fasta_cat, shell = True)

# Blast with the custom database provided by the user
    print("starting BLAST to custom database")
    start_time_bl = time.time()

    # Blast with the custom database provided by the user or taxon.fasta
    if args.user_db:
        blast_db = args.user_db
    else:
        # Use taxon.fasta as the database if user_db is not provided
        taxon_fasta = f"{args.taxon.replace(' ', '_')}.fasta" if ' ' in args.taxon else f"{args.taxon}.fasta"
        if not os.path.exists(taxon_fasta):
            print(f"Taxon fasta file {taxon_fasta} not found.")
            return
        blast_db = taxon_fasta

    blastn_nt_path = config["blastn"]["path"]
    blast_nt = blastn_nt_path + ' ' + "-db"  + ' '+ args.wd + "/" + blast_db +' ' + "-query"+ ' ' + args.wd + "/BLAST/" + sra_id + "_readscat.fasta" + ' ' + "-out" + ' '+ args.wd + "/BLAST_RESULT/" + sra_id + "_blastalignblast.txt" + ' '+ "-num_threads"+ ' ' + args.threads + ' ' + "-max_target_seqs 10 -evalue 10 -outfmt '6 qseqid sseqid sscinames sblastnames staxids pident length mismatch gapopen qstart qend sstart send stitle sskingdoms evalue bitscore'"
    subprocess.call(blast_nt, shell = True)

    end_time_bl = time.time()
    time_taken_bl = end_time_bl - start_time_bl
    time_taken_bl_minutes = int(time_taken_bl // 60)
    time_taken_bl_seconds = time_taken_bl % 60
    print("Custom BLAST for" + ' ' + sra_id + ' ' +":", time_taken_bl_minutes, "minutes and {:.2f} seconds".format(time_taken_bl_seconds))

# Extarcts unique entries based on the first field (column) from a BLAST result file.
# It uses the "sort" command with the options "-u" (unique) and "-k1,1" (sort based on the first field only).
    extract_top_blast = " sort -u -k1,1" +' ' + args.wd + "/BLAST_RESULT/" + sra_id + "_blastalignblast.txt > " + ' '+ args.wd + "/BLAST_RESULT/" + sra_id + "_alignblast.txt"
    subprocess.call(extract_top_blast, shell = True)

# Concatenate the contents of a BLAST result file,
# then filter the lines based on the sixth field (column) is greater than percentage identity.
# The filtered lines are then redirected to a new file.
    extract_top_blast1 = " cat" + ' ' + args.wd + "/BLAST_RESULT/"+ sra_id + "_alignblast.txt | awk -F '\t' '{if($6>" + args.pi + "){print$0}}' > " + ' ' + args.wd + "/BLAST_RESULT/" + sra_id + "_pi_alignblast.txt"
    subprocess.call(extract_top_blast1, shell = True)
    print("Extracted the hits greater than" + ' ' + args.pi + ' ' + "% into" + ' ' + sra_id + "_pi_alignblast.txt")

# Concatenate the contents of a BLAST result file,
# then filter the lines based on the seventh field (column) is greater than base-pair length size.
# The filtered lines are then redirected to a new file.
    extract_top_blast_bit =  "cat" +' '+ args.wd + "/BLAST_RESULT/" + sra_id  +  "_pi_alignblast.txt | awk -F '\t' '{if($7>" + args.len + "){print$0}}' > " + ' '+ args.wd + "/BLAST_RESULT/" + sra_id + "_bp_alignblast.txt"
    subprocess.call(extract_top_blast_bit, shell = True)
    print("Extracted the hits greater than" + ' ' + args.len + ' ' + "base pair into" + ' ' + sra_id + "_bp_alignblast.txt")

#  Making files to run sam2.sh
# Making different folders for each sra_id to run sam2.sh
    mkdir_key ="mkdir" + ' ' + sra_id +"_sam"
    subprocess.call(mkdir_key, shell = True)

# Copying the extracted BLAST result based on percentage identity/ base-pair length file in a specific sam folder for the sra_id
    keyword_txt = "cp" + ' ' + args.wd + "/BLAST_RESULT/" + sra_id + "_bp_alignblast.txt" + ' ' + sra_id +"_sam/"+"keyword.txt"
    subprocess.call(keyword_txt, shell = True)

# Extract the first column from tab-seperated file for specific sra run.
    good_contigs = " awk -F '\t' '{print $1}' "  + ' '+ args.wd + "/BLAST_RESULT/" + sra_id + "_bp_alignblast.txt > " + ' ' + args.wd + "/BLAST_RESULT/" + sra_id + "_align_blast_ID.txt"
    subprocess.call(good_contigs, shell = True)

# Copying the file with 1st column extracted in a specific sam folder for the sra_id
    cp_keywordID =  "cp " + ' '+ args.wd + "/BLAST_RESULT/" + sra_id + "_align_blast_ID.txt" + ' ' + sra_id +"_sam/"+"keyword_IDlist.txt"
    subprocess.call(cp_keywordID, shell = True)

# Extracting fasta file with unique identifiers
    grab_id = r"perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' " + ' '+ args.wd + "/BLAST_RESULT/"  + sra_id + "_align_blast_ID.txt"+ ' '+ args.wd + "/BLAST/" +  sra_id + "_readscat.fasta > " +' '+ args.wd + "/BLAST_RESULT/"  + sra_id + "_align.fa"
    subprocess.call(grab_id, shell = True)

# Copying the file made with extracting unique identifiers in a specific sam folder for the sra_id
    keyword_contigs =  "cp " + ' '+ args.wd + "/BLAST_RESULT/"  + sra_id + "_align.fa" + ' ' + sra_id+"_sam/"+"keywords_contigs.fa"
    subprocess.call(keyword_contigs, shell = True)

# Remove the space from the keyword.txt
    remove_spaces = "sed -i 's/ /_/g'" + ' ' + sra_id+"_sam/"+"keyword.txt"
    subprocess.call(remove_spaces, shell = True)

# Copying the sam3.sh into all the sam folders for every sra_id
    cp_sam3 = "cp" + ' '+"sam3.sh"+ ' ' + sra_id+"_sam/sam3.sh"
    subprocess.call(cp_sam3, shell = True)

# Changing file permissions for the files required for the execution of sam2.sh
    perm = "chmod 777" + ' '+  sra_id +"_sam/"+"keyword.txt | chmod 777" +' ' + sra_id +"_sam/keyword_IDlist.txt | chmod 777" + ' '+ sra_id+"_sam/keywords_contigs.fa | chmod 777" + ' ' + sra_id+"_sam/sam3.sh"
    subprocess.call(perm, shell =True)

# Change the current working directory to args.wd+"/"+sra_id +"_sam"
    os.chdir(args.wd+"/"+sra_id +"_sam")

# Run the bash script 'sam3.sh' using subprocess
    subprocess.run(['bash','sam3.sh'])

# Change the current working directory to the specified directory (args.wd)
    os.chdir(args.wd)

# Construct the awk command to remove duplicate entries from a FASTA file
    rm_duplicates = "awk '/^>/{f=!d[$1];d[$1]=1}f' " +  ' ' + sra_id+"_sam/sam_duplicates.fa >" + ' '+args.wd +"/BLAST_RESULT/" + sra_id + "_samalign.fa"
    subprocess.call(rm_duplicates, shell = True)

# Blast with the nucleotide database provided by the user

    print("Starting blast to nt database")
    start_time_bl2 = time.time()
# Changing permissions to read, write and execute
    perm1 = "chmod 777" + ' '+  args.wd + "/BLAST_RESULT/" + sra_id + "_samalign.fa"
    subprocess.call(perm1, shell = True)

    blastn_nt_path = config["blastn"]["path"]

    blastn_nt_path_nt = config["nt"]["path"]
# Threshold size in bytes
    threshold_size_bytes = 50 * 1024  # 150 kilobytes
# Check the size of the input fasta file
    input_fasta_file = os.path.join(args.wd, "BLAST_RESULT", sra_id + "_samalign.fa")
    input_file_size = os.path.getsize(input_fasta_file)

# Determine whether to split the fasta file or not
    if input_file_size < threshold_size_bytes:
# Initialize blast_results_paths as an empty list
        blast_results_paths = []
# Skip splitting and directly perform BLAST
        blast_output_file = os.path.join(args.wd, 'BLAST_NT_RESULTS', f"{sra_id}_samalign_BLASTtont.txt")
        blast_command = (
            f"{blastn_nt_path} -db {blastn_nt_path_nt} -query {input_fasta_file} "
            f"-out {blast_output_file} "
            f"-num_threads {args.threads} "
            "-max_target_seqs 10 -evalue 10 -outfmt '6 qseqid sseqid sscinames sblastnames staxids "
            "pident length mismatch gapopen qstart qend sstart send stitle sskingdoms evalue bitscore'"
        )
        subprocess.call(blast_command, shell=True)
        blast_results_paths.append(blast_output_file)

    else:
# Split fasta file

# Initialize blast_results_paths as an empty list
        blast_results_paths = []

    # Split fasta file
        fasta_splitter_sam_command = (
            "./fasta-splitter.pl --n-parts 200 --measure count " +
            f"{input_fasta_file} --out-dir {os.path.join(args.wd, 'BLAST_SPLIT_NT/')}"
        )
        subprocess.call(fasta_splitter_sam_command, shell=True)

        # Get the number of split files
        split_files = [f for f in os.listdir(os.path.join(args.wd, 'BLAST_SPLIT_NT')) if f.endswith('.fa')]
        num_splits = len(split_files)

    # Perform BLAST on each split fasta file
        for i in range(1, num_splits + 1):
            blast_output_file = os.path.join(args.wd, 'BLAST_NT_RESULTS', f"{sra_id}_samalign_BLASTtont_{i:03}.txt")
            blast_command = (
                f"{blastn_nt_path} -db {blastn_nt_path_nt} -query "
                f"{os.path.join(args.wd, 'BLAST_SPLIT_NT', sra_id + '_samalign.part-')}{'{:03}'.format(i)}.fa "
                f"-out {blast_output_file} "
                f"-num_threads {args.threads} "
                "-max_target_seqs 10 -evalue 10 -outfmt '6 qseqid sseqid sscinames sblastnames staxids "
                "pident length mismatch gapopen qstart qend sstart send stitle sskingdoms evalue bitscore'"
            )
            subprocess.call(blast_command, shell=True)
# Combine BLAST results into one file
# Append the path of the BLAST result file to the list
        blast_results_paths.append(blast_output_file)

# Combine all BLAST result files into a single file
    output_combined_file = os.path.join(args.wd, 'BLAST_NT_RESULTS', f"{sra_id}_samalign_BLASTtont_combined.txt")
    with open(output_combined_file, 'w') as combined_file:
        for result_file in blast_results_paths:
            with open(result_file, 'r') as f:
                combined_file.write(f.read() + '\n')

    for result_file in blast_results_paths:
        os.remove(result_file)

    print("blast to nucleotide database finished ")

    end_time_bl2 = time.time()
    time_taken_bl2 = end_time_bl2 - start_time_bl2
    time_taken_bl2_minutes = int(time_taken_bl2 // 60)
    time_taken_bl2_seconds = time_taken_bl2 % 60
    print("NT BLAST for" + ' ' + sra_id + ' ' +":", time_taken_bl2_minutes, "minutes and {:.2f} seconds".format(time_taken_bl2_seconds))


# Extarcts unique entries based on the first field (column) from a BLAST result file.
# It uses the "sort" command with the options "-u" (unique) and "-k1,1" (sort based on the first field only).

    grab_top = "sort -u -k1,1" + ' ' + args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_samalign_BLASTtont_combined.txt > " + ' '+ args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_aligntoNT.txt"
    subprocess.call(grab_top, shell = True)

# Grep keyword
    grep_keyword = shlex.quote(args.grep_keyword)
    grab_taxon = f"grep -E {grep_keyword}" + ' ' + args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_aligntoNT.txt > "+ ' ' + args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_grep_align_toNT.txt"
    subprocess.call(grab_taxon, shell = True)

# Extract the first column from tab-seperated file for specific sra run
    grab_id_taxon = "awk -F '\t' '{print $1}'" +' ' + args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_grep_align_toNT.txt > " +' '+ args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_grepID_align_toNT.txt"
    subprocess.call(grab_id_taxon, shell = True)

# this command seems to be extracting the first field from a specific file and saving the output to another file
    grab_cut =  " cut -f 1 -d :" +' '+args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_grepID_align_toNT.txt > "+ ' '+ args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_grepIDb_align_toNT.txt"
    subprocess.call(grab_cut, shell = True)

# Extracting fasta file with unique identifiers
    final = r"perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV'" + ' ' + args.wd + "/BLAST_NT_RESULTS/" + sra_id + "_grepIDb_align_toNT.txt" +' ' + args.wd + "/BLAST/" + sra_id + "_readscat.fasta > " +' '+ args.wd + "/BLAST_NT_RESULTS/"  + sra_id + "_final.fa"
    subprocess.call(final, shell = True)

# Results
    move_final_contigs = "cp" + ' ' + args.wd +"/BLAST_NT_RESULTS/" + sra_id + "_final.fa" + ' ' + args.wd + "/Potential_positive_SRA_Runs/"
    subprocess.call(move_final_contigs, shell = True)
    delete_empty = "find" + ' ' + args.wd +"/Potential_positive_SRA_Runs/ -type f -empty -print -delete "
    subprocess.call(delete_empty, shell = True)
    move_blast_table = "cp" + ' ' + args.wd +  "/BLAST_NT_RESULTS/"  + sra_id + "_final.fa" + ' ' + args.wd + "/blast_table/"
    subprocess.call(move_blast_table, shell = True)
    print("BLAST TABLE containing SRA Runs positive for" + ' ' + args.grep_keyword)

# Delete all the unnecessary files and folders
    files_to_remove = [
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_blastalignblast.txt"),
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_alignblast.txt"),
    os.path.join(args.wd, "BLAST", f"{sra_id}_readscat.fasta"),
    os.path.join(args.wd, "BLAST", f"{sra_id}_readscat.fastq"),
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_pi_alignblast.txt"),
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_bp_alignblast.txt"),
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_align_blast_ID.txt"),
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_align.fa"),
    os.path.join(args.wd, "BLAST_RESULT", f"{sra_id}_samalign.fa"),
    os.path.join(args.wd, "BLAST_NT_RESULTS", f"{sra_id}_aligntoNT.txt"),
    os.path.join(args.wd, "BLAST_NT_RESULTS", f"{sra_id}_grep_align_toNT.txt"),
    os.path.join(args.wd, "BLAST_NT_RESULTS", f"{sra_id}_grepID_align_toNT.txt"),
    os.path.join(args.wd, "BLAST_NT_RESULTS", f"{sra_id}_grepIDb_align_toNT.txt")
    ]

# Iterate over files to remove
    for file_to_remove in files_to_remove:
        if os.path.exists(file_to_remove):
            os.remove(file_to_remove)
            print(f"Deleted: {file_to_remove}")
        else:
            print(f"The file {file_to_remove} does not exist.")

# Delete SAM files
    sam_folder = os.path.join(args.wd, sra_id + "_sam")
    shutil.rmtree(sam_folder)


# Construct the pattern to match files
    pattern = os.path.join(args.wd, "BLAST_SPLIT_NT", f"{sra_id}_samalign*.fa")

# Get a list of files matching the pattern
    files_to_remove_1 = glob.glob(pattern)

# Iterate over the matched files and remove them
    for file_to_remove_1 in files_to_remove_1:
        os.remove(file_to_remove_1)
        print(f"Deleted: {file_to_remove_1}")

# Construct the pattern to match files in BLAST_NT_RESULTS directory
    pattern_2 = os.path.join(args.wd, "BLAST_NT_RESULTS", f"{sra_id}_samalign_BLASTtont_*.txt")

# Get a list of files matching the pattern in BLAST_NT_RESULTS directory
    files_to_remove_2 = glob.glob(pattern_2)

# Iterate over the matched files in BLAST_NT_RESULTS directory and remove them
    for file_to_remove_2 in files_to_remove_2:
        os.remove(file_to_remove_2)
        print(f"Deleted: {file_to_remove_2}")

    end_time = time.time()
    time_taken = end_time - start_time
    minutes, seconds = divmod(time_taken, 60)
    print(f"Processing {sra_id} on core {rank}")
    print(f"1st 2-step BLAST SRA ID: {sra_id}, Started: {time.ctime(start_time)}, Ended: {time.ctime(end_time)}, Time taken: {int(minutes)}m {seconds:.2f}s")


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    start_time = time.time()
    max_runtime_seconds = 48 * 3600
    stop_before_seconds = max_runtime_seconds - 1800

    # File to track progress
    progress_file = os.path.join(args.wd, "processed_sra_ids.txt")

    if rank == 0:
        for directory in directories:
            pathExist(directory)

        id_list = get_ncbi_ids(args.mail, 'sra', args.keyword, args.field, args.maximum_returned_items, args.threads)
        sra_list, sra_final_list = get_sra_runs(id_list, args.db)

        if args.user_db:
            create_blast_db(args.user_db)
        elif args.taxon:
            config = configparser.ConfigParser()
            config.read("config.ini")
            datasets_path = config.get("datasets", "path")
            download_fasta_by_taxon(args.taxon, datasets_path)
            create_blast_db()
        else:
            print("Please provide either --user_db or a taxon name.")

        # Read already processed SRA IDs if resuming
        if os.path.exists(progress_file):
            with open(progress_file) as f:
                processed_sra_ids = set(line.strip() for line in f if line.strip())
            sra_final_list = [x for x in sra_final_list if x not in processed_sra_ids]
            print(f"Resuming: {len(sra_final_list)} remaining SRA IDs.")
        else:
            processed_sra_ids = set()

    else:
        sra_final_list = None

    # Broadcast SRA list to all ranks
    sra_final_list = comm.bcast(sra_final_list, root=0)

    # Distribute IDs among ranks
    chunk_size = len(sra_final_list) // size
    remainder = len(sra_final_list) % size
    start_index = rank * chunk_size + min(rank, remainder)
    end_index = start_index + chunk_size + (1 if rank < remainder else 0)
    sra_final_list_chunk = sra_final_list[start_index:end_index]
    print(f"Rank {rank} processing {len(sra_final_list_chunk)} SRA IDs: {sra_final_list_chunk}")

    # Each rank keeps track of its processed IDs
    processed_local = []

    # Process each SRA ID
    for sra_id in sra_final_list_chunk:
        elapsed_time = time.time() - start_time
        if elapsed_time >= stop_before_seconds:
            print(f"Rank {rank}: Reached time limit (~47.5h). Saving progress and exiting gracefully.")
            break

        try:
            analysis_fun(sra_id)
            processed_local.append(sra_id)

        except Exception as e:
            print(f"Rank {rank}: Error processing {sra_id}: {e}")
            continue

    # Gather all processed IDs from all ranks
    all_processed = comm.gather(processed_local, root=0)

    # Rank 0 writes combined processed IDs to file
    if rank == 0:
        flat_list = [item for sublist in all_processed for item in sublist]
        with open(progress_file, "a") as f:
            for sra_id in flat_list:
                f.write(sra_id + "\n")
            f.flush()
            os.fsync(f.fileno())

        print(f"Updated processed SRA list written to {progress_file}")

    comm.Barrier()

    # Final memory and runtime logging (only rank 0)
    if rank == 0:
        end_time = time.time()
        total_runtime = end_time - start_time
        minutes, seconds = divmod(total_runtime, 60)

        process = psutil.Process(os.getpid())
        memory_gb = process.memory_info().rss / (1024 ** 3)

        print(f"Script finished. Runtime: {int(minutes)}m {seconds:.2f}s, Memory: {memory_gb:.2f} GB")

        # Log performance
        csv_file = os.path.join(args.wd, "script_performance_log.csv")
        file_exists = os.path.isfile(csv_file)
        with open(csv_file, mode='a', newline='') as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow(["Date", "Runtime_seconds", "Runtime_minutes", "Memory_GB"])
            writer.writerow([time.strftime("%Y-%m-%d %H:%M:%S"), total_runtime, total_runtime / 60, memory_gb])

        print(f"Performance log saved to {csv_file}")

    sys.stdout.flush()
    sys.stderr.flush()
    comm.Barrier()

