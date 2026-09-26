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
   SraMiner - version 0.3(03/25/26)-Part2

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
from prettytable import PrettyTable
from mpi4py import MPI

parser = argparse.ArgumentParser()
parser.add_argument('--wd', type=str)
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
    args.wd + "/BLAST_CT_AFTER_ASSEMBLY",
    args.wd + "/BLAST_AFTER_ASSEMBLY",
    args.wd + "/Positive_SRA",
    args.wd + "/ASSEMBLY",
    args.wd + "/BLAST_NT_AFTER_ASSEMBLY"
]

#Function-0
#Delete SRR folder
def delete_srr_folder():
    # Directory name
    directory = "SRR_FILES"
    path = os.path.join(args.wd, directory)
#Remove the Directory
    shutil.rmtree(path)

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
# Function: fetch_sra_metadata_for_Potential_positive_SRA_Runs
# Description: uses the Entrez.esummary to take parameters such as email, id_list and Positive folder location
# Parameters:
#   - email
#   - id_list
#   - positive_sra_folder
# Returns:
# Metadata

def fetch_sra_metadata_for_Potential_positive_SRA_Runs(email, id_list, positive_sra_folder):
    Entrez.email = email
    retmax = 10000
    combined_records = []
    metadata_dict = {}

    # Splitting id_list into batches of 10,000 IDs
    for i in range(0, len(id_list), 10000):
        id_batch = id_list[i:i+10000]
        id_str = ','.join(id_batch)

        try:
            # Fetching metadata for the current batch of IDs
            handle = Entrez.esummary(db="sra", id=id_str)
            records = Entrez.read(handle)
            handle.close()

            combined_batch_records = [{"Runs": record["Runs"], "ExpXml": record["ExpXml"]} for record in records]
            combined_records.extend(combined_batch_records)

            for record in combined_batch_records:
                sra_id = record['Runs'].split('acc="')[1].split('"')[0]
                if os.path.exists(os.path.join(positive_sra_folder, sra_id + "_target.fa")):
                    metadata = {
                        'Sra Id': sra_id,
                        'Title': '',
                        'Platform': '',
                        'Organism': '',
                        'Library Name': '',
                        'Library Strategy': '',
                        'Library Source': '',
                        'Bioproject': '',
                        'Biosample': ''
                    }
                    if '<Title>' in record['ExpXml']:
                        metadata['Title'] = record['ExpXml'].split('<Title>')[1].split('</Title>')[0]
                    if 'instrument_model="' in record['ExpXml']:
                        metadata['Platform'] = record['ExpXml'].split('instrument_model="')[1].split('"')[0]
                    if 'ScientificName="' in record['ExpXml']:
                        metadata['Organism'] = record['ExpXml'].split('ScientificName="')[1].split('"')[0]
                    if '<LIBRARY_NAME>' in record['ExpXml']:
                        metadata['Library Name'] = record['ExpXml'].split('<LIBRARY_NAME>')[1].split('</LIBRARY_NAME>')[0]
                    if '<LIBRARY_STRATEGY>' in record['ExpXml']:
                        metadata['Library Strategy'] = record['ExpXml'].split('<LIBRARY_STRATEGY>')[1].split('</LIBRARY_STRATEGY>')[0]
                    if '<LIBRARY_SOURCE>' in record['ExpXml']:
                        metadata['Library Source'] = record['ExpXml'].split('<LIBRARY_SOURCE>')[1].split('</LIBRARY_SOURCE>')[0]
                    if '<Bioproject>' in record['ExpXml']:
                        metadata['Bioproject'] = record['ExpXml'].split('<Bioproject>')[1].split('</Bioproject>')[0]
                    if '<Biosample>' in record['ExpXml']:
                        metadata['Biosample'] = record['ExpXml'].split('<Biosample>')[1].split('</Biosample>')[0]

                    metadata_dict[sra_id] = metadata

        except Exception as e:
            # Handle any exceptions or errors that may occur during the fetching process
            print(f"Error fetching batch {i}-{i + retmax}: {str(e)}")
            return None  # Return None if fetching fails

    return metadata_dict or {}  # Return an empty dictionary if no data fetched

#Function-3
# Function: fetch_sra_metadata_biosample
# Description: uses the Entrez.efetch to take parameters such as email, id_list
# Parameters:
#   - email
#   - id_list
# Returns:
# Biosample_metadata
def fetch_sra_metadata_biosample(email, biosample_ids):
    Entrez.email = email
    retmax = 10000
    biosample_metadata = {}

    # Splitting biosample_ids into batches of 10,000 IDs
    for i in range(0, len(biosample_ids), 10000):
        id_batch = biosample_ids[i:i+10000]
        id_str = ','.join(id_batch)

        try:
            handle = Entrez.efetch(db="biosample", id=id_str, retmode="xml")
            xml_data = handle.read()
            handle.close()
            root = ET.fromstring(xml_data)
#            print(ET.tostring(root, encoding='utf8').decode('utf8'))

            for bio_root in root.findall(".//BioSample"):
                biosample_acc = bio_root.attrib.get("accession")
                biosample_uid = bio_root.attrib.get("id")

                metadata = {
                    "Biosample": biosample_acc,
                    "BioSample_UID": biosample_uid,
                }

                org_name = bio_root.findtext(".//Description/Organism/OrganismName")
                if org_name:
                    metadata["Organism"] = org_name

                for id_tag in bio_root.findall(".//Ids/Id"):
                    key = id_tag.attrib.get("db") or id_tag.attrib.get("db_label")
                    if key and id_tag.text:
                        metadata[key] = id_tag.text


                for attr in bio_root.findall(".//Attributes/Attribute"):
                    key = attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or attr.attrib.get("display_name")
                    if key:
                        value = (attr.text or "").strip()
                        if value:
                            # If the same key appears multiple times, store as list
                            if key in metadata:
                                if isinstance(metadata[key], list):
                                    metadata[key].append(value)
                                else:
                                    metadata[key] = [metadata[key], value]
                            else:
                                metadata[key] = value

                if biosample_acc:
                    biosample_metadata[biosample_acc] = metadata

        except Exception as e:
            print(f"Error fetching batch {i}-{i + batch_size}: {str(e)}")
    print(biosample_metadata)
    return biosample_metadata



#Function-3
# Function: export_metadata_to_excel
# Description: Combine the biosample and the run metadata into a single txt file for a particular SRA ID
# Parameters:
#   - sra_metadata
#   - biosample_metadata
# Returns:
# Combined Metadata table
def export_metadata_to_excel(sra_metadata, biosample_metadata, filename="metadata_export.xlsx"):
    combined_list = []
    all_fields = set()

    # 1) Combine metadata
    for sra_id, sra_md in sra_metadata.items():
        biosample_acc = sra_md.get("Biosample")
        combined = dict(sra_md)
        combined.update(biosample_metadata.get(biosample_acc, {}))

        # Ensure Sra Id is inside the dictionary for the row
        combined['Sra Id'] = sra_id
        combined_list.append(combined)
        all_fields.update(combined.keys())

    # 2) Create DataFrame
    df = pd.DataFrame(combined_list)

    # 3) Define column order
    priority_fields = [
        'Sra Id', 'Title', 'Platform', 'Organism',
        'Library Strategy', 'Library Source', 'Bioproject', 'Biosample'
    ]

    # Get extra fields and sort them
    extra_fields = sorted([f for f in all_fields if f not in priority_fields])
    final_column_order = [f for f in (priority_fields + extra_fields) if f in df.columns]

    # Reorder columns
    df = df[final_column_order]

    # 4) Clean up lists (Excel doesn't like Python lists in cells)
    for col in df.columns:
        df[col] = df[col].apply(lambda x: "; ".join(map(str, x)) if isinstance(x, list) else x)

    df.to_excel(filename, index=False)
    print(f"Success! Metadata saved to {filename}")


def index(user_db):
    config = configparser.ConfigParser()
    config.read("config.ini")

    # Get BWA path
    bwa_path = config["bwa"]["path"]

    # Build and run command
    bwa_index = f"{bwa_path} index {user_db}"
    subprocess.call

#Function-3
# Function: assembly
# Description: Run Assembly followed by 2- step BLAST

def assembly(sra_id):
    start_time_p1 = time.time()

    config = configparser.ConfigParser() # Load configuration file
    config.read("config.ini")

    bwa_path = config["bwa"]["path"]
    samtools_path = config["samtools"]["path"]
    megahit_path = config["megahit"]["path"]

    bwa_se = f"{bwa_path}"  + ' ' + "mem -t" + ' ' + args.threads + ' '+ "-v 3 -B 1 -E 1 -O 1 -k 17" + ' ' + args.user_db + ' '+ args.wd + "/FASTQ_FILES/" + sra_id + ".fastq > " +args.wd + "/ASSEMBLY/"+ sra_id + "_to_targetSE.sam2"
    subprocess.call(bwa_se, shell = True)

    grab_map_se = f"{samtools_path}" + ' ' +"view -F4 -h -@ 46"+ ' '+args.wd + "/ASSEMBLY/"+  sra_id + "_to_targetSE.sam2 > " + ' '+args.wd + "/ASSEMBLY/"+  sra_id + "_to_targetSE_mapped.sam"
    subprocess.call(grab_map_se, shell = True)

    samfastqSE = f"{samtools_path}" + ' ' +"fastq" + ' ' + args.wd + "/ASSEMBLY/" + sra_id +  "_to_targetSE_mapped.sam > " + ' ' + args.wd + "/ASSEMBLY/" + sra_id + "_to_targetSE_single.fastq"
    subprocess.call(samfastqSE, shell = True)

    metaspades = f"{megahit_path}" + ' ' + "-r" + ' ' + args.wd + "/ASSEMBLY/" + sra_id + "_to_targetSE_single.fastq -o" + ' ' + args.wd + "/ASSEMBLY/" + sra_id + "_to_target_assembly"
    subprocess.call( metaspades, shell = True)

    end_time_p1 = time.time()
    time_taken_p1 = end_time_p1 - start_time_p1
    minutes, seconds = divmod(time_taken_p1, 60)

    print(f" Assembly for SRA ID: {sra_id}, Started: {time.ctime(start_time_p1)}, Ended: {time.ctime(end_time_p1)}, Time taken: {int(minutes)}m {seconds:.2f}s")

    start_time = time.time()  # Record the start time for this SRA ID
    config = configparser.ConfigParser() # Load configuration file
    config.read("config.ini")
    blastn_path = config["blastn"]["path"]

    print("starting 2nd BLAST to custom database after assembly")

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

    blast_nt = blastn_path + ' ' + "-db"  + ' '+ args.wd + "/" + blast_db +' ' + "-query" + ' ' + args.wd + "/ASSEMBLY/" + sra_id + "_to_target_assembly/final.contigs.fa"+ ' ' + "-out" + ' '+ args.wd + "/BLAST_AFTER_ASSEMBLY/" + sra_id + "_scaffolds.txt -num_threads" + ' ' + args.threads + ' ' + "-max_target_seqs 10 -evalue 10 -outfmt '6 qseqid sseqid sscinames sblastnames staxids pident length mismatch gapopen qstart qend sstart send stitle sskingdoms evalue bitscore'"
    subprocess.call(blast_nt, shell = True)

#extract the top blast hit to the custom datatbase
    extract_top_blast = " sort -u -k1,1" +' ' + args.wd + "/BLAST_AFTER_ASSEMBLY/" + sra_id +"_scaffolds.txt > " + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_scaffolds_align.txt"
    subprocess.call(extract_top_blast, shell = True)

# Concatenate the contents of a BLAST result file,
# then filter the lines based on the sixth field (column) is greater than percentage identity.
# The filtered lines are then redirected to a new file.
    extract_top_blast1 = " cat" + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_scaffolds_align.txt | awk -F '\t' '{if($6>" + args.pi + "){print$0}}' > " + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id + "_pi_alignblast.txt"
    subprocess.call(extract_top_blast1, shell = True)
    print("Extracted the hits greater than" + ' ' + args.pi + ' ' + "% into" + ' ' + sra_id + "_pi_alignblast.txt")

#filter out the short blast hit
    extract_top_blast1 = " cat" + ' '+  args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_pi_alignblast.txt | awk -F '\t' '{if($7>" + args.len + "){print$0}}' > " + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_grep_bp_scaffolds.txt"
    subprocess.call(extract_top_blast1, shell = True)

#make directory to save all the inputs required by the sam2.sh
    mkdir_key ="mkdir" + ' ' + sra_id +"_sam2"
    subprocess.call(mkdir_key, shell = True)

#filter to get keywords.txt
    key1 = "cp" + ' '+ args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_grep_bp_scaffolds.txt"+ ' ' + sra_id +"_sam2/" + "keywords.txt"
    subprocess.call(key1, shell = True)

#from the hit list extract just the ID list of good contigs
    extract_id = "awk -F '\t' '{print $1}'" + ' ' +args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_grep_bp_scaffolds.txt > " + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_grep_scaffolds_ID.txt"
    subprocess.call(extract_id, shell = True)

    key2 = "cp" + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_grep_scaffolds_ID.txt" + ' ' + sra_id +"_sam2/" + "keyword_IDlist.txt"
    subprocess.call(key2, shell = True)

#grab the good contigs"
    extract_contigs = r"perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV'" + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id +"_grep_scaffolds_ID.txt" +' ' + args.wd + "/ASSEMBLY/" + sra_id + "_to_target_assembly/" + "final.contigs.fa > " + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id + "_grep_scaffolds.fasta"
    subprocess.call(extract_contigs, shell = True)
    key3 = "cp" + ' '+ args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id + "_grep_scaffolds.fasta" + ' ' + sra_id +"_sam2/"+ "keywords_contigs.fa"
    subprocess.call(key3, shell = True)

    rm_key1 = " sed -i 's/ /_/g'" + ' ' + sra_id +"_sam2/"+ "keywords.txt"
    subprocess.call(rm_key1, shell = True)
#copy sam2.sh in each directory
    cp_sam2 = "cp" + ' '+"sam2.sh"+ ' ' + sra_id+"_sam2/sam2.sh"
    subprocess.call(cp_sam2, shell = True)
    os.chdir(args.wd+"/"+sra_id +"_sam2")
    subprocess.run(['bash','sam2.sh'])
    os.chdir(args.wd)
#remove duplicates in sam2.sh
    rm_sam = "awk '/^>/{f=!d[$1];d[$1]=1}f'" + ' ' + sra_id +"_sam2/sam_duplicates.fasta > " + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id + "_samcontigs.fasta"
    subprocess.call(rm_sam, shell = True)

    blastn_nt = config["nt"]["path"]

    blast_nt_ass = blastn_path + ' ' + "-db"  + ' '+ blastn_nt  +' ' + "-query" + ' ' + args.wd + "/BLAST_CT_AFTER_ASSEMBLY/" + sra_id + "_samcontigs.fasta -out" + ' '+ args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_samcontigs.txt -num_threads"+ ' ' + args.threads + ' ' + "-max_target_seqs 10 -evalue 10 -outfmt '6 qseqid sseqid sscinames sblastnames staxids pident length mismatch gapopen qstart qend sstart send stitle sskingdoms evalue bitscore'"
    subprocess.call(blast_nt_ass, shell = True)
    grep_keyword = shlex.quote(args.grep_keyword)

    grab_1 =  f"grep -E {grep_keyword}" + ' ' +  args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id +"_samcontigs.txt > " + ' '+   args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id +"_grep_samcontigs.txt"
    subprocess.call(grab_1, shell = True)
    grab_2 = "  awk -F '\t' '{print $1}'" + ' ' + args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id +"_grep_samcontigs.txt > " +  ' ' + args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_grep_samcontigsID.txt"
    subprocess.call(grab_2, shell = True)
    grab_3 = " cut -f 1 -d : " + ' ' + args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_grep_samcontigsID.txt > " + ' ' +  args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_grep_samcontigsIDb.txt"
    subprocess.call(grab_3, shell = True)
    grab_4 = r" perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV'" + ' ' + args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_grep_samcontigsIDb.txt"+ ' ' + args.wd + "/ASSEMBLY/" + sra_id + "_to_target_assembly/final.contigs.fa > " + ' ' + args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_target.fa"
    subprocess.call(grab_4, shell = True)

    move_final_contigs = "cp" + ' ' +  args.wd + "/BLAST_NT_AFTER_ASSEMBLY/" + sra_id + "_target.fa" + ' ' + args.wd + "/Positive_SRA/"
    subprocess.call(move_final_contigs, shell = True)
    delete_empty = "find" + ' ' + args.wd +"/Positive_SRA/ -type f -empty -print -delete "
    subprocess.call(delete_empty, shell = True)

    end_time = time.time()
    time_taken = end_time - start_time
    minutes, seconds = divmod(time_taken, 60)

# Files to remove
    files_to_remove = [
    os.path.join(args.wd, "BLAST_AFTER_ASSEMBLY", f"{sra_id}_scaffolds.txt"),
    os.path.join(args.wd, "BLAST_CT_AFTER_ASSEMBLY", f"{sra_id}_scaffolds_align.txt"),
 #   os.path.join(args.wd, "BLAST_CT_AFTER_ASSEMBLY", f"{sra_id}_grep_bp_scaffolds.txt"),
    os.path.join(args.wd, "BLAST_CT_AFTER_ASSEMBLY", f"{sra_id}_grep_scaffolds_ID.txt"),
    os.path.join(args.wd, "BLAST_CT_AFTER_ASSEMBLY", f"{sra_id}_grep_scaffolds.fasta"),
    os.path.join(args.wd, "BLAST_CT_AFTER_ASSEMBLY", f"{sra_id}_samcontigs.fasta"),
    os.path.join(args.wd, "BLAST_NT_AFTER_ASSEMBLY", f"{sra_id}_samcontigs.txt"),
#    os.path.join(args.wd, "BLAST_NT_AFTER_ASSEMBLY", f"{sra_id}_grep_samcontigs.txt"),
    os.path.join(args.wd, "BLAST_NT_AFTER_ASSEMBLY", f"{sra_id}_grep_samcontigsID.txt"),
    os.path.join(args.wd, "BLAST_NT_AFTER_ASSEMBLY", f"{sra_id}_grep_samcontigsIDb.txt"),
    os.path.join(args.wd, "FASTQ_FILES", f"{sra_id}_2.fastq"),
    os.path.join(args.wd, "FASTQ_FILES", f"{sra_id}_1.fastq"),
    os.path.join(args.wd, "FASTQ_FILES", f"{sra_id}.fastq")
    ]

# Iterate over files to remove
    for file_to_remove in files_to_remove:
        if os.path.exists(file_to_remove):
            os.remove(file_to_remove)
            print(f"Deleted: {file_to_remove}")
        else:
            print(f"The file {file_to_remove} does not exist.")
    # Delete SAM files
    sam_folder = os.path.join(args.wd, sra_id + "_sam2")
    shutil.rmtree(sam_folder)


    print(f"2nd 2-Step BLAST for SRA ID: {sra_id}, Started: {time.ctime(start_time)}, Ended: {time.ctime(end_time)}, Time taken: {int(minutes)}m {seconds:.2f}s")

#Function
#Delete FASTQ folder
def delete_fastqfiles_folder():
    # Directory name
    directory = "FASTQ_FILES"
    path = os.path.join(args.wd, directory)
#Remove the Directory
    shutil.rmtree(path)
#Function-0
#Delete assembly
def delete_assembly():
    # Directory name
    directory = "ASSEMBLY"
    path = os.path.join(args.wd, directory)
#Remove the Directory
    shutil.rmtree(path)


if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print("MPI ranks:", size)

    start_time = time.time()
    max_runtime_seconds = 48 * 3600       # 48-hour HPC limit
    stop_before_seconds = max_runtime_seconds - 1800  # Stop 30 mins early

    # File names specific to this script
    progress_file = os.path.join(args.wd, "processed_sra_assembly_ids.txt")

    log_file = os.path.join(args.wd, "assembly_performance_log.csv")
    metadata_file = os.path.join(args.wd, "metadata_table.txt")

    if rank == 0:
        delete_srr_folder()
        for directory in directories:
            pathExist(directory)

        id_list = get_ncbi_ids(args.mail, 'sra', args.keyword, args.field,
                               args.maximum_returned_items, args.threads)
        index(args.user_db)

        # Collect SRA run IDs from folder
        sra_final_list = []

        for root, dirs, files in os.walk('./Potential_positive_SRA_Runs'):
            sra_final_list = [filename.split('_')[0] for filename in files if filename.endswith('_final.fa')]

        print(f"SRA Final List ({len(sra_final_list)}):", sra_final_list)

        # Resume: read already processed SRA IDs
        if os.path.exists(progress_file):
            with open(progress_file) as f:
                processed_ids = set(line.strip() for line in f if line.strip())
            sra_final_list = [x for x in sra_final_list if x not in processed_ids]
            print(f"Resuming: {len(sra_final_list)} remaining SRA IDs.")
        else:
            processed_ids = set()
    else:
        sra_final_list = None
        id_list = None
        processed_ids = None

    # Broadcast data to all ranks
    sra_final_list = comm.bcast(sra_final_list, root=0)
    id_list = comm.bcast(id_list, root=0)
    processed_ids = comm.bcast(processed_ids, root=0)

    # Distribute SRA IDs
    chunk_size = len(sra_final_list) // size
    remainder = len(sra_final_list) % size
    start_index = rank * chunk_size + min(rank, remainder)
    end_index = start_index + chunk_size + (1 if rank < remainder else 0)
    sra_final_list_chunk = sra_final_list[start_index:end_index]

    print(f"Rank {rank} processing {len(sra_final_list_chunk)} SRA IDs")

    processed_this_round = []

    for sra_id in sra_final_list_chunk:
        elapsed_time = time.time() - start_time
        if elapsed_time >= stop_before_seconds:
            print(f"Rank {rank}: Near 48h limit. Saving progress and exiting safely.")
            break

        try:
            assembly(sra_id)
            processed_this_round.append(sra_id)

            # all ranks write progress
            with open(progress_file, "a") as f:
                f.write(sra_id + "\n")

        except Exception as e:
            print(f"Rank {rank}: Error processing {sra_id} — {e}")
            continue

    # Gather results from all ranks
    all_processed = comm.gather(processed_this_round, root=0)

    comm.Barrier()

    # Rank 0 writes all processed SRA IDs to progress file
    if rank == 0:
        with open(progress_file, "a") as f:
            for sublist in all_processed:
                for sra_id in sublist:
                    f.write(sra_id + "\n")
        print(f"Updated progress file: {progress_file}")

    comm.Barrier()

    # Rank 0 handles metadata and logging
    if rank == 0:
        elapsed_time = time.time() - start_time

        if elapsed_time < stop_before_seconds:
            sra_metadata = fetch_sra_metadata_for_Potential_positive_SRA_Runs(
                args.mail, id_list, args.wd + "/Positive_SRA/"
            )

            if sra_metadata:
                biosample_ids = [m['Biosample'] for m in sra_metadata.values()]
                biosample_metadata = fetch_sra_metadata_biosample(args.mail, biosample_ids)

    # Check if sra_metadata is None
                if biosample_metadata is None:
                   print("Error: Failed to fetch SRA metadata.")
                else:
        # Extract Biosample IDs from SRA metadata
                    biosample_ids = [metadata['Biosample'] for metadata in sra_metadata.values()]

        # Fetch Biosample metadata
                    biosample_metadata = fetch_sra_metadata_biosample(args.mail, biosample_ids)

        # Check if biosample_metadata is None
                    if biosample_metadata is None:
                        print("Error: Failed to fetch Biosample metadata.")
                    else:
            # Create combined table
                        export_metadata_to_excel(sra_metadata, biosample_metadata)

        # Performance logging
        end_time = time.time()
        total_runtime = end_time - start_time
        minutes, seconds = divmod(total_runtime, 60)
        memory_gb = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 3)

        file_exists = os.path.isfile(log_file)
        with open(log_file, mode='a', newline='') as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow(["Date", "Runtime_seconds", "Runtime_minutes", "Memory_GB"])
            writer.writerow([time.strftime("%Y-%m-%d %H:%M:%S"), total_runtime, total_runtime / 60, memory_gb])

        print(f"Runtime: {int(minutes)}m {seconds:.2f}s | Memory: {memory_gb:.2f} GB")
        print(f"Performance log saved to {log_file}")

    # Clean up safely (if enough time left)
    elapsed_time = time.time() - start_time
    if rank == 0 and elapsed_time < stop_before_seconds:
        delete_fastqfiles_folder()
        delete_assembly()

    sys.stdout.flush()
    sys.stderr.flush()
    comm.Barrier()
