#!/usr/bin/env python

from ast import arg
from genericpath import exists
from operator import is_
from pathlib import Path
import sys
import argparse
import subprocess
import readline
import os
import glob
from xmlrpc.client import boolean 
import yaml 
import re
from Bio import SeqIO
import gzip

# Define command-line arguments
parser = argparse.ArgumentParser("\ncreate configfile and run all:\n\tpython porA_sub.py,\nrun with current configfile:\n\tpython porA_sub.py --config <configfile>")
parser.add_argument('--config', type=str, help='Run with previous configfile, enter path to configfile')
args = parser.parse_args()

def is_number(nr):
    try:
        int(nr)
        return True
    except ValueError:
        return False

# Utility function to replace '/dpssi' with '/users'
def ensure_users_path(path):
        return path.replace('/dpssi', '/users')

while True:
    cores = input("Enter how many cores you want to use for the snakemake run in total: ")
    if is_number(cores):
        break
    else:
        print("Input is not a number, please input number of cores for the run: ")

if not args.config:
    ##get path tab complete with input from prompt
    def complete(text, state):
            """
            Tab autocomplete for prompts.

            Args:
                text (str): current user input
                state (int): index from 0 to n until the function returns a non-string value

            Returns:
                str: matching file at the current state index
            """
            matches = glob.glob(text+'*')
            return matches[state] + os.sep if os.path.isdir(matches[state]) else matches[state]

    readline.set_completer_delims('')
    readline.parse_and_bind("tab: complete")
    readline.set_completer(complete)

    # Prompt the user for input and check path exists 

    while True:
        results_folder = input("Enter the path to the results folder (. for current directory): ")
        if results_folder and os.path.isdir(results_folder):
            res= ensure_users_path(os.path.abspath(results_folder))
            break
        elif results_folder: 
            res= ensure_users_path(os.path.abspath(results_folder))
            os.makedirs(res)
            break
        else:
            print("Path cannot be empty. Please enter the path to the results folder (. for current directory): ")

    while True:
        sample_folder = input("Enter the path to where the links to the samples (reads) are (. for current directory): ")
        if os.path.isdir(sample_folder):
            sample= ensure_users_path(os.path.abspath(sample_folder))
            break
        else: 
            print(f"The path '{sample_folder}' is either not correct or not a directory, please enter the path to the folder where the sample read files (fastq.gz) are: ")
   
    runid = input("Enter full name of runid (name of output folders): ")
    
    while True:
        partition = input("Enter which partition to use (priority in the queue, standard, surveillance, outbreak): ")
        if partition in {"standard", "surveillance", "outbreak"}:
            break
        else: 
            print("Not a accepted partition please enter either, standard, surveillance or outbreak: ")

    while True:
        Dotrunc = input("Do you want to add a truncation length for the filtering? [yes/y or no/n]: ")
        if Dotrunc.lower() in {"yes", "y"}:
            truncs = input("Do you want to change the truncation length for the filtering? [yes/y or no/n] (Default: R1 280bp, R2 260bp, only change if quality plots and tables includes/removes too many reads): ")
            if truncs.lower() in {"yes", "y"}:
                trunc1 = input("Input truncation length for read1 in bp: ")
                trunc2 = input("Input truncation length for read2 in bp: ")
                if is_number(trunc1) and is_number(trunc2):
                    break
                else:
                    print("Input is not a number, please input a number for both read 1 (fwd) and read 2 (rev)")
            elif (truncs.lower() in {"no", "n"}):
                trunc1=280
                trunc2=260
                break
            else:
                print("Please only enter either 'yes' or 'y' if you want to change the default truncation lengths of R1 280bp and R2 260bp otherwise enter either 'n' or 'no'")
        elif (Dotrunc.lower() in {"no", "n"}):
            trunc1=0
            trunc2=0
            break
        else:
            print("Please only enter either 'yes' or 'y' if you want to add truncation lengths, otherwise enter either 'n' or 'no'")


    print("Creating configfile, please wait...")

    ##get sample name and paths 
    def extract_info(files_list, sampleFolder):
        extracted_info = []
        for file in files_list:
            if re.search(r'_S[0-9].*_R.*', file):
                sample_name = re.sub(r'_S[0-9].*_R.*', '', file)
                reads=f'"{sample_name}": {sampleFolder}/{file}'
                extracted_info.append(reads)
            elif re.search(r'_R[12]', file):
                sample_name = re.sub(r'_R[12].*', '', file)
                reads=f'"{sample_name}": {sampleFolder}/{file}'
                extracted_info.append(reads)
            else:
                print("input files does not appear to be read files, not R1 or R2 present in filename, please control and possibly rename files")
                sys.exit(1)
        return extracted_info
 

    # Get all files in the directory
    all_files = os.listdir(sample)

    # Filter R1 files and retrieve relevant information
    R1_files = [file for file in all_files if re.search(r'_R1.*\.gz', file)]
    formatted_output_R1 = extract_info(R1_files, sample)

    # Filter R2 files and retrieve relevant information
    R2_files = [file for file in all_files if re.search(r'_R2.*\.gz', file)]
    formatted_output_R2 = extract_info(R2_files, sample)

    with open(f'{res}/Resources/configs/config_{runid}.yml', 'w') as f:
        f.write(f'folder: {res}\n')
        f.write(f'partition: {partition}\n')
        f.write(f'runid: {runid}\n')
        f.write(f'trunc1: {trunc1}\n')
        f.write(f'trunc2: {trunc2}\n')
        #f.write(f'min_reads: {nrReads}\n')
        f.write(f'samplesR1:\n')
        for line in formatted_output_R1:
            f.write(f'  {line}\n')
        f.write("samplesR2:\n")
        for line in formatted_output_R2:
            f.write(f'  {line}\n')

 
    #make snakemake dryrun and actual run
    print("Configfile created, now running dryrun of snakemake, please wait...")
    dryrun = f"""snakemake --snakefile snakefile_PORA_dada2.smk --configfile {res}/Resources/configs/config_{runid}.yml --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 -np"""
    dryrun_sub = subprocess.run(dryrun, shell=True, executable="/bin/bash")

    snakemake_errors = input("Did the command run without any errors (red)? [yes/y or no/n] Only green and yellow lines should be there: ")
    if snakemake_errors.lower() in {"yes", "y"}:
        snakemake_jobs = input("Does the number of jobs listed concur with the number of samples? [yes/y or no/n]: ")
        if snakemake_jobs.lower() in {"yes", "y"}:
            runs=f"""snakemake --snakefile snakefile_PORA_dada2.smk --configfile {res}/Resources/configs/config_{runid}.yml --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 --cores {cores} --latency-wait 60"""
            runs_sub = subprocess.run(runs, shell=True, executable="/bin/bash")
        else:
            print("check errors in configfile or snakefile and try again")
            sys.exit(1)
    else:
        print("check errors in configfile or snakefile and try again")
        sys.exit(1)

    print("Checking number of reads in files with primers removed, will keep files with more than 100 reads, please wait...")

    def test_reads(rmprimer_file):
        count=0
        with gzip.open(rmprimer_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                count += 1
        if count > int(post_nrReads): 
            return True
        else:
            return False


    with open(f'{res}/Resources/configs/config_{runid}.yml', 'r') as file_conf:
        config = yaml.safe_load(file_conf)

    #post_nrReads=nrReads/2
    #post_nrReads=1500 ##1K threshold
    post_nrReads=100 ##test miseq


    conf_keys=list(config["samplesR1"].keys())

    for key in conf_keys:
        if test_reads(f'{res}/Resources/data/primers_removed/{runid}/{key}_rmprimer_R1.fastq.gz') == False:
            del config["samplesR1"][key]
            del config["samplesR2"][key]

    out_file=f'{res}/Resources/configs/config2_{runid}.yml'
    with open(out_file, 'w') as file:
            yaml.dump(config, file, default_flow_style=False)


    print("Now running dryrun of snakemake for dada2, please wait...")
    dryrun = f"""snakemake --snakefile snakefile_PORA_dada2_part2.smk --configfile {res}/Resources/configs/config2_{runid}.yml --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 -np"""
    dryrun_sub = subprocess.run(dryrun, shell=True, executable="/bin/bash")

    snakemake_errors = input("Did the command run without any errors (red)? [yes/y or no/n] Only green and yellow lines should be there: ")
    if snakemake_errors.lower() in {"yes", "y"}:
        runs=f"""snakemake --snakefile snakefile_PORA_dada2_part2.smk --configfile {res}/Resources/configs/config2_{runid}.yml --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 --cores {cores} --latency-wait 60"""
        runs_sub = subprocess.run(runs, shell=True, executable="/bin/bash")
    else:
        print("check errors in configfile or snakefile and try again")
        sys.exit(1)


elif args.config:
    print("Now running dryrun of snakemake, please wait...")
    dryrun = f"""snakemake --snakefile snakefile_PORA_dada2.smk --configfile {args.config} --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 -np"""
    dryrun_sub = subprocess.run(dryrun, shell=True, executable="/bin/bash")

    snakemake_errors = input("Did the command run without any errors (red)? (yes/y or no/n) Only green and yellow lines should be there: ")
    if snakemake_errors.lower() in {"yes", "y"}:
        snakemake_jobs = input("Does the number of jobs listed concur with the number of samples? [yes/y or no/n]: ")
        if snakemake_jobs.lower() in {"yes", "y"}:
            runs=f"""snakemake --snakefile snakefile_PORA_dada2.smk --configfile {args.config} --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 --cores {cores} --latency-wait 60"""
            runs_sub = subprocess.run(runs, shell=True, executable="/bin/bash")
        else:
            print("check errors in configfile or snakefile and try again")
            sys.exit(1)
    else:
        print("check errors in configfile or snakefile and try again")
        sys.exit(1)

    
    with open(args.config, 'r') as file_conf:
        config = yaml.safe_load(file_conf)
    
    conf_folder=config["folder"]
    conf_runid=config["runid"]

    out_file=os.path.abspath(f'{conf_folder}/Resources/configs/config2_{conf_runid}.yml')

    if os.path.isfile(out_file):
        print("Now running dryrun of snakemake for dada2, please wait...")
        dryrun = f"""snakemake --snakefile snakefile_PORA_dada2_part2.smk --configfile {out_file} --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 -np"""
        dryrun_sub = subprocess.run(dryrun, shell=True, executable="/bin/bash")

        snakemake_errors = input("Did the command run without any errors (red)? (yes/y or no/n) Only green and yellow lines should be there: ")
        if snakemake_errors.lower() in {"yes", "y"}:
            runs=f"""snakemake --snakefile snakefile_PORA_dada2_part2.smk --configfile {out_file} --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 --cores {cores} --latency-wait 60"""
            runs_sub = subprocess.run(runs, shell=True, executable="/bin/bash")
        else:
            print("check errors in configfile or snakefile and try again")
            sys.exit(1)

    else:
        print("Checking number of reads in files with primers removed, will keep files with than 1500 reads, please wait...")
        #conf_nrReads=config["min_reads"]
        #post_nrReads=conf_nrReads/2
        post_nrReads=1500

        def test_reads(rmprimer_file):
            count=0
            with gzip.open(rmprimer_file, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    count += 1
            if count > int(post_nrReads): 
                return True
            else:
                return False

        sample_primersrm= os.path.abspath(f'{conf_folder}/data/primers_removed/{conf_runid}')
        conf_path=os.path.abspath(args.config)

        conf_keys=list(config["samplesR1"].keys())

        for key in conf_keys:
            if test_reads(f'{sample_primersrm}/{key}_rmprimer_R1.fastq.gz') == False:
                del config["samplesR1"][key]
                del config["samplesR2"][key]

        with open(out_file, 'w') as file:
                yaml.dump(config, file, default_flow_style=False)

        print("Now running dryrun of snakemake for dada2, please wait...")
        dryrun = f"""snakemake --snakefile snakefile_PORA_dada2_part2.smk --configfile {out_file} --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 -np"""
        dryrun_sub = subprocess.run(dryrun, shell=True, executable="/bin/bash")

        snakemake_errors = input("Did the command run without any errors (red)? (yes/y or no/n) Only green and yellow lines should be there: ")
        if snakemake_errors.lower() in {"yes", "y"}:
            runs=f"""snakemake --snakefile snakefile_PORA_dada2_part2.smk --configfile {out_file} --use-conda --conda-base-path /users/data/Tools/Conda/Miniconda3-py311_23.5.2-0-Linux-x86_64 --cores {cores} --latency-wait 60"""
            runs_sub = subprocess.run(runs, shell=True, executable="/bin/bash")
        else:
            print("check errors in configfile or snakefile and try again")
            sys.exit(1)

else:
    print("Run either without any options or with config")




