#!/usr/bin/env python3
"""# marit.hetland@outlook.com
# github marithetland
# October 2020
# Process fast5-files from ONT sequencing
# Basecalling, demultiplexing, assembly
# This is based on Ryan Wick (rrwick)'s
# MinION-Desktop repo on his github page.
"""

#import modules
import logging, time
import glob
import datetime
import itertools
import os, sys, re
import shutil
import csv
import pandas as pd
from argparse import ArgumentParser
import subprocess
from Bio import SeqIO
import pathlib
import collections
import dateutil.parser
import h5py
import random
import statistics
import tempfile
from subprocess import call
from subprocess import check_output, CalledProcessError, STDOUT

#Options for basecalling and barcoding ## From Ryan
BASECALLING = collections.OrderedDict([
    ('r9.4_fast', ['--config dna_r9.4.1_450bps_fast.cfg ']),
    ('r9.4_hac',  ['--config dna_r9.4.1_450bps_hac.cfg ']),
    ('r9.5',  ['--config dna_r9.5_450bps.cfg ']),
    ('r10_fast',  ['--config dna_r10_450bps_fast.cfg ']),
    ('r10_hac',   ['--config dna_r10_450bps_hac.cfg ']),
    ('r10.4_sup',   ['--config dna_r10.4_e8.1_sup.cfg ']),
])

BARCODING = collections.OrderedDict([
    ('native_1-12',  ['--barcode_kits "EXP-NBD104" ']),
    ('native_13-24', ['--barcode_kits "EXP-NBD114" ']),
    ('native_1-24',  ['--barcode_kits "EXP-NBD104 EXP-NBD114" ']),
    ('native_1-96',  ['--barcode_kits "EXP-NBD196" ']),
    ('rapid_1-12',   ['--barcode_kits "SQK-RBK004" ']),
    ('none',         ['--disable_trim_barcodes'])
])

#Defs
def parse_args():
    #Version
    parser = ArgumentParser(description='Basecall, demultiplex and assemble reads from ONT sequencing')
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + 'v.1.0.0')

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    output_args = parser.add_argument_group('Output options (required)')
    optional_args = parser.add_argument_group('Optional flags')
    advanced_args = parser.add_argument_group('Advanced options')

    #Input
    input_args.add_argument('-i', '--input_dir', type=pathlib.Path, required=True, help='Input directory, which will be recursively searched for fast5-files.')
    input_args.add_argument('-b', '--basecalling_model', type=str, required=True, choices=["r9.4_fast","r9.4_hac","r9.5","r10_fast","r10_hac","r10.4_sup"], help='Indicate which basecalling mode to use. In most cases you probably want to use a HAC option.')
    input_args.add_argument('-k', '--barcode_kit', type=str, required=True, choices=["none","native_1-12","native_13-24","native_1-24","rapid_1-12"], help='Indicate which barcode-kits were used, if any.')

    #Output - currently writes to same dir as input
    output_args.add_argument('-o', '--outdir', type=pathlib.Path, required=False, default='.', help='Output directory for all output files.')


    #Options
    optional_args.add_argument('-f', '--filtlong', type=str, choices=["on","off"], required=False, default="on", help='Subsample fastq-files with Filtlong? Default: on.')
    optional_args.add_argument('-a', '--assemble', action='store_true', required=False, help='Assemble the fastq-files from ONT only with unicycler. Default: off.')
    optional_args.add_argument('--resume', action='store_true', required=False, help='Use this flag if your first run was interrupted and you want to resume. Default: off.')
    optional_args.add_argument('--cpu', action='store_true', required=False, help='If GPU is busy, use CPU with this flag. Will use 4 threads and 6 callers. Default: GPU.')

    #Advanced options
    advanced_args.add_argument('--chunks_per_runner', type=str, required=False, help='Advanced option. Change chunks per runner. Default = 300')
    
    return parser.parse_args()


def run_command(command, **kwargs): #yekwahs
    command_str = ''.join(command)
    #logging.info('Running shell command: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"Error:": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"Error:": message})

def check_guppy_version():
    run_command(['guppy_basecaller --version > guppy_basecaller_version.txt'], shell=True)     ##Check for empty results, skip
    pass

def check_filtlong_version():
    run_command(['filtlong --version > filtlong_version.txt'], shell=True)     ##Check for empty results, skip
    pass

def check_python_version():
    try:
        assert sys.version_info >= (3, 5)
    except AssertionError:
        sys.exit('Error: Python 3.5 or greater is required')

#def merge_fastq():


def check_arguments(args):
    '''Check that the arguments are satisfied.'''
    barcode_choices = list(BARCODING.keys())
    args.barcode_kit = args.barcode_kit.lower()
    if args.barcode_kit not in barcode_choices:
        sys.exit('Error: valid --barcodes choices are: {}'.format(join_with_or(barcode_choices)))

    model_choices = list(BASECALLING.keys())
    args.basecalling_model = args.basecalling_model.lower()
    if args.basecalling_model not in model_choices:
        sys.exit('Error: valid --model choices are: {}'.format(join_with_or(model_choices)))

    if not args.input_dir.is_dir():
        sys.exit('Error: {} is not a directory'.format(args.in_dir))

    if args.outdir.is_file():
        sys.exit('Error: {} is a file (must be a directory)'.format(args.outdir))
    elif not args.outdir.is_file():
        directory=args.outdir
        if not os.path.exists(directory):
             os.makedirs(directory)
             print('Created output directory: '+ str(args.outdir))
        elif os.path.exists(args.outdir) and os.listdir(directory): 
            print("Warning: Specified output directory is not empty.")
        elif os.path.exists(args.outdir) and not os.listdir(directory):
            print("Specified output directory is empty.")
    
    return directory


def get_guppy_command(input_dir, save_dir, barcode_kits, basecalling_model, resume, cpu):
    guppy_command = ['guppy_basecaller',
                     '--input_path ', input_dir, 
                     '--recursive ',
                     '--save_path ', save_dir,
                     '--compress_fastq ']
    guppy_command += BASECALLING[basecalling_model]
    guppy_command += BARCODING[barcode_kits]
    if cpu:
        guppy_command += ['--num_callers 6'] #change to auto / add option to use CPU or GPU with options
    else:
        guppy_command += ['--device ', 'cuda:all:100%'] #change to auto / add option to use CPU or GPU with options

    if chunks:
        guppy_command += ['--chunks_per_runner ', chunks] #change to auto / add option to use CPU or GPU with options
    else:
        guppy_command += ['--chunks_per_runner ', '300'] #change to auto / add option to use CPU or GPU

    if resume:
         guppy_command += ['--resume']  #add option to resume
    
    print(guppy_command)
    return guppy_command

def listToString(s):  
    # initialize an empty string 
    str1 = " " 
    # return string   
    return (str1.join(s)) 

def main():
    # Set up log to stdout
    logfile= None
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running susamr-bascealling script v.1.0.0')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    #Check arguments, set variables and set up output directory

    logging.info("Hello! Just checking that all parameters and versions are in place before we start...")
    args = parse_args()
    check_python_version()
    check_guppy_version()
    check_filtlong_version()
    check_arguments(args)

    outdir=os.path.abspath(str(args.outdir))
    if outdir[-1] != '/':
        outdir = outdir + '/'
        print(outdir)


    raw_fast5s=os.path.abspath(str(args.input_dir))
    basecalled_fastq=(outdir+'002_basecalled_demultiplexed/')
    cat_fastq=(outdir+'003_fastq/')
    subsampled_fastq=(outdir+'004_subsampled_fastq/')
    #ont_only_assemblies=(outdir+'005_ont_only_assemblies/')
    barcode_kit=args.barcode_kit
    print("Specified barcode kit is: " + barcode_kit)
    basecaller_mode=args.basecalling_model
    print("Using basecaller mode: " + basecaller_mode)
    
    logging.info("We're good to go! Basecalling and demultiplexing with Guppy now")
    logging.info("PS. You can see which version of guppy was used in the file: guppy_basecaller_version.txt")

    ##Part 1: Run Guppy    
    #Get Guppy command, run guppy  (create def basecall_reads():)
    guppy_command=(get_guppy_command(raw_fast5s, basecalled_fastq,barcode_kit, basecaller_mode, args.resume, args.cpu, args.chunks_per_runner))
    run_command([listToString(guppy_command)], shell=True)
    
    logging.info("Guppy is done, now the fastq-files for each genome are being concatenated - bear with me")

    ##Part 2: Merge FASTQs for each genome
    if barcode_kit == 'none':
        run_command(['mkdir ',cat_fastq,' ; cat ',basecalled_fastq,'*fastq.gz >> ',cat_fastq,'reads.fastq.gz'], shell=True)
    else:
        run_command(['mkdir ',cat_fastq,' ; for dir in $(ls -d ',basecalled_fastq,'barcode[0-9][0-9] ',basecalled_fastq,'unclassified) ; do cd ${dir} ; base=$(basename $dir) ; cat *gz >> ',cat_fastq,'${base}.fastq.gz ; cd .. ; done' ], shell=True)


    ##Part 3: Run Filtlong
    if args.filtlong:
        logging.info("Less is more, so subsampling the reads now with FiltLong")
        if barcode_kit == 'none':
            run_command(['mkdir ',subsampled_fastq,' ; filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ',cat_fastq,'reads.fastq.gz| gzip > ',subsampled_fastq,'reads_subsampled.fastq.gz'], shell=True)
        else:
            run_command(['mkdir ',subsampled_fastq,' ; for f in $(ls ',cat_fastq,'*.fastq.gz | cut -d"." -f1) ; do base=$(basename ${f}) ;  filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ${f}.fastq.gz| gzip > ',subsampled_fastq,'${base}_subsampled.fastq.gz ; done'], shell=True)
        logging.info("The FASTQ files are now ready for further analysis :) You'll find them in: " + subsampled_fastq)
    else:
        logging.info("The FASTQ files are now ready for further analysis :) You'll find them in: " + cat_fastq)

    ## Part 4: Count number of reads
    logging.info("Counting number of reads in each raw FASTQ file and printing to: counts_raw_fastq_tmp.txt")
    run_command(['fast_count > counts_raw_fastq_tmp.txt ; sed -i "s/seq_count/Number_of_reads/g" counts_raw_fastq_tmp.txt  ; sed -i "s/total_length/Number_of_bases_in_reads/g" counts_raw_fastq_tmp.txt  ;  ',cat_fastq,'*fastq.gz >> counts_raw_fastq_tmp.txt ; cat counts_raw_fastq_tmp.txt | rev | cut -d"/" -f-1 | rev >> counts_raw_fastq.txt ; rm counts_raw_fastq_tmp.txt' ], shell=True)
    if args.filtlong:
        logging.info("Counting number of reads in each subsampled FASTQ file and printing to: counts_subsampled_fastq.txt")
        run_command(['fast_count > counts_subsampled_fastq_tmp.txt ;sed -i "s/seq_count/Number_of_reads/g" counts_subsampled_fastq_tmp.txt  ; sed -i "s/total_length/Number_of_bases_in_reads/g" counts_subsampled_fastq_tmp.txt   ; fast_count ',subsampled_fastq,'*fastq.gz  >> counts_subsampled_fastq_tmp.txt  ; cat counts_subsampled_fastq_tmp.txt | rev | cut -d"/" -f-1 | rev >> counts_subsampled_fastq.txt ; rm counts_subsampled_fastq_tmp.txt ' ], shell=True)


if __name__ == '__main__':
    main()
#EOF

#TODO: Add option to specify path to filtlong, fastq_count and guppy-basecaller
#TODO: Make in proper python code...
