# ont-basecalling

## Basecalling of ONT reads

This is a script to perform ONT-basecalling on SUS-AMR's Linux computers. It takes as input the FAST5 files from an ONT run (in directories or subdirectories) and outputs the raw fastq and fastq that have been subsampled (with filtlong).

Please note: This script is based on Ryan Wick's script basecall.py from https://github.com/rrwick/MinION-desktop and adapted for the use at the SUS AMR lab.


## Table of Contents

[Requirements](#Requirements)  
[Example Usage](#Basic-usage)  
[Usage](#Usage)  


## Requirements
These need to be installed and in path for the entire pipeline to work. Other versions of these tools will possibly work too, but these are the ones I have tested.

* guppy_basecaller v6.4.2 (Available for customers at https://community.nanoporetech.com/)
* filtlong v0.2.1 (https://github.com/rrwick/Filtlong)
* fast_count (https://github.com/rrwick/MinION-desktop)

## Example usage
``` 
python ont-basecalling.py --input_dir . --basecalling_model r9.4.1_sup --barcode_kit native_1-96 --key_file barcode_sample_key.csv
```

## Usage:

```
usage: ont-basecalling.py [-h] [-v] -i INPUT_DIR -b
                          {r9.4.1_fast,r9.4.1_hac,r9.4.1_sup,r9.5,r10_fast,r10_hac,r10.4_sup}
                          -k
                          {none,native_1-12,native_13-24,native_1-24,native_1-24_r10,native_1-96,rapid_1-12}
                          [-o OUTDIR] [--key_file KEY_FILE]
                          [--filtlong {on,off}] [--fast_count {on,off}]
                          [--resume] [--cpu]
                          [--chunks_per_runner CHUNKS_PER_RUNNER]

Basecall, demultiplex and assemble reads from ONT sequencing

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input options (required):
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Input directory, which will be recursively searched
                        for fast5-files.
  -b {r9.4.1_fast,r9.4.1_hac,r9.4.1_sup,r9.5,r10_fast,r10_hac,r10.4_sup,r10.4.1_sup}, --basecalling_model {r9.4.1_fast,r9.4.1_hac,r9.4.1_sup,r9.5,r10_fast,r10_hac,r10.4_sup,r10.4.1_sup}
                        Indicate which basecalling mode to use. In most cases
                        you probably want to use a HAC option.
  -k {none,native_1-12,native_13-24,native_1-24,native_1-24_r104,native_1-96,rapid_1-12,native_1-24_r1041}, --barcode_kit {none,native_1-12,native_13-24,native_1-24,native_1-24_r104,native_1-96,rapid_1-12,native_1-24_r1041}
                        Indicate which barcode-kits were used, if any.

Output options:
  -o OUTDIR, --outdir OUTDIR
                        Output directory for all output files. Default is
                        current working directory.

Optional flags:
  --key_file KEY_FILE   Provide a csv file with barcode,sample_name (one per
                        line) to rename the files.
  --filtlong {on,off}   Subsample fastq-files with Filtlong. Default: on.
  --fast_count {on,off}
                        QC of the raw ONT fastq files (003_fastq). Default:
                        on.
  --resume              Use this flag if your first run was interrupted and
                        you want to resume. Default: off.
  --cpu                 If GPU is busy, use CPU with this flag. Will use 4
                        threads and 6 callers. Default: GPU.

Advanced options:
  --chunks_per_runner CHUNKS_PER_RUNNER
                        Advanced option. Change chunks per runner. Default =
                        300

```
