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

* guppy_basecaller (Available for customers at https://community.nanoporetech.com/)
* filtlong (https://github.com/rrwick/Filtlong)
* fast_count (https://github.com/rrwick/MinION-desktop)

## Example usage
``` 
python ont-basecalling.py -i . -b r9.4_hac -k native_1-24
```

## Usage:

```
usage: ont-basecalling.py [-h] [-v] -i INPUT_DIR -b
                                             {r9.4_fast,r9.4_hac,r9.5,r10_fast,r10_hac}
                                             -k
                                             {none,native_1-12,native_13-24,native_1-24,rapid_1-12}
                                             [-o OUTDIR] [-f {on,off}] [-a]
                                             [--resume] [--cpu]

Basecall, demultiplex and assemble reads from ONT sequencing

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input options (required):
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Input directory, which will be recursively searched
                        for fast5-files.
  -b {r9.4_fast,r9.4_hac,r10_fast,r10_hac}, --basecalling_model {r9.4_fast,r9.4_hac,r10_fast,r10_hac}
                        Indicate which basecalling mode to use. In most cases
                        you probably want to use a HAC option.
  -k {none,native_1-12,native_13-24,native_1-24,rapid_1-12}, --barcode_kit {none,native_1-12,native_13-24,native_1-24,rapid_1-12}
                        Indicate which barcode-kits were used, if any.

Output options (required):
  -o OUTDIR, --outdir OUTDIR
                        Output directory for all output files.

Optional flags:
  -f {on,off}, --filtlong {on,off}
                        Subsample fastq-files with Filtlong? Default: on.
  -a, --assemble        Assemble the fastq-files from ONT only with unicycler.
                        Default: off.
  --resume              Use this flag if your first run was interrupted and
                        you want to resume. Default: off.
  --cpu                 If GPU is busy, use CPU with this flag. Will use 4
                        threads and 6 callers. Default: GPU.
```


## ToDo-list
* Add option to provide barcode-filename key to rename files in script
