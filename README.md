# Demutliplexer
Demultiplexer for Primers or Barcodes

This pipeline is associated with the research paper: [Analysis of clinically relevant large tandem repeats using nanopore sequencing](https://www.nature.com/articles/s41598-025-30441-3)

## Introduction

This Script can be fed with FASTQ formatted sequencing 
information to detect and sort reads based on one or more barcodes
or primers. This script is universal and works with all sequencing
technologies.

## Table of Contents
- [Installation](#Installation)
- [Options](#options)
- [Usage](#usage)
- [Output](#output)
- [License](#license)


## Installation

This script require python >= 3.8

The necessary packages:

```
pip install regex biopython
```

## Options

```
#################################
Short	Long	Type	Default	Description
-i	--input	string	Required	Path to the input FASTQ file
-b	--barcodes	string	Required	Path to the barcode file
-o	--output	string	Required	Path to the output folder
-m	--mismatches	int	2	Maximum number of allowed mismatches
-w	--window	int	60	Search window size (bp from end of read to search)
-r	--require_both	flag	False	Require a valid barcode at both ends of the sequence
-p	--processes	int	4	Number of parallel processes to use
-v	--verbose	flag	False	Enable verbose output for detailed logging
-s	--split	flag	False	Allows barcodes to be split if two barcode sequences are given
#################################

```


## Usage

```
python demultiplex.py -i <input_fastq> -b <barcodes_file> -o <output_folder> [options]
```

If --split is used --require_both is not needed to find barcodes at both ends and --require_both will be ignored.
If --require_both is used and not --split only the first barcode will be read from the sample sheet.

--mismatches uses regex re logic to allow this max number of errors, these can be subsititions, insertions or deletions.

There is no limit to the size of the sequences used or the number of mismatches allowed, however increasing both will lead
to much longer script running times. For normal use 2 or fewer mismatches should be used, and if primer/barcode sequences
are lengthened the search --window should be increased accordingly.

--window do not forget to include linker and adapter sequence lengths that remain in the fastq file if trimming has not been
perfomed.

--processes this number should not exceed the number of cores available for the process as processing times will not be shortened.

--verbose is only at a preliminary level of development.

## Output

**Files in \*output Folder (one folder per sample)**
| **_File name_**                        | **_Description_**                                                                                                                                                                                                                                                                                                        |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| \*.fastq           | Reads split by barcode/primers

* Name given to barcode or primer found in sample sheet column one

## Sample Sheet

Barcodes/Primers and Names must be kept at the fifth line of the csv or they will be incorrectly read.
If --split barcodes is used as an option a second primer/barcode sequence must be provided
If --require_both is used as an option, only the first barcode will be read, but only if --split is not used.

Specificity for barcodes and there primer sequences can be increased by adding unique sequences to the primer/barcode.
  -Use this if you would like to seperate barcodes based on primer direction.



## License

MIT License

Copyright (c) 2025 MedGenMedUniWien

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
