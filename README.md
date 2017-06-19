MPD - Multiplex PCR Design
============================

by Thomas Wingo and David Cutler

## Citation

Please cite our [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1453-3) if you use MPD in your work. Thanks.

## Description

MPD is a program designed to automate creation of multiplex primer design written in C. The `mpd_lessGreedy` and `mpd_moreGreedy` binaries differ in which primer pool they choose to start with for pool creation. Either binary can be used as stand-alone or in conjunction with the [MPD perl package](http://github.com/wingolab-org/mpd-perl).

## Installation
- Clone the repository
- Make with `make all`
- Binaries will be compiled and saved to the `build` directory

## Required files
- You will need a hashed copy of the genome to run the primer software.
- Instructions below show how you can create one yourself. A prebuild hg38 genome with flat dbSnp files is available from [this repository](https://bitbucket.org/wingolab/mpd-dat/). It may be cloned like so, `git clone https://bitbucket.org/wingolab/mpd-dat.git`.

### Build Hashed Genome
- Download the genome of interest as a fasta file
- Use `bin/run_index.pl`, which creates a sh script to run `index_genome`

### Flat dbSnp
- These can be obtained from this [this repository](https://bitbucket.org/wingolab/mpd-dat/), which were prepared from dbSNP version 140.
- To create your own flat snp file set based on criteria of your own devising, each line should contain tab-delimited fields of the following:
```
name numberOfReporters chrom position MinorAlleleFrequency allele1/allele2
```
- The `numberOfReporters` field is no longer used but retained for backwards compatibility. 
- Prepare a `sdx` file that contains the number of chromsome files to include as the 1st line and then a list of the names of all chromosome files. On the command line you might try: `ls -1 *.line | wc -l > db_flat.sdx; ls -1 >> db_flat.sdx`. Note that the sdx should be in the same order that the chromsomes are in for the indexed genome. See the genome's sdx file (e.g., `cat hg38.d14.sdx`) to see the order.

##  Run mpd
- The easiest way of using MPD is to use the [Perl pacakge MPD](http://github.com/wingolab-org/mpd-perl), but either `mpd_lessGreedy` and `mpd_moreGreedy` binaries may be executed from the command line interactively.

