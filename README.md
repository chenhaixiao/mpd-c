MPD - Multiplex PCR Design
============================

by Thomas Wingo and David Cutler

## Description

This package assists in the automation of multiplex primer design. It is under 
development. This code can be used as stand-alone software or in conjunction 
with the [mpd perl package](http://github.com/wingolab-org/mpd-perl).

## Installation

```bash
make 
```

## Build genome index

- Prepare or download a file with a copy of the whole genome in fasta format.
- run `index_genome` either interactively or by redirecting from a file 
this:

````
d
1000
hg38.fasta
hg38
````

## Obtain dbSNP files

- TODO...

## Usage

- Prepare a bed file of unique targets.
- run `mpd` either interactively or by redirecting from a file


````
d
TODO...

````
