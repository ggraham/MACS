# Randomly downsample the data

Although in the MACS3 main function [`callpeak`](/.callpeak.md), by
default, when comparing treatment sample with control sample, the
pileup signals from the larger sample will be downscaled to match the
size of the smaller one, over the years, many users ask if they can do
'random sample' instead of 'scaling' and have more controls on this
behavior. Therefore, we made this `randsample` function as a
pre-processing tool to downsample an alignment file (as long as the
format is accepted by MACS3) to any number of reads or any percentage
of reads, then write the output in BED format. Since the output is BED
format, *this tool can also be used to convert any other formats of
alignment file into BED format or BEDPE format (if the input is BAMPE
format)*.

To see a full description of command-line options for `randsample`,
type `macs3 randsample -h`.

## Options

### `-i IFILE [IFILE ...]` or `--ifile IFILE [IFILE ...]`

Required option. Alignment file. If multiple files are given as '-t A
B C', then they will all be read and combined. 

### `-p PERCENTAGE` or `--percentage PERCENTAGE`

Percentage of tags you want to keep. Input 80.0 for 80%. This option
can't be used at the same time with -n/--num. The default is to keep
all, so it can be used to convert any formatof alignment file into BED
format.
  
### `-n NUMBER` or `--number NUMBER`

Number of reads (or fragments if `-f BAMPE` or `-f BEDPE`) you want to
keep. Input 8000000 or 8e+6 for 8 million. This option can't be used
at the same time with -p/--percent. Please note that the number of
reads/fragments in the output may not be exactly the number specified
in this option. MACS will first use this number to calculate the
percentage ( or probability to keep each read ), so it may lose some
precision.

### `--seed SEED`
Set the random seed while down sampling data. Must be a non-negative
integer in order to be effective. Set this if you want to have
reproducible results or plan to produce some pseudo-technical
replicates. 
  
###  `-o OUTPUTFILE` or  `--ofile OUTPUTFILE`

Output BED file name. If not specified, will write to standard
output. Note, if the input format is BAMPE or BEDPE, the output will
be in BEDPE format.  DEFAULT: stdout

### `--outdir OUTDIR`

If specified all output files will be written to that
directory. Default: the current working directory
  
### `-f` or `--format`

Format of tag file. It can be "BED" or "ELAND" or "ELANDMULTI" or
"ELANDEXPORT" or "SAM" or "BAM" or "BOWTIE" or "BAMPE" or "BEDPE". Or
you can set it as "AUTO". The default "AUTO" setting will let `macs3
randsample` decide which format the file is. Please check the
definition in [`callpeak` function](./callpeak.md) if you choose
ELAND/ELANDMULTI/ELANDEXPORT/SAM/BAM/BOWTIE or BAMPE/BEDPE. DEFAULT:
"AUTO".

### `--buffer-size BUFFER_SIZE`
Buffer size for incrementally increasing internal array size to store
reads alignment information. In most cases, you don't have to change
this parameter.  However, if there are large number of
chromosomes/contigs/scaffolds in your alignment, it's recommended to
specify a smaller buffer size in order to decrease memory usage (but
it will take longer time to read alignment files). Minimum memory
requested for reading an alignment file is about # of CHROMOSOME *
BUFFER_SIZE * 8 Bytes. DEFAULT: 100000

## Example of usage

## Downsample the data generated for a species with a small genome

Often, when we are analyzing the ChIP-seq or ATAC-seq data from a
species with a small genome, such as microbes, we will encounter a
problem that, since the sequencing depth (usually a full lane from
illumina machine) is too big, MACS will behave strangely. Either that
the background signal is too high since almost entire the genome can
be covered by a supposely target-enriched DNA assay; or the p-value is
overestimated since the infamous large number problem of p-value
calculation. I suggest that the large sample should be downsampled
into multiple smaller technical replicates. As for a normal TF
ChIP-seq on human genome, a 1x sequencing depth (= number of reads *
read length / mappable genome size) should be more than enough --
about 60million 50bp reads. So if you have a ChIP-seq on E. coli of
which the genome is about 4.6million bps, even 1 million reads are
already enough. So you have 100million reads from a single lane of
sequencing, you can do this to generate three tech replicates from
your alignment file:

```
$ macs3 randsample -i original.bam -n 1000000 -f BAM --seed 100 -o tech_rep1.bed
$ macs3 randsample -i original.bam -n 1000000 -f BAM --seed 200 -o tech_rep2.bed
$ macs3 randsample -i original.bam -n 1000000 -f BAM --seed 300 -o tech_rep3.bed
```

## Convert data into BED or BEDPE format

Although less informative than BAM format, the file in BED format can
be stored using much smaller space. This `randsample` command can be
used to *convert* any format that can be recognized by MACS3 into BED
format. It is more convinient if you are dealing with pair-end data
since with `-f BAMPE`, MACS3 will generate the `BEDPE` format output
instead which only contains three columns -- chromosome, leftmost
position of mapped pair, and rightmost position of mapped pair. And
this `BEDPE` can be further `gzip`ped into much smaller file, and can
be used in other commands of MACS3 directly. Please note that our
`BEDPE` format only contains three columns (or actually it's a
3-column BED format) and it's not the same as the BEDPE format used by
`bedtools`. To do the conversion:

```
$ macs3 randsample -i original.bam -f BAM -o converted.bed
$ macs3 randsample -i original_pe.bam -f BAMPE -o converted_pe.bedpe
```
