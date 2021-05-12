# Call peaks

In this page, we will discuss usage of the main function in MACS3 --
`callpeak` . It can be invoked by `macs3 callpeak` . The goal of this
function is to find enrichment of a target-enriched DNA assay such as
ChIP-seq, ATAC-seq, Cut&Tag and etc. The `callpeak` function contains
steps from preprocess, callling peaks, to some downstream
anslysis. MACS3 also provides a series of subcommands that you can use
to build a pipeline to implement the function of `callpeak`. You can
find the information
[here](https://github.com/macs3-project/MACS/wiki/Advanced%3A-Call-peaks-using-MACS2-subcommands).

To see a full description of command-line options for `callpeak`, type
`macs3 callpeak -h`. Here we list the essential options.

## Essential Options

### `-t TFILE [TFILE ...]`/`--treatment TFILE [TFILE ...]`

This is the only REQUIRED parameter for MACS. The file can be in any
supported format, although we highly recommend BAM format -- see
detail in the `--format` option. If you have more than one alignment
file, you can specify them as `-t A B C`. MACS will pool up all these
files together. Regarding how to properly process data with
replicates, please see [here](./Replicates.md).

### `-c [CFILE ...]`/`--control [CFILE ...]`

The control, genomic input or mock IP data file. Please follow the
same direction as for `-t`/`--treatment`. Some experiment doesn't have
matching controls, for example, the ATAC-seq. In these cases, please
ignore this option. Some users also use this option to *compare* two
treatment files from two conditions. However, please keep in mind that
MACS processes treatment and control file differently so in order to
do a fair comparison, you'd better consider [`bdgcmp`](./bdgcmp.md)
and [`bdgdiff`](./bdgdiff.md) command.

### `-n NAME`/`--name NAME`

The name of the experiment. MACS will use this NAME to create output
files including `NAME_peaks.xls`, `NAME_negative_peaks.xls`,
`NAME_peaks.bed` , `NAME_summits.bed`, `NAME_model.r` and so on. So
please avoid any confliction between these filenames and your existing
files. Some users use a side effect of this option to save files in a
subdirectory. For example, `-n A/B` will make MACS to save output
files in a subdirectory named `A`. However, we suggest the users to
use `--outdir` option for this purpose, since the `NAME` will also to
be used in the peak names in the output files.

### `--outdir OUTDIR`

MACS3 will save all output files into a specified folder for this
option. A new folder will be created if necessary.

### `-f FORMAT`/`--format FORMAT`

Format of tag file can be `ELAND`, `BED`, `ELANDMULTI`, `ELANDEXPORT`,
`SAM`, `BAM`, `BOWTIE`, `BAMPE`, or `BEDPE`. Default is `AUTO` which
will allow MACS to decide the format automatically. `AUTO` is also
useful when you combine different formats of files. Note that MACS
can't detect `BAMPE` or `BEDPE` format with `AUTO`, and you have to
implicitly specify the format for `BAMPE` and `BEDPE`.

Nowadays, the most common formats are `BED` or `BAM` (including
`BEDPE` and `BAMPE`). Our recommendation is to convert your data to
`BED` or `BAM` first. As a recommendation, if you have a plain text
`SAM` file, please convert it into the binary counterpart -- `BAM`
file, since through the converting process we can make sure the binary
format file is valid. 

Also, MACS3 can detect and read gzipped file. For example, `.bed.gz`
file can be directly used with `--format BED` option, without being
uncompressed first.

Here are detailed explanation of the recommanded formats:

#### `BED`

The BED format can be found at [UCSC genome browser
website](http://genome.ucsc.edu/FAQ/FAQformat#format1).

The essential columns in BED format input are the 1st column
`chromosome name`, the 2nd `start position`, the 3rd `end position`,
and the 6th, `strand`.

Note that, for `BED` format, the 6th column of strand information is
required by MACS. And please pay attention that the coordinates in BED
format are zero-based and half-open. See more detail at
[UCSC site](http://genome.ucsc.edu/FAQ/FAQtracks#tracks1).

The `BEDPE` format is in fact a `BED` format with only three columns,
which define the location of the entire fragment. Please refer to the
corresponding section below.

#### `BAM`/`SAM`

If the format is `BAM`/`SAM`, please check the definition at
[samtools project site](https://samtools.github.io/hts-specs/).  **If
the `BAM` file is generated for paired-end data, MACS will only keep
the left mate(5' end) tag of properly paired reads.** However, when
format `BAMPE` is specified, MACS will use the real fragments inferred
from alignment results for reads pileup.

MACS uses the `bwflag` in BAM/SAM to filter the alignment
results. Check
[Picard site](https://broadinstitute.github.io/picard/explain-flags.html)
for a nice explanation on the `bwflag`, also the
[wikipedia page](https://en.wikipedia.org/wiki/SAM_(file_format)) . If
the `bwflag` matches `2820`, i.e., 'read unmapped' or 'not primary
alignment' or 'read fails QC' or 'supplementary', the alignment record
will be discarded. If `bwflag` indicates (`1`) that the read is from a
read pair, then: if the read is not in a proper pair (not `2`), it
will be discarded; if the mate (other read in the pair) is not mapped,
it will be discarded; if the this read is the second read ( or 3'
read) in the pair, it will also be discarded. Please note that the
read is the second read or a 3' read doesn't mean that the read is the
'rightmost' read according to mapping coordination.

#### `BEDPE` or `BAMPE`

A special mode will be triggered while the format is specified as
`BAMPE` or `BEDPE`. In this way, MACS3 will process the `BAM` or `BED`
files as paired-end data. Instead of building a bimodal distribution
of plus and minus strand reads to predict fragment size, MACS3 will
use actual insert sizes of pairs of reads to build fragment pileup.

The `BAMPE` format is just a `BAM` format containing paired-end
alignment information, such as those from `BWA` or `BOWTIE`.

The `BEDPE` format is a simplified and more flexible `BED` format,
which only contains the first three columns defining the chromosome
name, left and right position of the fragment from Paired-end
sequencing. Please note, this is NOT the same format used by
`BEDTOOLS`, and the `BEDTOOLS` version of `BEDPE` is actually not in a
standard `BED` format. **You can use MACS3 subcommand
[`randsample`](./randsample.md) to convert a `BAM` file containing
paired-end information to a `BEDPE` format file**:

```
macs3 randsample -i the_BAMPE_file.bam -f BAMPE -p 100 -o the_BEDPE_file.bed
```

The trick is to keep 100% of reads with [`randsample`](./randsample.md) .

### `-g GSIZE`/`--gsize GSIZE`

It's the mappable genome size or effective genome size which is
defined as the genome size which can be sequenced. Because of the
repetitive features on the chromosomes, the actual mappable genome
size will be smaller than the original size, about 90% or 70% of the
genome size. The default *hs* -- 2.7e9 is recommended for human
genome. Here are all precompiled parameters for effective genome size:

 * hs: 2.7e9
 * mm: 1.87e9
 * ce: 9e7
 * dm: 1.2e8

You can use the word `hs`, `mm`, `ce`, or `dm` to access the
precompiled numbers, or give an actual number of the mappable genome
size. 

Users may want to use k-mer tools to simulate mapping of Xbps long
reads to target genome, and to find the ideal effective genome
size. However, usually by taking away the simple repeats and Ns from
the total genome, one can get an approximate number of effective or
mappable genome size. A slight difference in the number won't cause a
big difference of peak calls, because this number is used to estimate
a genome-wide noise level which is usually the least significant one
compared with the *local biases* modeled by MACS.

### `-s TSIZE`/`--tsize TSIZE`

The size of sequencing tags/reads. If you don't specify it, MACS will
try to use the first 10 sequences from your input treatment file to
determine the tag size or the read length. Specifying it will override
the value automatically detected. Check the section on `--min-gap` on
a side effect of read length value.

### `-q QVALUE`/`--qvalue QVALUE`

The q-value (minimum FDR) cutoff to call significant regions. Default
is 0.05. Q-values are calculated from p-values using the
[Benjamini-Hochberg
procedure](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure).

### `-p PVALUE`/`--pvalue PVALUE`

The p-value cutoff. The p-value is calculated using Poisson CDF
function where the lambda is the local bias MACS3 estimated from
control or from a genome background, and the observed value is the
pileup signal from treatment. If `-p` is specified, MACS3 will use
p-value instead of q-value. In another word, `-p` will override `-q`
cutoff. In almost all cases (if there is no bug in the codes), given
the same value, `-p` provides a more relexed cutoff comparing to
`-q`. Regarding how to find "proper" cutoff, please check the
discussion on [cutoff analysis](./cutoff-analysis.md).

### `--min-length MINLEN` and `--max-gap MAXGAP`

These two options can be used to fine-tune the peak calling behavior
by specifying the minimum length of a called peak and the maximum
allowed a gap between two nearby regions to be merged. In other words,
a called peak has to be longer than `min-length`, and if the distance
between two nearby peaks is smaller than `max-gap` then they will be
merged as one. If they are not set, MACS3 will set the DEFAULT value
for `min-length` as the predicted fragment size `d`, and the DEFAULT
value for `max-gap` as the detected read length. Note, if you set a
`min-length` value smaller than the fragment size, it may have NO
effect on the result. For broad peak calling with `--broad` option
set, the DEFAULT `max-gap` for merging nearby stronger peaks will be
the same as narrow peak calling, and 4 times of the `max-gap` will be
used to merge nearby weaker (broad) peaks. You can also use
[`--cutoff-analysis`](./cutoff-analysis.md) option with the default
setting, and check the column `avelpeak` under different cutoff values
to decide a reasonable `min-length` value.

### `--nolambda`

With this flag on, MACS will use the whole genome background as local
lambda. This means MACS will not consider the local (1kbps, 10kbps)
bias at peak candidate regions. This option is recommended while
processing assay with no matching control, for example, ATAC-seq.

### `--slocal SMALLLOCAL` and `--llocal LARGELOCAL`

These two parameters control which two levels of regions will be
checked around the peak regions to calculate the maximum lambda as
local lambda. By default, MACS considers 1000bp for small local
region(`--slocal`), and 10000bps for large local region(`--llocal`)
which captures the bias from a long-range effect like an open
chromatin domain. You can tweak these according to your
project. Remember that if the region is set too small, a sharp spike
in the input data may kill a significant peak. In most cases, you
don't have to tweak these two parameters.

### `--nomodel`

While on, MACS will bypass building the read shifting model in order
to determine the 'fragment size' or 'd'. In PE mode, this option is
ignored since the actual insertion size of each read pair will be used
as the fragment size.

### `--extsize EXTSIZE`

While `--nomodel` is set, MACS uses this parameter to extend reads in
5'->3' direction to fix-sized fragments. For example, if the size of
the binding region for your transcription factor is 200 bp, and you
want to bypass the model building by MACS, this parameter can be set
as 200. This option is only valid when `--nomodel` is set or when MACS
fails to build model and `--fix-bimodal` is on. In most cases, if you
plan to control the behavior of read extension, please use `--nomodel
--extsize XXX`.

### `--shift SHIFT`

Note, this is NOT the legacy `--shiftsize` option which is replaced by
`--extsize`! You can set an arbitrary shift in bp here. Please Use
discretion while setting it other than the default value (0). When
`--nomodel` is set, MACS will use this value to move cutting ends (5')
then apply `--extsize` from 5' to 3' direction to extend them to
fragments. When this value is negative, ends will be moved toward
3'->5' direction, otherwise 5'->3' direction. Recommended to keep it
as default 0 for ChIP-Seq datasets, or -1 * half of *EXTSIZE* together
with `--extsize` option for detecting enriched cutting loci such as
certain DNAseI-Seq or ATAC-seq datasets. Note, you can't set values
other than 0 if the format is BAMPE or BEDPE for paired-end data. The
default is 0.

Here are some examples for combining `--shift` and `--extsize`:

1. To find enriched cutting sites such as some DNAse-Seq or ATAC-seq
datasets, all 5' ends of sequenced reads should be extended in both
directions to smooth the pileup signals. If the wanted smoothing
window is 200bps, then use `--nomodel --shift -100 --extsize 200`.

2. For certain nucleosome-seq data, we need to pile up the centers of
nucleosomes using a half-nucleosome size for wavelet analysis
(e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about
147bps, this option can be used: `--nomodel --shift 37 --extsize 73`.

### `--keep-dup KEEPDUPLICATES`

It controls the MACS behavior towards duplicate tags ( or pairs in PE
mode) at the exact same location -- the same coordination and the same
strand. The default `auto` option makes MACS to calculate the maximum
tags at the exact same location based on binomial distribution using
1e-5 as p-value cutoff; and the `all` option keeps every tag.  If an
integer is given, at most this number of tags will be kept at the same
location. The default is to keep one tag at the same
location. However, in many cases, users want to use their own
preprocessing tool to remove duplicates. If so, please make sure that
the duplicates have been actually *removed* from the BAM or BED file
instead of only simply marking them (e.g. `Picard` `MarkDuplicates`),
then use `--keep-dup all` for `callpeak`. Default: 1

### `--broad`

When this flag is on, MACS will enter the 'broad peak calling mode'
and try to composite broad regions in BED12 ( a gene-model-like format
) by putting nearby highly enriched regions into a broad region with
loose cutoff. The highly enriched regions will be controled by `-p` or
`-q` option, and the broad region is controlled by another cutoff
through `--broad-cutoff`. Please note that, the `max-gap` value for
merging nearby weaker/broad peaks is 4 times of `max-gap` for merging
nearby stronger peaks. The later one can be controlled by `--max-gap`
option, and by default it is the average fragment/insertion length in
the PE data. If your purpose is only to detect a broad region and you
have no interest in knowing how the highly enriched regions are
organized in the broad region, you can simply use a loose cutoff `-p`
or `-q` (of e.g. 0.1), bigger `--max-gap` and `--min-length` settings,
*without* using `--broad`. DEFAULT: False

### `--broad-cutoff`

Cutoff for the broad region. This option is not available unless
`--broad` is set. If `-p` is set, this is a p-value cutoff, otherwise,
it's a q-value cutoff.  DEFAULT: 0.1

### `--scale-to <large|small>`

MACS will decide which, treatment of control, has to be scaled in
order to balance the total number of reads. When set to `large`,
linearly scale the smaller dataset to the same depth as the larger
dataset. By default or being set as `small`, the larger dataset will
be scaled towards the smaller dataset. We recommend the default
behavior since to scale up small data would cause more false
positives. Please beware that the scaling will affect the resulting
'bedGraph' output of fragment pileup. For example, if the consequence
is that the treatment has to be scaled (either up or down-scaled), the
pileup value in the `NAME_treat_pileup.bdg` will not be integer
number. And the same is for the control pileup
`NAME_control_lambda.bdg` if the control has to be scaled.

### `-B`/`--bdg`

If this flag is on, MACS will store the fragment pileup, control
lambda in bedGraph files. The bedGraph files will be stored in the
current directory named `NAME_treat_pileup.bdg` for treatment data,
`NAME_control_lambda.bdg` for local lambda values from control.

### `--SPMR`

If this flag is on, MACS will store the 'Signal Per Million Reads' in
the bedGraph files mentioned for the `-B`/`--bdg` option. The values
in the 4th column of bedGraph file equal to the actual fragment pileup
divided by the number of million reads in the library. As for PE data,
the pileup will be divided by the number of million pairs. If you plan
to generate signals file that can be loaded into genome browser and
you plan to visually compare the pileup from different samples, this
`--SPMR` option is recommended.

### `--call-summits`

When this option is on, MACS will reanalyze the shape of signal
profile (p or q-score depending on the cutoff setting) to deconvolve
subpeaks within each peak called from the general procedure. It's
highly recommended to detect adjacent binding events. While used, the
output subpeaks of a big peak region will have the same peak
boundaries, and different scores and peak summit positions. Please
note that, if a peak has multiple summits, it will be written in
multiple lines in the output peak file with the same content in the
first three columns, which are used to describe the location of the
whole peak region. 

### `--buffer-size`

MACS uses a buffer size for incrementally increasing internal array
size to store reads alignment information for each chromosome or
contig. To increase the buffer size, MACS can run faster but will
waste more memory if certain chromosome/contig only has very few
reads. In most cases, the default value 100000 works fine. However, if
there are a large number of chromosomes/contigs in your alignment and
reads per chromosome/contigs are few, it's recommended to specify a
smaller buffer size in order to decrease memory usage (but it will
take longer time to read alignment files). Minimum memory requested
for reading an alignment file is about # of CHROMOSOME * BUFFER_SIZE *
8 Bytes. DEFAULT: 100000

## Output files

1. `NAME_peaks.xls` is a tabular file which contains information about
   called peaks. You can open it in excel and sort/filter using excel
   functions. Information include:
   
    - chromosome name
    - start position of peak
    - end position of peak
    - length of peak region
    - absolute peak summit position
    - pileup height at peak summit
    - -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then
      this value should be 10)
    - fold enrichment for this peak summit against random Poisson
      distribution with local lambda,
    - -log10(qvalue) at peak summit
   
   Coordinates in XLS is 1-based which is different from BED
   format. When `--broad` is enabled for broad peak calling, the
   pileup, p-value, q-value, and fold change in the XLS file will be
   the mean value across the entire peak region, since peak summit
   won't be called in broad peak calling mode.

2. `NAME_peaks.narrowPeak` is BED6+4 format file which contains the
   peak locations together with peak summit, p-value, and q-value. You
   can load it to the UCSC genome browser. Definition of some specific
   columns are:
   
   - 5th: integer score for display. It's calculated as
     `int(-10*log10pvalue)` or `int(-10*log10qvalue)` depending on
     whether `-p` (pvalue) or `-q` (qvalue) is used as score
     cutoff. Please note that currently this value might be out of the
     [0-1000] range defined in [UCSC ENCODE narrowPeak
     format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12). You
     can let the value saturated at 1000 (i.e. p/q-value = 10^-100) by
     using the following 1-liner awk: `awk -v OFS="\t"
     '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak`
   - 7th: fold-change at peak summit
   - 8th: -log10pvalue at peak summit
   - 9th: -log10qvalue at peak summit
   - 10th: relative summit position to peak start
   
   The file can be loaded directly to the UCSC genome browser. Remove
   the beginning track line if you want to analyze it by other tools.

3. `NAME_summits.bed` is in BED format, which contains the peak
   summits locations for every peak. The 5th column in this file is
   the same as what is in the `narrowPeak` file. If you want to find
   the motifs at the binding sites, this file is recommended. The file
   can be loaded directly to the UCSC genome browser. Remove the
   beginning track line if you want to analyze it by other tools.

4. `NAME_peaks.broadPeak` is in BED6+3 format which is similar to
   `narrowPeak` file, except for missing the 10th column for
   annotating peak summits. This file and the `gappedPeak` file will
   only be available when `--broad` is enabled. Since in the broad
   peak calling mode, the peak summit won't be called, the values in
   the 5th, and 7-9th columns are the mean value across all positions
   in the peak region. Refer to `narrowPeak` if you want to fix the
   value issue in the 5th column.

5. `NAME_peaks.gappedPeak` is in BED12+3 format which contains both
   the broad region and narrow peaks. The 5th column is the score for
   showing grey levels on the UCSC browser as in `narrowPeak`. The 7th
   is the start of the first narrow peak in the region, and the 8th
   column is the end. The 9th column should be RGB color key, however,
   we keep 0 here to use the default color, so change it if you
   want. The 10th column tells how many blocks including the starting
   1bp and ending 1bp of broad regions. The 11th column shows the
   length of each block and 12th for the start of each block. 13th:
   fold-change, 14th: *-log10pvalue*, 15th: *-log10qvalue*. The file can
   be loaded directly to the UCSC genome browser. Refer to
   `narrowPeak` if you want to fix the value issue in the 5th column.

6. `NAME_model.r` is an R script which you can use to produce a PDF
   image of the model based on your data. Load it to R by:

   `$ Rscript NAME_model.r`

   Then a pdf file `NAME_model.pdf` will be generated in your current
   directory. Note, R is required to draw this figure.

7. The `NAME_treat_pileup.bdg` and `NAME_control_lambda.bdg` files are
   in bedGraph format which can be imported to the UCSC genome browser
   or be converted into even smaller bigWig files. The
   `NAME_treat_pielup.bdg` contains the pileup signals (normalized
   according to `--scale-to` option) from ChIP/treatment sample. The
   `NAME_control_lambda.bdg` contains local biases estimated for each
   genomic location from the control sample, or from treatment sample
   when the control sample is absent. The subcommand `bdgcmp` can be
   used to compare these two files and make a bedGraph file of scores
   such as p-value, q-value, log-likelihood, and log fold changes.
 
