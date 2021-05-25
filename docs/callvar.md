# Call variants

This function in MACS3 can be used to **directly** call genetics
variants and infer allele specific events from the BAM files inside
the ChIP-seq or ATAC-seq peaks, without the need of deep whole genome
sequencing. It can be invoked by `macs3 callvar` . The method was
previously implemented in a separate project named
[`SAPPER`](https://github.com/taoliu/SAPPER). The MACS3 port of SAPPER
has been modified to directly work on the BAM files people usually use
for peak calling. However, since `callvar` will randomly access chunks
of alignments within given peaks, the BAM files as input to`callvar`
have to be sorted by coordinates and indexed through `samtools`. In
another word, the `.bai` files are required together with the `.bam`
files.

## Algorithm

Targeted DNA assays such as ChIP-Seq or ATAC-seq provide high coverage
at functional regulatory elements where proteins interact with DNA, we
built this `callvar` algorithm based on the hypothesis that the DNA
sequences at the targeted regulatory elements can be piled up or
assembled ( into ‘unitigs’ ) for predicting SNV and short INDELs and
infering allele specific events.

TBA...

## Workflow of `callvar`

<img src="./callvar_workflow.jpg" width="800" />

## Commandline options

If you type this command with `-h`, you will see a full description of command-line options.

``` $macs3 callvar -h usage: macs3 callvar [-h] -b PEAKBED -t TFILE
[-c CFILE] [--outdir OUTDIR] -o OFILE [--verbose VERBOSE]
[-g GQCUTOFFHETERO] [-G GQCUTOFFHOMO] [-Q Q] [-D MAXDUPLICATE]
[-F FERMI] [--fermi-overlap FERMIMINOVERLAP]
[--top2alleles-mratio TOP2ALLELESMINRATIO]
[--altallele-count ALTALLELEMINCOUNT] [--max-ar MAXAR] [-m NP] ```

```

### Common usage

In general, user has to provide 1 BAM file from targeted DNA
sequencing assay (ChIP-seq or ATAC-seq) `TREAT_sorted.bam`, and 1 BED
file for regions of interest (i.e. peak regions)
`peaks.bed`. Optionally, user can also provide a BAM file from a
matching unbiased control sample, i.e. genomic input, or igg control
`CTRL_sorted.bam`. The idea is that the BAM file from targeted DNA
assay may have allele specific biases, whereas the control sample
doesn't.  The typical command is

```
$ macs3 callvar -b peaks.bed -t TREAT_sorted.bam -c CTRL_sorted.bam -o peaks.vcf
```

The result of variants in peaks, and allele specific events can be found in the `peaks.vcf` file.

### Prepare Input Files

The input files in BAM format should be sorted and indexed. You can use `samtools`:

```
$ samtools sort -o sorted.bam original.bam
$ samtools index -b sorted.bam
```

### Essential options

#### Region of interest in BED format `-b`

You can specify the regions you plan to call variants from. Usually,
it should be the peak file (narrowPeak or any BED format) from
MACS. REQUIRED.

#### Treatment file in sorted BAM `-t`

ChIP-seq/ATAC-seq treatment file in BAM format, sorted by
coordinates. Make sure the `.bai` file is avaiable in the same
directory. REQUIRED.

#### Optional `-c CFILE`

Optional control file in BAM format, sorted by coordinates. Make sure
the `.bai` file is avaiable in the same directory. When the control
file is provided, MACS3 will consider that the allelic bias in control
should be close to 1 (i.e. unbiased) whereas the allelic enrichment in
treatment can be biased. If available, a low-depth (or even
high-depth) whole genome sequencing libary can be used as control (or
genomic input).

#### Output directory and output filename `--outdir` and required `-o`

You can let MACS3 write all output files in a specific directory by
using `--outdir`. You can also specify the output file name through
`-o`. The `-o` option is REQUIRED.

#### Cutoff for calling heterozygous allele `-g`

You can also use `--gq-hetero`. This option specifies the Genotype
Quality (GQ) score `(-10log10((L00+L11)/(L01+L00+L11)))` cutoff for
Heterozygous allele type. Default:0, or there is no cutoff on GQ. You
can set a high value to keep more confident Heterozygous Alleles,
however, you can also keep the default and filter the result later.

#### Cutoff for calling homozygous allele (different with reference allele) `-G`

You can also use `--gq-homo`. This option specifies the Genotype
Quality (GQ) score (-10log10(L01/(L01+L00+L11))) cutoff for Homozygous
allele (not the same as reference) type. Default:0, or ther is no
cutoff on GQ. You can set a high value to keep more confident
Homozygous Alleles, however, you can also keep the default and filter
the result later.

#### Controling when to use assembler `fermi-lite` `-F`



#### The option for `fermi-lite` assembler `--fermi-overlap`

#### Quality control threshold for calling variants `--top2alleles-mratio` `--altallele-count` and `--max-ar`

### Interpret Output file
