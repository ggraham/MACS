# Call variants

This function in MACS3 can be used to **directly** call genetics variants and infer allele specific events from the BAM files inside the ChIP-seq or ATAC-seq peaks, without the need for deep whole genome sequencing. It can be invoked by `macs3 callvar` . The method was previously implemented in a separate project named [`SAPPER`](https://github.com/taoliu/SAPPER). The MACS3 port of SAPPER has been modified to directly work on the BAM files people usually use for peak calling. However, since `callvar` will randomly access chunks of alignments within given peaks, the BAM files as input to`callvar` have to be sorted by coordinates and indexed through `samtools`. In another word, the `.bai` files are required.

## Algorithm

Targeted DNA assays such as ChIP-Seq or ATAC-seq provide high coverage at functional regulatory elements where proteins interact with DNA, we built this `callvar` algorithm based on the hypothesis that the DNA sequences at the targeted regulatory elements can be piled up or assembled ( into ‘unitigs’ ) for predicting SNV and short INDELs and infering allele specific events.

TBA...

## Workflow of `callvar`

<img src="./callvar_workflow.jpg" width="800" />

## Commandline options

If you type this command with `-h`, you will see a full description of command-line options.

```
$macs3 callvar -h
usage: macs3 callvar [-h] -b PEAKBED -t TFILE [-c CFILE] [--outdir OUTDIR] -o OFILE [--verbose VERBOSE] [-g GQCUTOFFHETERO] [-G GQCUTOFFHOMO] [-Q Q]
                     [-D MAXDUPLICATE] [-F FERMI] [--fermi-overlap FERMIMINOVERLAP] [--top2alleles-mratio TOP2ALLELESMINRATIO]
                     [--altallele-count ALTALLELEMINCOUNT] [--max-ar MAXAR] [-m NP]

optional arguments:
  -h, --help            show this help message and exit

Input files arguments:
  -b PEAKBED, --peak PEAKBED
                        Peak regions in BED format, sorted by coordinates. REQUIRED.
  -t TFILE, --treatment TFILE
                        ChIP-seq/ATAC-seq treatment file in BAM format, sorted by coordinates. Make sure the .bai file is avaiable in the same directory.
                        REQUIRED.
  -c CFILE, --control CFILE
                        Optional control file in BAM format, sorted by coordinates. Make sure the .bai file is avaiable in the same directory.

Output arguments:
  --outdir OUTDIR       If specified all output files will be written to that directory. Default: the current working directory
  -o OFILE, --ofile OFILE
                        Output VCF file name.
  --verbose VERBOSE     Set verbose level of runtime message. 0: only show critical message, 1: show additional warning message, 2: show process information, 3:
                        show debug messages. DEFAULT:2

Variant calling arguments:
  -g GQCUTOFFHETERO, --gq-hetero GQCUTOFFHETERO
                        Genotype Quality score (-10log10((L00+L11)/(L01+L00+L11))) cutoff for Heterozygous allele type. Default:0, or there is no cutoff on GQ.
  -G GQCUTOFFHOMO, --gq-homo GQCUTOFFHOMO
                        Genotype Quality score (-10log10((L00+L01)/(L01+L00+L11))) cutoff for Homozygous allele (not the same as reference) type. Default:0, or
                        ther is no cutoff on GQ.
  -Q Q                  Only consider bases with quality score greater than this value. Default: 20, which means Q20 or 0.01 error rate.
  -D MAXDUPLICATE       Maximum duplicated reads allowed per mapping position, mapping strand and the same CIGAR code. Default: 1. When sequencing depth is
                        high, to set a higher value might help evaluate the correct allele ratio.
  -F FERMI, --fermi FERMI
                        Option to control when to apply local assembly through Fermi. By default (set as 'auto'), while SAPPER detects any INDEL variant in a
                        peak region, it will utilize Fermi to recover the actual DNA sequences to refine the read alignments. If set as 'on', Fermi will be
                        always invoked. It can increase specificity however sensivity and speed will be significantly lower. If set as 'off', Fermi won't be
                        invoked at all. If so, speed and sensitivity can be higher but specificity will be significantly lower. Default: auto
  --fermi-overlap FERMIMINOVERLAP
                        The minimal overlap for fermi to initially assemble two reads. Must be between 1 and read length. A longer fermiMinOverlap is needed
                        while read length is small (e.g. 30 for 36bp read, but 33 for 100bp read may work). Default:30
  --top2alleles-mratio TOP2ALLELESMINRATIO
                        The reads for the top 2 most frequent alleles (e.g. a ref allele and an alternative allele) at a loci shouldn't be too few comparing to
                        total reads mapped. The minimum ratio is set by this optoin. Must be a float between 0.5 and 1. Default:0.8 which means at least 80% of
                        reads contain the top 2 alleles.
  --altallele-count ALTALLELEMINCOUNT
                        The count of the alternative (non-reference) allele at a loci shouldn't be too few. By default, we require at least two reads support
                        the alternative allele. Default:2
  --max-ar MAXAR        The maximum Allele-Ratio allowed while calculating likelihood for allele-specific binding. If we allow higher maxAR, we may mistakenly
                        assign some homozygous loci as heterozygous. Default:0.95

Misc arguments:
  -m NP, --multiple-processing NP
                        CPU used for mutliple processing. Please note that, assigning more CPUs does not guarantee the process being faster. Creating too many
                        parrallel processes need memory operations and may negate benefit from multi processing. Default: 1
```
	

### Common usage

In general, user has to provide 1 BAM file from targeted DNA sequencing assay (ChIP-seq or ATAC-seq) `TREAT_sorted.bam`, and 1 BED file for regions of interest (i.e. peak regions) `peaks.bed`. Optionally, user can also provide a BAM file from a matching unbiased control sample, i.e. genomic input, or igg control `CTRL_sorted.bam`. The idea is that the BAM file from targeted DNA assay may have allele specific biases, whereas the control sample doesn't.  The typical command is

```
$ macs3 callvar -b peaks.bed -t TREAT_sorted.bam -c CTRL_sorted.bam -o peaks.vcf
```

The result of variants in peaks, and allele specific events can be found in the `peaks.vcf` file.

### Essential options

