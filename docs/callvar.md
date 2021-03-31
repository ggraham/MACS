# Call variants

This function in MACS3 can be used to **directly** call genetics variants and infer allele specific events from the BAM files inside the ChIP-seq or ATAC-seq peaks, without the need for deep whole genome sequencing. It can be invoked by `macs3 callvar` . The method was previously implemented in a separate project named [`SAPPER`](https://github.com/taoliu/SAPPER). The MACS3 port of SAPPER has been modified to directly work on the BAM files people usually use for peak calling. However, since `callvar` will randomly access chunks of alignments within given peaks, the BAM files as input to`callvar` have to be sorted by coordinates and indexed through `samtools`. In another word, the `.bai` files are required.

## Algorithm

Targeted DNA assays such as ChIP-Seq or ATAC-seq provide high coverage at functional regulatory elements where proteins interact with DNA, we built this `callvar` algorithm based on the hypothesis that the DNA sequences at the targeted regulatory elements can be piled up or assembled ( into ‘unitigs’ ) for predicting SNV and short INDELs and infering allele specific events.

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

In general, user has to provide 1 BAM file from targeted DNA sequencing assay (ChIP-seq or ATAC-seq) `TREAT_sorted.bam`, and 1 BED file for regions of interest (i.e. peak regions) `peaks.bed`. Optionally, user can also provide a BAM file from a matching unbiased control sample, i.e. genomic input, or igg control `CTRL_sorted.bam`. The idea is that the BAM file from targeted DNA assay may have allele specific biases, whereas the control sample doesn't.  The typical command is

```
$ macs3 callvar -b peaks.bed -t TREAT_sorted.bam -c CTRL_sorted.bam -o peaks.vcf
```

The result of variants in peaks, and allele specific events can be found in the `peaks.vcf` file.

### Prepare Input Files

### Essential options

#### Region of interest in BED format `-b`

#### Treatment file in sorted BAM `-t`

#### Optional `-c CFILE`

#### Output directory and output filename `--outdir` and `-o`

#### Cutoff for calling heterozygous allele `-g`

#### Cutoff for calling homozygous allele (different with reference allele) `-G`

#### Controling when to use assembler `fermi-lite` `-F`

#### The option for `fermi-lite` assembler `--fermi-overlap`

#### Quality control threshold for calling variants `--top2alleles-mratio` `--altallele-count` and `--max-ar`

### Interpret Output file
