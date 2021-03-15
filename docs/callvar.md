# Call variants

This function in MACS3 can be used to **directly** call genetics variants from the BAM files inside the ChIP-seq or ATAC-seq peaks. It can be invoked by `macs3 callvar` . The method was previously implemented in a separate project named [`SAPPER`](https://github.com/taoliu/SAPPER). The MACS3 port of SAPPER has been modified to directly work on the BAM files people usually use for peak calling. However, since `callvar` will randomly access chunks of alignments within given peaks, the BAM files as input to`callvar` have to be sorted by coordinates and indexed through `samtools`. In another word, the `.bai` files are required.

## Algorithm

## Workflow of `callvar`

<img src="./callvar_workflow.jpg" width="600" />

## Commandline options

If you type this command with `-h`, you will see a full description of command-line options. Here we only list the essentials.
