
# DaPars2 {Method name as specified in Algorithm table.}

## Rulegraph

![rulegraph](rulegraph.DaPars2.png)

## Input & pre-processing
> * Note: the chromosome name conatin no 'chr' prefix in the files below, this workflow will add chr to each file.
>   * Bam files
>   * Bed files (sample)
>   * Gtf annotation
>   * Gff3 annotation
>   * Bed file of the genome

{Describe input files and how they will be processed in order for the method to work. Describe how sample tables have to look like, and any other input that is needed (e.g. genome).}

## Params
>   * Coverage Threshold
{Describe parameters needed by your METHOD.}

## Output & post-processing
> * Output directory chr1-22, chrX, chrY for each sample. Each directory contains DaPars2 results.
> * 01_DaPars2.bed

{Describe output files and postprocessing steps if necessary.}

## Notes
> * This workflow uses `Dapars2_Multi_Sample.py` where one can assign chromosome name as a command line argument, whereas `DaPars2_Multi_Sample_Multi_Chr.py` is hardcoded for standard human chromosomes with the 'chr' prefix. Both scripts otherwise produce identical output with test data.
> * In [DaPars2 documentation](http://bioinfo.szbl.ac.cn/DaPars2/DaPars2.html), the input files are in wiggle format; however, the testing data that it provides is in bedgraph format.*
> * The [DaPars2 documentation](http://bioinfo.szbl.ac.cn/DaPars2/DaPars2.html) states a dependency on R but no scripts appears to use R. The Docker image does not install R and Dapars2 runs successfully without error on test data.

{Notes about the METHOD.
e.g. Did you have to adjust the method's soure code?
}
