# CSI-UTR

https://github.com/UofLBioinformatics/CSI-UTR

## Input & pre-processing
Input file in bed reformat
CSI annotation file --> probably thats why it does not work yet

I make a sample table .txt from the samplesheet via an R file

## Params

-genome=<genome>                    (default: Rn6 -- Other options are Mm10, Hg38, Rn6_extended)
 -r=<read_length>                    (default: 75)
 -bed=<CSI_bed_file>                 (default: ./data/locations/Rn6.CSIs.bed)
 -annot=<CSI_annotation_file>        (default: ./data/annotations/Rn6.CSIs.annot.bed)
 -out=<output directory>             (default: ./CSI_OUT/)
 -coverage_cut=<coverage cutoff>     (default: 0.08)
 -usage_cut=<usage cutoff>           (default: 1.0)
 -p=<p value significance cutoff)    (default: 0.05)
 -q=<FDR significance cutoff)        (default: 0.10)

## Output & post-processing

Not done yet

## Notes

Zero output at the moment
