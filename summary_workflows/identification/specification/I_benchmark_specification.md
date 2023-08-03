# Identification benchmarks specification

## Synopsis

Benchmarks to assess the performance of RNAseq-based poly(A) site (PAS) identification compared with ground truth (poly(A)-sites identified from orthogonal 3'end seq data or existing PAS databases).

### Input data:

Comparison of PAS predicted from RNA-Seq data with 3'end sequencing data

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. Ground truth: Poly(A) sites identified based on orthogonal 3'end seq dataset

Comparison of PAS predicted from simulated RNA-Seq data with dataset used for simulation

1. Poly(A) sites identified based on simulated RNAseq data using the benchmarked tool
2. Ground truth: Poly(A) sites of transcripts used to simulate RNAseq reads

Calculating proportion of identified PAS assigned to different genomic features:

1. Poly(A) sites identified based on RNAseq data using the benchmarked tool
2. GTF file with genome annotation

### Metrics:

Based on the input data the following metrics are computed:

1. Sensitivity (TPR, True Positive Rate) = TP/(TP+FN)
2. Precision = TP/(TP+FP)
3. Poly(A) sites matched to multiple ground truth sites as proportion of all identified sites
4. Percentage of genes with correctly identified number of PAS
5. Poly(A) sites assigned to different genomic features: terminal exons, exons (excluding terminal exons), introns, intergenic regions (less and more than 1kb downstream of 3'-terminal exon)

TP - true positives - PAS identified by the tool and present in the orthogonal dataset  
FP - false positives - PAS identified by the tool and not present in the orthogonal dataset  
FN - false negatives - PAS not identified by the tool but present in the orthogonal dataset

The metrics should be computed for different distance thresholds between PAS identified by the tool from RNAseq dataset and PAS identified from the orthogonal 3'end seq dataset, i.e. the PAS identified by the tool should be within X nucleotides from the PAS identified from the orthogonal dataset for the prediction to be considered true.

### OpenEBench challenges

The metrics are visualised using 2D scatter plots and barplots, as described in _Plots_ section, which are then used for ranking the participating tools.
The plots for a given ground truth dataset belong to one benchmarking challenge as understood in OpenEBench schema.
Not all plots have to be prepared for each dataset, as described in _Plots_ section.
All identification challenges belong to the same benchmarking event.

## Inputs

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :----------|
| 1 | BED | [Specification][spec-bed] | [Link][in1] | BED file with PAS identification by the benchmarked tool |
| 2 | BED | [Specification][spec-bed] | [Link][in2] | BED file with PAS identification from ground truth dataset |
| 3 | GTF | [Specification][spec-gtf] | [Link][in3] | GTF file with genome annotation including 3' UTRs |
  
### Additional info
  
#### Input 1

This BED file contains single-nucleotide position of poly(A) sites identified by the benchmarked method.  
Fields:

- **chrom** - the name of the chromosome
- **chromStart** - the starting position of the feature in the chromosome; this corresponds to the last nucleotide just upstream of the cleavage and polyadenylation reaction; the starting position is 0-based, i.e. the first base on the chromosome is numbered 0
- **chromEnd** - the ending position of the feature in the chromosome; as identified PAS are single-nucleotide, the ending position is equal to `chromStart + 1`
- **name** - defines the name of the identified poly(A) site
- **score** - not used, leave as "."
- **strand** - defines the strand; either "." (=no strand) or "+" or "-".

#### Input 2

This BED file contains single-nucleotide position of poly(A) sites identified from the orthogonal 3'end-seq dataset.  
Fields are the same as in Input 1.

#### Input 3

This GTF file contains reference genome annotation for the RNAseq dataset including untranslated regions.
	
## Plots

The results of this benchmark will be visualised in OpenEBench using the following plots:

### 1. 2D scatter plot visualizing Sensitivity and Precision of poly(A) site identification. 

**Plot type**: 2D scatter plot

**Metric X**: Sensitivity(d)  
**Metric Y**: Precision(d)  
where _d - distance threshold_

**Ranking**: The best performing tool is the one with the highest Sensitivity combined with highest Precision (top right part of the plot) and the worst performing tool is the one with the lowest Sensitivity combined with lowest Precision (bottom left part of the plot).
The plot should be divided into diagonal quartiles based on the distance from optimal performance.

Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### 2. Bar plot visualizing Jaccard index for each tool.

**Plot type**: Bar plot

**Metric**: Jaccard index

**Ranking**: The best performing tool is the one with the highest Jaccard index.

Input data:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### 3. Bar plot visualizing percentage of genes with correctly identified number of PAS

**Plot type**: Bar plot

**Metric**: Percentage of genes with correctly identified number of PAS

**Ranking**: The best performing tool is the one with the highest ercentage of genes with correctly identified number of PAS.

Input data:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### Optional plots

This section includes plots that can be visualised outside of OpenEBench.

### 1. **2D line plot** visualising Sensitivity (Recall) and Precision in the form of Precision-Recall curve. 

**Plot type**: 2D line plot in which Sensitivity and Precision calculated for a range of distance thresholds and plotted together.

**Metric X**: Sensitivity(d)  
**Metric Y**: Precision(d)  
where _d - distance threshold_ in range 0-200 nt with 10 nt increments.

**Ranking**: Best performance is top-right.

Input datasets:

- RNA-Seq data compared with 3'end sequencing data

Note: 2D line plot is not supported in OpenEBench yet. If it's not implemented, the data should be visualised outside of OpenEBench.

### 2. Grouped bar plot visualizing proportions of PAS assigned to different genomic features in one plot.

**Plot type**: Grouped bar plot ([example][barplot_example]) in which bars corresponding to multiple featues are plotted using different colors and grouped by benchmarked tool.

**Metric**: Proportions of PAS assigned to different genomic features

**Ranking**: None

Features included in the plot:

- 3'-terminal exons
- exons (excluding terminal exons)
- introns
- intergenic region within 1 kb downstream of 3'-terminal exon
- intergenic region more than 1kb downstream of 3'-terminal exon

Input datasets:

- RNA-Seq data compared with 3'end sequencing data

Note: Grouped bar plot is not supported on OEB. Metrics should be visualized only outside of OpenEBench.
This plot should include results for benchmarked tools as well as for ground truth.

### 3. Bar plot visualizing proportion of PAS matched to multiple ground truth sites.

**Plot type**: Bar plot

**Metric**: Proportion of PAS matched to multiple ground truth sites

**Ranking**: The best performing tool is the one with the lowest proportion of PAS matched to multiple ground truth sites

Separate plots should be prepared for different values of distance threshold:

- 10 nt
- 50 nt
- 100 nt

Input data:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

## Outputs

Calculated metrics are saved in JSON file adhering to OpenEBench schema. 
Assessment output is generated for each tool separately and contains values of calculated metrics for a given input dataset.
Consolidation output contains summarized data from all benchmarked tools within one challenge in a format suitable for plotting, e.g. single values for barplot or X,Y value pairs for 2D scatter plot.

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :-------- |
| 1 | JSON | [Specification][spec-json] | [Link][assessment_out] | Assessment output JSON |
| 2 | JSON | [Specification][spec-json] | [Link][aggregation_out] | Aggregation output JSON
| 3 | JSON | [Specification][spec-json] | N/A | Consolidation output JSON |
  
### Additional info
 
#### Output 1

The OpenEBench assessment file contains the following attributes:

- **\_id** - follows the format: [COMMUNITY_ID]:[CHALLENGE_ID]\_[PARTICIPANT_ID]\_[METRIC_ID]\_[WINDOW_SIZE]nt
- **challenge_id** - unique identifier, corresponding to ground truth file name without extension. Must be of the format `[sample_name].([additional_info].)[genome_version]` where the part inside `()` is optional
- **participant_id** - benchmarked tool
- **metrics**:
	- **value** - metric value
	- **metric_id** - metric name; metric names used in this benchmark are specified in the table below
 
The following tables list the metric names, value types and units, and a description:

| Metric_id | Type | Unit | Description |
| :--- | :--- | :--- | :----------------- |
| `Precision_10nt` | `float` | N/A | Precision of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for 10 nt distance threshold |
| `Precision_50nt` | `float` | N/A | Precision of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for 50 nt distance threshold |
| `Precision_100nt` | `float` | N/A | Precision of PAS identification compared with orthogonal dataset; Precision = TP/(TP+FP); calculated for 100 nt distance threshold |
| `Sensitivity_10nt` | `float` | N/A | Sensitivity of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for 10 nt distance threshold |
| `Sensitivity_50nt` | `float` | N/A | Sensitivity of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for 50 nt distance threshold |
| `Sensitivity_100nt` | `float` | N/A | Sensitivity of PAS identification compared with orthogonal dataset; Sensitivity = (TP/(TP+FN)); calculated for 100 nt distance threshold |
| `Multi-matched_10nt` | `float` | % | Proportion of PAS identified from RNA-seq that were matched to multiple sites in grount truth; calculated for 10 nt distance threshold |
| `Multi-matched_50nt` | `float` | % | Proportion of PAS identified from RNA-seq that were matched to multiple sites in grount truth; calculated for 50 nt distance threshold |
| `Multi-matched_100nt` | `float` | % | Proportion of PAS identified from RNA-seq that were matched to multiple sites in grount truth; calculated for 100 nt distance threshold |
| `Genes_correct_PAS` | `float` | % | Percentage of genes with correctly identified number of PAS |

#### Output 2

The OpenEBench aggregation file contains information required to produce one plot.

#### Output 3

The OpenEBench consolidation file contains all the information about the new benchmarking run and specifies visualization types (2D scatter plot, barplot), descriptions of metrics used for X and Y axis and (X,Y) value pairs for each challenge and challenge participants.

[//]: # (References)
  
[in1]: ./example_files/input1.bed
[in2]: ./example_files/input2.bed
[in3]: ./example_files/input3.gtf
[assessment_out]: ./example_files/assessment_out.json
[aggregation_out]: ./example_files/aggregation_out.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[spec-bed]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
[spec-gtf]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format4>
[barplot_example]: <https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2_files/figure-html/thecode-1.png>
  
