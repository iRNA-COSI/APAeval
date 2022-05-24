# Differential usage benchmark specification

## Synopsis

Benchmark to test True Positive Rate (TPR) and False Positive Rate (FPR) of RNAseq-based identification of genes with differentially used poly(A) sites compared with ground truth (genes identified from orthogonal 3'end seq data).

The metrics are computed for two types of datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

### Input data:

Comparison of differential PAS usage predicted from RNA-Seq data with 3'end sequencing data

1. List of genes with information on differentially used poly(A) sites identified from RNA-Seq data using the benchmarked tool
2. List of genes with information on differentially used poly(A) sites identified from orthogonal 3'end seq dataset

Comparison of differential PAS usage predicted from simulated RNA-Seq data with dataset used for simulation

1. List of genes with information on differentially used poly(A) sites identified from simulated RNA-Seq data using the benchmarked tool
2. List of genes with information on differentially used poly(A) sites identified from the dataset used for simulation

### Metrics

Based on the input data the following metrics are computed:

1. TPR (True Positive Rate) = TP/(TP+FN)
2. FPR (False Positive Rate) = FP/(FP+TN)

TP - true positives - genes identified by the tool and present in the orthogonal dataset  
FP - false positives - genes identified by the tool and not present in the orthogonal dataset  
FN - false negatives - genes not identified by the tool but present in the orthogonal dataset  
TN - true negatives - genes not identified by the tool and not present in the orthogonal dataset

### OpenEBench challenges

The metrics are visualised using 2D scatter plots and barplots, as described in _Plots_ section, which are then used for ranking the participating tools.
The plots for a given input dataset belong to one benchmarking challenge as understood in OpenEBench schema.
Not all plots have to be prepared for each dataset, as described in _Plots_ section.
All differential usage challenges belong to the same benchmarking event.

## Inputs

| # | Format | Link | Example data |
  | --- | --- | --- | --- |
  | 1 | TSV | [Wikipedia][wiki-tsv] | [Link][in1] |
  | 2 | TSV | [Wikipedia][wiki-tsv] | [Link][in2] |

### Additional info inputs
  
#### Input 1

This TSV file contains information on differential usage of poly(A) sites identified from RNA-Seq data by the benchmarked method.  
Columns:

1. gene ID
2. significance of differential PAS usage

Column names are not added to the file.

#### Input 2

This TSV file contains information on differential usage of poly(A) sites identified from the orthogonal 3'end-seq dataset.  
Fields are the same as in Input 1.

## Plots

The results of this benchmark will be visualised in OpenEBench using the following plots:

### 1. **2D scatter plot** visualizing **TPR and FPR** of identification of genes with differentially used PAS.

**Plot type**: 2D scatter plot

**Metric X**: FPR  
**Metric Y**: TPR

**Ranking**: The best performing tool is the one with the highest TPR combined with lowest FPR (top left part of the plot) and the worst performing tool is the one with the lowest TPR combined with highest FPR (bottom right part of the plot). The plot should be divided into diagonal quartiles based on the distance from optimal performance.

Input datasets:

- RNA-Seq data compared with 3'end sequencing data
- Simulated RNA-Seq data compared with dataset used for simulation

## Outputs

Calculated metrics are saved in JSON file adhering to OpenEBench schema. 
Assessment output is generated for each tool separately and contains values of calculated metrics for a given input dataset.
Aggregation output contains summarized metrics values from all benchmarked tools for all plots of the given challenge.
Consolidation output combines assessment data for the participant currently benchmarked with aggregation of previously benchmarked participants and will be uploaded to the OEB DB.

| # | Format | Link | Example data | Description |
| :-- | :--- | :--- | :--- | :-------- |
| 1 | JSON | [Specification][spec-json] | [Link][assessment_out] | Assessment output JSON |
| 2 | JSON | [Specification][spec-json] | [Link][aggregation_out] | Aggregation output JSON
| 3 | JSON | [Specification][spec-json] | N/A | Consolidation output JSON |

### Additional info
 
#### Output 1

The OpenEBench assessment file contains the following attributes:

- **\_id** - follows the format: [COMMUNITY_ID]:[CHALLENGE_ID]\_[PARTICIPANT_ID]\_[METRIC_ID]
- **challenge_id** - challenge name
- **participant_id** - benchmarked tool
- **metrics**
	- **value** - metric value
	- **metric_id** - metric name; metric names used in this benchmark are specified in the table below
 
The following tables list the metric names, value types and units, and a description:

| Metric_id | Type | Unit | Description |
| :--- | :--- | :--- | :----------------- |
| `TPR` | `float` | N/A | True Positive Rate (TPR) of detection of genes with differential PAS usage compared with orthogonal dataset; TPR = (TP/(TP+FN)) |
| `FPR` | `float` | N/A | False Positive Rate of detection of genes with differential PAS usage compared with orthogonal dataset; FPR = FP/(FP+TN) |
 
#### Output 2

The OpenEBench aggregation file contains information required to produce the plots.

#### Output 3

The OpenEBench consolidation file contains all the information about the new benchmarking run and specifies visualization types (2D scatter plot, barplot), descriptions of metrics used for X and Y axis and (X,Y) value pairs for each challenge and challenge participants.
 

[//]: # (References)
  
[in1]: ./example_files/input1.tsv
[in2]: ./example_files/input2.tsv
[assessment_out]: ./example_files/assessment_out.json
[aggregation_out]: ./example_files/aggregation_out.json
[spec-json]: <https://www.ecma-international.org/publications-and-standards/standards/ecma-404/>
[wiki-tsv]: <https://en.wikipedia.org/wiki/Tab-separated_values>
  
