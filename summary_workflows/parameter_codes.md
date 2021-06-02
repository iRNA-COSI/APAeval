# Parameter codes

Parameter codes uniquely identify what combinations of input files and parameters were used to run an execution workflow. 
These codes are used in summary workflows to determine appropriate comparison files.

## Formatting

Input to summary workflows must be in the following format:

```
PARAMCODE/METHOD/PARAMCODE_METHOD_OUTCODE.ext
```

This document specifies the `PARAMCODE` component of this string. For more information on the full naming scheme refer to the [execution workflow README][ex_readme]. 


# Working Parameter code table


| PARAMCODE | dataset/samples                                                                                    | file name                         | matching ground truth                                                                              | run type | genome/annotation | explicit param\_settings | processing directive                             | benchmarks     |
| --------- | -------------------------------------------------------------------------------------------------- | --------------------------------- | -------------------------------------------------------------------------------------------------- | -------- | ----------------- | ------------------------ | ------------------------------------------------ | -------------- |
| AA        | SRR6795720                                                                                         | Mayr\_NB\_R1                      | SRR1005606<br>                                                                                     | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AB        | SRR6795721                                                                                         | Mayr\_NB\_R2                      | SRR1005607<br>                                                                                     | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AC        | SRR6795723                                                                                         | Mayr\_NB\_R3                      | SRR6795688<br>                                                                                     | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AD        | SRR6795724                                                                                         | Mayr\_NB\_R4                      | SRR6795689                                                                                         | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AE        | SRR6795718                                                                                         | Mayr\_CD5B\_R3                    | SRR6795684                                                                                         | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AF        | SRR6795719                                                                                         | Mayr\_CD5B\_R4                    | SRR6795685                                                                                         | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AG        | SRR6795726                                                                                         | Mayr\_M\_R2                       | SRR6795691                                                                                         | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AH        | SRR6795713                                                                                         | Mayr\_GC\_R2                      | SRR6795693                                                                                         | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AI        | SRR6795715                                                                                         | Mayr\_GC\_R1                      | SRR6795692                                                                                         | indiv    | GRCh38            | 4 cores                  | individual sample                                | I2,I3,I4,Q2,Q3 |
| AJ        | SRR6795720, SRR6795721, SRR6795723, SRR6795724, SRR6795718, SRR6795719, SRR6795726, SRR6795713, SRR6795715 | Mayr\_allBcell\_group             | SRR1005606, SRR1005607, SRR6795688, SRR6795689, SRR6795684, SRR6795685, SRR6795691, SRR6795693, SRR6795692 | group    | GRCh38            | 4 cores                  | merge all samples                                | I2,I3,I4       |
| AK        | SRR6795720, SRR6795721, SRR6795723, SRR6795724                                                        | Mayr\_NB\_group                   | SRR1005606, SRR1005607, SRR6795688, SRR6795689                                                        | group    | GRCh38            | 4 cores                  | merge all samples                                | Q2,Q3          |
| AL        | SRR6795718, SRR6795719                                                                              | Mayr\_CD5B\_group                 | SRR6795684, SRR6795685                                                                              | group    | GRCh38            | 4 cores                  | merge all samples                                | Q2,Q3          |
| AM        | SRR6795713, SRR6795715                                                                              | Mayr\_GC\_group                   | SRR6795693, SRR6795692                                                                              | group    | GRCh38            | 4 cores                  | merge all samples                                | Q2,Q3          |
| AN        | SRR6795720, SRR6795721, SRR6795723, SRR6795724; SRR6795718, SRR6795719                                  | Mayr\_NB\_group; Mayr\_CD5B\_group | SRR1005606, SRR1005607, SRR6795688, SRR6795689; SRR6795684, SRR6795685                                  | group    | GRCh38            | 4 cores                  | compare first set against second (";" delimited) | D2             |
| AO        | SRR6795720, SRR6795721, SRR6795723 ,SRR6795724; SRR6795713, SRR6795715                                  | Mayr\_NB\_group; Mayr\_GC\_group   | SRR1005606, SRR1005607, SRR6795688, SRR6795689; SRR6795693, SRR6795692                                  | group    | GRCh38            | 4 cores                  | compare first set against second (";" delimited) | D2             |
| AP        | SRR6795718, SRR6795719; SRR6795713, SRR6795715                                                        | Mayr\_CD5B\_group; Mayr\_GC\_group | SRR6795684, SRR6795685; SRR6795693, SRR6795692                                                        | group    | GRCh38            | 4 cores                  | compare first set against second (";" delimited) | D2             |

**Note:** The semi-colon (`;`) is used to delimit sets of samples that will be used for differential testing. 

# Examples

## MISO Identification on Grouped Samples
The output for PAS identification (output code "01") from running MISO on the all Mayr B cell replicates (Mayr_allBcell_group) would be

```
AJ/MISO/AJ_MISO_01.bed
```

## QAPA Quantification 
The output for quantification (output code "02") from running QAPA on the first naive B cell replicate (Mayr_NB_R1) would be 

```
AA/QAPA/AA_QAPA_02.bed
```

[//]: # (References)
  
[ex_readme]: /execution_workflows/README.md
