# Summary workflows
This is where the summary workflows for APAeval live.

## Overview
![apaeval-swfs][apaeval-swfs]

APAeval consists of three benchmarking events to evaluate the performance of different tasks that the methods of interest (=participants) might be able to perform: poly(A) site identification, quantification, and assessment of their differential usage. A method can participate in one, two or all three events, depending on its functions.   

Within a benchmarking event, one or more challenges will be performed. A challenge is primarily defined by the input dataset used for performance assessment. A challenge is evaluated within a summary workflow, which is run on the OEB infrastructure. The summary workflow will compute all metrics relevant for the challenge. Summary workflows can be re-used between challenges, however, depending on the input dataset, different metrics might be calculated, and summary workflows might thus be adapted to individual challenges (Example here: in challenge Ix, metrics I1 and I2 are computed, whereas in challenge Iy, an additional metric I5 is assessed. Apart from calculating metric I5 however, the summary workflows for the challenges Ix and Iy are the same.)    

In order to compare the performance of participants within a challenge, the respective summary workflow will be run on output files from all eligible participant execution workflows. In a first step the provided files are validated. Subsequently, all required metrics (scripts In for Identification, Qn for quantification, Dn for differential usage) are computed, using the matched ground truth files, if applicable. Finally, the results will be gathered in OEB specific .json files per participant.
Based on the created .json files, OEB will visualize all results per challenge, such that performance of participants can be compared for each metric.


## HOW TO
For an example of a summary workflow and further instructions, refer to the [Q2 summary workflow][q2-swf].

[//]: # (References)
[apaeval-swfs]: ../images/SWFs.png
[q2-swf]: quantification/Q2/README.md