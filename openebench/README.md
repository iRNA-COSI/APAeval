# OpenEBench requirements

The following schema groups needed to be provided in order to register benchamrks
with [OpenEBench][oeb]:

* **Tools** and **References** (to tools): for [participants][methods] benchmarked in APAeval.

* **Contacts**: Contact persons representing the APAeval community.

To see how a complete case looks like, including the
directory structure, here is a [full set of examples][example-full].

* **Challenges**/**Metrics**: Have been created for APAeval in close collaboration with OEB team members and are not published here.

When sending/uploading APAeval result json files to OEB, filtering has to be done in order to comply with the APAeval schemas registered in OEB. Use the [filter jsons util][filter-jsons] from the APAeval utils directory to **remove metrics that are NOT listed below**.

### OEB registered metrics for APAeval Identification
- F1_score:[10/50/100]nt
- Jaccard_index:[10/50/100]nt
- Precision:[10/50/100]nt
- Sensitivity:[10/50/100]nt
- percentage_genes_w_correct_nPAS

### OEB registered metrics for APAeval absolute Quantification
- Sum_FP_TPM:[10/50/100]nt
- Percent_FP_TPM:[10/50/100]nt
- Sensitivity:[10/50/100]nt
- Precision:[10/50/100]nt
- F1_score:[10/50/100]nt
- Jaccard_index:[10/50/100]nt
- Pearson_r:intersection:[10/50/100]nt

### Example: Filtering for 2023 OEB upload
```
python filter_jsons.py --file-list [PARTICIPANT]_consolidated.[EVENT].json --out_prefix "filtered_" --b_metrics "union all_GT Spearman relative 25nt multi_matched FDR" --b_challenges "TE"
```
[example-full]: <https://github.com/inab/benchmarking-data-model/tree/master/prototype-data/1.0.x/QfO>
[filter-jsons]: ../utils/README.md
[methods]: <https://docs.google.com/spreadsheets/d/1jvZOpF8iKzFRltrg99dgmkcHNXQVvaUla2JEjt5vyCU/edit>
[oeb]: <https://openebench.bsc.es/>
