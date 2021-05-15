# {Descriptive title of benchmark}

## Notes for writing benchmark specifications

{Carefully read this section, but do NOT include it in your final
specifications}

> * Specifications MUST be written in [GitHub markdown][spec-markdown]
>   following this template.
> * You SHOULD provide example data for each input and output data type. You MUST
>   provide example data for custom data types (i.e, types that are not well
>   known in the community).
> * You SHOULD avoid specifying custom file formats. If custom file formats
>   cannot be avoided, you SHOULD use the simplest option possible, ideally
>   building on existing, well-known file formats wherever possible (e.g.,
>   extend BED format by one or more custom columns). For custom file formats,
>   you MUST strive to describe them without ambiguities and with sufficient
>   detail, clarity and precision that a skilled reader can implement a simple
>   parser for them.
> * If referencing a standard or conventional file format, you SHOULD provide
>   a reference to its official specification. If not available, you MUST
>   provide a link to a detailed description of the file format. For custom
>   file formats, you MUST provide enough detailin enough detail in clarity
>   that it can be parsed.
> * If describing atomic data types, you MUST use the corresponding Python data
>   types `int`, `float`, `str` and `bool`; note that in these cases, values
>   MUST be provided in a form that allows casting them to the indicated types
>   without errors or further conversions, e.g., values specified as type `int`
>   MUST only include numerical characters.
> * For vectors and matrices, use types `vector` and `matrix` respectively.
> * Descriptions and units of metrics SHOULD be specified in a way that it will
>   be possible to create meaningful axis labels when visualizing metrics.
> * Compute performance benchmarks (run time, memory usage) MUST be computed on
>   infrastructure that guarantees comparable performance. Reach out to
>   technical support for further information.

## Synopsis

{Describe the benchmark in a few sentences, focusing on its rational}

## General info

* **Challenge:** {One of "Identification", "Quantification", and "Differential
  Usage"}
* **Identifier:** {Benchmark identifier as listed in the spreadsheet of
  bechmarks, e.g., "Q1"}

## Inputs

| # | Format | Link | Example data |
| --- | --- | --- | --- |
| 1 | {Format} | {URL to specification} | {link to example data in repo} |
| 2 | {Format} | {URL to specification} | {link to example data in repo} |
| ... | ... | ... | ... |
| n | {Format} | {URL to specification} | {link to example data in repo} |

### Additional info inputs

{Use this section to describe any custom formats or add any other additional
information. Use subeading to structure the section, e.g., when further
describing formats, create a subheading referring to the object you would like
to describe, e.g., "#### Format 1"}

## Outputs

| # | Format | Link | Example data |
| --- | --- | --- | --- |
| 1 | {Format} | {URL to specification} | {link to example data in repo} |
| 2 | {Format} | {URL to specification} | {link to example data in repo} |
| ... | ... | ... | ... |
| n | {Format} | {URL to specification} | {link to example data in repo} |

### Additional info outputs

{Use this section to describe any custom formats or add any other additional
information. Use subeading to structure the section, e.g., when further
describing formats, create a subheading referring to the object you would like
to describe, e.g., "#### Format 1"}

## Metrics

| # | Description | Unit | Compute from | Transformations | Type after transformations | Additional comments |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | {Description} | {Unit} | {Reference to output object} | {Describe in detail how the output sould be transformed to obtain the metric} | {Data type of metric after gathering and transformation, typically one of "vector" and "matrix"} | {Anything else you would like to add} |
| 2 | {Description} | {Unit} | {Reference to output object} | {Describe in detail how the output sould be transformed to obtain the metric} | {Data type of metric after gathering and transformation, typically one of "vector" and "matrix"} | {Anything else you would like to add} |
| ... | ... | ... | ... | ... | ... | ... |
| n | {Description} | {Unit} | {Reference to output object} | {Describe in detail how the output sould be transformed to obtain the metric} | {Data type of metric after gathering and transformation, typically one of "vector" and "matrix"} | {Anything else you would like to add} |

### Additional info metrics

{Use this section to describe any custom formats or add any other additional
information. Use subeading to structure the section, e.g., when further
describing formats, create a subheading referring to the object you would like
to describe, e.g., "#### Format 1"}

[//]: # (References)

[short-hand-ref]: <https://my-url-target.edu>
