# APAeval
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->

[![GitHub license](https://img.shields.io/github/license/iRNA-COSI/APAeval?color=orange)](https://github.com/iRNA-COSI/APAeval/blob/main/LICENSE)
[![All Contributors](https://img.shields.io/badge/all_contributors-26-orange.svg)](#contributors-)
[![DOI:zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.8290348-informational)](https://doi.org/10.5281/zenodo.8290348)
[![DOI:biorxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2023.06.23.546284-informational)](https://doi.org/10.1101/2023.06.23.546284)

<!-- ALL-CONTRIBUTORS-BADGE:END -->

Welcome to the [APAeval][apa-eval] GitHub repository.

**Quick links**
- [Use a benchmarked method on your own RNA-seq data](#use-a-benchmarked-method-on-your-own-rna-seq-data)
- [Benchmark a new method](#benchmark-a-new-method)
- [Extend APAeval's benchmarks](#extend-apaevals-benchmarks)

**APAeval** is a community effort that was born as the **APAeval hackathon at the RNA 2021 Conference**. We are aiming to evaluate computational methods for the detection and quantification of poly(A) sites from RNA-seq samples in an open, reproducible and extensible manner.




[![logo][apa-eval-logo]][apa-eval]

<!-- TOC -->

- [Overview of APAeval benchmarking](#overview-of-apaeval-benchmarking)
- [What can you do?](#what-can-you-do)
    - [Use a benchmarked method on your own RNA-seq data](#use-a-benchmarked-method-on-your-own-rna-seq-data)
    - [Benchmark a new method](#benchmark-a-new-method)
    - [Extend APAeval's benchmarks](#extend-apaevals-benchmarks)
- [Some technical stuff](#some-technical-stuff)
    - [OpenEBench](#openebench)
    - [APAeval conda environment](#apaeval-conda-environment)
    - [Tutorials](#tutorials)
- [Code of Conduct](#code-of-conduct)
- [Open Science, licenses & attribution](#open-science-licenses--attribution)
- [Get in touch](#get-in-touch)
- [Contributors ✨](#contributors-)

<!-- /TOC -->

## Overview of APAeval benchmarking

APAeval currently consists of three **benchmarking events**, each consisting of a set of **challenges** for bioinformatics methods (=**participants**) that use RNA-seq data to:

1. **Identify** polyadenylation sites
2. Report poly(A) site expression as **absolute quantification** in TPM
3. Report **relative expression** of poly(A) sites within transcripts
 
> We'd still like to set up a fourth event to evaluate tools that calculate **differential usage** of polyadenylation sites. If you'd like to contribute, continue reading [below](#extend-apaevals-benchmarks).

![schema][apa-eval-overview]

1. As described above, APAeval consists of three benchmarking events to evaluate the performance of different tasks that the methods of interest (=participants) might be able to perform: PAS identification, absolute quantification, and relative quantification. A method can participate in one, two or all three events, depending on its functions.
2. Raw data: For challenges within the benchmarking events, APAeval is using data from several different selected publications. Generally, one dataset (consisting of one or more samples) corresponds to one challenge (here, datasets for challenges x and y are depicted). All raw RNA-seq data is processed with nf-core/rna-seq for quality control and mapping. For each dataset we provide a matching ground truth file, created from 3’ end seq data from the same publications as the raw RNA-seq data, that will be used in the challenges to assess the performance of participants. You can find an overview of RNA-seq and matching ground truth samples in [the APAeval Zenodo snapshot][apaeval-zenodo].
3. Sanctioned input files: The processed input data is made available in .bam format. Additionally, for each dataset a gencode annotation in .gtf format, as well as a reference PAS atlas in .bed format for participants that depend on pre-defined PAS (not shown), are provided. 
4. In order to evaluate each participant in different challenges, a re-usable [“method workflow”][apaeval-mwf-readme] has to be written in either [Snakemake][snakemake] or [Nextflow][nf]. Within this workflow, all necessary pre- and post-processing steps that are needed to get from the input formats provided by APAeval (see 3.), to the output specified by APAeval in their metrics specifications (see 5.) have to be performed. 
5. To ensure compatibility with the workflows of the benchmarking events, [specifications for file formats][apaeval-specs] (output of method workflows = input for benchmarking workflows) are provided by APAeval. 
6. Within a benchmarking event, one or more challenges will be performed. A challenge is primarily defined by the input dataset used for performance assessment. Results of a challenge (metrics) are computed for each participant within a ["benchmarking workflow"][apaeval-bwfs]. 
7. In order to compare the performance of participants, results for each participant are uploaded to the [OEB database](#openebench), where metrics for all participants are visualized per challenge.

## What can you do?

### Use a benchmarked method on your own RNA-seq data
Firstly, you might want to check our [manuscript][manuscript] or our [OpenEBench site][apaeval-oeb] to find the method that would perform best for your use case. If you have decided on a method to use, head over to the [method workflows section in this repo][apaeval-mwf-readme] and follow the instructions in the `README.md` of the method of your choice. All our method workflows are built in either [Snakemake][snakemake] or [Nextflow][nf], and use [containers][docker] for individual steps to ensure reproducibility and reusability. For instructions on how to set up a [conda environment][conda] for running APAeval workflows [see here](#apaeval-conda-environment).

> You'll need to have your RNA-seq data ready in `.bam` format. No idea how to get there? You could check out the [nf-core][nf-core] [RNA-Seq analysis pipeline][nf-core-rna-seq] or other tools such as [ZARP][zarp].


### Benchmark a new method
Have you developed a new computational method for investigating APA from RNA-seq data? Or are you interested in one of the tools we haven't managed to include in APAeval yet? We'd be very happy if you decided to contribute to APAeval!

In order to ensure reproducibility of the benchmarks, as well as reusability and shareability of the benchmarked method, you'd start by writing an APAeval style [method workflow][apaeval-mwf-readme]. That workflow will take `.bam` files as an input, and create `.bed` files compatible with the [specification for the respective APAeval benchmarking event][apaeval-specs]. Create a PR (pull request; please ask in our [Github discussions board][discussions] to be added to APAeval as a collaborator, or create the PR from a fork) in this repo and wait for your request to be approved. You can then run the workflow on the [data for all APAeval challenges][apaeval-zenodo] and use the resulting `.bed` files in the corresponding [APAeval benchmarking workflow][apaeval-bwfs] in order to compare the performance of your tool to the [APAeval ground truths][apaeval-zenodo]. Finally you can submit your metrics `.json` files to us and we'll take care of including them in our [OEB site][apaeval-oeb]. 

### Extend APAeval's benchmarks
One of the main goals of APAeval is to provide *extensible* benchmarking, such that new tools, new challenges or new metrics can be added at any time. Therefore we warmly welcome any contribution to the project. A good starting point would be to visit our [issue][issues] and [discussion][discussions] boards. The latter one is also the place where you can reach out to us and request we add you to the repo as a collaborator (alternatively, create your PRs from a fork). You can then take on an existing task, suggest a new one, or start a discussion. 



## Some technical stuff
### OpenEBench

We are partnering with [OpenEBench][oeb], a benchmarking and technical
monitoring platform for bioinformatics tools. OpenEBench development,
maintenance and operation is coordinated by [Barcelona Supercomputing Center
(BSC)][bsc] together with partners from the European Life Science
infrastructure initiative [ELIXIR][elixir].

OpenEBench tooling will facilitate the computation and visualization of
benchmarking results and store the results of all benchmarking events and challenges in their databases, making it easy for others to explore results. This should
also make it easy to add additional participants to existing benchmarking events
later on. OpenEBench developers are also advising us on creating benchmarks
that are compatible with good practices in the wider community of
bioinformatics challenges.


### APAeval conda environment

For reproducible execution of our workflows (both method and benchmarking workflows) we're using a conda environment with fixed versions of Snakemake, Nextflow, some python packages, and Singularity. Make sure you have [conda][conda] installed and from the root directory of this repo create the APAeval environment with

```bash
conda env create -f apaeval_env.yaml
```

You can then activate it with:

```
conda activate apaeval
```

> NOTE: If you're working on Windows or Mac, you might have to google about setting up a virtual machine for running Singularity. 

> ANOTHER NOTE: If you run into problems regarding root access & Singularity with the described setup, try removing Singularity installation from the `apaeval_env.yaml` and [install it independently][singularity].


You can now execute the workflows!

### Tutorials
Here are some pointers and tutorials for the main software tools that we are using at APAeval:

Conda: [tutorial][tutorial-conda]   
Docker: [tutorial][tutorial-docker]   
Git: [tutorial][tutorial-git]   
GitHub: [general tutorial][tutorial-gh] / [GitHub flow tutorial][tutorial-gh-flow]  
Nextflow: [tutorial][tutorial-nextflow]  
Singularity: [tutorial][tutorial-singularity]  
Snakemake: [tutorial][tutorial-snakemake]  


## Code of Conduct

Please be kind to one another and mind the [Contributor Covenant's Code of
Conduct][coc-original] for all interactions with the community. A copy of the
Code of Conduct is also [shipped with this repository][coc-local]. Please
report any violations to the Code of Conduct to [apaeval@irnacosi.org][contact].

## Open Science, licenses & attribution

Following best practices for writing software and sharing data and code is
important to us, and therefore we want to apply, as much as possible, [FAIR
Principles][fair] to data and software alike. This includes publishing all
code open source, under permissive [licenses approved][osi-licenses] by the
[Open Source Initiative][osi] and all data by a permissive [Creative
Commons][cc] license.

In particular, we publish all code under the [MIT license][license-mit] and all
data under the [CC0 license][license-cc0]. An exception are all _benchmarking
workflows_, which are published under the [GPLv3 license][license-gplv3], as
the provided template is derived from an [OpenEBench][oeb] [example
workflow][oeb-example-workflow] that is itself licensed under GPLv3. A copy of
the MIT license is also [shipped with this repository][license].

We also believe that attribution, provenance and transparency are crucial for
an open and fair work environment in the sciences, especially in a community
effort like APAeval. Therefore, we would like to make clear from the beginning
that in all publications deriving from APAeval (journal manuscript, data and
code repositories), any non-trivial contributions will be acknowledged by
authorship. 

We expect that all contributors accept the license and attribution policies
outlined above.

## Get in touch
If you would like to contribute to APAeval or have any questions, we'd be happy to hear from you via our [Github Discussions board][discussions]. If you already have a specific issue in mind, feel free to add it to our [issues board][issues]. You can also reach out to [apaeval@irnacosi.org][contact].

## How to cite APAeval
If APAeval was useful for you in your work, please cite our [manuscript][manuscript]:

**Extensible benchmarking of methods that identify and quantify polyadenylation sites from RNA-seq data**  
Sam Bryce-Smith, Dominik Burri, Matthew R. Gazzara, Christina J. Herrmann, Weronika Danecka, Christina M. Fitzsimmons, Yuk Kei Wan, Farica Zhuang, Mervin M. Fansler, José M. Fernández, Meritxell Ferret, Asier Gonzalez-Uriarte, Samuel Haynes, Chelsea Herdman, Alexander Kanitz, Maria Katsantoni, Federico Marini, Euan McDonnel, Ben Nicolet, Chi-Lam Poon, Gregor Rot, Leonard Schärfen, Pin-Jou Wu, Yoseop Yoon, Yoseph Barash, Mihaela Zavolan  
*bioRxiv 2023.06.23.546284*; doi: https://doi.org/10.1101/2023.06.23.546284 

[apa-eval]: <https://irnacosi.org/2021/01/04/rna-society-2021-apaeval-challenge/>
[apa-eval-logo]: images/logo.png
[apaeval-oeb]: <https://openebench.bsc.es/benchmarking/OEBC007?event=OEBE0070000003>
[apa-eval-overview]: images/overview.png
[apaeval-mwf-readme]: ./method_workflows/README.md
[apaeval-specs]: ./method_workflows/method_workflow_file_specifications.md
[apaeval-bwfs]: ./benchmarking_workflows/README.md
[apaeval-zenodo]: <https://zenodo.org/record/8290348>
[bsc]: <https://www.bsc.es/>
[cc]: <https://creativecommons.org/>
[coc-local]: CODE_OF_CONDUCT.md
[coc-original]: <https://www.contributor-covenant.org/>
[conda]: <https://docs.conda.io/en/latest/>
[contact]: <mailto:apaeval@irnacosi.org>
[discussions]: <https://github.com/iRNA-COSI/APAeval/discussions>
[docker]: <https://www.docker.com/>
[elixir]: <https://elixir-europe.org/>
[fair]: <https://www.go-fair.org/fair-principles/>
[issues]: <https://github.com/iRNA-COSI/APAeval/issues>
[license]: LICENSE
[license-mit]: <https://opensource.org/licenses/MIT>
[license-cc0]: <https://creativecommons.org/publicdomain/zero/1.0/>
[license-gplv3]: <https://www.gnu.org/licenses/gpl-3.0.en.html>
[manuscript]: <https://www.biorxiv.org/content/10.1101/2023.06.23.546284v1>
[nf]: <https://www.nextflow.io/>
[nf-core]: <https://nf-co.re/>
[nf-core-rna-seq]: <https://nf-co.re/rnaseq>
[oeb]: <https://openebench.bsc.es/>
[oeb-example-workflow]: <https://github.com/inab/TCGA_benchmarking_dockers>
[osi]: <https://opensource.org/>
[osi-licenses]: <https://opensource.org/licenses>
[singularity]: <https://sylabs.io/singularity/>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[tutorial-conda]: <https://conda.io/projects/conda/en/latest/user-guide/getting-started.html>
[tutorial-docker]: <https://docs.docker.com/get-started/>
[tutorial-git]: <https://git-scm.com/docs/gittutorial>
[tutorial-gh]: <https://guides.github.com/activities/hello-world/>
[tutorial-gh-flow]: <https://www.youtube.com/watch?v=GgjIvUrOpmg>
[tutorial-nextflow]: <https://www.nextflow.io/blog/2020/learning-nextflow-in-2020.html>
[tutorial-singularity]: <https://singularity-tutorial.github.io/>
[tutorial-snakemake]: <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>
[zarp]: <https://github.com/zavolanlab/zarp>

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="http://chelseaherdman.com"><img src="https://avatars.githubusercontent.com/u/50838086?v=4?s=100" width="100px;" alt="Chelsea Herdman"/><br /><sub><b>Chelsea Herdman</b></sub></a><br /><a href="#projectManagement-chelseaherdman" title="Project Management">📆</a> <a href="#eventOrganizing-chelseaherdman" title="Event Organizing">📋</a> <a href="#ideas-chelseaherdman" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Achelseaherdman" title="Reviewed Pull Requests">👀</a> <a href="#talk-chelseaherdman" title="Talks">📢</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=chelseaherdman" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/ninsch3000"><img src="https://avatars.githubusercontent.com/u/36634279?v=4?s=100" width="100px;" alt="ninsch3000"/><br /><sub><b>ninsch3000</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=ninsch3000" title="Code">💻</a> <a href="#data-ninsch3000" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=ninsch3000" title="Documentation">📖</a> <a href="#design-ninsch3000" title="Design">🎨</a> <a href="#eventOrganizing-ninsch3000" title="Event Organizing">📋</a> <a href="#mentoring-ninsch3000" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-ninsch3000" title="Project Management">📆</a> <a href="#question-ninsch3000" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Aninsch3000" title="Reviewed Pull Requests">👀</a> <a href="#talk-ninsch3000" title="Talks">📢</a> <a href="#ideas-ninsch3000" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Aninsch3000" title="Bug reports">🐛</a> <a href="#tutorial-ninsch3000" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/EuancRNA"><img src="https://avatars.githubusercontent.com/u/46812323?v=4?s=100" width="100px;" alt="Euan McDonnell"/><br /><sub><b>Euan McDonnell</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=EuancRNA" title="Code">💻</a> <a href="#ideas-EuancRNA" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-EuancRNA" title="Mentoring">🧑‍🏫</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://git.scicore.unibas.ch/kanitz"><img src="https://avatars.githubusercontent.com/u/10855418?v=4?s=100" width="100px;" alt="Alex Kanitz"/><br /><sub><b>Alex Kanitz</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Auniqueg" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=uniqueg" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=uniqueg" title="Documentation">📖</a> <a href="#example-uniqueg" title="Examples">💡</a> <a href="#eventOrganizing-uniqueg" title="Event Organizing">📋</a> <a href="#ideas-uniqueg" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-uniqueg" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-uniqueg" title="Maintenance">🚧</a> <a href="#mentoring-uniqueg" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-uniqueg" title="Project Management">📆</a> <a href="#question-uniqueg" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Auniqueg" title="Reviewed Pull Requests">👀</a> <a href="#talk-uniqueg" title="Talks">📢</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://yuukiiwa.github.io/"><img src="https://avatars.githubusercontent.com/u/41866052?v=4?s=100" width="100px;" alt="Yuk Kei Wan"/><br /><sub><b>Yuk Kei Wan</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Ayuukiiwa" title="Bug reports">🐛</a> <a href="#blog-yuukiiwa" title="Blogposts">📝</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=yuukiiwa" title="Code">💻</a> <a href="#data-yuukiiwa" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=yuukiiwa" title="Documentation">📖</a> <a href="#example-yuukiiwa" title="Examples">💡</a> <a href="#eventOrganizing-yuukiiwa" title="Event Organizing">📋</a> <a href="#ideas-yuukiiwa" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-yuukiiwa" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-yuukiiwa" title="Project Management">📆</a> <a href="#question-yuukiiwa" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=yuukiiwa" title="Tests">⚠️</a> <a href="#tutorial-yuukiiwa" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/BenNicolet"><img src="https://avatars.githubusercontent.com/u/59935248?v=4?s=100" width="100px;" alt="Ben"/><br /><sub><b>Ben</b></sub></a><br /><a href="#data-BenNicolet" title="Data">🔣</a> <a href="#ideas-BenNicolet" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-BenNicolet" title="Project Management">📆</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/pjewell-biociphers"><img src="https://avatars.githubusercontent.com/u/63738218?v=4?s=100" width="100px;" alt="pjewell-biociphers"/><br /><sub><b>pjewell-biociphers</b></sub></a><br /><a href="#maintenance-pjewell-biociphers" title="Maintenance">🚧</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mzavolan"><img src="https://avatars.githubusercontent.com/u/44894507?v=4?s=100" width="100px;" alt="mzavolan"/><br /><sub><b>mzavolan</b></sub></a><br /><a href="#data-mzavolan" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mzavolan" title="Documentation">📖</a> <a href="#eventOrganizing-mzavolan" title="Event Organizing">📋</a> <a href="#financial-mzavolan" title="Financial">💵</a> <a href="#ideas-mzavolan" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-mzavolan" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-mzavolan" title="Project Management">📆</a> <a href="#question-mzavolan" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Amzavolan" title="Reviewed Pull Requests">👀</a> <a href="#talk-mzavolan" title="Talks">📢</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mfansler"><img src="https://avatars.githubusercontent.com/u/1182216?v=4?s=100" width="100px;" alt="Mervin Fansler"/><br /><sub><b>Mervin Fansler</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Amfansler" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mfansler" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mfansler" title="Documentation">📖</a> <a href="#eventOrganizing-mfansler" title="Event Organizing">📋</a> <a href="#ideas-mfansler" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-mfansler" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-mfansler" title="Project Management">📆</a> <a href="#question-mfansler" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Amfansler" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mkatsanto"><img src="https://avatars.githubusercontent.com/u/31883096?v=4?s=100" width="100px;" alt="Maria Katsantoni"/><br /><sub><b>Maria Katsantoni</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=mkatsanto" title="Code">💻</a> <a href="#ideas-mkatsanto" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-mkatsanto" title="Mentoring">🧑‍🏫</a> <a href="#question-mkatsanto" title="Answering Questions">💬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/daneckaw"><img src="https://avatars.githubusercontent.com/u/30384499?v=4?s=100" width="100px;" alt="daneckaw"/><br /><sub><b>daneckaw</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=daneckaw" title="Code">💻</a> <a href="#data-daneckaw" title="Data">🔣</a> <a href="#eventOrganizing-daneckaw" title="Event Organizing">📋</a> <a href="#ideas-daneckaw" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-daneckaw" title="Project Management">📆</a> <a href="#tutorial-daneckaw" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/dominikburri"><img src="https://avatars.githubusercontent.com/u/7873536?v=4?s=100" width="100px;" alt="Dominik Burri"/><br /><sub><b>Dominik Burri</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Adominikburri" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=dominikburri" title="Code">💻</a> <a href="#data-dominikburri" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=dominikburri" title="Documentation">📖</a> <a href="#example-dominikburri" title="Examples">💡</a> <a href="#eventOrganizing-dominikburri" title="Event Organizing">📋</a> <a href="#ideas-dominikburri" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-dominikburri" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#mentoring-dominikburri" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-dominikburri" title="Project Management">📆</a> <a href="#question-dominikburri" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=dominikburri" title="Tests">⚠️</a> <a href="#tutorial-dominikburri" title="Tutorials">✅</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mrgazzara"><img src="https://avatars.githubusercontent.com/u/58347523?v=4?s=100" width="100px;" alt="mrgazzara"/><br /><sub><b>mrgazzara</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=mrgazzara" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mrgazzara" title="Documentation">📖</a> <a href="#data-mrgazzara" title="Data">🔣</a> <a href="#eventOrganizing-mrgazzara" title="Event Organizing">📋</a> <a href="#ideas-mrgazzara" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-mrgazzara" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-mrgazzara" title="Maintenance">🚧</a> <a href="#projectManagement-mrgazzara" title="Project Management">📆</a> <a href="#mentoring-mrgazzara" title="Mentoring">🧑‍🏫</a> <a href="#talk-mrgazzara" title="Talks">📢</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/FitzsimmonsCM"><img src="https://avatars.githubusercontent.com/u/33811247?v=4?s=100" width="100px;" alt="Christina Fitzsimmons"/><br /><sub><b>Christina Fitzsimmons</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=FitzsimmonsCM" title="Documentation">📖</a> <a href="#eventOrganizing-FitzsimmonsCM" title="Event Organizing">📋</a> <a href="#ideas-FitzsimmonsCM" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-FitzsimmonsCM" title="Project Management">📆</a> <a href="#talk-FitzsimmonsCM" title="Talks">📢</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/lschaerfen"><img src="https://avatars.githubusercontent.com/u/76001100?v=4?s=100" width="100px;" alt="Leo Schärfen"/><br /><sub><b>Leo Schärfen</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=lschaerfen" title="Code">💻</a> <a href="#ideas-lschaerfen" title="Ideas, Planning, & Feedback">🤔</a> <a href="#talk-lschaerfen" title="Talks">📢</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://chilampoon.github.io"><img src="https://avatars.githubusercontent.com/u/43943114?v=4?s=100" width="100px;" alt="poonchilam"/><br /><sub><b>poonchilam</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=chilampoon" title="Code">💻</a> <a href="#ideas-chilampoon" title="Ideas, Planning, & Feedback">🤔</a> <a href="#question-chilampoon" title="Answering Questions">💬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/dseyres"><img src="https://avatars.githubusercontent.com/u/85220637?v=4?s=100" width="100px;" alt="dseyres"/><br /><sub><b>dseyres</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=dseyres" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=dseyres" title="Documentation">📖</a> <a href="#ideas-dseyres" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/plger"><img src="https://avatars.githubusercontent.com/u/9786697?v=4?s=100" width="100px;" alt="Pierre-Luc"/><br /><sub><b>Pierre-Luc</b></sub></a><br /><a href="#data-plger" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=plger" title="Documentation">📖</a> <a href="#eventOrganizing-plger" title="Event Organizing">📋</a> <a href="#ideas-plger" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-plger" title="Project Management">📆</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/SamBryce-Smith"><img src="https://avatars.githubusercontent.com/u/49978382?v=4?s=100" width="100px;" alt="SamBryce-Smith"/><br /><sub><b>SamBryce-Smith</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=SamBryce-Smith" title="Code">💻</a> <a href="#ideas-SamBryce-Smith" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3ASamBryce-Smith" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=SamBryce-Smith" title="Documentation">📖</a> <a href="#maintenance-SamBryce-Smith" title="Maintenance">🚧</a> <a href="#mentoring-SamBryce-Smith" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-SamBryce-Smith" title="Project Management">📆</a> <a href="#question-SamBryce-Smith" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3ASamBryce-Smith" title="Reviewed Pull Requests">👀</a> <a href="#tutorial-SamBryce-Smith" title="Tutorials">✅</a> <a href="#talk-SamBryce-Smith" title="Talks">📢</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/pinjouwu9325"><img src="https://avatars.githubusercontent.com/u/41280353?v=4?s=100" width="100px;" alt="Pin-Jou Wu"/><br /><sub><b>Pin-Jou Wu</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=pinjouwu9325" title="Code">💻</a> <a href="#ideas-pinjouwu9325" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/yoseopyoon"><img src="https://avatars.githubusercontent.com/u/84806078?v=4?s=100" width="100px;" alt="yoseopyoon"/><br /><sub><b>yoseopyoon</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=yoseopyoon" title="Code">💻</a> <a href="#ideas-yoseopyoon" title="Ideas, Planning, & Feedback">🤔</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/faricazjj"><img src="https://avatars.githubusercontent.com/u/25573986?v=4?s=100" width="100px;" alt="Farica Zhuang"/><br /><sub><b>Farica Zhuang</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Afaricazjj" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=faricazjj" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=faricazjj" title="Documentation">📖</a> <a href="#ideas-faricazjj" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-faricazjj" title="Maintenance">🚧</a> <a href="#projectManagement-faricazjj" title="Project Management">📆</a> <a href="#question-faricazjj" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Afaricazjj" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/AsierGonzalez"><img src="https://avatars.githubusercontent.com/u/11410715?v=4?s=100" width="100px;" alt="Asier Gonzalez"/><br /><sub><b>Asier Gonzalez</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3AAsierGonzalez" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=AsierGonzalez" title="Code">💻</a> <a href="#example-AsierGonzalez" title="Examples">💡</a> <a href="#ideas-AsierGonzalez" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-AsierGonzalez" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#mentoring-AsierGonzalez" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-AsierGonzalez" title="Project Management">📆</a> <a href="#question-AsierGonzalez" title="Answering Questions">💬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/txellferret"><img src="https://avatars.githubusercontent.com/u/63742994?v=4?s=100" width="100px;" alt="txellferret"/><br /><sub><b>txellferret</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=txellferret" title="Code">💻</a> <a href="#example-txellferret" title="Examples">💡</a> <a href="#ideas-txellferret" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-txellferret" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#mentoring-txellferret" title="Mentoring">🧑‍🏫</a> <a href="#question-txellferret" title="Answering Questions">💬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://grexor.github.io"><img src="https://avatars.githubusercontent.com/u/1798838?v=4?s=100" width="100px;" alt="Gregor Rot"/><br /><sub><b>Gregor Rot</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Agrexor" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=grexor" title="Code">💻</a> <a href="#ideas-grexor" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-grexor" title="Maintenance">🚧</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Agrexor" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/jmfernandez"><img src="https://avatars.githubusercontent.com/u/2072150?v=4?s=100" width="100px;" alt="José María Fernández"/><br /><sub><b>José María Fernández</b></sub></a><br /><a href="#ideas-jmfernandez" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-jmfernandez" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#mentoring-jmfernandez" title="Mentoring">🧑‍🏫</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
