# APAeval
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-15-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

Welcome to the [APAeval][apa-eval] GitHub repository.

**APAeval** is a community effort to evaluate computational methods for the
detection and quantification of poly(A) sites and the estimation of their
differential usage across RNA-seq samples.

[![logo][apa-eval-logo]][apa-eval]

## What is APAeval?

APAeval consists of three challenges, each consisting of a set of benchmarks
for tools that:

1. **Identify** polyadenylation sites
2. **Quantify** polyadenylation sites
3. Calculate **differential usage** of polyadenylation sites

For more info, please refer to our [landing page][apa-eval].

## How to get involved?

If you would like to contribute to APAeval, the first things you would need to
do are:

- [Register][apa-eval] for one or more of the APAeval teams
- Please use [this form][form-service-accounts] to provide us with your user
  names/handles of the various service accounts we are using (see
  [below](#how-do-we-work)) so that we can add you to the corresponding
  repositories/organizations
- Wait for us to reach out to you with an invitation to our Slack space, then
  use the invitation to sign up and post a short intro message in the
  [`#general`][slack-general] channel

## What is there to do?

The bulk of the work falls into roughly two tasks, writing method _execution
workflows_ and benchmark _summary workflows_ as highlighted in the figure
below:

![schema][apa-eval-schema]

### Execution workflows

_Execution workflows_ contain all steps that need to be run _per method_:

1. **Pre-processing:** Convert the input files the APAeval team has prepared
  into the input files a given method consumes, if applicable.
2. **Method execution:** Execute the method in any way necessary to compute the
  output files for _all_ challenges (may require more than one run of the tool
  if, e.g., run in different execution modes).
3. **Post-processing:** Convert the output files of the method into the formats
  consumed by the _summary workflows_ as specified by the APAeval team, if
  applicable.

_Execution workflows_ should be implemented in either [Nexflow][nf] or
[Snakemake][snakemake], and individual steps should be isolated through the
use of either [Conda][conda] virtual environments (deprecated; to run on AWS we need containerized workflows) or
[Docker][docker]/[Singularity][singularity] containers.

### Summary workflows

**_Summary workflows_** contain all steps that need to be run _per benchmark_,
using outputs of the invididual method _execution workflows_ as inputs. They
follow the [OpenEBench][oeb] workflow model, described
[here][oeb-example-workflow], implemented in Nextflow. OpenEBench workflows
consist of the following 3 steps:

1. **Validation:** Output data of the various _execution workflows_ is
  validated against the provided specifications to data consistency.
2. **Metrics computation:** Actual benchmarking metrics are computed as
  specified, e.g., by comparisons to ground truth/gold standard data sets.
3. **Consolidation:** Data is consolidated for consistency with other
  benchmarking efforts, based on the [OpenEBench/ELIXIR benchmarking data
  model][oeb-data-model].

Following the OpenEBench workflows model also ensures that result
visualizations are automatically generated, as long as the desired graph types
are supported by OpenEBench.

### Miscellaneous

Apart from writing _execution_ and _summary workflows_, there are various other
smaller jobs that participants can work on, including, e.g.:

- Pre-processing RNA-Seq input data via the [nf-core][nf-core] [RNA-Seq
  analysis pipeline][nf-core-rna-seq]
- Writing additional benchmark specifications
- Housekeeping jobs (reviewing code, helping to keep the repository clean,
  enforce good coding practices, etc.)
- Work with our partner OpenEBench on their open issues, e.g., by extending
  their portfolio of [supported visualization][oeb-open-issues]

If you do not know where to start, simply ask us!

## How do we work?

To account for everyone's different agendas and time zones, the hackathon is
organized such that contributors can work, as much as possible, in their own
time.

### Open Science, licenses & attribution

Following best practices for writing software and sharing data and code is
important to us, and therefore we want to apply, as much as possible, [FAIR
Principles][fair] to data and software alike. This includes publishing all
code open source, under permissive [licenses approved][osi-licenses] by the
[Open Source Initiative][osi] and all data by a permissive [Creative
Commons][cc] license.

In particular, we publish all code under the [MIT license][license-mit] and all
data under the [CC0 license][license-cc0]. An exception are all _summary
workflows_, which are published under the [GPLv3 license][license-gplv3], as
the provided template is derived from an [OpenEBench][oeb] [example
workflow][oeb-example-workflow] that is itself licensed under GPLv3. A copy of
the MIT license is also [shipped with this repository][license].

We also believe that attribution, provenance and transparency are crucial for
an open and fair work environment in the sciences, especially in a community
effort like APAeval. Therefore, we would like to make clear from the beginning
that in all publications deriving from APAeval (journal manuscript, data and
code repositories), any non-trivial contributions will be acknowledged by
authorship. All authors will be strictly listed alphabetically, by last name,
with no exceptions, wherever possible under the name of **The APAeval Team**
and accompanied by a more detailed description of how people contributed.

We expect that all contributors accept the license and attribution policies
outlined above.

### Communication

#### Chat

We are making use of [**Slack**][slack] (see [above](#how-to-get-involved) to
see how you can become a member) for asynchronous communication. Please use the
most appropriate channels for discussions/questions etc.:

- [`#general`][slack-general]: Introduce yourself and find pointers to get you
  started. All APAeval-wide announcements will be put here!
- [`#admin`][slack-admin]: Ask questions about the general organization of
  APAeval.
- [`#tech-support`][slack-tech-support]: Ask questions about the
  technical infrastructure and relevant software, e.g., AWS, GitHub, Nextflow,
  Snakemake.
- [`#challenge-differential-usage`][slack-diff-usage]: Discuss the poly(A) site
  differential usage challenge.
- [`#challenge-identification`][slack-identification]: Discuss the poly(A) site
  identification challenge.
- [`#challenge-quantification`][slack-quantification]: Discuss the poly(A) site
  quantification challenge.
- `#tool-{METHOD_NAME}`: Discussions channels for each indidvidual method.
- [`#random`][slack-random]: Post anything that doesn't fit into any of the
  other channels.

#### Video calls

Despite the event taking place mostly asynchrounously, we do have a few video
calls to increase the feeling of collaboration. In particular, we hae a kickoff
and a wrap-up meeting, three reporting calls and the presentation of
APAeval during the RNA Meeting. Please join meetings if you can, we have tried
to shuffle the times a bit, so that most people should be able to attend at
least some of the calls, irrespective of your location.

This calendar contains all video call events, including the necessary login
info, and we would like to kindly ask you to subscribed to it:

- Calendar ID: `59bboug9agv30v32r6bvaofdo4@group.calendar.google.com`
- [Public address][calendar-url]

> Please do not download the ICS file and then import it, as any updates to the
> calendar will not be synced. Instead, copy the calendar ID or public address
> and paste it in the appropriate field of your calendar application. Refer to
> your calendar application's help pages if you do not know how to subscribeh
> to a calendar.

With the exception of presenting APAeval during the RNA Meeting, all video
calls will take place in the following [**Zoom**][zoom] room:

- [Direct link][vc-direct-link]
- Meeting ID: `656 9429 1427`
- Passcode: `APAeval`

There is also a [meeting agenda][vc-agenda].

> For more lively meetings, participants are encouraged to switch on their
> cameras. But please mute your microphones if you are not currently speaking.

### Social coding

We are making extensive use of [**GitHub**][gh]'s project management resources
to allow contributors to work independently on individual, largely
self-contained, issues. There are several Kanban [project boards][gh-projects],
listing relevant issues for different kinds of tasks, such as drafting
benchmarking specifications and implementing/running method execution
workflows.

The idea is that people assign themselves to open issues (i.e., those issues
that are not yet assigned to someone else). Note that in order to be able to
do so, you will need to be a member of this GitHub repository (see
[above](#how-to-get-involved) to see how you can become a member). Once you
have assigned yourself, you can move/drag the issue from the **To do** to the
**In progress** column of the Kanban board.

When working on an issue, please start by cloning (preferred) or forking the
repository. Then create a feature branch and implement your code/changes. Once
you are happy with them, please create a pull request against the `main`
branch, making sure to fill in the provided template (in particular, please
refer to the original issue you are addressing with this pull request) and to
assign two reviewers. This workflow ensures collaborative coding and is
sometimes referred to as [GitHub flow][gh-flow]. If you are not familiar with
Git, GitHub or the GitHub flow, there are many useful tutorials online, e.g.,
those [listed below](#software).

### Cloud infrastructure

[AWS][aws] kindly sponsored credits for their compute and storage
infrastructure that we can use to run any heavy duty computations in the cloud
(e.g., RNA-Seq data pre-processing by or method execution workflows).

This also includes credits to run [Seqara Lab][seqera-labs]'s
[Nextflow Tower][nf-tower], a convenient web-based platform to run
[Nextflow][nf] workflows, such as the [nf-core][nf-core]
[RNA-Seq analysis workflow][nf-core-rna-seq] we are using for pre-processing
RNA-Seq data. Seqera Labs has kindly offered to hold a workshop on Nextflow
and Nextflow Tower and will be providing technical support during the
hackathon.

Setting up the AWS organization and infrastructure is still ongoing, and we
will update this section with more information as soon as that is done.

### OpenEBench

We are partnering with [OpenEBench][oeb], a benchmarking and technical
monitoring platform for bioinformatics tools. OpenEBench development,
maintenance and operation is coordinated by [Barcelona Supercomputing Center
(BSC)][bsc] together with partners from the European Life Science
infrastructure initiative [ELIXIR][elixir].

OpenEBench tooling will facilitate the computation and visualization of
benchmarking results and store the results of all challenges and benchmarks
in their databases, making it easy for others to explore results. This should
also make it easy to add additional methods to existing benchmarking challenges
later on. OpenEBench developers are also advising us on creating benchmarks
that are compatible with good practices in the wider community of
bioinformatics challenges.

The OpenEBench team will give a presentation to APAeval participants on how the
service is structured and what it offers to its visitors and the APAeval
effort.

### Software

Here are some pointers and tutorials for the main software tools that we are
going to use throughout the hackathon:

- [Conda][conda]: [tutorial][tutorial-conda]
- [Docker][docker]: [tutorial][tutorial-docker]
- [Git][git]: [tutorial][tutorial-git]
- [GitHub][gh]: [general tutorial][tutorial-gh] / [GitHub flow
  tutorial][tutorial-gh-flow]
- [Nextflow][nf]: [tutorial][tutorial-nextflow]
- [Singularity][singularity]: [tutorial][tutorial-singularity]
- [Snakemake][snakemake]: [tutorial][tutorial-snakemake]

Note that you don't need to know about all of these, e.g., one of Conda (deprecated; to run on AWS we need containerized workflows), Docker
and/or Singularity will typically be enough. [See
below](#nextflow-or-snakemake), for a discussion of the supported workflow
languages/management systems. Again, working with one will be enough for most
issues.

In addition to these, basic programming/scripting skills may be required for
most, but not for all issues. For those that do, you are generally free to
choose your preferred language, although for those people who have experience
with Python, we recommend you to go with that. It just makes it easier for
others to review your code, and it typically integrates better with our
templates and the general bioinformatics ecosystem/community.

Note that even if you don't have experience with any of these tools/languages,
and you don't feel like or have no time learning them, there is surely still
something that you can help us with. Just ask us and we will try to put your
individual skills to good use! :muscle:

#### Nextflow or Snakemake?

As mentioned [further above](#what-is-there-to-do), we would like _execution
workflows_ to be written in one of two "workflow languages": [Nextflow][nf] or
[Snakemake][snakemake]. Specifying workflows in such a language rather than,
say, by stringing together Bash commands, is considered good practice, as it
increases reusability and reproducibility, and so is in line with our goal of
adhering to [FAIR][fair] software principles.

But why only Nextflow and Snakemake, and not, e.g., the [Common Workflow
Language][cwl], the [Workflow Definition Language][wdl] or [Galaxy][galaxy]?
There are no particular reasons other than that APAeval organizers have
experience with these workflows languages and are thus able to provide
technical support. If you are an experienced workflow developer and prefer
another workflow language, you are welcome to use that one instead, but note
that we have no templates available and will not be able to help you much in
case you encounter issues during development or execution.

As for _summary workflows_, we are bound to implement these in Nextflow, as
they are executed on [OpenEBench][oeb], which currently only accepts Nextflow
workflows.

For this reason, as well as the fact that we will provide [Nextflow
Tower][nf-tower] for convient execution of Nextflow workflows on [AWS][aws]
cloud infrastructure ([see above](#cloud-infrastructure)) and use a [Nextflow
analysis pipeline][nf-core-rna-seq] for pre-processing RNA-Seq data sets, we
recommend novices without any other considerations (e.g., colleagues already
working with Snakemake) to use Nextflow.

#### Conda environment file

In order to execute scripts with either Nextflow or Snakemake in a reproducible
manner, we need to ensure the versions of these software are specified. In order 
to do that, we created a Conda environment file that contains specific versions 
Nextflow, Snakemake and some core libraries. To use this environment, you first
need to create it by using:

```bash
conda env create -f apaeval_env.yaml`
```

You then need to activate the environment with:

```
conda activate apaeval_execution_workflows
```

You can now execute the workflows!

### Code of Conduct

Please be kind to one another and mind the [Contributor Covenant's Code of
Conduct][coc-original] for all interactions with the community. A copy of the
Code of Conduct is also [shipped with this repository][coc-local]. Please
report any violations to the Code of Conduct to either or both of
[Christina][coc-contact-christina] and [Alex][coc-contact-alex] via Slack.

[apa-eval]: <https://irnacosi.org/2021/01/04/rna-society-2021-apaeval-challenge/>
[apa-eval-logo]: images/logo.png
[apa-eval-members]: <https://docs.google.com/document/d/1G7u-WQ6C-I_sXZ-15CIBw2iNgw6jkTNo7hnRTjci_b4/edit#heading=h.tarrapa8v8n6>
[apa-eval-schema]: images/schema.png
[aws]: <http://aws.amazon.com/>
[bsc]: <https://www.bsc.es/>
[calendar-url]: <https://calendar.google.com/calendar/ical/59bboug9agv30v32r6bvaofdo4%40group.calendar.google.com/public/basic.ics>
[cc]: <https://creativecommons.org/>
[coc-contact-alex]: <https://app.slack.com/client/T01PW9SAN7K/D01PP4WK7TL/user_profile/U01PEJ5TW4V>
[coc-contact-christina]: <https://app.slack.com/client/T01PW9SAN7K/D01PP4WK7TL/user_profile/U01PV9T8V9A>
[coc-local]: CODE_OF_CONDUCT.md
[coc-original]: <https://www.contributor-covenant.org/>
[conda]: <https://docs.conda.io/en/latest/>
[contact]: <mailto:apaeval@irnacosi.org>
[cwl]: <https://www.commonwl.org/>
[docker]: <https://www.docker.com/>
[elixir]: <https://elixir-europe.org/>
[fair]: <https://www.go-fair.org/fair-principles/>
[form-service-accounts]: <https://forms.gle/eKCHe5GWtvGrriek8>
[galaxy]: <https://usegalaxy.org/>
[gh]: <http://github.com/>
[gh-flow]: <https://guides.github.com/introduction/flow/>
[gh-join]: <https://github.com/join>
[gh-projects]: <https://github.com/iRNA-COSI/APAeval/projects/>
[git]: <https://git-scm.com/>
[license]: LICENSE
[license-mit]: <https://opensource.org/licenses/MIT>
[license-cc0]: <https://creativecommons.org/publicdomain/zero/1.0/>
[license-gplv3]: <https://www.gnu.org/licenses/gpl-3.0.en.html>
[nf]: <https://www.nextflow.io/>
[nf-core]: <https://nf-co.re/>
[nf-core-rna-seq]: <https://nf-co.re/rnaseq>
[nf-tower]: <https://tower.nf/>
[oeb]: <https://openebench.bsc.es/>
[oeb-data-model]: <https://github.com/inab/benchmarking-data-model>
[oeb-example-workflow]: <https://github.com/inab/TCGA_benchmarking_dockers>
[oeb-open-issues]: <https://github.com/inab/OpenEBench_scientific_visualizer/issues>
[osi]: <https://opensource.org/>
[osi-licenses]: <https://opensource.org/licenses>
[seqera-labs]: <https://seqera.io/>
[singularity]: <https://sylabs.io/singularity/>
[slack]: <http://slack.com/>
[slack-admin]: <https://apaeval.slack.com/archives/C01PEJQEUMT>
[slack-diff-usage]: <https://apaeval.slack.com/archives/C022DBAFSE7>
[slack-identification]: <https://apaeval.slack.com/archives/C0232UJQREU>
[slack-quantification]: <https://apaeval.slack.com/archives/C022A2PM4GM>
[slack-general]: <https://apaeval.slack.com/archives/C01PHLQKNH0>
[slack-random]: <https://apaeval.slack.com/archives/C01Q7FMRJ3A>
[slack-tech-support]: <https://apaeval.slack.com/archives/C022RNSAUV7>
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[tutorial-conda]: <https://conda.io/projects/conda/en/latest/user-guide/getting-started.html>
[tutorial-docker]: <https://docs.docker.com/get-started/>
[tutorial-git]: <https://git-scm.com/docs/gittutorial>
[tutorial-gh]: <https://guides.github.com/activities/hello-world/>
[tutorial-gh-flow]: <https://www.youtube.com/watch?v=GgjIvUrOpmg>
[tutorial-nextflow]: <https://www.nextflow.io/blog/2020/learning-nextflow-in-2020.html>
[tutorial-singularity]: <https://singularity-tutorial.github.io/>
[tutorial-snakemake]: <https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>
[vc-agenda]: <https://docs.google.com/document/d/1Cl3xq7_uwApAYxUbzeVSBsRfGUmtRc0jSnZ3yrWM3ks/edit#>
[vc-direct-link]: <https://unibas.zoom.us/j/65694291427?pwd=QUMyMjQ2SSt2eS9iZW50YVZCOC8wQT09>
[wdl]: <https://github.com/openwdl/wdl>
[zoom]: <https://zoom.us/>

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="http://chelseaherdman.com"><img src="https://avatars.githubusercontent.com/u/50838086?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Chelsea Herdman</b></sub></a><br /><a href="#projectManagement-chelseaherdman" title="Project Management">📆</a> <a href="#eventOrganizing-chelseaherdman" title="Event Organizing">📋</a> <a href="#ideas-chelseaherdman" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Achelseaherdman" title="Reviewed Pull Requests">👀</a> <a href="#talk-chelseaherdman" title="Talks">📢</a></td>
    <td align="center"><a href="https://github.com/ninsch3000"><img src="https://avatars.githubusercontent.com/u/36634279?v=4?s=100" width="100px;" alt=""/><br /><sub><b>ninsch3000</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=ninsch3000" title="Code">💻</a> <a href="#data-ninsch3000" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=ninsch3000" title="Documentation">📖</a> <a href="#design-ninsch3000" title="Design">🎨</a> <a href="#eventOrganizing-ninsch3000" title="Event Organizing">📋</a> <a href="#mentoring-ninsch3000" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-ninsch3000" title="Project Management">📆</a> <a href="#question-ninsch3000" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Aninsch3000" title="Reviewed Pull Requests">👀</a> <a href="#talk-ninsch3000" title="Talks">📢</a></td>
    <td align="center"><a href="https://github.com/EuancRNA"><img src="https://avatars.githubusercontent.com/u/46812323?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Euan McDonnell</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=EuancRNA" title="Code">💻</a> <a href="#ideas-EuancRNA" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-EuancRNA" title="Mentoring">🧑‍🏫</a></td>
    <td align="center"><a href="https://git.scicore.unibas.ch/kanitz"><img src="https://avatars.githubusercontent.com/u/10855418?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Alex Kanitz</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Auniqueg" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=uniqueg" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=uniqueg" title="Documentation">📖</a> <a href="#example-uniqueg" title="Examples">💡</a> <a href="#eventOrganizing-uniqueg" title="Event Organizing">📋</a> <a href="#ideas-uniqueg" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-uniqueg" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-uniqueg" title="Maintenance">🚧</a> <a href="#mentoring-uniqueg" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-uniqueg" title="Project Management">📆</a> <a href="#question-uniqueg" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Auniqueg" title="Reviewed Pull Requests">👀</a> <a href="#talk-uniqueg" title="Talks">📢</a></td>
    <td align="center"><a href="https://yuukiiwa.github.io/"><img src="https://avatars.githubusercontent.com/u/41866052?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Yuk Kei Wan</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Ayuukiiwa" title="Bug reports">🐛</a> <a href="#blog-yuukiiwa" title="Blogposts">📝</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=yuukiiwa" title="Code">💻</a> <a href="#data-yuukiiwa" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=yuukiiwa" title="Documentation">📖</a> <a href="#example-yuukiiwa" title="Examples">💡</a> <a href="#eventOrganizing-yuukiiwa" title="Event Organizing">📋</a> <a href="#ideas-yuukiiwa" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-yuukiiwa" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-yuukiiwa" title="Project Management">📆</a> <a href="#question-yuukiiwa" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=yuukiiwa" title="Tests">⚠️</a> <a href="#tutorial-yuukiiwa" title="Tutorials">✅</a></td>
    <td align="center"><a href="https://github.com/BenNicolet"><img src="https://avatars.githubusercontent.com/u/59935248?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ben</b></sub></a><br /><a href="#data-BenNicolet" title="Data">🔣</a> <a href="#ideas-BenNicolet" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-BenNicolet" title="Project Management">📆</a></td>
    <td align="center"><a href="https://github.com/pjewell-biociphers"><img src="https://avatars.githubusercontent.com/u/63738218?v=4?s=100" width="100px;" alt=""/><br /><sub><b>pjewell-biociphers</b></sub></a><br /><a href="#maintenance-pjewell-biociphers" title="Maintenance">🚧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/mzavolan"><img src="https://avatars.githubusercontent.com/u/44894507?v=4?s=100" width="100px;" alt=""/><br /><sub><b>mzavolan</b></sub></a><br /><a href="#data-mzavolan" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mzavolan" title="Documentation">📖</a> <a href="#eventOrganizing-mzavolan" title="Event Organizing">📋</a> <a href="#financial-mzavolan" title="Financial">💵</a> <a href="#ideas-mzavolan" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-mzavolan" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-mzavolan" title="Project Management">📆</a> <a href="#question-mzavolan" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Amzavolan" title="Reviewed Pull Requests">👀</a> <a href="#talk-mzavolan" title="Talks">📢</a></td>
    <td align="center"><a href="https://github.com/mfansler"><img src="https://avatars.githubusercontent.com/u/1182216?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mervin Fansler</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Amfansler" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mfansler" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mfansler" title="Documentation">📖</a> <a href="#eventOrganizing-mfansler" title="Event Organizing">📋</a> <a href="#ideas-mfansler" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-mfansler" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-mfansler" title="Project Management">📆</a> <a href="#question-mfansler" title="Answering Questions">💬</a> <a href="https://github.com/iRNA-COSI/APAeval/pulls?q=is%3Apr+reviewed-by%3Amfansler" title="Reviewed Pull Requests">👀</a></td>
    <td align="center"><a href="https://github.com/mkatsanto"><img src="https://avatars.githubusercontent.com/u/31883096?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Maria Katsantoni</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=mkatsanto" title="Code">💻</a> <a href="#ideas-mkatsanto" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-mkatsanto" title="Mentoring">🧑‍🏫</a> <a href="#question-mkatsanto" title="Answering Questions">💬</a></td>
    <td align="center"><a href="https://github.com/daneckaw"><img src="https://avatars.githubusercontent.com/u/30384499?v=4?s=100" width="100px;" alt=""/><br /><sub><b>daneckaw</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=daneckaw" title="Code">💻</a> <a href="#data-daneckaw" title="Data">🔣</a> <a href="#eventOrganizing-daneckaw" title="Event Organizing">📋</a> <a href="#ideas-daneckaw" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-daneckaw" title="Project Management">📆</a> <a href="#tutorial-daneckaw" title="Tutorials">✅</a></td>
    <td align="center"><a href="https://github.com/dominikburri"><img src="https://avatars.githubusercontent.com/u/7873536?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dominik Burri</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/issues?q=author%3Adominikburri" title="Bug reports">🐛</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=dominikburri" title="Code">💻</a> <a href="#data-dominikburri" title="Data">🔣</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=dominikburri" title="Documentation">📖</a> <a href="#example-dominikburri" title="Examples">💡</a> <a href="#eventOrganizing-dominikburri" title="Event Organizing">📋</a> <a href="#ideas-dominikburri" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-dominikburri" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#mentoring-dominikburri" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-dominikburri" title="Project Management">📆</a> <a href="#question-dominikburri" title="Answering Questions">💬</a></td>
    <td align="center"><a href="https://github.com/mrgazzara"><img src="https://avatars.githubusercontent.com/u/58347523?v=4?s=100" width="100px;" alt=""/><br /><sub><b>mrgazzara</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=mrgazzara" title="Code">💻</a> <a href="https://github.com/iRNA-COSI/APAeval/commits?author=mrgazzara" title="Documentation">📖</a> <a href="#data-mrgazzara" title="Data">🔣</a> <a href="#eventOrganizing-mrgazzara" title="Event Organizing">📋</a> <a href="#ideas-mrgazzara" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-mrgazzara" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-mrgazzara" title="Maintenance">🚧</a> <a href="#projectManagement-mrgazzara" title="Project Management">📆</a> <a href="#mentoring-mrgazzara" title="Mentoring">🧑‍🏫</a> <a href="#talk-mrgazzara" title="Talks">📢</a></td>
    <td align="center"><a href="https://github.com/FitzsimmonsCM"><img src="https://avatars.githubusercontent.com/u/33811247?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Christina Fitzsimmons</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=FitzsimmonsCM" title="Documentation">📖</a> <a href="#eventOrganizing-FitzsimmonsCM" title="Event Organizing">📋</a> <a href="#ideas-FitzsimmonsCM" title="Ideas, Planning, & Feedback">🤔</a> <a href="#projectManagement-FitzsimmonsCM" title="Project Management">📆</a> <a href="#talk-FitzsimmonsCM" title="Talks">📢</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/lschaerfen"><img src="https://avatars.githubusercontent.com/u/76001100?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Leo Schärfen</b></sub></a><br /><a href="https://github.com/iRNA-COSI/APAeval/commits?author=lschaerfen" title="Code">💻</a> <a href="#ideas-lschaerfen" title="Ideas, Planning, & Feedback">🤔</a> <a href="#talk-lschaerfen" title="Talks">📢</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
