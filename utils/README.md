# APAeval utils

This directory provides Dockerfiles and corresponding 'recipes' for scripting/processing tasks that are common across workflows.

Contents:
- [Available utilities](#Available-utilities)
- [Contributing instructions](#Contributing-instructions)

# Available utilities

## Convert GTF to BED12 (aka 'gene model') file

**Description:** Convert a [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) annotation file to [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (also known as a 'gene model' file) format.

**subdirectory name**: **{update this}**

**Dockerhub link:** https://hub.docker.com/r/apaeval/gtftobed

**Dockerhub url (for pipeline)**: `docker.io/apaeval/gtftobed:1.0`

**Usage (shell command to run the tool)**:
```
python3 /app/gtfTobed.py --gtf <path to GTF file>
```


# Contributing instructions

Have code for a pre-processing / file wrangling task that you think may be general to other tools? To add a tool to the directory, we need you to create a pull request with the following tasks completed:

- Create a suitably named subdirectory under `utils` to house the subsequent files (e.g. gtf_to_bed12)
- Add a Dockerfile to the subdirectory to provide necessary dependencies and scripts. This should then be added to the `apaeval` team Dockerhub (if you don't have access, you can request as part of the pull request)
- Add description of task and how to run it to the README. This should follow the format below:


\## Header description of tool (placed under Available utilities)

(remove the escape slash when you copy and paste to interpret the hashes as titles)

Description: short 1/2 sentence description of the task the utility achieves

Dockerhub link: `https://hub.docker.com/r/apaeval/<tool_name>`

Docker.io url (to use in pipeline): <docker.io/apaeval/<tool_name>:<tag>

Usage (shell command to run the tool):

```
python utility.py -i <input> -o <output>
```
