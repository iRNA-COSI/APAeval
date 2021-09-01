# APAeval utils directory

This directory provides Dockerfiles and 'code recipes' for scripting/processing tasks that are common across workflows. This resource should simplify development for incoming execution workflow writers and ease  
review by using universal code.

Contents:
- Available utilities
- Contributing instructions

# Available utilities

## Convert GTF to BED12 (aka 'gene model') file

Description: Convert a [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) annotation file to [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (also known as a 'gene model' file) format.

Dockerhub link: https://hub.docker.com/r/apaeval/gtftobed

Dockerhub url (for pipeline): `docker.io/apaeval/gtftobed:1.0`

Usage (shell command to run the tool):

```
python3 /app/gtfTobed.py --gtf <path to GTF file>
```



# Contributing instructions

Have code for a pre-processing / file wrangling task that you think may be general to other tools? To add a tool to the directory, we need the following:

- Create a suitably named subdirectory under `utils` to house the subsequent files
- Dockerfile to provide necessary package requirements. This should be added to the `apaeval` team Dockerhub
- Add Description of task and how to run it to the README. This should follow the format below:


\## Header description of tool (placed under Available utilities, remove the escape slash)

Description: short 1/2 sentence description of the task the utility achieves

Dockerhub link: `https://hub.docker.com/r/apaeval/<tool_name>`

Docker.io url (to use in pipeline): <docker.io/apaeval/<tool_name>:<tag>

Usage (shell command to run the tool):

```
python utility.py -i <input> -o <output>
```
