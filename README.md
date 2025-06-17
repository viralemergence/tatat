<img src="https://github.com/user-attachments/assets/cebd70c4-f279-4f25-b1b0-6f83cf9a0b35" width=50%><br/>
# Transcriptome Assembly, Thinning, and Annotation Tool
TATAT is a Dockerized tool for generating an annotated transcriptome from raw RNA-seq data. It leverages several open-source bioinformatic tools, that are well established for transcriptome generation, and provides python and bash scripts to coordinate these tools and parse their output. Consequently, it can be used to create reproducible workflows that run in high-performance computing (HPC) environments, reducing transcriptome generation runtime from days to hours. It also has been shown to work well with non-model organisms.
<br><br>
### Acquiring TATAT
The Docker image for TATAT is available here in the [GitHub registry](https://github.com/viralemergence/tatat/pkgs/container/transcriptome-assembly-thinning-annotation-tool). Once downloaded it is immediately ready to use with Docker. However, it can be used by Singularity as well if converted to a ".sif", which is done automatically if downloaded via the "singularity pull" command:
```
singularity pull tatat.sif docker://ghcr.io/viralemergence/transcriptome-assembly-thinning-annotation-tool:latest
```
If a user wishes to generate the Docker image locally, the github repo should be cloned and within the repo folder the following command run:
```
sudo docker build --platform=linux/amd64 -t tatat .
```
However, it is recommended that the publically available TATAT Docker image be used. Lastly, there is no versioning of the TATAT image, though this may be implemented later.
<br><br>
### Hardware Requirements
*De novo* assembly remains a resource-intensive computation. We observed that to complete the rousettus assemblies described in the [rousettus_tutorial](rousettus_tutorial) folder, 10 CPUs and 50 GB RAM were required for the most demanding single assembly generation. Similarly, the thinning step required 10 CPUs and 60 GB RAM, and the annotation step required 20 CPUs and 60 GB RAM. The RAM values listed are definitively required, but the CPUs are not; they just make the overall runtime bearable.
<br><br>
Furthermore, each of these steps generally requires hours to run, with the assemblies taking the longest. By running the rousettus assemblies in parallel, and using 100 CPUs and 500 GB RAM on a HPC, they were able to finish in ~3.5 hours. To run sequentially could become at least 1 day of runtime.
<br><br>
Consequently, while it is theoretically possible to run TATAT locally on a computer, it is not recommended. It would be better to run on a HPC, or cloud-based environment like AWS or GCP.
<br><br>
### Workflow Overview
The following diagram shows an overview of the intended workflow for TATAT and which tools are used to facilitate each step. <br>
NOTE: It is not required to pull sequencing data from the SRA; it is a feature that can be used.
<br><br>
<p align=center>
  <img src="https://github.com/user-attachments/assets/602fcf0e-0a07-4572-babb-d311a5a47711" width=60%>
</p>

### Benchmarking
The 20 samples detailed in the [rousettus_tutorial](rousettus_tutorial) folder totalled up to 690 GB of uncompressed sequencing data. The assemblies were processed in parallel (as described above), but all other stages were run as single, unified jobs. The entire TATAT coding gene workflow was run for the tutorial 3 times and the median taken from the runs to determine overall runtime and runtime per stage, using the same resources. Times were rounded to the nearest 10 minutes.

| Stage | Runtime | CPUs | RAM | Parallel |
| :--- | :---: | :---: | :---: | :---: |
| Overall | 8 hr 10 min | | | |
| Assembly | 3 hr 40 min | 10 | 50 GB | 10 jobs |
| Thinning | 1 hr 50 min | 10 | 60 GB | No |
| Annotation | 2 hr 30 min | 20 | 60 GB | No |

In summary, 20 non-model organism samples of raw RNA-seq data, totalling 690 GB were converted to a comprehensive coding transcriptome in ~8 hours (a single workday). The assembly stage required the most resources at 100 CPUs and 500 GB of RAM, but this step could be run sequentially rather than in parallel, and each assembly then only requires 10 CPUs and 50 GB RAM.

Notes:
- The individual stage runtimes do not perfectly add up to the "Overall" value because there was some additional time used for requesting and waiting on resource allocation between stages, and the rounding to 10 minutes per stage slightly under reports the individual stage runtimes.
- To reiterate, the resource usage for Assembly is *per* job, but 10 jobs were run in parallel, totalling 100 CPUs and 500 GB RAM.
- For some of these stages the full RAM requested was not used, but often was close. E.g. Thinning in one run used 58 Gb of RAM, 50 Gb in another, but never the full 60 GB. The extra RAM requested allows wiggle room to prevent an Out Of Memory (OOM) error.
- This was just for the coding genes; the non-coding genes took much longer and are described in the tutorial README.md.
<br><br>
### Getting Started
A detailed tutorial is available in the [rousettus_tutorial](rousettus_tutorial) folder and new users are encouraged to work through it. However, some general notes are included here:
- All the bash scripts provided require the .env file to be filled out correctly. If you wish to use them, fill out the .env file before running anything, otherwise the scripts will throw errors.
- The annotation step requires a BLAST database. The vertebrate database used in the initial publication is available at: [vertebrata core nt BLAST db](https://zenodo.org/records/15685806). However, if a different or updated database is desired, follow the steps outlined in the rousettus tutorial.
- The default assembly tools include FASTP and rnaSPAdes, and the default thinning tool is EvidentialGene. TATAT has been designed to allow other tools to be used at these steps and then fed into the workflow. However, currently no other tools have been tested and we cannot guarantee the workflow will function as expected. Therefore, users are encouraged to use the bioinformatic tools already provided in the TATAT Docker image.
<br><br>
### Limitations
- Currently TATAT is not designed to capture all isoforms in the final transcriptome, and is instead designed to select the longest potential isoform to represent the gene. This means that analyses, such as Differential Gene Expression (DGE) analysis, using TATAT transcriptomes are able to show overall changes in gene expression, but not per isoform.
- TATAT is primarily designed to generate coding transcriptomes and has been optimized and validated on coding sequences. TATAT has the capability to annotate non-coding RNAs (ncRNA) and add them to the transciptomes generated, but this feature takes significantly longer (up to 5-10x longer!) and has not been extensively validated.
<br><br>
### Citing
If you use TATAT in your work, please cite the publication listed below:
<br><br>
*Publication pending*
