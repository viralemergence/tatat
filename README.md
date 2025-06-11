<img src="https://github.com/user-attachments/assets/c743b2a9-bb2c-4fac-802e-0809aa08ebfd"><br/>
# Transcriptome Assembly, Thinning, and Annotation Tool
TATAT is a Dockerized tool for generating an annotated transcriptome from raw RNA-seq data. It leverages several open-source bioinformatic tools, that are well established for transcriptome generation, and provides python and bash scripts to coordinate these tools and parse their output. Consequently, it can be used to create reproducible workflows that run in high-performance computing (HPC) environments, reducing transcriptome generation runtime from days to hours. It also has been shown to work well with non-model organisms.
<br><br>
### Acquiring TATAT
The Docker image for TATAT is available at: url_address. Once downloaded it is immediately ready to use with Docker. However, it can be used by Singularity as well if converted to a ".sif", which is done automatically if downloaded via the "singularity pull" command:
```
singularity pull tatat.sif docker://alexanderbrown1313/transcriptome-assembly-thinning-annotation-tool:latest
```
If a user wishes to generate the Docker image locally, the github repo should be cloned and within the repo folder the following command run:
```
sudo docker build --platform=linux/amd64 -t tatat .
```
However, it is recommended that the publically available TATAT Docker image be used. Lastly, there is no versioning of the TATAT image, though this may be implemented later.
<br><br>
### Hardware Requirements
*De novo* assembly remains a resource-intensive computation. We observed that to complete the rousettus assemblies described in the [rousettus_tutorial](rousettus_tutorial) folder, 10 CPUs and 50 Gb RAM were required for the most demanding assembly generation. Similarly, the thinning step required 10 CPUs and 60 Gb RAM, and the annotation step required 20 CPUs and 60 Gb RAM. The RAM values listed are definitively required, but the CPUs are not; they just make the overall runtime bearable.
<br><br>
Furthermore, each of these steps generally requires hours to run, with the assemblies taking the longest. By running the rousettus assemblies in parallel, and using 100 CPUs and 500 Gb RAM on a HPC, they were able to finish in ~3.5 hours. To run sequentially could become at least 1 day of runtime.
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

### Getting Started
A detailed tutorial is available in the [rousettus_tutorial](rousettus_tutorial) folder and new users are encouraged to work through it. However, some general notes are included here:
- All the bash scripts provided require the .env file to be filled out correctly. If you wish to use them, fill out the .env file before running anything, otherwise the scripts will throw errors.
- The annotation step requires a BLAST database. The vertebrate database used in the initial publication is available at: url_address. However, if a different or updated database is desired, follow the steps outlined in the rousettus_tutorial.
- The default assembly tools include fastp and rnaSPAdes, and the default thinning tool is EvidentialGenes. TATAT has been designed to allow other tools to be used at these steps and then fed into the workflow. However, currently no other tools have been tested and we cannot guarantee the workflow will function as expected. Therefore, users are encouraged to use the bioinformatic tools already provided in the TATAT Docker image.
<br><br>
### Limitations
- Currently TATAT is not designed to capture all isoforms in the final transcriptome, and is instead designed to select the longest potential isoform to represent the gene. This means that Differential Gene Expression (DGE) analysis using TATAT transcriptomes is able to show overall changes in gene expression, but not per isoform.
- TATAT is primarily designed to generate coding transcriptomes and has been optimized and validated on coding sequences. TATAT has the capability to annotate non-coding RNAs (ncRNA) and add them to the transciptomes generated, but this feature takes significantly longer (up to 5-10x longer!) and has not been extensively validated.
<br><br>
### Citing
If you use TATAT in your work, please cite the publication listed below:
<br><br>
*Publication still pending*
