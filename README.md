<img src="https://github.com/user-attachments/assets/c743b2a9-bb2c-4fac-802e-0809aa08ebfd"><br/>
# Transcriptome Assembly, Thinning, and Annotation Tool
TATAT is a Dockerized tool for generating an annotated transcriptome from raw RNA-seq data. It leverages several open-source bioinformatic tools, that are well established for transcriptome generation, and provides python and bash scripts to coordinate these tools and parse their output. Consequently, it can be used to create reproducible workflows that run in high-performance computing (HPC) environments, reducing transcriptome generation runtime from days to hours. It also has been shown to work well with non-model organisms.
<br><br>
### Acquiring TATAT
The Docker image for TATAT is available at: url.address. Once downloaded it is immediately ready to use with Docker. However, it can be used by Singularity as well if converted to a ".sif", which is done automatically if downloaded via the "singularity pull" command:
```
singularity pull tatat.sif docker://alexanderbrown1313/transcriptome-assembly-thinning-annotation-tool:latest
```
If a user wishes to generate the Docker image locally, the github repo should be cloned and within the repo folder the following command run:
```
sudo docker build --platform=linux/amd64 -t tatat .
```
However, it is recommended that the publically available TATAT Docker image be used. Lastly, there is no versioning of the TATAT image, though this may be implemented later.
<br><br><br><br><br>
Additional thoughts:
Consequently, the entire workflow may be run with a single bash script. However, the processes are also decoupled, so that a user may choose to only run part of the workflow or introduce their own preferenced tools at different steps (e.g. using Trinity instead of rnaSPAdes for de novo assembly).
