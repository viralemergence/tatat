# TATAT Rousettus Tutorial
The following is a tutorial for using TATAT to process raw RNA-seq data from the *Rousettus aegyptiacus* bat and make an annotated transcriptome.
The data includes 20 samples from 11 tissue types, with half the samples from an adult male bat and the other half from an adult female bat.
The details can be found in the [sample_metadata.csv](rousettus_tutorial/sample_metadata.csv) file.

### A Brief Overview
The main function of TATAT is to generate a transcriptome of coding genes, and so that will be covered first.
However, TATAT also has the ability to add non-coding RNA (ncRNA) genes to the transcriptome. It requires that the coding transcriptome already
be generated, and we suspect will generally be of less interest to users, so that will be covered second. For coding genes, the workflow can be depicted by
the following diagram:
<br>
<p align=center>
  <img src="https://github.com/user-attachments/assets/602fcf0e-0a07-4572-babb-d311a5a47711" width=60%>
</p>

In this tutorial, TATAT will first be used to download the sequencing data from the Sequencing Read Archive (SRA), remove sequencing errors and adaptors using FASTP, assemble the reads into contigs via *de novo* assembly with rnaSPAdes, thin out the excess contigs, then annotate the remaining contigs and select which best represent each coding gene. Then the initial set of contigs will be screened for matches to the coding genes, and the remaining contigs annotated for ncRNA genes.

### TATAT "Core" Coding Genes
It should be noted that TATAT works best either in HPC or cloud based environments. Unless your local machine meets the [hardware requirements](https://github.com/viralemergence/tatat/blob/readme/README.md#hardware-requirements), do not attempt to run locally. Once the Docker image (or sif) has been downloaded or generated (see [Acquiring TATAT](https://github.com/viralemergence/tatat/blob/readme/README.md#acquiring-tatat)), it can be run via bash scripts or on the command line. We will briefly discuss some of the preparation necessary to run TATAT:

#### The .env File
Since TATAT runs in a virtual container, it is necessary to explicitly grant access to directories and files on the host system. When running TATAT via Docker or Singularity, those file paths would have to be typed out repeatedly in the commands. In order to simplify providing that information, the tutorial uses an .env file which stores all the directory and file paths as environmental variables.
