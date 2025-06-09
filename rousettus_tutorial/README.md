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

### TATAT "Core" Coding Genes: Preparation
It should be noted that TATAT works best either in HPC or cloud based environments. Unless your local machine meets the [hardware requirements](https://github.com/viralemergence/tatat/blob/readme/README.md#hardware-requirements), do not attempt to run locally. Once the Docker image (or sif) has been downloaded or generated (see [Acquiring TATAT](https://github.com/viralemergence/tatat/blob/readme/README.md#acquiring-tatat)), it can be run via bash scripts or on the command line. We will briefly discuss some of the additional preparation necessary to run TATAT as follows:

#### The .env File
Since TATAT runs in a virtual container, it is necessary to explicitly grant access to directories and files on the host system. When running TATAT via Docker or Singularity, those file paths would have to be typed out repeatedly in the commands. In order to simplify providing that information, the tutorial uses an .env file which stores all the directory and file paths as environmental variables. An example file can be found at [.env-example](../.env-example). Likewise, there are some additional variables that are convenient for storing potentially sensitive information, like the NCBI API key. E.g.:
```
SINGULARITY_IMAGE="/path/to/image.sif"
APP_DIR="/path/to/dir"
DATA_DIR="/path/to/dir"
SQLITE_DB_DIR=$DATA_DIR"/path/to/dir"
TRANSCRIPTOME_DATA_DIR=$DATA_DIR"/path/to/dir"
```
Take the time to fill out these variables with the appropriate information. TATAT can be used without this feature, but for this tutorial it will be expected that
this file is completed. Even if the directories don't exist, the bash scripts in the tutorial are written to generate most of the directories, with the information provided. The main exception to this are for the SINGULARITY_IMAGE, APP_DIR, DATA_DIR, SQLITE_DB_DIR, TRANSCRIPTOME_DATA_DIR, and SCRATCH_DIR. These must be generated manually.

#### Job Scheduling
Many of the stages in TATAT run for hours and attempting to run the code with "live" sessions via tools like ssh is not generally recommended, as loss of the connection may terminate the process or make it difficult to tell when the process has completed. Likewise, generating assemblies sequentially via a ssh connection could take days, whereas we observed with the rousettus data performing assembly in parallel completed in ~3.5 hours.
<br><br>
Consequently, we strongly recommend using a job scheduler to launch the TATAT bash scripts and use the job scheduler SLURM in this tutorial. We will briefly cover aspects of SLURM, but an in-depth discussion is beyond the scope of this tutorial, and may be irrelevant if a different job scheduler is used. However, if a user has limited resources or feels comfortable running TATAT just through the command line, it is not necessary to use a job scheduler.

#### The BLAST Database
The final annotation stage of TATAT requires a BLAST database. For this tutorial the NCBI's "core_nt" database was downloaded and all vertebrata sequences extracted and built into a new BLAST database. The user may follow the code in [make_vertebrata_core_nt_blastdb.sh](make_vertebrata_core_nt_blastdb.sh) to generate their own version, however, our version is available at url_address. If a user wishes to generate their own version for this tutorial, note that the NCBI constantly updates their database and so final results may differ slightly. If using the .env file, the path to this database will need to be stored in:
```
TATAT_BLASTDB_DIR="/path/to/dir"
```
Lastly, outside of this tutorial the same methodology could be used to generate BLAST databases for different clades, organisms, sequences, etc.

### TATAT "Core" Coding Genes: Running The Whole Thing
If the .env has been correctly filled out, the BLAST database is prepared, the job scheduler used is SLURM, and all the bash scripts have been correctly modified to use whichever partition is avialble to the user, at this point the entire TATAT workflow can be run using [tatat.sh](tatat.sh). However, it is likely that there will be some minor errors. Consequently, for the purposes of this tutorial, we will walk through the individual stages of TATAT, which comprise of Assembly, Thinning, and Annotation.
