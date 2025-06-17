# TATAT Rousettus Tutorial
The following is a tutorial for using TATAT to process raw RNA-seq data from the *Rousettus aegyptiacus* bat and make an annotated transcriptome.
The data includes 20 samples from 11 tissue types, with half the samples from an adult male bat and the other half from an adult female bat.
The details can be found in the [sample_metadata.csv](sample_metadata.csv) file.

### A Brief Overview
The main function of TATAT is to generate a transcriptome of coding genes, and so that will be covered first.
However, TATAT also has the ability to add non-coding RNA (ncRNA) genes to the transcriptome. It requires that the coding transcriptome already
be generated, and we suspect will generally be of less interest to users, so that will be covered second. For coding genes, the workflow can be depicted by
the following diagram:
<br>
<p align=center>
  <img src="https://github.com/user-attachments/assets/602fcf0e-0a07-4572-babb-d311a5a47711" width=60%>
</p>

In this tutorial, TATAT will first be used to download the sequencing data from the Sequencing Read Archive (SRA), remove sequencing errors and adapters using FASTP, assemble the reads into contigs via *de novo* assembly with rnaSPAdes, thin out the excess contigs, then annotate the remaining contigs and select which best represent each coding gene.

### TATAT Coding Genes: Preparation
It should be noted that TATAT works best either in HPC or cloud based environments. Unless your local machine meets the [hardware requirements](https://github.com/viralemergence/tatat/blob/readme/README.md#hardware-requirements), do not attempt to run locally. Once the Docker image (or sif) has been downloaded or generated (see [Acquiring TATAT](https://github.com/viralemergence/tatat/blob/readme/README.md#acquiring-tatat)), it can be run via bash scripts or on the command line. We will briefly discuss some of the additional preparation necessary to run TATAT as follows:

#### The .env File
Since TATAT runs in a virtual container, it is necessary to explicitly grant access to directories and files on the host system. When running TATAT via Docker or Singularity, those file paths would have to be typed out repeatedly in the commands. In order to simplify providing that information, the tutorial uses an .env file which stores all the directory and file paths as environmental variables. An example file can be found at [.env-example](../.env-example). Likewise, there are some additional variables that are convenient for storing potentially sensitive information, like the NCBI API key. E.g.:
```
NCBI_API_KEY=ncbikey
SINGULARITY_IMAGE="/path/to/image.sif"
DATA_DIR="/path/to/dir"
SQLITE_DB_DIR=$DATA_DIR"/path/to/dir"
TRANSCRIPTOME_DATA_DIR=$DATA_DIR"/path/to/dir"
```
Take the time to fill out these variables with the appropriate information. TATAT can be used without this feature, but for this tutorial it will be expected that
this file is completed. Even if the directories don't exist, the bash scripts in the tutorial are written to generate the directories with the information provided.

#### Job Scheduling
Many of the stages in TATAT run for hours and attempting to run the code with "live" sessions via tools like ssh is not generally recommended, as loss of the connection may terminate the process or make it difficult to tell when the process has completed. Likewise, generating assemblies sequentially via a ssh connection could take days, whereas we observed with the rousettus data performing assembly in parallel completed in ~3.5 hours.
<br><br>
Consequently, we strongly recommend using a job scheduler to launch the TATAT bash scripts and we use the job scheduler SLURM in this tutorial. We will briefly cover aspects of SLURM, but an in-depth discussion is beyond the scope of this tutorial, and may be irrelevant if a different job scheduler is used. However, if a user has limited resources or feels comfortable running TATAT just through the command line, it is not necessary to use a job scheduler.

#### The BLAST Database
The final annotation stage of TATAT requires a BLAST database. For this tutorial the NCBI's "core_nt" database was downloaded and all vertebrata sequences extracted and built into a new BLAST database. The user may follow the code in [make_vertebrata_core_nt_blastdb.sh](make_vertebrata_core_nt_blastdb.sh) to generate their own version, however, our version is available at [vertebrata core nt BLAST db](https://zenodo.org/records/15685806). If a user wishes to generate their own version for this tutorial, note that the NCBI constantly updates their database and so final results may differ slightly. If using the .env file, the path to this database will need to be stored in:
```
TATAT_BLASTDB_DIR="/path/to/dir"
```
Lastly, outside of this tutorial the same methodology could be used to generate BLAST databases for different clades, organisms, sequences, etc.

#### Docker vs Singularity
It is this software engineer's opinion that Docker is better than Singularity, and running TATAT via Docker would be more ideal. However, in many HPC environments Singularity is preferentially chosen as it does not require root permissions to run. Consequently, for this tutorial we use Singularity.

#### Sample Metadata
Lastly, TATAT requires a metadata file on the samples to be analyzed. The data is stored in a sqlite database and subsequently expedites performing assemblies in parallel, running TATAT on different transcriptomes (e.g. species) in parallel, and aids in certain other functions such as looking for consensus sequences during contig thinning. An example file is [sample_metadata.csv](sample_metadata.csv). Any number of columns may be included, and they will be imported into the database, but the following columns are required: uid, transcriptome, r1_reads, r2_reads.
- uid: A Unique IDentifier. For SRA data, it must be the SRA number. For personal samples, it may be any string of alphanumeric values without spaces, as long as no other sample has the same uid (otherwise it is not unique).
- transcriptome: The transcriptome to which the sample belongs. This could be a species, or age group, or a random string. But the final transcriptome(s) will be generated with this information.
- r1_reads & r2_reads: The paths for the R1 and R2 reads. For the tutorial and SRA data, this does not initially need to be supplied, as TATAT will autopopulate the R1 and R2 read paths when downloading the SRA samples. However, if personal samples are being used, this information must be provided in the sample_metadata.csv file.

### TATAT Coding Genes: Running The Whole Thing
If the .env has been correctly filled out, the BLAST database is prepared, the job scheduler used is SLURM, and all the bash scripts have been correctly modified to use whichever partition is available to the user, at this point the entire TATAT workflow can be run using [tatat.sh](tatat.sh). The command would look like:
```
sbatch tatat.sh path/to/.env
```
However, it is likely that there will be some minor errors. Consequently, for the purposes of this tutorial, we will generate the initial sqlite database, then walk through the individual stages of TATAT, which comprise of Assembly, Thinning, and Annotation.

#### Sqlite Database
TATAT uses a sqlite database to store metadata about the samples, contigs, and annotations. To generate it, the code in [tatat.sh](tatat.sh) can be run manually, or the bash script could be run with subsequent lines deleted or commented out. Running through the command line would be:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -sample_metadata /src/sqlite_db/sample_metadata.csv
```
Once generated, TATAT will use the same database for all the following stages.
<br><br>
As a brief aside on Singularity and the args used here:
- "exec" creates a container from the image, runs whatever code is in the command, then shuts the container down.
- "--pwd /src" sets the working directory of the container to "/src".
- "--no-home" keeps other directories on the host system from secretly being mounted to the container.
- "--bind" gives the container access to directories on the host system, using the format of "--bind host_directory:container_directory".

### TATAT Coding Genes: Assembly
The [assembly_array.sh](tatat_main/assembly_array.sh) script contains all the steps necessary to generate *de novo* assemblies from the samples. The script itself can be submitted to slurm and will request 20 jobs be processed, with a maximum of only 10 at a time:
```
#SBATCH --array=0-19:1%10
```

The first part of the script uses the array number generated by slurm to retrieve the sample uid, which for this tutorial is the SRA number:
```
SRA_NUMBER=$(singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sample_metadata_extraction.py \
    -sqlite_db /src/sqlite_db/tatat.db \
    -array_index $SLURM_ARRAY_TASK_ID \
    -return_uid)
```

The SRA number is then used to download the sequencing reads, and update the R1 and R2 paths in the sqlite database:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SRA_DOWNLOAD_DIR:/src/data/download_dir \
    --bind $SRA_COLLATE_DIR:/src/data/collate_dir \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/sra_read_download.py \
    -sra_number $SRA_NUMBER \
    -download_dir /src/data/download_dir \
    -collate_dir /src/data/collate_dir \
    -sqlite_db /src/sqlite_db/tatat.db
```

Once downloaded, the reads are processed with FASTP, removing low quality reads, fixing sequences errors when possible, and removing the adapter sequences:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SRA_COLLATE_DIR:/src/data/fastq_dir \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/fastp_orchestration.py \
    -fastq_dir /src/data/fastq_dir \
    -outdir /src/data/trimmed \
    -sqlite_db /src/sqlite_db/tatat.db \
    -uid $SRA_NUMBER -r1_adapter $R1_ADAPTER -r2_adapter $R2_ADAPTER -cpus $SLURM_CPUS_PER_TASK
```
It is worth noting that the code is written to automatically detect if the reads are paired. Some of these samples do not have paired reads. When this occurs, the R2 adapter sequence is not used and the FASTP command is altered under the hood.

Finally, the trimmed reads are passed to rnaSPAdes to generate contigs:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $FASTQ_TRIMMED_DIR:/src/data/trimmed \
    --bind $RNASPADES_ASSEMBLY_DIR:/src/data/assembly \
    --bind $RNASPADES_COLLATED_ASSEMBLY_DIR:/src/data/collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/assembly/rnaspades_orchestration.py \
    -fastq_dir /src/data/trimmed \
    -assembly_dir /src/data/assembly \
    -collated_dir /src/data/collated \
    -unique_identifier $SRA_NUMBER -cpus $SLURM_CPUS_PER_TASK -memory $MEM_IN_GB
```

Each of these steps is purposely designed to be carried out separately by a specific script so that a user has the option of using their own tools as desired. E.g. Cutadapt could be used instead of FASTP and Trinity could be used instead of rnaSPAdes. Hypothetically this should cause no issues with the rest of TATAT. However, this has not currently been tested and it is highly recommended to use the default tools, especially as they have been shown in many publications to be fast and require less resources than other tools. Regardless, once these commands have finished, all the final contigs should be in the folder indicated by $RNASPADES_COLLATED_ASSEMBLY_DIR.

**Troubleshooting:** If the assembly code was launched via slurm and failed, check the output error log. If run on the command line, errors will be printed to the terminal. While any number of errors could potentially be generated, the most common include:
- file/directory does not exist: This usually results from incorrectly assigning path variables in the .env file, and needs to be corrected.
- OOM error: This stands for Out Of Memory and should not occur for the tutorial, but if the RAM amount requested was lowered or the command was run locally on a machine lacking sufficient RAM, this error could occur.
- directory already exists: This usually occurs when re-running the script and is not a problem. The mkdir commands create this error if the directory has already been generated.

### TATAT Coding Genes: Thinning
After the assembly stage has completed, the thinning stage can be run with [thinning.sh](tatat_main/thinning.sh).
<br><br>
The first step generates a sqlite table for the *de novo* assemblies called "transcripts" and a table for candidate CDS called "cds":
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_transcripts_table \
    -create_cds_table
```
If for any reason the thinning stage needs to be redone, it is important to make sure this step is also re-run, as it will delete the old tables and generate new, empty tables.

Then all the contigs need to be merged into a single file and their metadata added to the "transcripts" table:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $RNASPADES_COLLATED_ASSEMBLY_DIR:/src/data/collated \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/merge_fastas_and_set_metadata.py \
    -assembly_fasta_dir /src/data/collated \
    -merged_path /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db
```
The single file, called "raw_transcriptome.fna" here, is a fasta file that contains all the *de novo* assemblies generated previously. This file is used repeatedly later on, so make sure not to delete it. Also, because this file tends to be large, it was decided it was better to leave it out of the sqlite database. However, information such as the sequence ids, sequence lengths, and which sample they come from is all stored in the sqlite "transcripts" table.

Querying either the sqlite database or raw_transcriptome.fna will reveal at this point ~14.5 million transcripts have been generated, which far exceed the expected number of genes. This is likely because these potential candidates are not all biologically relevant, as outlined in the diagram below:
<br>
<p align=center>
  <img src="https://github.com/user-attachments/assets/850ef67e-5531-43d6-b3d0-ec572caa5c0f" width=70%>
</p>

The transcripts could include assemblies generated from DNA contamination, degraded RNA, transcriptional noise, chimeric assemblies, and other anomalies, in addition to real mRNA transcripts and their isoforms. Also, when many samples are included for the same transcriptome there may be redundancies (e.g. we observed an apparent transcript for Ubiquitin generated separately for almost every sample, leading to essentially the same transcript being reported 14 times!). Consequently, it becomes necessary to "thin" these candidates by removing suspected anomalies and redundancies. TATAT accomplishes this by using EvidentialGene in the following command:
```
singularity exec \
    --pwd /src \
    --no-home \
    --env LC_ALL=C \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    --bind $EVIGENE_OUTPUT_DIR:/src/evigene_output \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/evigene_orchestration.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -outdir /src/evigene_output \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcriptome rousettus -prefix_column sample_uid \
    -run_evigene -cpus $SLURM_CPUS_PER_TASK -mem $SLURM_MEM_PER_NODE -phetero 2 -minaa 99 \
    -run_transcript_metadata_appender -run_cds_and_metadata
```
Additionally, for this tutorial we pass the "phetero" arg, as we expect there to be some sequence discrepencies due to heterozygosity in the samples, and "minaa", as mammals tend to have longer genes and this removes genes with fewer than 100 amino acids. For more details on optimizing these args with other organisms, see the EvidentialGene [homepage](http://arthropods.eugenes.org/EvidentialGene/evigene/).

This step generally takes a couple hours to run, but once completed it will have populated the "cds" table with candidate CDS ids, start and end positions derived from the raw transcripts, strand directionality, the parental transcript id, and other information. However, ideally the "transcripts" table entries will have direct connections to the "cds" table entries. To quickly add this, the following command is run:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/thinning/evigene_orchestration.py \
    -assembly_fasta holder -outdir holder \
    -sqlite_db /src/sqlite_db/tatat.db \
    -transcriptome holder -prefix_column holder \
    -update_transcript_cds_ids
```
This is much faster and provides each transcript row with the CDS id(s) it corresponds to, if any. This step is also run separately so that the EvidentialGene step can be run in parallel for multiple transcriptomes. The "cds" table now contains all the information necessary to begin the Annotation stage.

**Troubleshooting:** See the Assembly stage troubleshooting section for similar tips.

### TATAT Coding Genes: Annotation
After the thinning stage the [annotation.sh](tatat_main/annotation.sh) script may be run.
<br><br>
We begin the annotation process by performing a BLAST search where the candidate sequences in the "cds" table are the queries and the "vertebrata_core_nt" database (described previously) contains the subject sequences. The query sequences are extracted:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -sql_queries /src/app/example_sql_queries/cds_for_blastn_sql_queries.json \
    -cds_fasta /src/transcriptome_data/cds.fna
```
Then BLASTed:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TATAT_BLASTDB_DIR:/src/blastdb \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    $SINGULARITY_IMAGE \
    blastn -db /src/blastdb/vertebrata_core_nt \
    -query /src/transcriptome_data/cds.fna \
    -out /src/blast_hits/cds_hits.tsv \
    -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -mt_mode 1 \
    -outfmt "6 qseqid sacc qlen" -max_target_seqs 10
```
This step also takes a couple hours, and it is recommended to use as many CPUs as possible, but eventually a diminishing returns effect is observed. 10-20 CPUs seems to yield best results. Also it is recommended to use 5-10 max_target_seqs, as some of the subject sequences in the database do not have gene symbols associated with them in the NCBI database, and subsequent scripts try to identify sequences with gene symbols to use in the final annotation.

After the BLAST search finishes, the resulting file will contain the accession numbers for the subject sequences to which our queries mapped. These accession numbers can be used to obtain gene symbols. First we generate another sqlite table called "accession_numbers":
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/sqlite_db_prep.py \
    -sqlite_db_dir /src/sqlite_db \
    -create_acc_num_table
```
And then use the NCBI's "Datasets" tool to submit the BLAST result's accession numbers in batches:
```
singularity exec \
    --pwd /src \
    --no-home \
    --env NCBI_API_KEY=$NCBI_API_KEY \
    --bind /etc:/etc \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/make_accession_gene_symbol_mapping.py \
    -blast_results /src/blast_hits/cds_hits.tsv \
    -sqlite_db /src/sqlite_db/tatat.db \
    -table_name accession_numbers -rna_type coding
```
Once completed, the "accession_numbers" table will be populated with accession numbers and gene symbols.

Notes:
- The rna_type "coding" is used to skip any genes that may be non-coding.
- Here we bind the "/etc" directory, as this contains credentials necessary for interacting with web services. Depending on the host system's directory structure, this may need to be changed.

**Disclaimer**: This is the most touchy part of TATAT. Unfortunately, since the Datasets tool uses the NCBI servers, sometimes the script crashes if the servers are experiencing high demand, maintenance, database updates, or other factors not fully understood. For instance, it has been observed the NCBI servers seem to reject requests via Datasets around midnight. However, most of the time it runs correctly.
<br><br>
That being said, we also highly recommend obtaining an NCBI API key to use for this part of the script. Members of the scientific community have claimed that frequent requests to the NCBI servers without a key may result in ALL requests being rejected. This is supposedly a protective measure to prevent overloading the servers, either from novice scripters or malicious attacks. Please obtain an NCBI API key to prevent being blocked from the servers.
<br><br>
After this the final annotation runs quickly:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $BLAST_HITS_DIR:/src/blast_hits \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/annotation/assign_gene_annotations_to_cds.py \
    -blast_results /src/blast_hits/cds_hits.tsv \
    -sqlite_db /src/sqlite_db/tatat.db -transcriptome rousettus
```
This script functions by identifying all the CDS that found a hit via BLAST, and then if that CDS had multiple hits, assigning the most significant hit to it with a real gene symbol (e.g. a "LOC" gene with a higher e-score would be skipped for a real gene like UBB, even if the e-score was slightly lower). Then all the CDS with a shared gene symbol are clustered (e.g. all sequences that matched UBB), and the longest sequence is chosen to represent that cluster; the logic here is there may be many isoforms for that gene, and the longest CDS represents the longest isoform. These sequences are flagged in the "cds" table as being "core" genes, i.e. they represent the "core" transcriptome, but other isoforms likely exist, and definitely ncRNA exists.
<br><br>
Finally the "core" coding transcriptome sequences may be extracted:
```
singularity exec \
    --pwd /src \
    --no-home \
    --bind $TRANSCRIPTOME_DATA_DIR:/src/transcriptome_data \
    --bind $SQLITE_DB_DIR:/src/sqlite_db \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/evigene_cds_aa_extraction.py \
    -assembly_fasta /src/transcriptome_data/raw_transcriptome.fna \
    -sqlite_db /src/sqlite_db/tatat.db \
    -sql_queries /src/app/example_sql_queries/core_gene_sql_queries.json \
    -cds_fasta /src/transcriptome_data/rousettus_cds_core.fna \
    -add_gene_name -transcriptome rousettus
```
This will produce a fna file where each sequence has the unique cds_id assigned by TATAT and the gene symbol as the header, and the CDS under the header. This file may be used for subsequent analyses.

**Troubleshooting:** As described previously, the main errors likely to occur here involve the NCBI servers and Datasets tool. To address them, consider:
- Obtaining a NCBI API key.
- Not running this part of the code around midnight.
- Checking if the host system stores its SSL certificates in a directory other than "/etc".

### TATAT Coding Genes: Post QC Summary
Normally, after the annotation stage of TATAT is complete, there are no additional steps. However, in the interest of validating TATAT we generated a number of scripts to ensure the final coding transcriptome from this tutorial/dataset produced reliable sequences for subsequent analysis. We will not cover each individual step, but they can all be found in [post_qc.sh](post_qc/post_qc.sh)
<br><br>
Similarly, an in-depth discussion of the QC results from each script are in the publication. Please refer to the publication if you want to compare your outputs, especially Figure 2.
<br><br>
The key takeaways from our QC analyses were:
- TATAT was able to recover most of the NCBI genes for *Rousettus aegyptiacus*, as well as identify new genes.
- TATAT appears to have some difficulty when replicate count is low, or transcript levels are low (both of which are common problems to transcriptomics).
- The nucleotide sequences produced by TATAT have high sequence similarity to the NCBI genes, but tend to be a little shorter.
- The final transcriptome is suitable for downstream analyses such as Differential Gene Expression (DGE), Gene Set Enrichment (GSE), Gene Ontology Enrichment (GOE), and other analyses that look at changes in the counts of reads mapping to gene sequences.
- However, we make no guarantees as to whether the amino acid sequence is highly accurate or that the longest isoform was correctly identified.

### TATAT Non-Coding Genes: Summary
In addition to generating a coding transcriptome, TATAT is able to generate a non-coding transcriptome as described in [ncrna.sh](ncrna.sh)
<br><br>
The code logic is very similar to what was described in the "TATAT Coding Genes" part of the tutorial, so we will not describe it in great detail here. However, it is dependent on the coding transcriptome for it to work, and there are some other key takeaways:
- Since ncRNA is not well understood, there are not many tools for thinning it. Consequently, the best we could do was remove transcripts with high length, remove any sequences that mapped to the coding transcripts, then cluster the remaining transcripts by sequence identity. This still left ~5 million transcripts.
- Consequently, performing a BLAST search of the remaining ~5 million sequences took about 1 day to run.
- The subsequent annotation steps also take longer, and the final ncRNA transcriptome was ~70,000 sequences, which is fairly high.
- Lastly, we did not have a clear way to validate that these sequences were biologically relevant, and can only rely on the fact that the coding transcriptome QC showed TATAT worked well.

In summary, we suggest the ncRNA part of TATAT be used with some caution, but in principle it should work. However, the dramatically increased runtime makes it inconvenient to run, so unless there is a clear use case for the ncRNA transcriptome, we would not include it in a standard workflow.
