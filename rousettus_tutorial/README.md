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
In this tutorial, TATAT will be used to download the sequencing data from the Sequencing Read Archive (SRA), remove sequencing errors and adaptors using FASTP, assemble the reads into contigs via *de novo* assembly with rnaSPAdes, thin out the large number of contigs, then annotate the remaning contigs and select which best represent each gene.
