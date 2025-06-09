# TATAT Rousettus Tutorial

The following is a tutorial for using TATAT to process raw RNA-seq data from the *Rousettus aegyptiacus* bat and make an annotated transcriptome. 
The data includes 20 samples from 11 tissues, with half the samples from an adult male bat and the other half from an adult female bat. 
The exact details can be found in the [sample_metadata.csv](rousettus_tutorial/sample_metadata.csv) file.
<br><br>
The main function of TATAT is to generate a transcriptome of coding genes, and so that will be covered first.
However, TATAT also has the ability to add non-coding RNA (ncRNA) genes to the transcriptome. It requires that the coding transcriptome already 
be generated, and we suspect will generally be of less interest to users, so that will be covered second.
