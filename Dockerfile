FROM ubuntu:22.04

WORKDIR /src

# Copy the tools directory which contains the Exonerate compiled binary
COPY ./src/tools /src/tools

# Install tools for compiling code and install dependencies
# --CD-HIT requires build-essential
# --SPAdes requires python3
RUN apt-get update \
    && apt-get install -y wget \
    && apt-get install build-essential -y \
    && apt install python3 -y \
    && apt install python3-pip -y

# Install SRA toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -P /src/tools \
    && tar xzvf /src/tools/sratoolkit.3.1.1-ubuntu64.tar.gz -C /src/tools \
    && rm /src/tools/sratoolkit.3.1.1-ubuntu64.tar.gz

# Configure SRA toolkit
ENV PATH="/src/tools/sratoolkit.3.1.1-ubuntu64/bin:$PATH"

# Install fastp
RUN wget http://opengene.org/fastp/fastp.0.23.4 -P /src/tools/fastp \
    && mv /src/tools/fastp/fastp.0.23.4 /src/tools/fastp/fastp \
    && chmod a+x /src/tools/fastp/fastp

# Configure fastp
ENV PATH="/src/tools/fastp:$PATH"

# Install SPAdes
RUN wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz -P /src/tools \
    && tar xzvf /src/tools/SPAdes-4.0.0-Linux.tar.gz -C /src/tools \
    && rm /src/tools/SPAdes-4.0.0-Linux.tar.gz

# Configure SPAdes
ENV PATH="/src/tools/SPAdes-4.0.0-Linux/bin:$PATH"

# Install for kmer coverage QC
RUN pip install pandas==2.2.3
RUN pip install seaborn==0.13.2

# Install CD-HIT-EST
RUN wget https://github.com/weizhongli/cdhit/archive/refs/tags/V4.8.1.tar.gz -P /src/tools \
    && tar xzvf /src/tools/V4.8.1.tar.gz -C /src/tools \
    && rm /src/tools/V4.8.1.tar.gz \
    && cd /src/tools/cdhit-4.8.1 \
    && make openmp=yes

# Configure CD-HIT-EST
ENV PATH="/src/tools/cdhit-4.8.1:$PATH"

# Install exonerate dependency
RUN apt-get update \
    && apt-get install -y ibglib2.0-dev

# Install exonerate
RUN tar xzvf /src/tools/exonerate-2.2.0-x86_64.tar.gz -C /src/tools \
    && rm /src/tools/exonerate-2.2.0-x86_64.tar.gz

# Configure exonerate fastanrdb
ENV PATH="/src/tools/exonerate-2.2.0-x86_64/bin:$PATH"

# Install BLAST+ dependencies
RUN apt-get install -y curl \
    && apt-get install -y perl \
    && apt-get install libgomp1

# Install BLAST+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz -P /src/tools \
    && tar xzvf /src/tools/ncbi-blast-2.16.0+-x64-linux.tar.gz -C /src/tools \
    && rm /src/tools/ncbi-blast-2.16.0+-x64-linux.tar.gz

# Configure BLAST+
ENV PATH="/src/tools/ncbi-blast-2.16.0+/bin:$PATH"

# Configure BLAST+ database
RUN mkdir -p "blastdb"
ENV BLASTDB="/src/blastdb"

# Install EvidentialGene
RUN wget http://arthropods.eugenes.org/EvidentialGene/other/evigene_old/evigene.tar -P /src/tools \
    && tar xvf /src/tools/evigene.tar -C /src/tools \
    && rm /src/tools/evigene.tar

# Configure EvidentialGene
ENV EVIGENE="/src/tools/evigene"