FROM ubuntu:22.04

WORKDIR /src

# Copy the app directory which contains the scripts necessary to run TATAT
COPY ./src/app /src/app

# Copy the tools directory which contains the Exonerate compiled binary
COPY ./src/tools /src/tools

# Install tools for compiling code and install dependencies
# --CD-HIT requires build-essential
# --SPAdes requires python3
RUN apt-get update \
    && apt-get install -y wget \
    && apt-get install build-essential -y \
    && apt install python3 -y \
    && apt install python3-pip -y \
    $$ unzip

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

# Install Diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz -P /src/tools/diamond \
    && tar xzvf /src/tools/diamond/diamond-linux64.tar.gz -C /src/tools/diamond \
    && rm /src/tools/diamond/diamond-linux64.tar.gz

# Configure Diamond
ENV PATH="/src/tools/diamond:$PATH"

# Install NCBI Datasets & Dataformat
RUN wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets -P /src/tools/datasets \
    && wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat -P /src/tools/dataformat

# Configure Datasets & Dataformat
ENV PATH="/src/tools/datasets:$PATH"
ENV PATH="/src/tools/dataformat:$PATH"
RUN chmod +x /src/tools/datasets/datasets /src/tools/dataformat/dataformat

# Install matplotlib-venn for post annotation analyses
RUN pip3 install matplotlib-venn==1.1.2

# Install BUSCO gene set analysis dependencies HMMER, biopython, requests
RUN apt-get install -y hmmer=3.3.2+dfsg-1 \
    && pip3 install biopython \
    && pip3 install requests

# Install BUSCO, but limited to gene set analysis
RUN wget https://gitlab.com/ezlab/busco/-/archive/5.8.3/busco-5.8.3.tar.gz -P /src/tools \
    && tar xzvf /src/tools/busco-5.8.3.tar.gz -C /src/tools \
    && rm /src/tools/busco-5.8.3.tar.gz \
    && cd /src/tools/busco-5.8.3 \
    && python3 -m pip install .

# Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz -P /src/tools \
    && tar xzvf /src/tools/salmon-1.10.0_linux_x86_64.tar.gz -C /src/tools \
    && rm /src/tools/salmon-1.10.0_linux_x86_64.tar.gz

# Configure Salmon
ENV PATH="/src/tools/salmon-latest_linux_x86_64/bin:$PATH"

# Install sklearn for post analysis
RUN pip install scikit-learn==1.6.1

# Install sqlite for metadata storage
RUN apt install -y sqlite3 libsqlite3-dev