FROM ubuntu:18.04

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython \
                   python3-pip \
                   python3-pysam \
                   wget \
                   unzip && \
    pip3 install biopython \
                 ete3 \
                 pycosat \
                 PyYAML \
                 requests \
                 numpy \
                 pandas \
                 bokeh

# Install dependencies from conda 
RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc" && \
    conda install -c bioconda bowtie2 krakenuniq kallisto gmap snap-aligner openssl=1.0 samtools bedtools bwa mafft bcftools tabix && \
    conda clean -afy
# Install Picard 

RUN wget https://github.com/broadinstitute/picard/releases/download/2.18.15/picard.jar -P /usr/bin/
# ENV PATH="/picard.jar:${PATH}"
# ENV picard.jar /picard.jar

# Install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.0.11.0/gatk-4.0.11.0.zip && unzip gatk-4.0.11.0.zip
ENV PATH=/gatk-4.0.11.0/:$PATH
#ENV GATK gatk-4.0.11.0/gatk

# Install VarScan 
RUN wget --no-check-certificate https://sourceforge.net/projects/varscan/files/latest/download && mv download VarScan
RUN mv VarScan /usr/local/bin/

# Install gff3ToGenePred
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
RUN chmod +x gff3ToGenePred
RUN mv gff3ToGenePred /usr/local/bin



###########
# ANNOVAR #
###########


# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz


# INSTALL
RUN wget https://github.com/vpeddu/lava/raw/docker/annovar.zip && unzip annovar.zip 
RUN mv annovar/*.pl /usr/local/bin/


##########
# JAVA 8 #
##########

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME



CMD ["/bin/bash"]