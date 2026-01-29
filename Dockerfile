FROM rocker/tidyverse:4.4.0

LABEL description="IMPACT Trial peptide discovery pipeline environment"

# =============================================================================
# System Dependencies
# =============================================================================
RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
        bedtools \
        bzip2 \
        cpanminus \
        curl \
        libbz2-dev \
        libcurl4-openssl-dev \
        libgdal-dev \
        liblzma-dev \
        libmagick++-dev \
        libreadline-dev \
        libssl-dev \
        libudunits2-dev \
        libxml2-dev \
        libxt-dev \
        python3-dev \
        python3-pip \
        zlib1g && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Link python3 to python for compatibility
RUN ln -sf /usr/bin/python3 /usr/bin/python

# =============================================================================
# Python Packages
# =============================================================================
RUN python3 -m pip install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir \
        rpg==2.0.3 \
        biopython

# =============================================================================
# R Package Management
# =============================================================================
# Fix matrixStats version conflict (Bioconductor 3.20 requires >= 1.4.1)
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/matrixStats_1.5.0.tar.gz', \
    repos=NULL, type='source')"

# Configure Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)" && \
    R -e "BiocManager::install(version='3.20', ask=FALSE)"

# Install Bioconductor packages
RUN R -e 'BiocManager::install(c( \
    "biomaRt", \
    "Biostrings", \
    "BSgenome.Hsapiens.UCSC.hg38", \
    "data.table", \
    "DBI", \
    "GenomicFeatures", \
    "GenomicRanges", \
    "ggthemes", \
    "optparse", \
    "RColorBrewer", \
    "R.utils", \
    "RSQLite", \
    "rtracklayer", \
    "tidytext" \
), ask=FALSE)'

# Install CRAN Archive packages
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/glue_1.8.0.tar.gz', \
    repos=NULL, type='source')"

# Install GitHub packages
RUN R -e "remotes::install_github('clauswilke/colorblindr', \
    ref='1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies=TRUE)" && \
    R -e "remotes::install_github('thomasp85/patchwork', \
    ref='1cb732b129ed6a65774796dc1f618558c7498b66', dependencies=TRUE)" && \
    R -e "remotes::install_github('qsbase/qs2', \
    ref='6f835e54eb7d10123051f44894ee1001566094fb', dependencies=TRUE)"

# =============================================================================
# Bioinformatics Tools
# =============================================================================
# BLAST+ v2.17.0
RUN curl -fsSL https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz \
    | tar -xz && \
    mv ncbi-blast-2.17.0+/bin/* /usr/local/bin/ && \
    rm -rf ncbi-blast-2.17.0+

# SeqKit v2.10.1
RUN curl -fsSL https://github.com/shenwei356/seqkit/releases/download/v2.10.1/seqkit_linux_amd64.tar.gz \
    | tar -xz -C /usr/local/bin/ && \
    chmod +x /usr/local/bin/seqkit

# =============================================================================
# Install Pipeline Scripts
# =============================================================================
WORKDIR /opt/impact_trial

# Copy scripts (exclude input and results directories)
COPY scripts/create-gtex-matrices/*.R /opt/impact_trial/scripts/create-gtex-matrices/
COPY scripts/splice-junction-filtering/*.R /opt/impact_trial/scripts/splice-junction-filtering/
COPY scripts/translation/*.py /opt/impact_trial/scripts/translation/
COPY scripts/translation/*.R /opt/impact_trial/scripts/translation/
COPY scripts/translation/*.sh /opt/impact_trial/scripts/translation/
COPY scripts/translation/util/ /opt/impact_trial/scripts/translation/util/

# Make scripts executable and readable
RUN find /opt/impact_trial/scripts -type f \( -name "*.py" -o -name "*.R" -o -name "*.sh" \) -exec chmod 755 {} \; && \
    find /opt/impact_trial/scripts -type d -exec chmod 755 {} \; && \
    chown -R 1000:1000 /opt/impact_trial/scripts

# Add scripts to PATH
ENV PATH="/opt/impact_trial/scripts/create-gtex-matrices:/opt/impact_trial/scripts/splice-junction-filtering:/opt/impact_trial/scripts/translation:${PATH}"

# =============================================================================
# Finalize
# =============================================================================
WORKDIR /data

CMD ["/bin/bash"]
