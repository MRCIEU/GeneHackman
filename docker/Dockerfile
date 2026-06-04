# This recipe uses amd64-only assets (Posit R .deb, PLINK, Miniconda, libpng12 .deb). On ARM hosts,
# build with: docker build --platform linux/amd64 ...
FROM python:3.12-bookworm

RUN apt update && \
    apt install -y --no-install-recommends \
    build-essential software-properties-common dirmngr gnupg pandoc lmodern texlive-latex-base unzip \
    libpng-dev libtiff5-dev libjpeg-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    libxml2-dev libcurl4-openssl-dev libssl-dev libgmp-dev libnlopt-dev cmake libcairo2-dev libxt-dev \
    texlive-latex-recommended r-cran-sass r-cran-mime libudunits2-dev libgdal-dev \
    libfontconfig1-dev libuv1-dev libicu-dev libbz2-dev liblzma-dev libglib2.0-dev \
    sqlite3 wget vim ripgrep git && \
    wget https://cdn.posit.co/r/debian-12/pkgs/r-4.4.2_1_amd64.deb && \
    apt-get update && \
    apt-get install -y ./r-4.4.2_1_amd64.deb && \
    rm r-4.4.2_1_amd64.deb && \
    ln -s /opt/R/4.4.2/bin/R /usr/local/bin/R && \
    ln -s /opt/R/4.4.2/bin/Rscript /usr/local/bin/Rscript && \
    mkdir -p /usr/lib/R/site-library && chmod 777 /usr/lib/R/site-library && \
    mkdir -p /home/scripts && \
    rm -rf /var/lib/apt/lists/* && \
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O /tmp/plink.zip && \
    unzip /tmp/plink.zip -d /usr/local/bin/ && \
    rm /tmp/plink.zip && \
    ln -s /usr/local/bin/plink /usr/local/bin/plink1.9 && \
    wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -O /usr/local/bin/liftOver && \
    chmod +x /usr/local/bin/liftOver && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    chmod +x miniconda.sh && \
   ./miniconda.sh -b -f -p /bin && \
    wget http://snapshot.debian.org/archive/debian/20160413T160058Z/pool/main/libp/libpng/libpng12-0_1.2.54-6_amd64.deb && \
    dpkg -i libpng12-0_1.2.54-6_amd64.deb

# Miniconda 25+ requires explicit acceptance of Anaconda repo channel ToS before installs.
RUN /bin/condabin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /bin/condabin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    cd /home && \
    git clone https://github.com/MRCIEU/PHESANT.git && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone https://github.com/samtools/bcftools.git && \
    cd /home/bcftools && \
    make && \
    cd /home && \
    git clone https://github.com/bulik/ldsc.git && \
    #hack to get around ldsc not working when RSID and SNP fields are defined
    sed -i "s/'SNP': 'SNP'/'PNS': 'PNS'/g" /home/ldsc/munge_sumstats.py && \
    cd /home/ldsc && \
    /bin/condabin/conda env create --file environment.yml

COPY . /home/GeneHackman
RUN Rscript /home/GeneHackman/docker/requirements.R
ENV R_LIBS_SITE=/usr/lib/R/site-library
RUN pip install --no-cache-dir -r /home/GeneHackman/docker/requirements.txt

CMD ["/bin/bash"]
