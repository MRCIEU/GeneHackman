FROM python:3.11-slim

RUN apt update && \
    apt install -y --no-install-recommends r-base=4.2.2.20221110-2 r-base-dev=4.2.2.20221110-2 \
    build-essential software-properties-common dirmngr gnupg pandoc lmodern texlive-latex-base \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \ 
    libxml2-dev libcurl4-gnutls-dev libssl-dev libgmp-dev libnlopt-dev cmake libcairo2-dev libxt-dev \
    texlive-latex-recommended r-cran-sass r-cran-mime libudunits2-dev libgdal-dev \
    plink1.9 plink2 sqlite3 wget vim ripgrep git && \
    mkdir -p /home/scripts && \
    rm -rf /var/lib/apt/lists/* && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    chmod +x miniconda.sh && \
   ./miniconda.sh -b -f -p /bin && \
    wget http://snapshot.debian.org/archive/debian/20160413T160058Z/pool/main/libp/libpng/libpng12-0_1.2.54-6_amd64.deb && \
    dpkg -i libpng12-0_1.2.54-6_amd64.deb

RUN cd /home && \
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

#TODO: other things to think about installing: bcftools, finemap, some sort of meta-analysis tool (METAL, GWAMA, etc)
#RUN wget http://christianbenner.com/finemap_v1.4.2_x86_64.tgz && \
#    tar -xf finemap_v1.4.2_x86_64.tgz && \
#    mv finemap_v1.4.2_x86_64 /usr/local/bin/finemap

COPY docker/requirements.* docker/
RUN Rscript docker/requirements.R
RUN pip install -r docker/requirements.txt

CMD ["/bin/bash"]
