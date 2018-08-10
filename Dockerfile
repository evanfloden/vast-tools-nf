FROM debian:jessie
MAINTAINER Evan Floden <evanfloden@gmail.com>

# Dockerfile with for VAST-TOOLS (Vertebrate Alternative Splicing and Transcription Tools)
# See https://github.com/vastgroup/vast-tools for more information on the software

# Install basic requirements
RUN apt-get update \ 
	&& apt-get install -y --no-install-recommends \
        ed \
        less \
	locales \
	vim-tiny \
	wget \
	ca-certificates \
        build-essential \
        cmake \
        curl \
        libcurl4-gnutls-dev \
        libssl-dev \
        perl \
	&& rm -rf /var/lib/apt/lists/*

# install R  (required for the R packages below)
RUN echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list &&\
    #apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' &&\
    apt-get update --fix-missing && \
    apt-get -y --force-yes install r-base
		
# install R packages
RUN R -e 'install.packages(c("optparse", \
                             "RColorBrewer", \
                             "reshape2", \
                             "ggplot2", \
                             "devtools"), repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_github("kcha/psiplot")'

# install Bowtie
RUN apt-get install -y python && \
    wget -q https://github.com/BenLangmead/bowtie/releases/download/v1.2.1.1/bowtie-1.2.1.1-linux-x86_64.zip && \
    unzip bowtie-1.2.1.1-linux-x86_64.zip && \
    rm bowtie-1.2.1.1-linux-x86_64.zip && \
    ln -s /bowtie-1.2.1.1/bowtie /usr/bin/bowtie

# install VAST-TOOLS
RUN apt-get install -y git && \
    git clone https://github.com/vastgroup/vast-tools.git && \
    cd vast-tools && \
    git checkout v2.0.2 && \
    chmod +x vast-tools && \
    ln -s /vast-tools/vast-tools /usr/bin/vast-tools 
