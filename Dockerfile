FROM debian:jessie
MAINTAINER Evan Floden <evanfloden@gmail.com>

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


# 
# install R 
#
RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list &&\
 apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 &&\
 apt-get update --fix-missing && \
 apt-get -y install r-base
		
#
# Install R libraries 
# https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#
RUN R -e 'install.packages(c("optparse", "RColorBrewer", "reshape2", "ggplot2", "devtools"), repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_github("kcha/psiplot")'

#
# Install Bowtie v1
# 
RUN wget -q https://github.com/BenLangmead/bowtie/releases/download/v1.2.1.1/bowtie-1.2.1.1-linux-x86_64.zip && \
    unzip bowtie-1.2.1.1-linux-x86_64.zip && \
    ln -s bowtie-1.2.1.1/bowtie /usr/bin/bowtie 
