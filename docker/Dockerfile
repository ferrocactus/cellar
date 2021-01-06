FROM library/archlinux:latest
MAINTAINER Euxhen Hasanaj ehasanaj@cs.cmu.edu

### Install apps
RUN pacman -Syyu --noconfirm

RUN pacman -S --noconfirm\
        base-devel\
        r\
        git\
        python\
        python-pip\
        python-virtualenv\
        wget\
        metis\
        xorg-server-xvfb\
        nss\
        gtk2\
        gtk3\
        libxss\
        lib32-glibc\
        gcc-fortran\
        pandoc\
        llvm\
        hdf5\
        libgit2

#### Install orca
RUN wget https://github.com/plotly/orca/releases/download/v1.1.1/orca-1.1.1-x86_64.AppImage -P /opt
RUN chmod 777 /opt/orca-1.1.1-x86_64.AppImage
RUN cd /opt && /opt/orca-1.1.1-x86_64.AppImage --appimage-extract
RUN printf '#!/bin/bash \nxvfb-run --auto-servernum --server-args "-screen 0 640x480x24" /opt/squashfs-root/app/orca "$@"' > /usr/bin/orca
RUN chmod 777 /usr/bin/orca
RUN cd /opt && wget https://archive.archlinux.org/packages/g/gconf/gconf-3.2.6%2B11%2Bg07808097-5-x86_64.pkg.tar.xz
RUN cd /opt && pacman -U --noconfirm gconf-3.2.6+11+g07808097-5-x86_64.pkg.tar.xz

EXPOSE 23123

COPY .Rprofile /root/.Rprofile
RUN mkdir ~/.R
RUN echo "MAKEFLAGS = -j8" > ~/.R/Makevars
RUN Rscript -e 'x <- file.path(R.home("doc"), "html"); if (!file.exists(x)) {dir.create(x, recursive=TRUE); file.copy(system.file("html/R.css", package="stats"), x)}'

#### Install R libraries
RUN Rscript -e 'install.packages(c("shiny", "shinydashboard", "shinyjs"))'
RUN Rscript -e 'install.packages(c("shinyBS", "shinyjqui", "shinyWidgets"))'
RUN Rscript -e 'install.packages(c("reticulate"))'
RUN Rscript -e 'install.packages(c("ggplot2", "gplots", "plotly"))'
RUN Rscript -e 'install.packages(c("rjson", "DT", "bsplus"))'
RUN Rscript -e 'install.packages("htmlTable")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'BiocManager::install("SingleR")'
RUN Rscript -e 'devtools::install_github("aertslab/cisTopic")'

RUN git clone https://github.com/ferrocactus/cellar /home/cellar
RUN cd /home/cellar && git config pull.rebase false

RUN mkdir -p /home/cellar/datasets/user_uploaded
RUN mkdir -p /home/cellar/datasets/server

# Install python3.7
RUN useradd -m nonroot
RUN echo "nonroot ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/nonroot
USER nonroot
RUN cd /home/nonroot && git clone https://aur.archlinux.org/python37.git
RUN cd /home/nonroot/python37 && makepkg -scri --noconfirm --skippgpcheck

RUN virtualenv --python=/usr/bin/python3.7 /home/nonroot/py3.7

# Breakpoint 1
ARG VER=unknown
RUN /home/nonroot/py3.7/bin/pip3 install numba matplotlib numpy pandas\
        scikit-learn sklearn umap-learn kneed statsmodels joblib h5py\
        anndata scanpy leidenalg psutil plotly bidict pydiffmap\
        diffxpy scipy BinToGene

RUN /home/nonroot/py3.7/bin/pip3 install tifffile

### Install Cluster Ensembles
RUN /home/nonroot/py3.7/bin/pip3 install git+https://github.com/ferrocactus/Cluster_Ensembles
USER root
WORKDIR /home/cellar

COPY gpmetis /home/nonroot/py3.7/lib/python3.7/site-packages/Cluster_Ensembles/Hypergraph_Partitioning/hmetis-1.5-linux/gpmetis
ENV PATH="/home/nonroot/py3.7/bin:$PATH:/home/nonroot/py3.7/lib/python3.7/site-packages/Cluster_Ensembles/Hypergraph_Partitioning/hmetis-1.5-linux"

# Breakpoint 2
ARG VER2=unknown
RUN pip install BinToGene --upgrade
RUN git pull

ENTRYPOINT ["Rscript", "app.R"]
