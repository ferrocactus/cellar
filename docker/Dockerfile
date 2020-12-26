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
        wget\
        metis\
        xorg-server-xvfb\
        nss\
        gtk2\
        gtk3\
        libxss\
        lib32-glibc\
        gcc-fortran\
        pandoc

RUN cd /opt && wget https://archive.archlinux.org/packages/l/llvm9-libs/llvm9-libs-9.0.1-1-x86_64.pkg.tar.zst
RUN cd /opt && pacman -U --noconfirm llvm9-libs-9.0.1-1-x86_64.pkg.tar.zst
RUN cd /opt && wget https://archive.archlinux.org/packages/l/llvm9/llvm9-9.0.1-1-x86_64.pkg.tar.zst
RUN cd /opt && pacman -U --noconfirm llvm9-9.0.1-1-x86_64.pkg.tar.zst

#### Install orca
RUN wget https://github.com/plotly/orca/releases/download/v1.1.1/orca-1.1.1-x86_64.AppImage -P /opt
RUN chmod 777 /opt/orca-1.1.1-x86_64.AppImage
RUN cd /opt && /opt/orca-1.1.1-x86_64.AppImage --appimage-extract
RUN printf '#!/bin/bash \nxvfb-run --auto-servernum --server-args "-screen 0 640x480x24" /opt/squashfs-root/app/orca "$@"' > /usr/bin/orca
RUN chmod 777 /usr/bin/orca
RUN cd /opt && wget https://archive.archlinux.org/packages/g/gconf/gconf-3.2.6%2B11%2Bg07808097-5-x86_64.pkg.tar.xz
RUN cd /opt && pacman -U --noconfirm gconf-3.2.6+11+g07808097-5-x86_64.pkg.tar.xz

### Install Cluster Ensembles
RUN pip3 install git+https://github.com/ferrocactus/Cluster_Ensembles
RUN echo "export PATH=$PATH:/usr/lib/python3.8/site-packages/Cluster_Ensembles/Hypergraph_Partitioning/hmetis-1.5-linux/" >> /root/.bashrc

EXPOSE 23123

COPY .Rprofile /root/.Rprofile
RUN mkdir ~/.R
RUN echo "MAKEFLAGS = -j8" > ~/.R/Makevars

#### Install R libraries
RUN Rscript -e 'install.packages(c("shiny", "shinydashboard", "shinyjs"))'
RUN Rscript -e 'install.packages(c("shinyBS", "shinyjqui", "shinyWidgets"))'
RUN Rscript -e 'install.packages(c("reticulate"))'
RUN Rscript -e 'install.packages(c("ggplot2", "gplots", "plotly"))'
RUN Rscript -e 'install.packages(c("rjson", "DT", "bsplus"))'
RUN Rscript -e 'install.packages(c("BiocManager"))'
RUN Rscript -e 'BiocManager::install("SingleR")'
RUN Rscript -e 'install.packages(c("htmlTable", "devtools"))'

#### Install python3 dependencies
# Breakpoint 1
ARG VER=unknown

RUN git clone https://github.com/ferrocactus/cellar /home/cellar
RUN cd /home/cellar && git config pull.rebase false

RUN mkdir -p /home/cellar/datasets/user_uploaded
RUN mkdir -p /home/cellar/datasets/server
WORKDIR /home/cellar

RUN pacman -S --noconfirm hdf5
RUN pip3 install wheel
RUN pip3 install numba==0.50.1 matplotlib numpy pandas scikit-learn sklearn umap-learn==0.3.10 kneed statsmodels joblib h5py anndata scanpy leidenalg psutil plotly bidict pydiffmap diffxpy scipy BinToGene

# Breakpoint 2
ARG VER2=unknown
RUN git pull

ENTRYPOINT ["Rscript", "app.R"]

