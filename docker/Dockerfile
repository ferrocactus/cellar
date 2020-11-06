FROM library/archlinux:latest
MAINTAINER Euxhen Hasanaj ehasanaj@cs.cmu.edu

### Install apps
RUN pacman -Syyu --noconfirm

RUN pacman -S --noconfirm\
        base-devel\
        python3\
        python-pip\
        r\
        git\
        wget\
        xorg-server-xvfb\
        gtk2\
        metis\
        lib32-glibc\
        gcc-fortran\
        libxtst\
        libxss\
        nss\
        alsa-lib

#### Install orca
RUN wget https://github.com/plotly/orca/releases/download/v1.1.1/orca-1.1.1-x86_64.AppImage -P /opt
RUN chmod 777 /opt/orca-1.1.1-x86_64.AppImage
RUN cd /opt && /opt/orca-1.1.1-x86_64.AppImage --appimage-extract
RUN printf '#!/bin/bash \nxvfb-run --auto-servernum --server-args "-screen 0 640x480x24" /opt/squashfs-root/app/orca "$@"' > /usr/bin/orca
RUN chmod 777 /usr/bin/orca
RUN cd /opt && wget https://archive.archlinux.org/packages/g/gconf/gconf-3.2.6%2B11%2Bg07808097-5-x86_64.pkg.tar.xz
RUN cd /opt && pacman -U --noconfirm gconf-3.2.6+11+g07808097-5-x86_64.pkg.tar.xz

#### Install Cluster Ensembles
RUN pip3 install git+https://github.com/ferrocactus/Cluster_Ensembles
RUN echo "export PATH=$PATH:/usr/lib/python3.8/site-packages/Cluster_Ensembles/Hypergraph_Partitioning/hmetis-1.5-linux/" >> /root/.bashrc

EXPOSE 23123

COPY .Rprofile /root/.Rprofile

#### Install R libraries
RUN Rscript -e 'install.packages("shiny")'
RUN Rscript -e 'install.packages("shinydashboard")'
RUN Rscript -e 'install.packages("shinyjs")'
RUN Rscript -e 'install.packages("shinyBS")'
RUN Rscript -e 'install.packages("reticulate")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("plotly")'
RUN Rscript -e 'install.packages("rjson")'
RUN Rscript -e 'install.packages("DT")'
RUN Rscript -e 'install.packages("gplots")'
RUN Rscript -e 'install.packages("bsplus")'
RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install("SingleR")'
RUN Rscript -e 'install.packages("htmlTable")'
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'install.packages("shinyWidgets")'
RUN Rscript -e 'install.packages("shinyjqui")'

#### Install python3 dependencies

ADD https://api.github.com/repos/ferrocactus/cellar/git/refs/heads/master version.json
RUN git clone https://github.com/ferrocactus/cellar /home/cellar
RUN cd /home/cellar && git config pull.rebase false

RUN mkdir -p /home/cellar/datasets/user_uploaded
RUN mkdir -p /home/cellar/datasets/server
WORKDIR /home/cellar

RUN pip3 install -r /home/cellar/py_requirements.txt

ADD https://api.github.com/repos/ferrocactus/cellar/git/refs/heads/master version.json
RUN git pull

ENTRYPOINT ["Rscript", "app.R"]
