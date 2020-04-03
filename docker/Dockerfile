FROM library/ubuntu:latest
MAINTAINER Euxhen Hasanaj ehasanaj@cs.cmu.edu

RUN apt-get update && apt-get install -y \
        sudo\
        vim\
        git\
        wget\
        apt-transport-https\
        software-properties-common\
        build-essential\
        python3-pip\
        curl\
        libssl-dev\
        libcurl4-openssl-dev\
        libv8-dev

RUN apt-key adv --keyserver keyserver.ubuntu.com\
                --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

ADD tzdata.sh /tzdata.sh
RUN /tzdata.sh
RUN apt-get install -y r-base

EXPOSE 23123

RUN adduser --system --group --shell /bin/bash -u 1000 cellar \
        && echo "ucellar ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/user \
        && chmod 0440 /etc/sudoers.d/user

USER cellar

WORKDIR /home/cellar
ADD bashrc /home/cellar/.bashrc

git clone https://github.com/ferrocactus/cellar /home/cellar/package/

RUN python3 /home/cellar/setup.py install

# Need to install reticulate shiny shinyjs plotly V8 stringr BiocManager limma
# from R
