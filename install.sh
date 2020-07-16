#!/bin/bash

# Script to install cellar
# Source https://github.com/ferrocactus/cellar
# Maintainer: Euxhen Hasanaj ehasanaj@cs.cmu.edu

# Verify python is installed
if ! command -v python3 &> /dev/null
then
    echo "python3 not found in path. Please install python first."
    exit
fi
py_version=`python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))'`
req_py_version="3.7"

if ! [ "$(printf '%s\n' "$req_py_version" "$py_version" | sort -V | head -n1)" = "$req_py_version" ]
then
    echo "Python version 3.7 or higher is required."
    exit
fi

# Verify R is installed
if ! command -v R &> /dev/null
then
    echo "R not found in path. Please install R first."
    exit
fi
R_version=`R --version | grep -oP '(?<=version )[0-9]+'`
req_R_version="4"
if ! [ "$(printf '%s\n' "$req_R_version" "$R_version" | sort -V | head -n1)" = "$req_R_version" ]
then
    echo "R version 4 or higher is required."
    exit
fi

# Verify pip is installed
if ! command -v pip3 &> /dev/null
then
    echo "pip3 not found in path. Please install python first."
    exit
fi

# Install python packages first
echo "Installing Python dependencies..."
pip3 install -r py_requirements.txt
echo "Note: If you wish to use Cluster_Ensembles please follow the installation instructions in https://pypi.org/project/Cluster_Ensembles/"

echo "Installing R dependencies..."
# Install R packages
Rscript -e 'install.packages("shiny")'
Rscript -e 'install.packages("shinydashboard")'
Rscript -e 'install.packages("shinyjs")'
Rscript -e 'install.packages("shinyBS")'
Rscript -e 'install.packages("reticulate")'
Rscript -e 'install.packages("ggplot2")'
Rscript -e 'install.packages("plotly")'
Rscript -e 'install.packages("rjson")'
Rscript -e 'install.packages("DT")'
Rscript -e 'install.packages("gplots")'
Rscript -e 'install.packages("bsplus")'
Rscript -e 'install.packages("htmlTable")'
Rscript -e 'install.packages("BiocManager")'
Rscript -e 'BiocManager::install("SingleR")'
echo "Note: If you wish to convert plots to images, first install orca by following the instructions in https://github.com/plotly/orca"

echo "Compiling Cython modules..."
python3 src/methods/setup.py build_ext --inplace
python3 src/methods/setup.py clean

echo "cellar installed successfully. Run it with $ Rscript app.R"
