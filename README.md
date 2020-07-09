<h1>Cell Type Annotation Robot (Cellar)</h1>

Single-cell transcriptomics coupled with recent technological
developments in high-throughput RNA sequencing have greatly
expanded our knowledge of complex tissues and their cellular
composition. As the number of samples increases, there is a
growing need to correctly and automatically identify cell
types by examining the gene expression profiles of individual
cells. Several methods have been proposed by the Computational
Biology and Machine Learning community for each step of the
pipeline. Cellar is an attempt to put many of these methods
into a single and easy-to-use interface, while providing the user
with maximum flexibility when choosing each step of the pipeline.

Cellar is developed and maintained by <a href="http://www.sb.cs.cmu.edu/"
    target="_blank">Systems Biology Group</a> at
    <a href="https://www.cmu.edu/" target="_blank">Carnegie Mellon University</a>.

<h2>How to use?</h2>
Cellar is freely available as an online app at
<a href="https://data.test.hubmapconsortium.org/app/cellar" target="_blank">
https://data.test.hubmapconsortium.org/app/cellar</a>. For a local
installation see below.

<h2>Installation</h2>
If you wish to install Cellar locally, follow one of these methods:

<h3>Method 1: Docker</h3>
If you have docker installed, you can pull the image with

```bash
docker pull euxhen/cellar
```

and then run it

```bash
docker run -p 23123:23123 euxhen/cellar
```

After a short delay, the container will start listening to port 23123.
Visit `localhost:23123` in your web browser to start the app.

Alternatively, if you have a local folder with your datasets,
you can avoid copying files in and out of the container by
instead running

```bash
docker run -p 23123:23123 -v /path/to/your/datasets/folder:/home/cellar/datasets/server euxhen/cellar
```

You can find your datasets by checking `Server Datasets` under
the `Dataset` menu once the app starts.

<h3>Method 2: Manual Installation</h3>

Cellar requires `python (>=3.7)` and `R (>=4.0)`.

Cellar requires the following python packages that can be installed via `pip`:

* matplotlib
* numpy
* pandas
* scikit-learn
* sklearn
* umap-learn
* kneed
* statsmodels
* joblib
* anndata
* leidenalg
* scanpy
* psutil
* plotly
* bidict
* pydiffmap
* Cython
* Cluster_Ensembles

NOTE: Cluster_Ensembles package is only available for linux systems.
If using the latest version of sklearn, you can install via

```bash
pip install git+https://github.com/ferrocactus/Cluster_Ensembles
```
and also must have metis in your path. Follow the instructions in
<a href="https://pypi.org/project/Cluster_Ensembles/" target="_blank">
https://pypi.org/project/Cluster_Ensembles/</a>.

You will also need the following R packages

* shiny
* shinyjs
* shinydashboard
* shinyBS
* ggplot2
* plotly
* reticulate
* rjson
* DT
* gplots
* htmltools
* bsplus
* SingleR (installed via BiocManager, follow instructions in https://bioconductor.org/packages/release/bioc/html/SingleR.html)

Assuming all dependencies have been met, clone the github repo
and cd into it

```bash
git clone https://github.com/ferrocactus/cellar
cd cellar
```

Next, you will need to compile Cython modules. Run

```bash
python src/methods/setup.py build_ext --inplace
python src/methods/setup.py clean
```

(Optional) To be able to convert plots into images
you also need to hava orca installed. Follow instructions in
<a href="https://github.com/plotly/orca" target="_blank">
https://github.com/plotly/orca</a>.

Now you are all set. To start the app run

```bash
Rscript app.R
```

and open `localhost:23123` in your web browser.