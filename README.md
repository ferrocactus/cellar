<h1>Cell Type Annotation Robot (Cellar)</h1>

Single-cell transcriptomics coupled with recent technological
developments in high-throughput RNA sequencing have greatly
expanded our knowledge of complex tissues and their cellular
composition. As the number of samples increases, there is a
growing need to correctly and automatically identify cell
types by examining the gene expression profiles of individual
cells. Several methods have been proposed by the Computational
Biology and Machine Learning community for each step of the data
analysis pipeline. Cellar is an attempt to put many of these methods
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

<h3>Method 2: Use the install script</h3>

Assuming `python (>=3.7)` and `R (>=4.0)` are installed, run the `install.sh` script
and that will take care of installing the dependencies for you.

```bash
git clone https://github.com/ferrocactus/cellar
cd cellar
bash install.sh
```

Finally, run with
```bash
Rscript app.R
```

Note: Cluster_Ensembles (Optional) and orca (Optional) need to be installed manually.
Follow instructions in <a href="https://pypi.org/project/Cluster_Ensembles/" target="_blank">
https://pypi.org/project/Cluster_Ensembles/</a> and
<a href="https://github.com/plotly/orca" target="_blank">https://github.com/plotly/orca</a>,
respectively.

<h3>Method 3: Manual Installation</h3>

Cellar requires `python (>=3.7)` and `R (>=4.0)`.

Cellar requires the following python packages that can be installed via `pip`:

* matplotlib (>=3.1.2)
* numpy (>=1.18.2)
* pandas (>=0.25.3)
* scikit-learn (>=0.22)
* umap-learn (>=0.3.10)
* kneed (>=0.5.1)
* statsmodels (>=0.11)
* joblib (>=0.14.1)
* anndata (>=0.7.1)
* leidenalg (>=0.8.0)
* scanpy (>=1.5.1)
* psutil (>=5.7.0)
* plotly (>=4.8.1)
* bidict (>=0.19.0)
* pydiffmap (>=0.2.0.1)
* scipy (>=1.4.0)
* Cluster_Ensembles (>=1.16)
* diffxpy (>=0.7.4)

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
* SingleR (installed via BiocManager)
* htmlTable
* doSNOW
* GenomicRanges (installed via BiocManager)
* plot3D
* SnapATAC (installed via devtools)
* shinyWidgets

Assuming all dependencies have been met, clone the github repo
and cd into it

```bash
git clone https://github.com/ferrocactus/cellar
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

<h2>Documentation</h2>

You can find the documentation <a href="doc/cellar_guide.md" target="_blank">here</a>.

Alternatively, you can find the documentation to the cellar python library <a href="doc/cellar_python_guide.rst" target="_blank">here</a>.
