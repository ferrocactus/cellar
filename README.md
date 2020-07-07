<h1><strong>Cell</strong> Type <strong>A</strong>nnotation
<strong>R</strong>obot (Cellar)</h1>

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
maximum flexibility when choosing each step of the pipeline.

Cellar is developed and maintained by <a href="http://www.sb.cs.cmu.edu/">
    Systems Biology Group</a> at <a href="https://www.cmu.edu/">
    Carnegie Mellon University</a>.

<h2>How to use?</h2>
Cellar is freely available as an online app at
<a href="https://data.test.hubmapconsortium.org/app/cellar">
https://data.test.hubmapconsortium.org/app/cellar</a>. For a local
installation see below.

<h2>Installation</h2>
If you wish to install Cellar locally, follow one of these methods:

<h3>Method 1: Docker</h3>
If you have docker installed, you can pull the image with

```console
docker pull euxhen/cellar:v1
```

and then run it

```console
docker run -p 23123:23123 euxhen/cellar:v1
```

After a short delay, the container will start listening to port 23123.
Visit `localhost:23123` in your web browser to start the app.

Alternatively, if you have a local folder where you have stored
your datasets, you can avoid having to copy files in and out of
the container by instead running

```console
docker run -p 23123:23123 -v /path/to/your/datasets/folder:/home/cellar/datasets/server euxhen/cellar:v1
```

You can find your datasets by checking `Server Datasets` under
the `Dataset` menu once the app starts.

<h3>Method 2: Manual Installation</h3>

Cellar requires `python>=3.7` and `R>=4.0`. Assuming these are
installed, clone the github repo into a local folder

```console
git clone https://github.com/ferrocactus/cellar
```