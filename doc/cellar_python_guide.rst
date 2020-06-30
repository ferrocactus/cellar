*******************************************
Cell Annotation Robot (cellar) Python guide
*******************************************
In addition to the interface provided in
`<https://data.test.hubmapconsortium.org/app/cellar/>`_,
cellar can also be used as a standalone python library,
although its functionality is more limited.
Here we provide a quick tour of the main components
and describe how to extend existing methods.

Installation
____________

To get started, you can install cellar using pip by running
the following command (requires Python >=3.7)

.. code:: bash

    $ pip install git+https://github.com/ferrocactus/cellar

This will install any required dependencies including popular
frameworks such as `numpy`_, `scikit-learn`_, `pandas`_ and more.
For a full list see `py_requirements`_.

You can test the installation by running the following simple snippet

.. code:: python

    import cellar as cl

    adata = cl.load_file('default')

    cl.reduce_dim(adata)
    cl.cluster(adata)
    cl.reduce_dim_vis(adata)
    cl.plot(adata)

This simple snippet also demonstrates the most basic pipeline
that you can run. Essentially, after loading the data (normalized)
into an `AnnData`_ object, we reduce the dimensionality
of the data (default: `PCA`_), followed by a clustering algorithm
(default: `Leiden`_) and a dimensionality reduction method for visualization
(default: `UMAP`_). Finally, we plot the results using `Plotly`_ (this will
open a new web page in your default browser).

Available Units
_______________

The following units are available in the python version of cellar:

* ``cl.reduce_dim``: depending on the chosen method, it performs a linear or a non-linear mapping of the data to a lower dimensional space.
* ``cl.cluster``: clusters the reduced data into a number of clusters as specified by the user, or as automatically detected using an evaluation method.
* ``cl.reduce_dim_vis``: further reduces the dimensionality of the data into 2D embeddings used for visualization.
* ``cl.name_genes``: finds and stores synonyms for gene ids (e.g., genes provided in ensembl format get converted to gene names).
* ``cl.de``: runs differential expression for a given subset of the data vs the rest or another subset.
* ``cl.ss_cluster``: runs a semi-supervised clustering algorithm. Useful when clusters have been modified by the user, and the user wishes to refine them.
* ``cl.transfer_labels``: given a reference labeled dataset, use that to transfer the labels into another dataset using the method specified.
* ``cl.plot``: plots the 2D embeddings colored by cluster, or colored by the expression value of a gene.

.. _numpy: https://numpy.org/
.. _scikit-learn: https://scikit-learn.org/stable/
.. _pandas: https://pandas.pydata.org/
.. _py_requirements: https://github.com/ferrocactus/cellar/blob/master/py_requirements.txt
.. _AnnData: https://anndata.readthedocs.io/en/stable/anndata.AnnData.html
.. _PCA: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
.. _Leiden: https://github.com/vtraag/leidenalg
.. _UMAP: https://umap-learn.readthedocs.io/en/latest/
.. _Plotly: https://plotly.com/
.. _cellar: https://data.test.hubmapconsortium.org/app/cellar
