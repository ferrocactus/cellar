*******************************************
Cell Annotation Robot (cellar) Python guide
*******************************************
In addition to the cellar interface provided in
`<https://data.test.hubmapconsortium.org/app/cellar/>`_,
cellar can also be used as a standalone python library,
although its functionality is more limited.
Here we provide a quick tour of the main components
and describe how to extend its functionality.

To get started, you can install cellar using pip by runnin
the following command (requires Python >3.7):

    $ pip install git+https://github.com/ferrocactus/cellar

This will install any required dependencies including popular
frameworks such as `numpy`_, `scikit-learn`_, `pandas`_ and more.
For a full list see `py_requirements`_.

You can test the installation by running the following simple snippet::

    import cellar as cl

    adata = cl.load_file('default')

    cl.reduce_dim(adata)
    cl.cluster(adata)
    cl.reduce_dim_vis(adata)
    cl.plot('adata')



.. _numpy: https://numpy.org/
.. _scikit-learn: https://scikit-learn.org/stable/
.. _pandas: https://pandas.pydata.org/
.. _py_requirements: https://github.com/ferrocactus/cellar/blob/master/_py_requirements.txt