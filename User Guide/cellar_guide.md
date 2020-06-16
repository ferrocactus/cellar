## README

Documentation for the single cell analysis app hosted by the Bar-Joseph [Systems Biology group](http://www.sb.cs.cmu.edu/) at Carnegie Mellon University accessible [here](https://data.test.hubmapconsortium.org/app/cellar).

---

### Choose Dataset

<img src="datasetMenu.png" alt="datasetMenu"
	title="Choose Dataset" width="300" height="300" align=left/>

In this tab, users can either upload their own single cell gene expression data
as a _csv_ or _h5ad_ file, with rows corresponding to cells and columns corresponding
to genes. For gene ids, please make sure they are either in _HGNC_, _entrez ID_, or _ENSEMBL_ format.

Clicking the **Load Dataset** button loads the dataset into the app for use.

---

### Clustering
<img src="clusteringMenu.png" alt="clusrterMenu"
	title="Choose Dataset" width="300" height="600" />

The **Run with current configuration** button is to be clicked after the user has selected the desired clustering parameters below, or is satisfied the provided defaults.

- **Dimensionality reduction** : Enables user to choose the method to reduce dimensions prior to clustering.

- **Number of components** : Enables the user to either manually choose the number of components after dimensionality reduction that clustering is carried on, or choose automatic to let the program choose optimal number of components based on explained variance.

- **Clustering** : Enables the user to choose the clustering algorithm.

- **Number of clusters** : Enables the user to choose the number of clusters they wish the data to correspond to.

- **Evaluation** :

- **Visualization Method** :

The following methods are for

- **Constrained**


---

### Label Transfer

---

### Selection and Labelling

---

### Analysis

---

### Appearance

---

### Saving and Loading sessions
