## README

Documentation for the single cell analysis app hosted by the Bar-Joseph [Systems Biology group](http://www.sb.cs.cmu.edu/) at Carnegie Mellon University accessible [here](https://data.test.hubmapconsortium.org/app/cellar).

---
### The Overall Layout


<img src="UI.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

The above graph show the layout of our single cell analysis app. 

On the left is the sidebar, with toggle panels enabling users to do various analysis as well as configuration settings, which will be introduced in more details later. 

<img src="UI_expression.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

On the right is the body, which shows the results of the analysis. The top of the body is the main plot visualizing the cells being analyzed. The colors can represent the clusters of cells  or the gene expression level of any selected gene. For example the expression level of the gene IL7R is shown in the above plot. When the expression level is shown, shapes of the points represent the clusters they belong to. 

The button "view cluster names" controls the collapsible table showing the user assigned names for each cluster and subset of cells. Detailed illustration of assigning names operation is in the "selection & labeling" section.

<img src="UI_analysis.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

The bottom of the body shows the results of various analysis that users can do, which will be illustrated in the "analysis" section.

The subsequent sections will introduce the functions in the sidebar panels.

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
<img src="selection.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

This panel helps users to accomplish 3 tasks: 

1. Selecting a group of cells and make them a new subset
2. Adding a new label for cells
3. Changing the label (cell type) of any defined subsets

For 1:

- **New Subset** : Input the name for the new subset.

- **Add Subset** : After selecting a group of cells in the main plot and entering the new subset name, clicking this button will define a new subset.

For 2:

- **New label** : Input the name of new label.

- **Add Label** : Clicking after entering the new label name will add the label name into the choice list.

For 3:

- **Select tissue** : Select the tissue that the cell subset belongs to. User defined labels are in the "user defined" option.

- **Select cell type** : The cell types in the selected tissue will be in the choice list. 

- **Choose subset** : Select the cell subset you want to update.

- **Update Subset Labels** : Update the selected subset with the select cell type.


---

### Analysis

<img src="analysis.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

Users can do various analysis including differentially expressed genes (DE genes) analysis in this panel.

- **View gene expression** : Let the colors in the plot show the expression level of the selected gene. When "Clusters" is selected, the colors show different clusters.

- **Select number of genes** : Select the number of genes used in the analysis.

- **alpha** : 

- **Correction** : 

- **Choose subset 1&2** : Choose the subsets you want to do DE gene analysis. Selecting "None" means selecting all the cells but those in the other subset.

- **Run DE analysis** : Calculate DE genes of Subset1 (vs. Subset2)

- **Search Gene card** : Enter the gene you want to search.

- **Search Card** : Clicking the button will open the gene card website of the entered gene.

After finishing DE analysis, the DE panel of the bottom of the UI body will show the results. Other panels will also be activated for subsequent analysis.

- **DE** :

  <img src="DE.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **heatmap** :

  <img src="heatmap.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **Gene Ontology** :

  <img src="onto.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **KEGG** :

  <img src="KEGG.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **MsigDB C2** :

  <img src="msigdb.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **Cell Type** :

  <img src="celltype.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **Disease** :

  <img src="disease.png" alt="datasetMenu"
	title="UI" width="600" height="300" align=left/>

- **User Markers** :

  <img src="usermarker.png" alt="datasetMenu"
	title="UI" width="600" height="100" align=left/>
---

### Appearance
<img src="appearance.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

This panel enable the users to adjust the appearance of the main plot.

- **Select theme** : Change the UI body to Light/Dark mode. 

- **Select dot size** : Select the size of the dots in the main plot.

- **Choose subset** : Select whether showing the cluster names or only IDs in the plot.

- **Choose subset** : Select the height of the plot.

---

### Export/Import

<img src="export&import.png" alt="datasetMenu"
	title="UI" width="300" height="300" align=left/>

This panel enables users to do the following:
1. Exporting the current session and import later for subsequent analysis or sharing the analysis results
2. Importing session after loading the dataset
3. Downloading the main plot.
4. Downloading a csv file including cell IDs with corresponding labels.

For 1:
- **Export Session** : Download the current session.

For 2:
- **Import Session** : Upload a session and load it.

For 3:
- **Select format** : Select which format of the plot to download.
- **Download Plot** : Download the plot.

For 4:
- **Input Subset IDs** : Enter the subset IDs you want to download. For example 1,2,5 means downloading subsets 1,2 and 5.
- **Download Selected Subsets** : Download subsets entered above together with their labels as a csv file.

