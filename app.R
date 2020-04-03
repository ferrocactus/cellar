reticulate::virtualenv_create(envname = "renv", python = "/usr/bin/python3")
reticulate::virtualenv_install('renv', packages=c(
       'matplotlib',
       'pandas',
       'numpy',
       'seaborn',
       'scikit-learn',
       'umap-learn',
       'tqdm',
       'kneed',
       'gtfparse',
       'statsmodels',
       'joblib',
       'gseapy',
       'anndata'))

source("gui/ui.R")
source("gui/server.R")

shinyApp(ui, server)
