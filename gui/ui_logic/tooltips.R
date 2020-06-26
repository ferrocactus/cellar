tooltips <- function(id, label='tooltips') {ns = NS(id); list(
    # Dataset Menu
    bsTooltip(ns("folder"),
              paste0("Select the datasets on the server or the datasets uploaded by the useR")),
    bsTooltip(ns("uploaded_dataset"),
              paste0("Choose a specific dataset")),
    bsTooltip(ns("file1"),
              paste0("Upload a dataset")),
    # Clustering Menu
    bsTooltip(ns("runconfigbtn"),
              paste0("Do clustering with current configuration")),
    bsTooltip(ns("dim_method"),
              paste0("Choose a dimension reduction method")),
    bsTooltip(ns("clu_method"),
              paste0("Choose a clustering method")),
    bsTooltip(ns("clu_method"),
              paste0("Choose a clustering method")),
    bsTooltip(ns("dim_options"),
              paste0("If automatic, will try and detect the optimal ",
              "number of components using the ankle heuristic of finding the ",
              "ankle in the explained variance graph.")),
    bsTooltip(ns("dim_n_components"),
              paste0("Integer value.")),
    bsTooltip(ns("clu_n_clusters"),
              paste0("Can be a single integer, ",
              "a list of comma separated values, ",
              "or a tuple specifying a range, e.g., (3, 9, 2) will start ",
              "at 3 and finish at 9 on increments of 2.")),
    bsTooltip(ns("eval_method"),
              paste0("Used to determine the best number of clusters ",
                     "if number of clusters is a list.")),
    bsTooltip(ns("ssc_method"),
              paste0("Semi-supervised clustering takes into account the ",
                     "existing labels and tries to refine the clusters.")),
    bsTooltip(ns("saved_clusters"),
              paste0("Clusters whose label should not be changed by ",
                     "the algorithm. Can be integers ",
                     "or a range specified by a dash, e.g., 1-5.")),
    bsTooltip(ns("clusters_to_merge"),
              paste0("Integer or range specified by a dash, e.g., 1-5.")),

    # Label Transfer
    bsTooltip(ns("uploaded_dataset_align"),
              paste0("Choose the reference dataset. Should be the same ",
                     "as the one on the session file.")),
    bsTooltip(ns("hubmap_dataset_align"),
              paste0("Choose the reference dataset. Should be the same ",
                     "as the one on the session file.")),

    # Selection & Labeling
    bsTooltip(ns("newsubset"),
              paste0("Store the selected points into a subset with this name.")),
    bsTooltip(ns("newlabelbox"),
              paste0("Enter the new label name you want to add")),
    bsTooltip(ns("tissue"),
              paste0("Select the tissue of the cell you want to label")),
    bsTooltip(ns("newlabels"),
              paste0("Select the cell type you want to label")),
    bsTooltip(ns("subset1_upd"),
              paste0("Update the name of this subset to the selected cell ")),

    # Analysis
    bsTooltip(ns("subset1"),
              paste0("Subset to run the analysis for against Subset 2.")),
    bsTooltip(ns("subset2"),
              paste0("If None, will consider all cells not in Subset 1.")),
    bsTooltip(ns("searchgene"),
             paste0("Gene name to search for.")),



    # Appearance

    bsTooltip(ns("dot_size"), paste0("Select the size of the dots in the plot")),
    bsTooltip(ns("show_names"), paste0("show cluster names in the plot")),

    # Export/Import Menu
    bsTooltip(ns("download_sess"),
              paste0("Download the current session for sharing or working later")),
    bsTooltip(ns("upload_sess"),
              paste0("Load a local session")),
    bsTooltip(ns("plot_download_format"),
              paste0("Select the format of the plot to be downloaded")),
    bsTooltip(ns("cell_subset_download"),
              paste0("Integer, list, or range specified by a dash, e.g., 1-7."))
)}
