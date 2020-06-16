tooltips <- function(id, label='tooltips') {ns = NS(id); list(
    # Dataset Menu

    # Clustering Menu
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
    #bsTooltip(ns("eval_method"),
    #          paste0("Used to determine the best number of clusters ",
    #                 "if number of clusters is a list.")),
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
    bsTooltip(ns("subset1_upd"),
              paste0("Update the name of this subset to the selected cell ",
                     "type.")),

    # Analysis
    #bsTooltip(ns("subset1"),
    #          paste0("Subset to run the analysis for against Subset 2.")),
    #bsTooltip(ns("subset2"),
    #          paste0("If None, will consider all cells not in Subset 1.")),
    #bsTooltip(ns("searchgene"),
    #          paste0("Gene name to search for.")),

    # Appearance
    #bsTooltip(ns(""), paste0()),

    # Export/Import Menu
    bsTooltip(ns("cell_subset_download"),
              paste0("Integer, list, or range specified by a dash, e.g., 1-7."))
)}
