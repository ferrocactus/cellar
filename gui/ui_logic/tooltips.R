tooltips <- function(id, label='tooltips') {ns = NS(id); list(
    # Dataset Menu

    # Clustering Menu
    bsTooltip(ns("runconfigbtn"),
              paste0("Do clustering with current configuration")),


    # Label Transfer
    bsTooltip(ns("uploaded_dataset_align"),
              paste0("Choose the reference dataset. Should be the same ",
                     "as the one on the session file.")),
    bsTooltip(ns("server_dataset_align"),
              paste0("Choose the reference dataset. Should be the same ",
                     "as the one on the session file.")),

    # Selection & Labeling
    bsTooltip(ns("subset1_upd"),
              paste0("Select the subset you want to update")),

    # Analysis

    bsTooltip(ns("subset1"),
                "If None, will consider all cells not in Subset 1.",
                placement= "bottom"
    ),
    bsTooltip(ns("subset2"),
                "Subset to run the analysis for against Subset 2.",
                placement= "bottom"
    ),

    # Appearance




    # Export/Import Menu
    bsTooltip(ns("download_sess"),
              paste0("Download the current session for sharing or working later")),
    bsTooltip(ns("download_plot"),
              paste0("Download the main plot in the selected format")),
    bsTooltip(ns("download_cells"),
              paste0("Download the selected subsets as a csv file."))

)}
