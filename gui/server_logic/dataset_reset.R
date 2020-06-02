dataset_reset <- function(input, output, session, reset, setNames, setPts,
                          labelList, deButtons, deGenes, pipe, newLabels) {
    observe({ if (reset() > 0) {
        print("Resetting")
        #converter = Con()
        updateSelectInput(
            session = session,
            inputId = "color",
            choices = c("Clusters", as.character(as.character(pipe()$col_names))),
            selected = "Clusters")
      

        setNames(c("None"))
        setPts(c(NA))
        labelList(c())
        newLabels(NULL)

        # Update sets
        for (i in 1:length(pipe()$n_clusters)) {
            setNames(c(setNames(),
                       (paste("Cluster_", as.character(i - 1), sep = ""))))
            setPts(c(setPts(), list(which(pipe()$labels == (i - 1)))))
        }

        n_clusters = pipe()$n_clusters
        updateSelectInput(
            session,
            "newlabels",
            choices = n_clusters
        )

        # Clear all analysis tabs
        if (length(deButtons()) > 0)
            for (i in 1:length(deButtons()))
                deButtons()[[i]]$destroy()
        deButtons(c())
        deGenes(c())
        output$DEbuttons = NULL
        output$genes = NULL
        output$GeneOntology = NULL
        output$KEGG = NULL
        output$Markers = NULL
        output$Msigdb = NULL

        reset(0)
    }})
}
