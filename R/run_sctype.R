#' @title Run sctype analysis on Seurat object
#' @description run an automated cell type annotation
#' @details This function is compatible with seurat versions 4 and 5 -- different methods(@,$) to retrieve counts,data,scale.data
#' @param seurat_object A Seurat object
#' @param known_tissue_type The tissue type of the input data (optional)
#' @param custom_marker_file Path to the custom marker file (optional)
#' @param assay e.g. RNA, SCT, integrated (default="RNA")
#' @param scaled logical indicating whether the matrix is scaled (default=TRUE)
#' @param plot_umap logical indicating whether to plot the UMAP (default= FALSE)
#' @param plot_bubble logical indicating whether to plot the bubble plot (default=FALSE)
#' @param name The name of the metadata column to store the scType results (default="sctype_classification")
#' @return A modified copy of the input Seurat object with a new metadata column
#'
#' @importFrom Seurat DimPlot
#' @importFrom dplyr %>% top_n
#' @examples
#' \dontrun{
#' seurat_object <- run_scType(seurat_object, "Immune system")
#' }
#' @author @IanevskiAleksandr
#' @export
run_sctype <- function(seurat_object, known_tissue_type = NULL, assay = "RNA", scaled = TRUE, custom_marker_file = NULL, plot_umap = FALSE, plot_bubble = FALSE, name = "sctype_classification") {
    # Check for missing arguments
    if (is.null(seurat_object)) {
        stop("Argument 'seurat_object' is missing")
    }
    if (!inherits(seurat_object, "Seurat")) {
        stop("Argument 'seurat_object' must be a Seurat object")
    }
    # Set default custom marker file
    if (is.null(custom_marker_file)) {
        data(ScTypeDB_full)
        custom_marker_file <- ScTypeDB_full
    }
    # Auto-detect tissue type if not provided
    if (is.null(known_tissue_type)) {
        print("Guessing tissue type: \n")
        tissue_type <- auto_detect_tissue_type(
            path_to_db_file = custom_marker_file,
            seuratObject = seurat_object,
            scaled = scaled, assay = assay
        )
        rownames(tissue_type) <- NULL
        tissue_type <- tissue_type$tissue[1]
    } else {
        tissue_type <- known_tissue_type
    }

    if (is.null(custom_marker_file)) {
        data(ScTypeDB_full)
        custom_marker_file <- ScTypeDB_full
    } else if (!is.data.frame(custom_marker_file) & (is.character(custom_marker_file))) {
        if (file.exists(custom_marker_file)) {
            custom_marker_file <- openxlsx::read.xlsx(custom_marker_file)
        } else {
            stop("`custom_marker_file` has not a valid value")
        }
    }

    # Prepare gene sets
    gs_list <- gene_sets_prepare(custom_marker_file, tissue_type)

    data_type <- if (scaled) "scale.data" else "counts"
    package_type <- data_type %in% names(attributes(seurat_object[[assay]]))

    # Calculate scType scores
    if (package_type) {
        print("Using Seurat v4 object")
        es.max <- sctype_score(
            scRNAseqData = slot(seurat_object[[assay]], data_type),
            scaled = TRUE, gs = gs_list$gs_positive,
            gs2 = gs_list$gs_negative
        )
    } else {
        print("Using Seurat v5 object")

        if (data_type == "scale.data") {
            scRNAseqData <- seurat_object[[assay]]$scale.data
        } else {
            scRNAseqData <- seurat_object[[assay]]$counts
        }

        es.max <- sctype_score(
            scRNAseqData = as.matrix(scRNAseqData),
            scaled = TRUE, gs = gs_list$gs_positive,
            gs2 = gs_list$gs_negative
        )
    }

    # Extract top cell types for each cluster
    cL_results <- do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl) {
        es.max.cl <- sort(rowSums(es.max[, rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters == cl)), 10)
    }))
    sctype_scores <- cL_results %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n = 1, wt = scores)

    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"

    seurat_object_res <- seurat_object
    seurat_object_res@meta.data[name] <- ""
    for (j in unique(sctype_scores$cluster)) {
        cl_type <- sctype_scores[sctype_scores$cluster == j, ]
        seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j, name] <- as.character(cl_type$type[1])
    }
    if (plot_umap) {
        plot_ <- Seurat::DimPlot(seurat_object_res, reduction = "umap", group.by = name)
        print(plot_)
    }
    if (plot_bubble) {
        bubbleplot_ <- make_bubbleplot(cL_results, sctype_scores, custom_marker_file)
        print(bubbleplot_)
    }
    text_ <- paste("New metadata added: ", name)
    print(text_)
    return(seurat_object_res)
}
