# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# auto_detect_tissue_type: automatically detect a tissue type of the dataset

#' Autotomatically detect a tissue of the dataset
#' @param celltype_markers_db (dataframe) DB file with cell types
#' @param seuratObject Seurat Object from wich to extract the input scRNA-seq matrix (rownames - genes, column names - cells),
#' @param scaled logical indicating whether the matrix is scaled (default=TRUE)
#' @param cluster_column column in metadata containing the clusters (default='seurat_clusters')
#' @param assay e.g. RNA, SCT, integrated
#' @param plot logical indicating whether to plot (default=FALSE)
#' @return list with two keys containing the postive and negative markers resp.
#' @author @IanevskiAleksandr
#' @importFrom dplyr %>% top_n
#' @export
auto_detect_tissue_type <- function(celltype_markers_db, seuratObject, scaled, cluster_column = "seurat_clusters", assay = "RNA", plot = FALSE, ...) {
    if (!is.logical(scaled)) {
        stop("Argument 'scaled' needs to be a logical/boolean")
    }
    if (!is.logical(plot)) {
        stop("Argument 'plot' needs to be a logical/boolean")
    }
    # get all tissue types in DB
    tissues_ <- unique(celltype_markers_db$tissueType)
    result_ <- c()

    for (tissue in tissues_) {
        print(paste0("Checking...", tissue))

        # prepare gene sets
        gs_list <- gene_sets_prepare(celltype_markers_db, tissue)

        # check Seurat version
        package_type <- substr(packageVersion("Seurat"), 1, 1)
        obj <- if (package_type == "5") {
            as.matrix(seuratObject[[assay]]$data_type)
        } else {
            as.matrix(seuratObject[[assay]]@data_type)
        }

        es.max <- sctype_score(
            scRNAseqData = obj, scaled = scaled, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative,
            marker_sensitivity = gs_list$marker_sensitivity, verbose = !0
        )


        cL_results <- do.call("rbind", lapply(unique(seuratObject@meta.data[[cluster_column]]), function(cl) {
            es.max.cl <- sort(rowSums(es.max[, rownames(seuratObject@meta.data[seuratObject@meta.data[[cluster_column]] == cl, ])]), decreasing = !0)
            head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
        }))

        dt_out <- cL_results %>%
            dplyr::group_by(cluster) %>%
            dplyr::top_n(n = 1)

        # return mean score for tissue
        result_ <- rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
    }

    # order by mean score
    result_ <- result_[order(-result_$score), ]

    # TODO move this to a separate function
    if (plot) {
        barplot(
            height = result_$score, names = result_$tissue, col = rgb(0.8, 0.1, 0.1, 0.6),
            xlab = "Tissue", ylab = "Summary score", main = "The higher summary score, the more likely tissue type is"
        )
    }
    return(result_)
}
