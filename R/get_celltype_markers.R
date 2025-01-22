#' Retrieve markers
#' @param ref_markers dataframe with reference markers (results from DE)
#' @param tissue tissue type
#' @param marker_per_cluster no. marker genes per cluster to select
#' @param min_pct_diff min. percentage difference (default=0.3)
#' @export
retrieve_markers <- function(ref_markers, tissue, marker_per_cluster = NULL, min_pct_diff = 0.3) {
    ScTypeDB <- data.frame(
        tissueType = rep(tissue, length(unique(ref_markers$cluster))),
        cellName = unique(ref_markers$cluster),
        geneSymbolmore1 = "",
        geneSymbolmore2 = ""
    )

    rownames(ScTypeDB) <- ScTypeDB$cellName
    for (cl in as.character(unique(ref_markers$cluster))) {
        df_ <- get_pos_neg(ref_markers, cl, min_pct_diff = min_pct_diff)
        df_pos <- df_[[1]][order(df_[[1]]$p_val_adj, (df_[[1]]$avg_log2FC) * -1), ]
        df_neg <- df_[[2]][order(df_[[2]]$p_val_adj, (df_[[2]]$avg_log2FC)), ]

        if (is.null(marker_per_cluster)) {
            positive_markers <- paste0(df_pos$gene, collapse = ",")
            negative_markers <- paste0(df_neg$gene, collapse = ",")
        } else {
            positive_markers <- paste0(df_pos$gene[seq_len(marker_per_cluster)],
                collapse = ","
            )
            negative_markers <- paste0(df_neg$gene[seq_len(marker_per_cluster)],
                collapse = ","
            )
        }
        ScTypeDB[cl, "geneSymbolmore1"] <- positive_markers
        ScTypeDB[cl, "geneSymbolmore2"] <- negative_markers
    }
    rownames(ScTypeDB) <- NULL

    return(ScTypeDB)
}
