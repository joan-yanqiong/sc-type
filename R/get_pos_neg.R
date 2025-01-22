#' Extract postive and negative markers from DE results
#' @param ref_markers dataframe with reference markers
#' @param cl celltype
#' @param min_pct_diff numeric
#' @importFrom dplyr %>% filter
#' @return list of positive and negative markers
#' @export
get_pos_neg <- function(ref_markers, cl, min_pct_diff) {
    cluster_markers <- ref_markers %>% dplyr::filter(cluster %in% cl)
    cluster_markers <- cluster_markers[abs(cluster_markers$pct.1 - cluster_markers$pct.2) > min_pct_diff, ]
    data_EV <- EnhancedVolcano::EnhancedVolcano(cluster_markers, rownames(cluster_markers), x = "avg_log2FC", y = "p_val_adj")
    pos_logFC <- data_EV$data %>% dplyr::filter(Sig == "FC_P" & avg_log2FC > 0)
    neg_logFC <- data_EV$data %>% dplyr::filter(Sig == "FC_P" & avg_log2FC < 0)
    return(list(pos_logFC, neg_logFC))
}
