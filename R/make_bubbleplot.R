#' Make bubbleplot
#' @param cL_results dataframe
#' @param sctype_scores scType scores
#' @author @IanevskiAleksandr
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggplot2 theme_void
#' @import ggraph
#' @export
make_bubbleplot <- function(
    cL_results,
    sctype_scores,
    custom_marker_file) {
    # prepare edges
    cL_results <- cL_results[order(cL_results$cluster), ]
    edges <- cL_results
    edges$type <- paste0(edges$type, "_", edges$cluster)
    edges$cluster <- paste0("cluster ", edges$cluster)
    edges <- edges[, c("cluster", "type")]
    colnames(edges) <- c("from", "to")
    rownames(edges) <- NULL

    # prepare nodes
    nodes_lvl1 <- sctype_scores[, c("cluster", "ncells")]
    nodes_lvl1$cluster <- paste0("cluster ", nodes_lvl1$cluster)
    nodes_lvl1$Colour <- "#f1f1ef"
    nodes_lvl1$ord <- 1
    nodes_lvl1$realname <- nodes_lvl1$cluster
    nodes_lvl1 <- as.data.frame(nodes_lvl1)
    nodes_lvl2 <- c()

    color_palette <- c("#5f75ae", "#92bbb8", "#64a841", "#e5486e", "#de8e06", "#eccf5a", "#b5aa0f", "#e4b680", "#7ba39d", "#b15928", "#ffff99", "#6a3d9a", "#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")

    for (i in seq_len(length(unique(cL_results$cluster)))) {
        dt_tmp <- cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]
        nodes_lvl2 <- rbind(
            nodes_lvl2,
            data.frame(
                cluster = paste0(dt_tmp$type, "_", dt_tmp$cluster),
                ncells = dt_tmp$scores, Colour = color_palette[i],
                ord = 2, realname = dt_tmp$type
            )
        )
    }
    nodes <- rbind(nodes_lvl1, nodes_lvl2)
    nodes$ncells[nodes$ncells < 1] <- 1
    files_db <- custom_marker_file[, c("cellName", "shortName")]
    files_db <- unique(files_db)
    nodes <- merge(nodes,
        files_db,
        all.x = TRUE, all.y = FALSE,
        by.x = "realname",
        by.y = "cellName", sort = FALSE
    )
    nodes$shortName[is.na(nodes$shortName)] <- nodes$realname[is.na(nodes$shortName)]
    nodes <- nodes[, c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

    mygraph <- igraph::graph_from_data_frame(edges, vertices = nodes)

    # Make the graph
    gggr <- ggraph::ggraph(mygraph, layout = "circlepack", weight = I(ncells)) +
        ggraph::geom_node_circle(aes(filter = ord == 1, fill = I("#F5F5F5"), colour = I("#D3D3D3")), alpha = 0.9) +
        ggraph::geom_node_circle(aes(filter = ord == 2, fill = I(Colour), colour = I("#D3D3D3")), alpha = 0.9) +
        ggplot2::theme_void() +
        ggraph::geom_node_text(aes(
            filter = ord == 2,
            label = shortName, colour = I("#ffffff"), fill = "white", repel = !1, parse = TRUE, size = I(log(ncells, 25) * 1.5)
        )) +
        ggraph::geom_node_label(aes(
            filter = ord == 1, label = shortName, colour = I("#000000"),
            size = I(3), fill = "white", parse = TRUE
        ), repel = !0, segment.linetype = "dotted")

    return(gggr)
}
