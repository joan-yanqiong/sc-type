# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type excel file

#' Prepare gene sets
#' @param celltype_markers_db (dataframe) DB file with cell types
#' @param cell_type cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#' @return list with two keys containing the postive and negative markers resp.
#' @author @IanevskiAleksandr
#' @importFrom HGNChelper checkGeneSymbols
#' @export
gene_sets_prepare <- function(celltype_markers_db, cell_type) {
    celltype_markers_db <- celltype_markers_db[celltype_markers_db$tissueType == cell_type, ]
    celltype_markers_db$geneSymbolmore1 <- gsub(" ", "", celltype_markers_db$geneSymbolmore1)
    celltype_markers_db$geneSymbolmore2 <- gsub(" ", "", celltype_markers_db$geneSymbolmore2)

    # correct gene symbols from the given DB (up-genes)
    celltype_markers_db$geneSymbolmore1 <- sapply(seq_len(nrow(celltype_markers_db)), function(i) {
        markers_all <- gsub(" ", "", unlist(strsplit(celltype_markers_db$geneSymbolmore1[i], ",")))
        markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
        markers_all <- sort(markers_all)

        if (length(markers_all) > 0) {
            suppressMessages({
                markers_all <- unique(na.omit(HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol))
            })
            paste0(markers_all, collapse = ",")
        } else {
            ""
        }
    })

    # correct gene symbols from the given DB (down-genes)
    celltype_markers_db$geneSymbolmore2 <- sapply(seq_len(nrow(celltype_markers_db)), function(i) {
        markers_all <- gsub(" ", "", unlist(strsplit(celltype_markers_db$geneSymbolmore2[i], ",")))
        markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
        markers_all <- sort(markers_all)

        if (length(markers_all) > 0) {
            suppressMessages({
                markers_all <- unique(na.omit(HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol))
            })
            paste0(markers_all, collapse = ",")
        } else {
            ""
        }
    })

    celltype_markers_db$geneSymbolmore1 <- gsub("///", ",", celltype_markers_db$geneSymbolmore1)
    celltype_markers_db$geneSymbolmore1 <- gsub(" ", "", celltype_markers_db$geneSymbolmore1)
    celltype_markers_db$geneSymbolmore2 <- gsub("///", ",", celltype_markers_db$geneSymbolmore2)
    celltype_markers_db$geneSymbolmore2 <- gsub(" ", "", celltype_markers_db$geneSymbolmore2)

    gs <- lapply(seq_len(nrow(celltype_markers_db)), function(j) gsub(" ", "", unlist(strsplit(toString(celltype_markers_db$geneSymbolmore1[j]), ","))))
    names(gs) <- celltype_markers_db$cellName
    gs2 <- lapply(seq_len(nrow(celltype_markers_db)), function(j) gsub(" ", "", unlist(strsplit(toString(celltype_markers_db$geneSymbolmore2[j]), ","))))
    names(gs2) <- celltype_markers_db$cellName

    list(gs_positive = gs, gs_negative = gs2)
}
