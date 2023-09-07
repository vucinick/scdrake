## -- Common functions related to gene annotation.

#' @title Collapse a character vector by `","`.
#' @description Used as `multiVals` parameter in [AnnotationDbi::mapIds()].
#' @param x A character vector.
#' @return A character scalar.
#'
#' @concept sc_annotation
#' @export
collapse_ensembl_multivals <- function(x) {
  str_comma(x, space = FALSE)
}

#' @title Create a `dataframe` with annotation of SCE object's genes using an `AnnotationDbi` package.
#' @description
#' The following columns are retrieved from the `AnnotationDbi` package: `SYMBOL`, `GENENAME`, `ENTREZID`.
#' Multi-mapping entities are collapsed by "," (i.e. a single ENSEMBL ID having multiple symbols, etc.).
#' Unknown symbols are replaced by the corresponding ENSEMBL IDs.
#' @param sce A `SingleCellExperiment` object.
#' @param annotation_db_file A character scalar: path to SQLite file of annotation DB, e.g. `org.Hs.eg.db$conn@dbname`.
#' @return A `dataframe` with annotation.
#'
#' @concept sc_annotation
#' @export
make_gene_annotation <- function(sce, annotation_db_file, db_NCBI = TRUE, organism = NULL, genome = NULL) {

  if (db_NCBI){
    genome_ann <- AnnotationDbi::loadDb(annotation_db_file)

    gene_annotation <- data.frame(
      ENSEMBL = rownames(sce),
      row.names = rownames(sce),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        SYMBOL = AnnotationDbi::mapIds(
          genome_ann,
          keys = .data$ENSEMBL,
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = collapse_ensembl_multivals
        ),
        SYMBOL = dplyr::if_else(is.na(SYMBOL), ENSEMBL, SYMBOL),
        GENENAME = AnnotationDbi::mapIds(
          genome_ann,
          keys = .data$ENSEMBL,
          column = "GENENAME",
          keytype = "ENSEMBL",
          multiVals = collapse_ensembl_multivals
        ),
        ENTREZID = AnnotationDbi::mapIds(
          genome_ann,
          keys = .data$ENSEMBL,
          column = "ENTREZID",
          keytype = "ENSEMBL",
          multiVals = collapse_ensembl_multivals
        )
      )
  }else{
    organism = gsub("_", " ", organism)
    AnnHub <- AnnotationHub()
    db_len <- length(query(AnnHub, c("EnsDB", genome)))
    if (db_len == 0){
      possible_genomes <- unique(mcols(query(AnnHub, c("EnsDB",  organism)))[,"genome"])
      cli::cli_alert_warning("The genome doesn't exist or the AnnHub doesn't have it.") # nolint
      cli::cli_alert_info(cat("Try: ", paste0(head(possible_genomes, -1), ", "), tail(possible_genomes, 1)))
    }else{
      # Take the last one (the newest version)
      genome_ann <- query(AnnHub, c("EnsDB", genome = genome))[[db_len]]
    }

    gene_annotation <- data.frame(
      ENSEMBL = rownames(sce),
      row.names = rownames(sce),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        SYMBOL = AnnotationDbi::mapIds(
          genome_ann,
          keys = ENSEMBL,
          column = "SYMBOL",
          keytype = "GENEID",
          # custom functions aren't implemented yet for EnsDB mapIds
          multiVals = "first"
        ),
        SYMBOL = dplyr::if_else(is.na(SYMBOL), ENSEMBL, SYMBOL),
        GENENAME = AnnotationDbi::mapIds(
          genome_ann,
          keys = ENSEMBL,
          column = "GENENAME",
          keytype = "GENEID",
          multiVals = "first"
        ),
        ENTREZID = AnnotationDbi::mapIds(
          genome_ann,
          keys = ENSEMBL,
          column = "ENTREZID",
          keytype = "GENEID",
          multiVals = "first"
        )
      )
  }

  gene_annotation <- gene_annotation %>%
    dplyr::select(ENSEMBL, SYMBOL, ENTREZID, dplyr::everything()) %>%
    as.data.frame() %>%
    set_rownames(.$ENSEMBL)

  duplicated_gene_symbols <- duplicated(gene_annotation$SYMBOL) | duplicated(gene_annotation$SYMBOL, fromLast = TRUE)

  gene_annotation[duplicated_gene_symbols, "SYMBOL"] <- gluec(
    "{gene_annotation[duplicated_gene_symbols, 'SYMBOL']}_{rownames(gene_annotation[duplicated_gene_symbols, ])}"
  )

  return(gene_annotation)
}

#' @title Load a SQL database file and run a function from the `AnnotationDbi` package.
#' @description An `AnnotationDbi` object cannot be used in parallel, so this is a workaround.
#' @param annotation_db_file A character scalar: path to SQLite file of annotation DB, e.g. `org.Hs.eg.db$conn@dbname`.
#' @param dbi_fun A function from the `AnnotationDbi` package which accepts `AnnotationDbi` object as parameter `x`,
#'   e.g. [AnnotationDbi::mapIds()].
#' @param ... Parameters passed to `dbi_fun`.
#' @return An output from `dbi_fun`.
#'
#' @concept sc_annotation
#' @export
with_dbi <- function(annotation_db_file, dbi_fun, ...) {
  dbi <- AnnotationDbi::loadDb(annotation_db_file)
  rlang::exec(dbi_fun, x = dbi, !!!list(...))
}
