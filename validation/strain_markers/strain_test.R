# SwapSeq 1.2 test on your own computer: do strain-level fingerprints catch swaps that species tables miss?
#
# Runs on your Mac (it downloads public data from Bioconductor that the cloud session cannot reach).
# Data: curatedMetagenomicData, shotgun studies that sampled the same people more than once,
# including infants (where species tables fail).
#   species : MetaPhlAn relative abundance  -> what SwapSeq uses today
#   markers : MetaPhlAn marker presence     -> strain-level fingerprint (swapseqr 0.3.0, type = "markers")
#
# Time: 30-90 minutes the first time (downloads), much less after that. Needs about 8 GB of free memory.
# Output: a short table at the end. Copy only that table back to the chat.

options(timeout = 3600)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("curatedMetagenomicData", "SummarizedExperiment")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
}
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("swapseqr", quietly = TRUE) ||
    utils::packageVersion("swapseqr") < "0.3.0") {
  remotes::install_github("bymankirat/swapseqr", upgrade = "never")
}
suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(swapseqr)
})

wanted <- c("HMP_2019_ibdmdb", "VatanenT_2016", "KosticAD_2015", "YassourM_2016",
            "BackhedF_2015", "ShaoY_2019")
have <- unique(sampleMetadata$study_name)
studies <- intersect(wanted, have)
message("Studies found: ", paste(studies, collapse = ", "))

get_table <- function(study, type) {
  x <- curatedMetagenomicData(paste0(study, ".", type), dryrun = FALSE,
                              counts = (type == "relative_abundance"))[[1]]
  list(m = as.matrix(assay(x)), meta = as.data.frame(colData(x)))
}

one_line <- function(study, kind, r, secs) {
  p <- attr(r, "power")
  pct <- function(v) if (is.null(v) || is.na(v)) "NA" else sprintf("%.0f%%", 100 * v)
  data.frame(study = study, data = kind, samples = nrow(r),
             strong_flags = sum(r$level == "flag"), check_list = sum(r$level == "check"),
             swaps_caught = pct(p$swap), with_check_list = pct(p$swapWithCheck),
             false_flags_per_1000 = if (is.null(p$falseFlagsPer1000)) NA else round(p$falseFlagsPer1000, 2),
             seconds = round(secs))
}

out <- list()
for (study in studies) {
  message("\n=== ", study, " ===")
  sp <- tryCatch(get_table(study, "relative_abundance"), error = function(e) NULL)
  mk <- tryCatch(get_table(study, "marker_presence"), error = function(e) NULL)
  if (is.null(sp) || is.null(mk)) { message("  could not download, skipped"); next }
  ids <- intersect(colnames(sp$m), colnames(mk$m))
  subj <- sp$meta[ids, "subject_id"]
  keep <- ids[!is.na(subj) & subj %in% names(which(table(subj) >= 3))]
  # keep the run time sensible on very large studies: at most 1,200 samples, whole people only
  if (length(keep) > 1200) {
    set.seed(1)
    people <- unique(sp$meta[keep, "subject_id"])
    people <- sample(people)
    n_per <- table(sp$meta[keep, "subject_id"])[people]
    people <- people[cumsum(n_per) <= 1200]
    keep <- keep[sp$meta[keep, "subject_id"] %in% people]
  }
  if (length(keep) < 30) { message("  too few repeat samples, skipped"); next }
  subj <- sp$meta[keep, "subject_id"]
  message("  ", length(keep), " samples from ", length(unique(subj)), " people")

  species <- sp$m[, keep]
  species <- species[rowSums(species > 0) > 0, , drop = FALSE]
  markers <- mk$m[, keep] > 0
  prev <- rowSums(markers)
  markers <- markers[prev >= 2 & prev <= 0.5 * ncol(markers), , drop = FALSE]
  if (nrow(markers) > 40000) { set.seed(1); markers <- markers[sort(sample(nrow(markers), 40000)), ] }
  storage.mode(markers) <- "integer"
  message("  ", nrow(markers), " informative markers")

  t1 <- system.time(r_sp <- suppressMessages(check_swaps(species, subj)))["elapsed"]
  t2 <- system.time(r_mk <- suppressMessages(check_swaps(markers, subj, type = "markers")))["elapsed"]
  out[[length(out) + 1]] <- one_line(study, "species", r_sp, t1)
  out[[length(out) + 1]] <- one_line(study, "markers", r_mk, t2)
  rm(sp, mk, species, markers); gc()
}

cat("\n==== COPY FROM HERE ====\n")
print(do.call(rbind, out), row.names = FALSE)
cat("swapseqr", as.character(utils::packageVersion("swapseqr")), "\n")
cat("==== TO HERE ====\n")
