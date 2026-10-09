#!/usr/bin/env Rscript
# =============================================================================
# Community-level analysis of MAG annotations: rhizosphere vs bulk soil,
# biochar vs control, across three sampling time scales.
#
#   seasonal : fortnightly dates, 17 May - 14 Aug 2024 (7 dates)
#   daily    : intensive week, days 1, 2 and 4 at 04:00 and 14:00
#   diel     : intensive week, day 3 every 4 h (02, 06, 10, 14, 18, 22 h)
#
# Annotation sets (each analysed separately, own output folders):
#   microtrait : fraction of community carrying each microTrait trait (granularity 3),
#                + abundance-weighted log10 min. generation time and optimum temperature
#   portraits  : abundance-weighted mean predicted probability per tool|trait[|model]
#   amr        : fraction of community carrying >= 1 RGI hit per drug class (+ any ARG)
#   bacmet     : fraction of community carrying each BacMet gene (+ any BacMet hit)
#   bgc        : abundance-weighted mean number of complete antiSMASH regions per class
#                (+ all complete regions)
#   gcf        : BiG-SLiCE GCF composition (Bray-Curtis; PERMANOVA + PCoA only)
#
# Weighting: for each sample, value = sum_i(a_i * x_i) / sum_i(a_i) over MAGs with a
# non-missing x_i, where a_i is the CoverM relative abundance (covered fraction >= 0.1).
#
# Statistics per set and time scale:
#   - PERMANOVA (vegan::adonis2), permutations restricted within plot; the treatment
#     main effect is plot-level and cannot be tested this way (reported NA)
#   - per-feature linear mixed models (nlme::lme, random = ~1 | plot), BH-adjusted
#     across features within each set, scale and term; diel time as a 24 h cycle
#
# Outputs:
#   <results>/analysis_tables/<set>/   cwm_by_sample.tsv, permanova_<scale>.tsv, lmm_<scale>.tsv
#   <results>/analysis_tables/summary_terms.tsv   significant features per set/scale/term
#   <results>/figures/<set>/<scale>_{ordination,trajectories,compartment_effects,heatmap}.pdf
#   <results>/figures/summary/significant_terms.pdf
#   <results>/figures/community_analysis_report.pdf   report: methods, log, summary, all figures, full code
#
# Runs in the pipeline's microTrait conda env (tidyverse, vegan, nlme, pheatmap).
# Usage:
#   Rscript community_analysis.R --results-dir <annotations> --metadata <sheet.tsv>
#       [--sets microtrait,portraits,amr,bacmet,bgc,gcf] [--nperm 999] [--min-mags 20]
#       [--granularity 3] [--abund coverm/coverm_rel_abund_masked_cf0.1.tsv]
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(vegan); library(nlme); library(pheatmap); library(grid)
})

# ---------------------------------------------------------------- arguments
args <- commandArgs(trailingOnly = TRUE)
opt <- list(sets = "microtrait,portraits,amr,bacmet,bgc,gcf", nperm = "999", `min-mags` = "20",
            granularity = "3", abund = "coverm/coverm_rel_abund_masked_cf0.1.tsv",
            `results-dir` = NA, metadata = NA)
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
if (is.na(opt$`results-dir`) || is.na(opt$metadata))
  stop("Usage: Rscript community_analysis.R --results-dir <annotations> --metadata <sheet.tsv>")
R <- opt$`results-dir`
SETS <- strsplit(opt$sets, ",")[[1]]
NPERM <- as.integer(opt$nperm); MIN_MAGS <- as.integer(opt$`min-mags`)
TAB <- file.path(R, "analysis_tables"); FIG <- file.path(R, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
logf <- file(file.path(TAB, "community_analysis.log"), open = "wt")
say <- function(...) { m <- paste0(...); message(m); writeLines(m, logf); flush(logf) }
say("community analysis -- ", format(Sys.time()), " -- sets: ", paste(SETS, collapse = ", "))
REPORT <- list()                                 # figures collected for the PDF report, in order
inp <- function(rel) { f <- file.path(R, rel); if (file.exists(f)) f else { say("  missing input: ", f); NULL } }

# Known date error in the sequencing sheet: 07/12/2024 is 12 July 2024 (no December sampling)
DATE_FIXES <- c("07/12/2024" = "12/07/2024")
SCALES <- c("seasonal", "daily", "diel")

# ---------------------------------------------------------------- metadata
read_sheet <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "latin1")
  hdr <- gsub(" ", "", trimws(strsplit(lines[1], "\t", fixed = TRUE)[[1]]))
  rows <- strsplit(lines[-1][nzchar(trimws(lines[-1]))], "\t", fixed = TRUE)
  if (max(lengths(rows)) == length(hdr) + 1 && "dna_plate_numdna_plate_col" %in% hdr) {
    i <- match("dna_plate_numdna_plate_col", hdr)
    hdr <- append(hdr[-i], c("dna_plate_num", "dna_plate_col"), after = i - 1)
    say("metadata: split fused header 'dna_plate_numdna_plate_col'")
  }
  rows <- lapply(rows, function(r) { r <- c(r, rep("", length(hdr))); trimws(r[seq_along(hdr)]) })
  setNames(as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE), hdr)
}
meta <- read_sheet(opt$metadata) %>%
  filter(dna_tube_num != "") %>%
  mutate(date_raw = date,
         date = as.Date(ifelse(date %in% names(DATE_FIXES), DATE_FIXES[date], date), format = "%d/%m/%Y"),
         sample = dna_tube_num,
         compartment = factor(recode(bulk_rhizo_plant, bulk = "bulk", rhizo = "rhizosphere"),
                              levels = c("bulk", "rhizosphere")),
         treatment = factor(recode(treatment, bio = "biochar", con = "control"), levels = c("control", "biochar")),
         plot = factor(plot),
         day = suppressWarnings(as.integer(day)),
         hour = suppressWarnings(as.integer(time) %/% 100),
         scale = case_when(is.na(day) ~ "seasonal", day == 3 ~ "diel",
                           hour %in% c(4, 14) ~ "daily", TRUE ~ NA_character_))
n_fixed <- sum(meta$date_raw %in% names(DATE_FIXES))
if (n_fixed) say("metadata: corrected date 07/12/2024 -> 12/07/2024 for ", n_fixed, " samples")
say("metadata: ", nrow(meta), " DNA samples; ",
    paste(names(table(meta$scale)), table(meta$scale), sep = "=", collapse = ", "))

# ---------------------------------------------------------------- abundances
ab <- read.delim(file.path(R, opt$abund), check.names = FALSE, row.names = 1)
keep <- intersect(colnames(ab), meta$sample)
say("abundance: ", ncol(ab), " samples in matrix; kept ", length(keep), "; excluded (not in sequencing sheet): ",
    paste(setdiff(colnames(ab), meta$sample), collapse = ", "))
ab <- as.matrix(ab[, keep, drop = FALSE]); ab[is.na(ab)] <- 0
empty <- colnames(ab)[colSums(ab) <= 0]
if (length(empty)) { say("abundance: no MAG coverage, dropped: ", paste(empty, collapse = ", "))
                     ab <- ab[, !colnames(ab) %in% empty, drop = FALSE] }
missing <- setdiff(meta$sample, colnames(ab))
say("abundance: sequencing-sheet samples without usable CoverM data: ", length(missing),
    if (length(missing)) paste0(" (", paste(missing, collapse = ", "), ")") else "")
sample_meta <- meta %>% filter(sample %in% colnames(ab)) %>%
  select(sample, date, day, hour, plot, treatment, compartment, scale)

# ---------------------------------------------------------------- helpers
# abundance-weighted mean of a MAG x feature matrix, ignoring missing values per feature
weighted_profile <- function(Fm) {
  mags <- intersect(rownames(Fm), rownames(ab))
  Fm <- Fm[mags, , drop = FALSE]; A <- ab[mags, , drop = FALSE]
  ok <- !is.na(Fm); F0 <- Fm; F0[!ok] <- 0
  num <- t(A) %*% F0; den <- t(A) %*% (ok * 1)
  out <- num / den; out[den == 0] <- NA
  out
}
presence_by_group <- function(df, group_col) {   # long (Genome, group) -> MAG x group presence (all MAGs)
  df <- df %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] != "") %>% distinct(Genome, g = .data[[group_col]])
  m <- matrix(0, nrow(ab), length(unique(df$g)), dimnames = list(rownames(ab), sort(unique(df$g))))
  df <- df %>% filter(Genome %in% rownames(m))
  m[cbind(df$Genome, df$g)] <- 1
  m
}
split_rows <- function(df, col) {                 # "a;b" -> two rows
  df[[col]] <- strsplit(as.character(df[[col]]), ";\\s*")
  tidyr::unnest(df, all_of(col))
}
prevalent <- function(Fm, n = MIN_MAGS) {        # features carried / non-zero in >= n MAGs
  keep <- colSums(!is.na(Fm) & Fm > 0) >= n
  Fm[, keep, drop = FALSE]
}

# ---------------------------------------------------------------- feature sets
build <- list()
build$microtrait <- function() {
  f <- inp(file.path("microtrait", "tables", paste0("trait_matrixatgranularity", opt$granularity, ".tsv")))
  if (is.null(f)) return(NULL)
  mt <- read.delim(f, check.names = FALSE); rownames(mt) <- mt$id
  tr <- setdiff(colnames(mt), c("id", "mingentime", "optimumT"))
  P <- prevalent(as.matrix(mt[, tr] > 0) * 1)
  extra <- cbind(log10_mingentime = log10(suppressWarnings(as.numeric(mt$mingentime))),
                 optimumT = suppressWarnings(as.numeric(mt$optimumT)))
  rownames(extra) <- rownames(mt)
  list(values = cbind(weighted_profile(P), weighted_profile(extra)),
       ylab = "community fraction carrying trait", label = function(x) sub("^.*:", "", x))
}
build$portraits <- function() {
  f <- inp("portraits/portraits_results.tsv.gz"); if (is.null(f)) return(NULL)
  d <- read_tsv(f, show_col_types = FALSE, progress = FALSE) %>%
    group_by(tool, feature) %>% mutate(n_models = n_distinct(tool_feature)) %>% ungroup() %>%
    mutate(key = ifelse(n_models > 1, paste(tool, feature, tool_feature, sep = "|"), paste(tool, feature, sep = "|")))
  w <- d %>% select(genome, key, value_probability) %>%
    pivot_wider(names_from = key, values_from = value_probability, values_fn = mean)
  Fm <- as.matrix(w[, -1]); rownames(Fm) <- w$genome
  list(values = weighted_profile(Fm), ylab = "community-weighted mean probability", label = identity)
}
build$amr <- function() {
  f <- inp("amr/rgi_all.tsv"); if (is.null(f)) return(NULL)
  d <- read_tsv(f, show_col_types = FALSE, progress = FALSE, col_types = cols(.default = "c")) %>%
    filter(Cut_Off %in% c("Perfect", "Strict")) %>% select(Genome, `Drug Class`)
  P <- presence_by_group(split_rows(d, "Drug Class"), "Drug Class")
  P <- cbind(prevalent(P), any_ARG = presence_by_group(mutate(d, g = "any_ARG"), "g")[, 1])
  list(values = weighted_profile(P), ylab = "community fraction carrying >= 1 ARG of class", label = identity)
}
build$bacmet <- function() {
  f <- inp("amr/bacmet_all.tsv"); if (is.null(f)) return(NULL)
  d <- read_tsv(f, show_col_types = FALSE, progress = FALSE, col_types = cols(.default = "c"))
  P <- prevalent(presence_by_group(d, "bacmet_gene"))
  P <- cbind(P, any_BacMet = presence_by_group(mutate(d, g = "any_BacMet"), "g")[, 1])
  list(values = weighted_profile(P), ylab = "community fraction carrying gene", label = identity)
}
build$bgc <- function() {
  f <- inp("bgc/antismash_regions_all.tsv"); if (is.null(f)) return(NULL)
  d <- read_tsv(f, show_col_types = FALSE, progress = FALSE, col_types = cols(.default = "c")) %>%
    filter(tolower(contig_edge) != "true")                       # complete regions only
  long <- split_rows(d %>% select(Genome, products), "products")
  cnt <- long %>% count(Genome, products) %>% filter(Genome %in% rownames(ab))
  m <- matrix(0, nrow(ab), length(unique(cnt$products)), dimnames = list(rownames(ab), sort(unique(cnt$products))))
  m[cbind(cnt$Genome, cnt$products)] <- cnt$n
  tot <- d %>% count(Genome) %>% filter(Genome %in% rownames(ab))
  all <- setNames(rep(0, nrow(ab)), rownames(ab)); all[tot$Genome] <- tot$n
  list(values = weighted_profile(cbind(prevalent(m), all_complete_BGCs = all)),
       ylab = "mean complete BGCs per genome (community-weighted)", label = identity)
}
build$gcf <- function() {
  f <- inp("bgc/bigslice/tables/gcf_sample_rel_abund.tsv"); if (is.null(f)) return(NULL)
  g <- read.delim(f, check.names = FALSE, row.names = 1)
  g <- as.matrix(g[, intersect(colnames(g), colnames(ab)), drop = FALSE])
  g <- g[rowSums(g > 0) >= 0.1 * ncol(g), , drop = FALSE]    # GCFs in >= 10% of samples
  say("  gcf: ", nrow(g), " GCFs present in >= 10% of samples")
  list(values = t(g), composition = TRUE, ylab = "relative abundance", label = identity)
}

# ---------------------------------------------------------------- analysis
time_factor <- function(d, sc) {
  if (sc == "seasonal") factor(format(d$date, "%d %b"), levels = unique(format(sort(d$date), "%d %b")))
  else if (sc == "daily") factor(paste0("d", d$day, " ", sprintf("%02d:00", d$hour)))
  else factor(sprintf("%02d:00", d$hour), levels = sprintf("%02d:00", sort(unique(d$hour))))
}
fig <- function(set, name) { dir.create(file.path(FIG, set), recursive = TRUE, showWarnings = FALSE)
                             file.path(FIG, set, name) }
save_plot <- function(p, set, name, width, height) {      # ggplot -> own PDF + report
  ggsave(fig(set, name), p, width = width, height = height)
  REPORT[[length(REPORT) + 1]] <<- list(type = "gg", obj = p, section = set, file = file.path(set, name))
}
save_heatmap <- function(ph, set, name, width, height) {  # pheatmap (silent) -> own PDF + report
  pdf(fig(set, name), width = width, height = height); grid.newpage(); grid.draw(ph$gtable); dev.off()
  REPORT[[length(REPORT) + 1]] <<- list(type = "grob", obj = ph$gtable, section = set, file = file.path(set, name))
}

analyse <- function(set, B, sc) {
  d <- sample_meta %>% filter(scale == sc, sample %in% rownames(B$values))
  if (nrow(d) < 10) { say("  ", sc, ": too few samples, skipped"); return(NULL) }
  d$time <- time_factor(d, sc)
  X <- B$values[d$sample, , drop = FALSE]
  v <- apply(X, 2, sd, na.rm = TRUE)
  X <- X[, is.finite(v) & v > 0 & colSums(is.na(X)) == 0, drop = FALSE]
  if (ncol(X) < 2) { say("  ", sc, ": fewer than 2 variable features, skipped"); return(NULL) }
  comp <- isTRUE(B$composition)
  say("  ", sc, ": ", nrow(d), " samples, ", ncol(X), " variable features")

  # PERMANOVA (within-plot permutations)
  set.seed(1)
  dis <- if (comp) vegdist(X, "bray") else dist(scale(X))
  pm <- adonis2(dis ~ compartment * treatment * time, data = d,
                permutations = permute::how(blocks = d$plot, nperm = NPERM), by = "terms")
  pm <- data.frame(set = set, scale = sc, term = rownames(pm), as.data.frame(pm), row.names = NULL, check.names = FALSE)
  pm$`Pr(>F)`[pm$term == "treatment"] <- NA
  write_tsv(pm, file.path(TAB, set, paste0("permanova_", sc, ".tsv")))

  # ordination figure (PCA, or PCoA on Bray-Curtis for compositions)
  ord <- if (comp) { o <- cmdscale(dis, k = 2, eig = TRUE); ve <- round(100 * o$eig[1:2] / sum(o$eig[o$eig > 0]), 1)
                     list(x = o$points, lab = c("PCoA1", "PCoA2"), ve = ve) }
         else { o <- prcomp(X, scale. = TRUE); ve <- round(100 * o$sdev^2 / sum(o$sdev^2), 1)
                list(x = o$x[, 1:2], lab = c("PC1", "PC2"), ve = ve[1:2]) }
  p <- ggplot(cbind(d, A1 = ord$x[, 1], A2 = ord$x[, 2]), aes(A1, A2, colour = time, shape = compartment)) +
    geom_point(size = 2.4, alpha = 0.85) + facet_wrap(~ treatment) + theme_bw() +
    labs(title = paste0(set, " - ", sc, ": ", if (comp) "GCF composition (Bray-Curtis PCoA)" else "community-weighted profiles (PCA)"),
         x = paste0(ord$lab[1], " (", ord$ve[1], "%)"), y = paste0(ord$lab[2], " (", ord$ve[2], "%)"))
  save_plot(p, set, paste0(sc, "_ordination.pdf"), 11, 6)
  if (comp) return(list(pm = pm, lm = NULL))

  # per-feature mixed models
  d$sin_h <- sin(2 * pi * d$hour / 24); d$cos_h <- cos(2 * pi * d$hour / 24)
  rhs <- if (sc == "diel") "compartment * treatment * (sin_h + cos_h)" else "compartment * treatment * time"
  res <- lapply(colnames(X), function(ft) {
    dd <- d; dd$y <- X[, ft]
    fit <- tryCatch(lme(as.formula(paste("y ~", rhs)), random = ~ 1 | plot, data = dd,
                        control = lmeControl(returnObject = TRUE)), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    a <- anova(fit)
    data.frame(feature = ft, term = rownames(a), F = a$`F-value`, p = a$`p-value`,
               mean_bulk = mean(dd$y[dd$compartment == "bulk"]),
               mean_rhizosphere = mean(dd$y[dd$compartment == "rhizosphere"]),
               mean_control = mean(dd$y[dd$treatment == "control"]),
               mean_biochar = mean(dd$y[dd$treatment == "biochar"]))
  })
  lm_res <- bind_rows(res)
  if (!nrow(lm_res)) return(list(pm = pm, lm = NULL))
  lm_res <- lm_res %>% filter(term != "(Intercept)") %>% group_by(term) %>%
    mutate(p_adj = p.adjust(p, "BH")) %>% ungroup() %>%
    mutate(set = set, scale = sc, rhizo_minus_bulk = mean_rhizosphere - mean_bulk,
           biochar_minus_control = mean_biochar - mean_control) %>%
    relocate(set, scale) %>% arrange(term, p_adj)
  write_tsv(lm_res, file.path(TAB, set, paste0("lmm_", sc, ".tsv")))

  # trajectories of the features with the strongest compartment x time effects
  top <- lm_res %>% filter(grepl("compartment", term), grepl("time|sin_h|cos_h", term)) %>%
    group_by(feature) %>% summarise(p = min(p_adj), .groups = "drop") %>% arrange(p) %>% head(12) %>% pull(feature)
  if (!length(top)) top <- lm_res %>% arrange(p_adj) %>% distinct(feature) %>% head(12) %>% pull(feature)
  long <- data.frame(d[, c("time", "compartment", "treatment")], X[, top, drop = FALSE], check.names = FALSE) %>%
    pivot_longer(all_of(top), names_to = "feature", values_to = "v") %>%
    group_by(feature, time, compartment, treatment) %>%
    summarise(m = mean(v), se = sd(v) / sqrt(n()), .groups = "drop") %>% mutate(feature = B$label(feature))
  p <- ggplot(long, aes(time, m, colour = compartment, linetype = treatment, group = interaction(compartment, treatment))) +
    geom_line() + geom_point(size = 1) + geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.15) +
    facet_wrap(~ feature, scales = "free_y") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), strip.text = element_text(size = 7)) +
    labs(title = paste0(set, " - ", sc, ": strongest compartment x time effects"), y = B$ylab, x = NULL)
  save_plot(p, set, paste0(sc, "_trajectories.pdf"), 13, 9)

  # rhizosphere - bulk effect sizes for features with a significant compartment effect
  eff <- lm_res %>% filter(term == "compartment", p_adj < 0.05) %>% arrange(desc(abs(rhizo_minus_bulk))) %>% head(30) %>%
    mutate(lab = B$label(feature))
  if (nrow(eff)) {
    p <- ggplot(eff, aes(reorder(lab, rhizo_minus_bulk), rhizo_minus_bulk, fill = rhizo_minus_bulk > 0)) +
      geom_col() + coord_flip() + theme_bw() + guides(fill = "none") +
      labs(title = paste0(set, " - ", sc, ": rhizosphere minus bulk (BH p < 0.05, top 30)"), x = NULL, y = paste("difference in", B$ylab))
    save_plot(p, set, paste0(sc, "_compartment_effects.pdf"), 10, 8)
  }

  # heatmap of the 40 most responsive features
  hm_f <- lm_res %>% group_by(feature) %>% summarise(p = min(p_adj), .groups = "drop") %>% arrange(p) %>% head(40) %>% pull(feature)
  if (length(hm_f) > 1) {
    hm <- data.frame(d[, c("compartment", "treatment", "time")], X[, hm_f, drop = FALSE], check.names = FALSE) %>%
      group_by(compartment, treatment, time) %>% summarise(across(all_of(hm_f), mean), .groups = "drop")
    mat <- t(as.matrix(hm[, hm_f])); rownames(mat) <- make.unique(B$label(rownames(mat)))
    colnames(mat) <- paste(substr(hm$compartment, 1, 5), hm$treatment, hm$time)
    mat <- mat[apply(mat, 1, sd) > 0, , drop = FALSE]
    if (nrow(mat) > 1)
      save_heatmap(pheatmap(mat, scale = "row", cluster_cols = FALSE, fontsize_row = 6, fontsize_col = 6, silent = TRUE,
                            main = paste0(set, " - ", sc, ": mean values (row z-scores; 40 most responsive features)")),
                   set, paste0(sc, "_heatmap.pdf"), 12, 9)
  }
  list(pm = pm, lm = lm_res)
}

# ---------------------------------------------------------------- run
all_pm <- list(); all_lm <- list()
for (set in SETS) {
  if (is.null(build[[set]])) { say(set, ": unknown set, skipped"); next }
  say(set, ":")
  B <- tryCatch(build[[set]](), error = function(e) { say("  build failed: ", conditionMessage(e)); NULL })
  if (is.null(B)) next
  B$values <- B$values[, colSums(!is.na(B$values)) > 0, drop = FALSE]
  dir.create(file.path(TAB, set), recursive = TRUE, showWarnings = FALSE)
  say("  ", ncol(B$values), " features x ", nrow(B$values), " samples")
  write_tsv(sample_meta %>% inner_join(data.frame(sample = rownames(B$values), B$values, check.names = FALSE), by = "sample"),
            file.path(TAB, set, "cwm_by_sample.tsv"))
  for (sc in SCALES) {
    r <- tryCatch(analyse(set, B, sc), error = function(e) { say("  ", sc, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) { all_pm[[paste(set, sc)]] <- r$pm; all_lm[[paste(set, sc)]] <- r$lm }
  }
}

# ---------------------------------------------------------------- summary across sets
pm_all <- bind_rows(all_pm); lm_all <- bind_rows(all_lm)
if (nrow(lm_all)) {
  summ <- lm_all %>% group_by(set, scale, term) %>%
    summarise(n_features = n(), n_sig_BH05 = sum(p_adj < 0.05, na.rm = TRUE), .groups = "drop") %>%
    full_join(pm_all %>% filter(!term %in% c("Residual", "Total")) %>%
                select(set, scale, term, permanova_R2 = R2, permanova_p = `Pr(>F)`),
              by = c("set", "scale", "term")) %>%
    mutate(term = sub("time$", "time", term)) %>% arrange(set, scale, term)
  write_tsv(summ, file.path(TAB, "summary_terms.tsv"))
  p <- ggplot(summ %>% filter(!is.na(n_sig_BH05)),
              aes(term, set, fill = n_sig_BH05 / n_features, label = paste0(n_sig_BH05, "/", n_features))) +
    geom_tile(colour = "white") + geom_text(size = 2.6) + facet_wrap(~ scale, scales = "free_x") +
    scale_fill_gradient(low = "grey95", high = "firebrick", name = "fraction\nsignificant") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Features with a significant effect (BH p < 0.05) per annotation set, time scale and model term",
         x = NULL, y = NULL)
  save_plot(p, "summary", "significant_terms.pdf", 14, 6)
  say("summary: ", file.path(TAB, "summary_terms.tsv"))
}
# ---------------------------------------------------------------- PDF report
text_pages <- function(title, lines, mono = TRUE, per_page = 62, width = 140) {
  lines <- unlist(lapply(lines, function(l) {            # hard-wrap long lines, keep indentation
    if (nchar(l) <= width) return(l)
    substring(l, seq(1, nchar(l), width), pmin(seq(width, nchar(l) + width - 1, width), nchar(l))) }))
  pages <- split(lines, ceiling(seq_along(lines) / per_page))
  for (k in seq_along(pages)) {
    grid.newpage()
    grid.text(if (length(pages) > 1) paste0(title, "  (", k, "/", length(pages), ")") else title,
              x = 0.03, y = 0.97, just = c("left", "top"), gp = gpar(fontsize = 13, fontface = "bold"))
    grid.text(paste(pages[[k]], collapse = "\n"), x = 0.03, y = 0.92, just = c("left", "top"),
              gp = gpar(fontsize = 7, fontfamily = if (mono) "mono" else "sans", lineheight = 1.05))
  }
}
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_lines <- if (length(script_path) && file.exists(script_path[1])) readLines(script_path[1]) else "(script source not available)"
header <- script_lines[seq_len(max(1, which(!grepl("^#", script_lines))[1] - 1))]
close(logf)
rep_file <- file.path(FIG, "community_analysis_report.pdf")
pdf(rep_file, width = 11.69, height = 8.27)               # A4 landscape
grid.newpage()
grid.text("Skyline MAGs: community-level analysis of MAG annotations", x = 0.5, y = 0.62,
          gp = gpar(fontsize = 22, fontface = "bold"))
grid.text(paste0("Rhizosphere vs bulk soil, biochar vs control; seasonal, daily and diel time scales\n",
                 "sets: ", paste(SETS, collapse = ", "), "\n",
                 "generated ", format(Sys.time(), "%d %b %Y %H:%M"), " with ", R.version.string, "\n",
                 "results: ", R), x = 0.5, y = 0.45, gp = gpar(fontsize = 11, lineheight = 1.3))
text_pages("Methods (from the script header)", sub("^# ?", "", header), mono = FALSE)
text_pages("Run log", readLines(file.path(TAB, "community_analysis.log")))
if (exists("summ")) {
  st <- summ %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))
  op <- options(width = 190)
  text_pages("Summary: significant features (BH p < 0.05) and PERMANOVA per set, scale and term",
             capture.output(print(as.data.frame(st), row.names = FALSE)), width = 190)
  options(op)
}
order_sec <- c("summary", SETS)
for (sec in order_sec) {
  items <- Filter(function(x) x$section == sec, REPORT)
  if (!length(items)) next
  grid.newpage(); grid.text(sec, gp = gpar(fontsize = 28, fontface = "bold"))
  grid.text(paste0(length(items), " figures; files in figures/", sec, "/"), y = 0.42, gp = gpar(fontsize = 11))
  for (it in items) {
    grid.newpage()
    if (it$type == "gg") print(it$obj) else grid.draw(it$obj)
    grid.text(it$file, x = 0.99, y = 0.01, just = c("right", "bottom"), gp = gpar(fontsize = 6, col = "grey40"))
  }
}
text_pages("Appendix: analysis code (scripts/community_analysis.R)",
           sprintf("%4d  %s", seq_along(script_lines), script_lines))
invisible(dev.off())
message("report: ", rep_file)
message("done -- tables in ", TAB, "; figures in ", FIG)
