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
# Style: tidyverse throughout; colours from ggsci; heatmaps with pheatmap. The abundance
# weighting is a matrix product (one function): in long format it would be ~650 million rows.
#
# Runs in the pipeline's microTrait conda env (tidyverse, vegan, nlme, pheatmap, ggsci).
# Usage:
#   Rscript community_analysis.R --results-dir <annotations> --metadata <sheet.tsv>
#       [--sets microtrait,portraits,amr,bacmet,bgc,gcf] [--nperm 999] [--min-mags 20]
#       [--granularity 3] [--abund coverm/coverm_rel_abund_masked_cf0.1.tsv]
# =============================================================================

needed <- c("dplyr", "tidyr", "readr", "purrr", "stringr", "tibble", "forcats", "lubridate", "ggplot2",
            "vegan", "nlme", "pheatmap", "ggsci")
absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent) > 0)
  stop("missing R packages: ", paste(absent, collapse = ", "),
       " -- install into this env (see scripts/install_r_packages.R)", call. = FALSE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr)
  library(tibble); library(forcats); library(lubridate); library(ggplot2)
  library(vegan); library(nlme); library(pheatmap); library(ggsci); library(grid)
})

# ---------------------------------------------------------------- arguments
args <- commandArgs(trailingOnly = TRUE)
opt <- list(sets = "microtrait,portraits,amr,bacmet,bgc,gcf", nperm = "999", `min-mags` = "20",
            granularity = "3", abund = "coverm/coverm_rel_abund_masked_cf0.1.tsv",
            `results-dir` = NA, metadata = NA) |>
  modifyList(set_names(as.list(args[c(FALSE, TRUE)]), str_remove(args[c(TRUE, FALSE)], "^--")))
if (is.na(opt$`results-dir`) || is.na(opt$metadata))
  stop("Usage: Rscript community_analysis.R --results-dir <annotations> --metadata <sheet.tsv>")
R        <- opt$`results-dir`
SETS     <- str_split_1(opt$sets, ",")
NPERM    <- as.integer(opt$nperm)
MIN_MAGS <- as.integer(opt$`min-mags`)
TAB      <- file.path(R, "analysis_tables")
FIG      <- file.path(R, "figures")
SCALES   <- c("seasonal", "daily", "diel")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

logf <- file(file.path(TAB, "community_analysis.log"), open = "wt")
say <- function(...) { m <- str_c(...); message(m); writeLines(m, logf); flush(logf) }
say("community analysis -- ", format(Sys.time()), " -- sets: ", str_flatten_comma(SETS))
inp <- function(rel) {
  f <- file.path(R, rel)
  if (file.exists(f)) f else { say("  missing input: ", f); NULL }
}

# Known date error in the sequencing sheet: 07/12/2024 is 12 July 2024 (no December sampling)
DATE_FIXES <- c("07/12/2024" = "12/07/2024")

# ---------------------------------------------------------------- colours (ggsci)
PAL_COMP <- set_names(pal_npg()(2), c("bulk", "rhizosphere"))
PAL_TRT  <- set_names(pal_npg()(4)[3:4], c("control", "biochar"))
PAL_HEAT <- pal_gsea("default", n = 50)(50)                       # diverging, for pheatmap
theme_set(theme_bw(base_size = 10))

# ---------------------------------------------------------------- metadata
read_sheet <- function(path) {
  lines <- read_lines(path, locale = locale(encoding = "latin1"))
  hdr <- str_split_1(lines[1], "\t") |> str_trim() |> str_remove_all(" ")
  rows <- lines[-1] |> discard(\(l) str_trim(l) == "") |> str_split("\t")
  if (max(lengths(rows)) == length(hdr) + 1 && "dna_plate_numdna_plate_col" %in% hdr) {
    i <- match("dna_plate_numdna_plate_col", hdr)
    hdr <- append(hdr[-i], c("dna_plate_num", "dna_plate_col"), after = i - 1)
    say("metadata: split fused header 'dna_plate_numdna_plate_col'")
  }
  rows |>
    map(\(r) c(r, rep("", length(hdr)))[seq_along(hdr)] |> str_trim() |> set_names(hdr) |> as.list() |> as_tibble_row()) |>
    list_rbind()
}

meta <- read_sheet(opt$metadata) |>
  filter(dna_tube_num != "") |>
  mutate(date_raw    = date,
         date        = dmy(coalesce(DATE_FIXES[date], date)),
         sample      = dna_tube_num,
         compartment = fct_recode(factor(bulk_rhizo_plant), rhizosphere = "rhizo") |> fct_relevel("bulk"),
         treatment   = fct_recode(factor(treatment), biochar = "bio", control = "con") |> fct_relevel("control"),
         plot        = factor(plot),
         day         = suppressWarnings(as.integer(day)),
         hour        = suppressWarnings(as.integer(time) %/% 100),
         scale       = case_when(is.na(day)          ~ "seasonal",
                                 day == 3            ~ "diel",
                                 hour %in% c(4, 14)  ~ "daily"))
n_fixed <- sum(meta$date_raw %in% names(DATE_FIXES))
if (n_fixed > 0) say("metadata: corrected date 07/12/2024 -> 12/07/2024 for ", n_fixed, " samples")
say("metadata: ", nrow(meta), " DNA samples; ",
    meta |> count(scale) |> mutate(s = str_c(scale, "=", n)) |> pull(s) |> str_flatten_comma())

# ---------------------------------------------------------------- abundances
ab_tbl <- read_tsv(file.path(R, opt$abund), show_col_types = FALSE)
genome_col <- names(ab_tbl)[1]
excluded <- setdiff(names(ab_tbl)[-1], meta$sample)
ab_long_samples <- intersect(names(ab_tbl)[-1], meta$sample)
say("abundance: ", ncol(ab_tbl) - 1, " samples in matrix; kept ", length(ab_long_samples),
    "; excluded (not in sequencing sheet): ", str_flatten_comma(excluded))
ab <- ab_tbl |>
  select(all_of(c(genome_col, ab_long_samples))) |>
  mutate(across(-all_of(genome_col), \(x) replace_na(x, 0))) |>
  column_to_rownames(genome_col) |>
  as.matrix()
empty <- colnames(ab)[colSums(ab) <= 0]
if (length(empty) > 0) {
  say("abundance: no MAG coverage, dropped: ", str_flatten_comma(empty))
  ab <- ab[, !colnames(ab) %in% empty, drop = FALSE]
}
missing <- setdiff(meta$sample, colnames(ab))
say("abundance: sequencing-sheet samples without usable CoverM data: ", length(missing),
    if (length(missing) > 0) str_c(" (", str_flatten_comma(missing), ")") else "")
sample_meta <- meta |>
  filter(sample %in% colnames(ab)) |>
  select(sample, date, day, hour, plot, treatment, compartment, scale)
MAGS <- rownames(ab)

# ---------------------------------------------------------------- helpers
# abundance-weighted mean of a MAG x feature matrix, NA-aware per feature (matrix product: see header)
weighted_profile <- function(Fm) {
  mags <- intersect(rownames(Fm), MAGS)
  Fm <- Fm[mags, , drop = FALSE]
  A  <- ab[mags, , drop = FALSE]
  ok <- !is.na(Fm)
  F0 <- replace(Fm, !ok, 0)
  num <- t(A) %*% F0
  den <- t(A) %*% (ok * 1)
  replace(num / den, den == 0, NA)
}
# tibble (Genome + feature columns) -> MAG x feature matrix covering all MAGs (absent = 0)
to_mag_matrix <- function(tbl) {
  tibble(Genome = MAGS) |>
    left_join(tbl, by = "Genome") |>
    mutate(across(-Genome, \(x) replace_na(x, 0))) |>
    column_to_rownames("Genome") |>
    as.matrix()
}
# long (Genome, group) -> MAG x group presence matrix
presence <- function(df, group_col) {
  df |>
    filter(!is.na(.data[[group_col]]), .data[[group_col]] != "", Genome %in% MAGS) |>
    distinct(Genome, g = .data[[group_col]]) |>
    mutate(v = 1) |>
    pivot_wider(names_from = g, values_from = v, values_fill = 0) |>
    to_mag_matrix()
}
any_presence <- function(df, name) {
  df |> filter(Genome %in% MAGS) |> distinct(Genome) |> mutate("{name}" := 1) |> to_mag_matrix()
}
# keep features carried by / non-zero in >= MIN_MAGS MAGs
prevalent <- function(Fm) Fm[, colSums(Fm > 0, na.rm = TRUE) >= MIN_MAGS, drop = FALSE]
read_chr <- function(f) read_tsv(f, col_types = cols(.default = "c"), progress = FALSE)

# ---------------------------------------------------------------- feature sets
SET_INFO <- c(
  microtrait = "fraction of the community carrying each microTrait trait (granularity 3); abundance-weighted log10 minimum generation time (gRodon) and optimum temperature",
  portraits  = "abundance-weighted mean predicted probability per tool|trait[|model] (BacDive-AI, GenomeSPOT, Traitar, MICROPHERRET)",
  amr        = "fraction of the community carrying >= 1 RGI hit (Perfect + Strict) per drug class; any_ARG",
  bacmet     = "fraction of the community carrying each BacMet2 (EXP) gene; any_BacMet",
  bgc        = "abundance-weighted mean number of complete antiSMASH regions per genome, per class and in total",
  gcf        = "BiG-SLiCE GCF composition (GCFs in >= 10% of samples): Bray-Curtis PERMANOVA and PCoA")

build <- list(
  microtrait = function() {
    f <- inp(file.path("microtrait", "tables", str_c("trait_matrixatgranularity", opt$granularity, ".tsv")))
    if (is.null(f)) return(NULL)
    mt <- read_tsv(f, show_col_types = FALSE, progress = FALSE) |> rename(Genome = id)
    traits <- mt |>
      select(-any_of(c("mingentime", "optimumT"))) |>
      mutate(across(-Genome, \(x) as.integer(x > 0))) |>
      to_mag_matrix() |>
      prevalent()
    extra <- mt |>
      transmute(Genome, log10_mingentime = log10(as.numeric(mingentime)), optimumT = as.numeric(optimumT)) |>
      column_to_rownames("Genome") |>
      as.matrix()
    list(values = cbind(weighted_profile(traits), weighted_profile(extra)),
         ylab = "community fraction carrying trait", label = \(x) str_remove(x, "^.*:"))
  },
  portraits = function() {
    f <- inp("portraits/portraits_results.tsv.gz"); if (is.null(f)) return(NULL)
    Fm <- read_tsv(f, show_col_types = FALSE, progress = FALSE) |>
      group_by(tool, feature) |>
      mutate(key = if (n_distinct(tool_feature) > 1) str_c(tool, feature, tool_feature, sep = "|")
                   else str_c(tool, feature, sep = "|")) |>
      ungroup() |>
      select(Genome = genome, key, value_probability) |>
      pivot_wider(names_from = key, values_from = value_probability, values_fn = mean) |>
      column_to_rownames("Genome") |>
      as.matrix()
    list(values = weighted_profile(Fm), ylab = "community-weighted mean probability", label = identity)
  },
  amr = function() {
    f <- inp("amr/rgi_all.tsv"); if (is.null(f)) return(NULL)
    d <- read_chr(f) |> filter(Cut_Off %in% c("Perfect", "Strict")) |> select(Genome, drug_class = `Drug Class`)
    P <- d |> separate_rows(drug_class, sep = ";\\s*") |> presence("drug_class") |> prevalent()
    list(values = weighted_profile(cbind(P, any_presence(d, "any_ARG"))),
         ylab = "community fraction carrying >= 1 ARG of class", label = identity)
  },
  bacmet = function() {
    f <- inp("amr/bacmet_all.tsv"); if (is.null(f)) return(NULL)
    d <- read_chr(f)
    P <- d |> presence("bacmet_gene") |> prevalent()
    list(values = weighted_profile(cbind(P, any_presence(d, "any_BacMet"))),
         ylab = "community fraction carrying gene", label = identity)
  },
  bgc = function() {
    f <- inp("bgc/antismash_regions_all.tsv"); if (is.null(f)) return(NULL)
    d <- read_chr(f) |> filter(str_to_lower(contig_edge) != "true")          # complete regions only
    per_class <- d |>
      separate_rows(products, sep = ";") |>
      count(Genome, products) |>
      pivot_wider(names_from = products, values_from = n, values_fill = 0) |>
      to_mag_matrix() |>
      prevalent()
    total <- d |> count(Genome, name = "all_complete_BGCs") |> to_mag_matrix()
    list(values = weighted_profile(cbind(per_class, total)),
         ylab = "mean complete BGCs per genome (community-weighted)", label = identity)
  },
  gcf = function() {
    f <- inp("bgc/bigslice/tables/gcf_sample_rel_abund.tsv"); if (is.null(f)) return(NULL)
    g <- read_tsv(f, show_col_types = FALSE, progress = FALSE)
    keep_samples <- intersect(names(g)[-1], colnames(ab))
    g <- g |>
      select(gcf = 1, all_of(keep_samples)) |>
      filter(rowSums(across(-gcf) > 0) >= 0.1 * length(keep_samples))         # GCFs in >= 10% of samples
    say("  gcf: ", nrow(g), " GCFs present in >= 10% of samples")
    list(values = g |> column_to_rownames("gcf") |> as.matrix() |> t(),
         composition = TRUE, ylab = "relative abundance", label = identity)
  })

# ---------------------------------------------------------------- figures and report store
REPORT <- list()
save_plot <- function(p, set, name, width, height) {
  dir.create(file.path(FIG, set), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(FIG, set, name), p, width = width, height = height)
  REPORT[[length(REPORT) + 1]] <<- list(type = "gg", obj = p, section = set, file = file.path(set, name))
}
save_heatmap <- function(ph, set, name, width, height) {
  dir.create(file.path(FIG, set), recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(FIG, set, name), width = width, height = height); grid.draw(ph$gtable); invisible(dev.off())
  REPORT[[length(REPORT) + 1]] <<- list(type = "grob", obj = ph$gtable, section = set, file = file.path(set, name))
}

# ---------------------------------------------------------------- analysis
add_time <- function(d, sc) {
  d |> mutate(time = case_when(
    sc == "seasonal" ~ format(date, "%d %b"),
    sc == "daily"    ~ str_c("d", day, " ", sprintf("%02d:00", hour)),
    TRUE             ~ sprintf("%02d:00", hour)),
    time = fct_reorder(time, as.numeric(date) * 100 + coalesce(hour, 0L)),
    sin_h = sin(2 * pi * hour / 24), cos_h = cos(2 * pi * hour / 24))
}

fit_feature <- function(ft, d, X, rhs) {
  dd  <- d |> mutate(y = X[, ft])
  fit <- lme(as.formula(str_c("y ~ ", rhs)), random = ~ 1 | plot, data = dd,
             control = lmeControl(returnObject = TRUE))
  a <- anova(fit)
  tibble(feature = ft, term = rownames(a), F = a$`F-value`, p = a$`p-value`) |>
    bind_cols(dd |> summarise(mean_bulk        = mean(y[compartment == "bulk"]),
                              mean_rhizosphere = mean(y[compartment == "rhizosphere"]),
                              mean_control     = mean(y[treatment == "control"]),
                              mean_biochar     = mean(y[treatment == "biochar"])))
}

analyse <- function(set, B, sc) {
  d <- sample_meta |> filter(scale == sc, sample %in% rownames(B$values)) |> add_time(sc)
  if (nrow(d) < 10) { say("  ", sc, ": too few samples, skipped"); return(NULL) }
  X <- B$values[d$sample, , drop = FALSE]
  sds <- apply(X, 2, sd, na.rm = TRUE)
  X <- X[, is.finite(sds) & sds > 0 & colSums(is.na(X)) == 0, drop = FALSE]
  if (ncol(X) < 2) { say("  ", sc, ": fewer than 2 variable features, skipped"); return(NULL) }
  comp <- isTRUE(B$composition)
  say("  ", sc, ": ", nrow(d), " samples, ", ncol(X), " variable features")

  # PERMANOVA (permutations within plot)
  set.seed(1)
  dis <- if (comp) vegdist(X, "bray") else dist(scale(X))
  pm <- adonis2(dis ~ compartment * treatment * time, data = d,
                permutations = permute::how(blocks = d$plot, nperm = NPERM), by = "terms") |>
    as.data.frame() |>
    rownames_to_column("term") |>
    as_tibble() |>
    mutate(set = set, scale = sc, `Pr(>F)` = if_else(term == "treatment", NA_real_, `Pr(>F)`), .before = 1)
  write_tsv(pm, file.path(TAB, set, str_c("permanova_", sc, ".tsv")))

  # ordination: PCA, or PCoA on Bray-Curtis for compositions
  if (comp) {
    o  <- cmdscale(dis, k = 2, eig = TRUE)
    sc_xy <- as_tibble(o$points, .name_repair = \(x) c("A1", "A2"))
    axl <- str_c(c("PCoA1", "PCoA2"), " (", round(100 * o$eig[1:2] / sum(o$eig[o$eig > 0]), 1), "%)")
  } else {
    o  <- prcomp(X, scale. = TRUE)
    sc_xy <- tibble(A1 = o$x[, 1], A2 = o$x[, 2])
    axl <- str_c(c("PC1", "PC2"), " (", round(100 * o$sdev[1:2]^2 / sum(o$sdev^2), 1), "%)")
  }
  p <- bind_cols(d, sc_xy) |>
    ggplot(aes(A1, A2, colour = time, shape = compartment)) +
    geom_point(size = 2.4, alpha = 0.85) +
    facet_wrap(~ treatment) +
    scale_colour_npg() +
    labs(title = str_c(set, " - ", sc, ": ", if (comp) "GCF composition (Bray-Curtis PCoA)" else "community-weighted profiles (PCA)"),
         x = axl[1], y = axl[2])
  save_plot(p, set, str_c(sc, "_ordination.pdf"), 11, 6)
  if (comp) return(list(pm = pm, lm = NULL))

  # per-feature mixed models
  rhs <- if (sc == "diel") "compartment * treatment * (sin_h + cos_h)" else "compartment * treatment * time"
  lm_res <- map(colnames(X), possibly(fit_feature, otherwise = NULL), d = d, X = X, rhs = rhs) |> list_rbind()
  if (nrow(lm_res) == 0) return(list(pm = pm, lm = NULL))
  lm_res <- lm_res |>
    filter(term != "(Intercept)") |>
    mutate(p_adj = p.adjust(p, "BH"), .by = term) |>
    mutate(set = set, scale = sc,
           rhizo_minus_bulk      = mean_rhizosphere - mean_bulk,
           biochar_minus_control = mean_biochar - mean_control, .before = 1) |>
    arrange(term, p_adj)
  write_tsv(lm_res, file.path(TAB, set, str_c("lmm_", sc, ".tsv")))

  # trajectories of the features with the strongest compartment x time effects
  top <- lm_res |>
    filter(str_detect(term, "compartment"), str_detect(term, "time|sin_h|cos_h")) |>
    summarise(p = min(p_adj), .by = feature) |> slice_min(p, n = 12, with_ties = FALSE) |> pull(feature)
  if (length(top) == 0) top <- lm_res |> summarise(p = min(p_adj), .by = feature) |>
    slice_min(p, n = 12, with_ties = FALSE) |> pull(feature)
  p <- d |>
    select(time, compartment, treatment) |>
    bind_cols(as_tibble(X[, top, drop = FALSE])) |>
    pivot_longer(all_of(top), names_to = "feature", values_to = "v") |>
    summarise(m = mean(v), se = sd(v) / sqrt(n()), .by = c(feature, time, compartment, treatment)) |>
    mutate(feature = B$label(feature)) |>
    ggplot(aes(time, m, colour = compartment, linetype = treatment, group = interaction(compartment, treatment))) +
    geom_line() + geom_point(size = 1) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.15) +
    facet_wrap(~ feature, scales = "free_y") +
    scale_colour_manual(values = PAL_COMP) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), strip.text = element_text(size = 7)) +
    labs(title = str_c(set, " - ", sc, ": strongest compartment x time effects"), y = B$ylab, x = NULL)
  save_plot(p, set, str_c(sc, "_trajectories.pdf"), 13, 9)

  # rhizosphere - bulk differences for features with a significant compartment effect
  eff <- lm_res |>
    filter(term == "compartment", p_adj < 0.05) |>
    slice_max(abs(rhizo_minus_bulk), n = 30, with_ties = FALSE) |>
    mutate(lab = fct_reorder(B$label(feature), rhizo_minus_bulk),
           higher_in = if_else(rhizo_minus_bulk > 0, "rhizosphere", "bulk"))
  if (nrow(eff) > 0) {
    p <- ggplot(eff, aes(lab, rhizo_minus_bulk, fill = higher_in)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = PAL_COMP, name = "higher in") +
      labs(title = str_c(set, " - ", sc, ": rhizosphere minus bulk (BH p < 0.05, top 30)"),
           x = NULL, y = str_c("difference in ", B$ylab))
    save_plot(p, set, str_c(sc, "_compartment_effects.pdf"), 10, 8)
  }

  # heatmap (pheatmap, ggsci colours) of the 40 most responsive features
  hm_f <- lm_res |> summarise(p = min(p_adj), .by = feature) |> slice_min(p, n = 40, with_ties = FALSE) |> pull(feature)
  if (length(hm_f) > 1) {
    hm <- d |>
      select(compartment, treatment, time) |>
      bind_cols(as_tibble(X[, hm_f, drop = FALSE])) |>
      summarise(across(all_of(hm_f), mean), .by = c(compartment, treatment, time)) |>
      arrange(compartment, treatment, time) |>
      mutate(col_id = str_c(str_sub(compartment, 1, 5), treatment, time, sep = " "))
    mat <- hm |> select(col_id, all_of(hm_f)) |> column_to_rownames("col_id") |> as.matrix() |> t()
    rownames(mat) <- make.unique(B$label(rownames(mat)))
    mat <- mat[apply(mat, 1, sd) > 0, , drop = FALSE]
    if (nrow(mat) > 1) {
      ann <- hm |> select(col_id, compartment, treatment) |> column_to_rownames("col_id")
      ph <- pheatmap(mat, scale = "row", cluster_cols = FALSE, color = PAL_HEAT,
                     annotation_col = ann, annotation_colors = list(compartment = PAL_COMP, treatment = PAL_TRT),
                     fontsize_row = 6, fontsize_col = 6, silent = TRUE,
                     main = str_c(set, " - ", sc, ": mean values (row z-scores; 40 most responsive features)"))
      save_heatmap(ph, set, str_c(sc, "_heatmap.pdf"), 12, 9)
    }
  }
  list(pm = pm, lm = lm_res)
}

# ---------------------------------------------------------------- run all sets
run_set <- function(set) {
  if (is.null(build[[set]])) { say(set, ": unknown set, skipped"); return(NULL) }
  say(set, ":")
  B <- tryCatch(build[[set]](), error = \(e) { say("  build failed: ", conditionMessage(e)); NULL })
  if (is.null(B)) return(NULL)
  B$values <- B$values[, colSums(!is.na(B$values)) > 0, drop = FALSE]
  dir.create(file.path(TAB, set), recursive = TRUE, showWarnings = FALSE)
  say("  ", ncol(B$values), " features x ", nrow(B$values), " samples")
  sample_meta |>
    inner_join(as_tibble(B$values, rownames = "sample"), by = "sample") |>
    write_tsv(file.path(TAB, set, "cwm_by_sample.tsv"))
  map(set_names(SCALES), \(sc) tryCatch(analyse(set, B, sc),
                                        error = \(e) { say("  ", sc, " failed: ", conditionMessage(e)); NULL }))
}
results <- map(set_names(SETS), run_set)
pm_all  <- results |> map(\(r) map(r, "pm") |> list_rbind()) |> list_rbind()
lm_all  <- results |> map(\(r) map(r, "lm") |> list_rbind()) |> list_rbind()

# ---------------------------------------------------------------- summary across sets
if (nrow(lm_all) > 0) {
  summ <- lm_all |>
    summarise(n_features = n(), n_sig_BH05 = sum(p_adj < 0.05, na.rm = TRUE), .by = c(set, scale, term)) |>
    full_join(pm_all |> filter(!term %in% c("Residual", "Total")) |>
                select(set, scale, term, permanova_R2 = R2, permanova_p = `Pr(>F)`),
              by = c("set", "scale", "term")) |>
    arrange(set, scale, term)
  write_tsv(summ, file.path(TAB, "summary_terms.tsv"))
  p <- summ |>
    filter(!is.na(n_sig_BH05)) |>
    ggplot(aes(term, set, fill = n_sig_BH05 / n_features, label = str_c(n_sig_BH05, "/", n_features))) +
    geom_tile(colour = "white") + geom_text(size = 2.6) +
    facet_wrap(~ scale, scales = "free_x") +
    scale_fill_material("red", name = "fraction\nsignificant") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Features with a significant effect (BH p < 0.05) per annotation set, time scale and model term",
         x = NULL, y = NULL)
  save_plot(p, "summary", "significant_terms.pdf", 14, 6)
  say("summary: ", file.path(TAB, "summary_terms.tsv"))
}
say("done -- tables in ", TAB, "; figures in ", FIG)

# ---------------------------------------------------------------- PDF report
text_pages <- function(title, lines, mono = TRUE, per_page = 62, width = 140, fontsize = 7) {
  lines <- lines |> map(\(l) if (nchar(l) <= width) l else str_sub(l, seq(1, nchar(l), width), seq(width, nchar(l) + width - 1, width))) |>
    list_c()
  pages <- split(lines, ceiling(seq_along(lines) / per_page))
  iwalk(unname(pages), \(pg, k) {
    grid.newpage()
    grid.text(if (length(pages) > 1) str_c(title, "  (", k, "/", length(pages), ")") else title,
              x = 0.03, y = 0.97, just = c("left", "top"), gp = gpar(fontsize = 13, fontface = "bold"))
    grid.text(str_flatten(pg, "\n"), x = 0.03, y = 0.92, just = c("left", "top"),
              gp = gpar(fontsize = fontsize, fontfamily = if (mono) "mono" else "sans", lineheight = 1.05))
  })
}
script_path  <- commandArgs(FALSE) |> str_subset("^--file=") |> str_remove("^--file=")
script_lines <- if (length(script_path) > 0 && file.exists(script_path[1])) read_lines(script_path[1]) else "(script source not available)"
header       <- script_lines[seq_len(max(1, which(!str_starts(script_lines, "#"))[1] - 1))]
close(logf)

rep_file <- file.path(FIG, "community_analysis_report.pdf")
pdf(rep_file, width = 11.69, height = 8.27)                   # A4 landscape
grid.newpage()
grid.text("Skyline MAGs: community-level analysis of MAG annotations", x = 0.5, y = 0.62,
          gp = gpar(fontsize = 22, fontface = "bold"))
grid.text(str_c("Rhizosphere vs bulk soil, biochar vs control; seasonal, daily and diel time scales\n",
                "sets: ", str_flatten_comma(SETS), "\n",
                "generated ", format(Sys.time(), "%d %b %Y %H:%M"), " with ", R.version.string, "\n",
                "results: ", R), x = 0.5, y = 0.45, gp = gpar(fontsize = 11, lineheight = 1.3))
text_pages("Methods (from the script header)", str_remove(header, "^# ?"), mono = FALSE)
text_pages("Run log", read_lines(file.path(TAB, "community_analysis.log")))
if (exists("summ")) {
  op <- options(width = 190)
  text_pages("Summary: significant features (BH p < 0.05) and PERMANOVA per set, scale and term",
             summ |> mutate(across(where(is.numeric), \(x) signif(x, 3))) |> as.data.frame() |>
               print(row.names = FALSE) |> capture.output(), width = 190)
  options(op)
}
for (sec in c("summary", SETS)) {
  items <- keep(REPORT, \(x) x$section == sec)
  if (length(items) == 0) next
  grid.newpage()
  grid.text(sec, gp = gpar(fontsize = 28, fontface = "bold"))
  grid.text(str_c(if (!is.na(SET_INFO[sec])) str_c(str_wrap(SET_INFO[sec], 100), "\n\n") else "",
                  length(items), " figures; files in figures/", sec, "/"),
            y = 0.38, gp = gpar(fontsize = 11, lineheight = 1.2))
  walk(items, \(it) {
    grid.newpage()
    if (it$type == "gg") print(it$obj) else grid.draw(it$obj)
    grid.text(it$file, x = 0.99, y = 0.01, just = c("right", "bottom"), gp = gpar(fontsize = 6, col = "grey40"))
  })
}
text_pages("Appendix: analysis code (scripts/community_analysis.R)",
           sprintf("%4d  %s", seq_along(script_lines), script_lines))
text_pages("Session info", capture.output(sessionInfo()))
invisible(dev.off())
message("report: ", rep_file)
