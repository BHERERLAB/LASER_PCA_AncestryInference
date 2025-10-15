library(data.table)
library(ggplot2)

# --- Paths ---
args <- commandArgs(trailingOnly = TRUE)
pca_file     <- args[1]
meta_file    <- args[2]
cag_id_file  <- args[3]
k            <- if (length(args) >= 4) as.integer(args[4]) else 10L
n_pcs        <- if (length(args) >= 5) as.integer(args[5]) else 10L
threshold_N  <- if (length(args) >= 6) as.integer(args[6]) else 7L

if (any(is.na(c(k, n_pcs, threshold_N)))) {
  stop(sprintf("Bad numeric args. Got k=%s n_pcs=%s threshold_N=%s",
               deparse(args[4]), deparse(args[5]), deparse(args[6])))
}


print(paste0("pca_file: ", pca_file))
print(paste0("meta_file: ", meta_file))
print(paste0("cag_id_file: ", cag_id_file))
print(paste0("k: ", k))
print(paste0("n_pcs: ", n_pcs))
print(paste0("threshold_N: ", threshold_N))


# --- Metadata (REF labels + continents) ---
metadata <- fread(meta_file)
setnames(metadata, old = c("s", "hgdp_tgp_meta.Genetic.region"),
         new = c("REF_ID", "Region"))


# --- Read PCA (popID, indivID, PC1..PC20) ---
dt <- fread(pca_file, sep = "\t", header = TRUE, data.table = TRUE, blank.lines.skip = TRUE)
# Drop any weird footer like a lone "!"
dt <- dt[popID != "!" & !is.na(popID)]

# Ensure expected names exist
if (!all(c("popID","indivID") %in% names(dt))) setnames(dt, 1:2, c("popID","indivID"))

# If PC columns aren’t exactly PC1..PC20 for some reason, fix them
pc_cols <- grep("^PC\\d+$", names(dt), value = TRUE)
if (length(pc_cols) < 2) {
  setnames(dt, 3:22, paste0("PC", 1:20))
  pc_cols <- paste0("PC", 1:20)
}

# read CaG IDs
cag_ids <- fread(cag_id_file, header = FALSE, col.names = "IID")[, IID]

# keep only IIDs present in metadata$REF_ID or in CaG list
keep_ids <- unique(c(metadata$REF_ID, cag_ids))
dt <- dt[indivID %chin% keep_ids]


print(head(dt))

# --- Label dataset + continent ---
dt[, IID := indivID]
dt[, FID := popID]
dt[, Dataset := ifelse(IID %in% metadata$REF_ID, "HGDP+1KG", "CaG")]
dt[, Continent := ifelse(Dataset == "HGDP+1KG",
                         metadata$Region[match(IID, metadata$REF_ID)],
                         "CaG")]
dt <- dt[!is.na(Continent)]

# --- Colors (same palette as before) ---
continent_colors <- c(
  "AFR"="#43AA8B","EUR"="#90BE6D","CSA"="#F8961E","EAS"="#F3722C",
  "MID"="#F9C74F","AMR"="#577590","OCE"="#F94144","CaG"="black"
)
dt[, Continent := factor(Continent, levels = names(continent_colors))]

# --- Sanity checks (optional) ---
print(dt[, .N, by = Dataset][order(-N)])
print(dt[, .N, by = Continent][order(-N)])

# --- Plot PC1 vs PC2 overlay ---
ref_data <- dt[Dataset == "HGDP+1KG"]
cag_data <- dt[Dataset == "CaG"]



p <- ggplot() +
  geom_point(data = ref_data, aes(PC1, PC2, color = Continent), size = 2.5, alpha = 0.85) +
  geom_point(data = cag_data, aes(PC1, PC2, color = Continent), size = 2.5, alpha = 0.25) +
  scale_color_manual(values = continent_colors) +
  theme_classic(base_size = 22) +
  theme(legend.position = "bottom") +
  labs(title = "PCA: CaG + HGDP/1KG (LASER ProPC)",
       x = "PC1", y = "PC2", color = "Continent")


ggsave("plot_PCA/PCA_PC1_vs_PC2_from_final_ProPCcoord.png", p, width = 12, height = 10, dpi = 300)


# --- Facets: REF alone vs REF + CaG (replacement for your p_facet block) ---
ref <- dt[Dataset == "HGDP+1KG"]
cag <- dt[Dataset == "CaG"]

ref_only <- copy(ref)[, Facet := "REF only"]
ref_plus <- rbindlist(list(ref, cag), use.names = TRUE)[, Facet := "REF + CaG"]
df_plot  <- rbindlist(list(ref_only, ref_plus), use.names = TRUE)
df_plot[, Facet := factor(Facet, levels = c("REF only", "REF + CaG"))]

# lock axes across facets for apples-to-apples comparison
xlim <- range(dt$PC1, na.rm = TRUE)
ylim <- range(dt$PC2, na.rm = TRUE)


p_facets <- ggplot(df_plot, aes(PC1, PC2, color = Continent)) +
  # Left facet: REF only
  geom_point(data = df_plot[Facet == "REF only"], size = 2.4, alpha = 0.85) +
  # Right facet: REF + CaG (REF solid, CaG faint)
  geom_point(data = df_plot[Facet == "REF + CaG" & Dataset == "HGDP+1KG"], size = 2.4, alpha = 0.85) +
  geom_point(data = df_plot[Facet == "REF + CaG" & Dataset == "CaG"],       size = 2.4, alpha = 0.25) +
  scale_color_manual(values = continent_colors, drop = FALSE) +
  theme_classic(base_size = 20) +
  theme(legend.position = "bottom") +
  labs(title = "PCA facets: REF only vs REF + CaG",
       x = "PC1", y = "PC2", color = "Continent") +
  facet_wrap(~ Facet, nrow = 1) +
  coord_cartesian(xlim = xlim, ylim = ylim)

ggsave("plot_PCA/PCA_PC1_vs_PC2_facets_REF_and_REFplusCaG.png",
       p_facets, width = 16, height = 8, dpi = 300)






# ---------------------------- KNN: find EUR CaG samples (strict majority) ----------------------------
# Requirements: dt, continent_colors already defined above
# ----------------------- kNN ancestry with threshold & breakdown (for CaG) -----------------------
library(FNN)

# Set parameters for KNN model
standardize_pcs <- FALSE  # optional


# ---------------- kNN ancestry with variable #PCs + threshold + breakdown ----------------
# Expects: dt with Dataset/Continent + PC1.., and 'continent_colors' defined.

#suppressPackageStartupMessages(require(FNN))

# --- Parameters (override earlier in your script if you like) ---
standardize_pcs <- get0("standardize_pcs", ifnotfound = FALSE)  # z-score PCs by REF

# --- Choose PC columns (PC1..PC{n_pcs}, or fewer if unavailable) ---
all_pc_names <- grep("^PC\\d+$", names(dt), value = TRUE)
max_pc_idx   <- min(n_pcs, length(all_pc_names))
pc_cols      <- paste0("PC", seq_len(max_pc_idx))

ref <- dt[Dataset == "HGDP+1KG"]
cag <- dt[Dataset == "CaG"]
stopifnot(nrow(ref) > 0, nrow(cag) > 0)



# Sanity / clamping
#if (nrow(ref) < k) { warning(sprintf("REF has %d samples; k -> %d.", nrow(ref), nrow(ref))); k <- nrow(ref) }
#if (threshold_N > k) { warning(sprintf("threshold_N (%d) > k (%d); threshold_N -> %d.", threshold_N, k, k)); threshold_N <- k }
#if (max_pc_idx < 1L) stop("No PC columns found in dt (expected PC1..).")

ref_mat <- as.matrix(ref[, ..pc_cols])
cag_mat <- as.matrix(cag[, ..pc_cols])

# Optional: standardize PCs by REF stats so each PC contributes equally
if (standardize_pcs) {
  mu  <- colMeans(ref_mat)
  sdv <- pmax(apply(ref_mat, 2, sd), 1e-9)
  ref_mat <- sweep(sweep(ref_mat, 2, mu, "-"), 2, sdv, "/")
  cag_mat <- sweep(sweep(cag_mat, 2, mu, "-"), 2, sdv, "/")
}


# --- kNN in PC space ---
nn       <- FNN::get.knnx(data = ref_mat, query = cag_mat, k = k)
nbr_idx  <- nn$nn.index                   # [n_query x k]
nbr_dist <- nn$nn.dist



# Neighbor continent labels as [n_query x k] matrix
ref_cont   <- as.character(ref$Continent)
labels_mat <- matrix(ref_cont[as.vector(nbr_idx)],
                     nrow = nrow(nbr_idx), ncol = ncol(nbr_idx), byrow = FALSE)


# Stable REF continent order (exclude "CaG")
ref_levels <- setdiff(names(continent_colors), "CaG")

# Per-continent neighbor COUNTS
counts_mat <- sapply(ref_levels, function(cc) rowSums(labels_mat == cc, na.rm = TRUE))
counts_mat <- as.matrix(counts_mat); colnames(counts_mat) <- ref_levels



# Assignment rule: if any continent has count >= threshold_N -> that label; else "Undefined".
# If multiple satisfy (possible when threshold_N <= floor(k/2)), tie-break by inverse-distance weights.
eps <- 1e-9
w   <- 1 / (nbr_dist + eps)
weighted_mat <- sapply(ref_levels, function(cc) rowSums((labels_mat == cc) * w, na.rm = TRUE))
weighted_mat <- as.matrix(weighted_mat); colnames(weighted_mat) <- paste0(ref_levels, "_w")




pred_pop <- character(nrow(cag))
for (i in seq_len(nrow(cag))) {
  row_counts <- counts_mat[i, ]
  hits <- names(row_counts)[row_counts >= threshold_N]
  if (length(hits) == 0L) {
    pred_pop[i] <- "Undefined"
  } else if (length(hits) == 1L) {
    pred_pop[i] <- hits
  } else {
    wrow  <- weighted_mat[i, paste0(hits, "_w"), drop = TRUE]
    hits2 <- hits[wrow == max(wrow)]
    pred_pop[i] <- if (length(hits2) == 1L) hits2 else "Ambiguous"
  }
}

# Build output with counts + proportions
counts_dt <- as.data.table(counts_mat)
setnames(counts_dt, ref_levels, paste0(ref_levels, "_count"))
props_dt <- counts_dt[, lapply(.SD, function(x) x / k)]
setnames(props_dt, names(props_dt), sub("_count$", "_prop", names(props_dt)))

result <- cbind(
  cag[, .(FID, IID, PC1, PC2)],
  Predicted_Pop = pred_pop,
  counts_dt,
  props_dt
)

outfile <- sprintf("plot_PCA/CaG_KNN_breakdown_k%d_thresh%d_PC%d.tsv", k, threshold_N, max_pc_idx)
fwrite(result, outfile, sep = "\t")

# Quick summary
print(result[, .N, by = Predicted_Pop][order(-N)])



# ---- CaG-only PCA plot colored by inferred population ----
# Join predictions onto CaG rows
dt_cag <- merge(
  dt[Dataset == "CaG"],
  result[, .(IID, Predicted_Pop)],
  by = "IID", all.x = TRUE
)

# Any missing or Ambiguous prediction -> "Undefined"
dt_cag[is.na(Predicted_Pop), Predicted_Pop := "Undefined"]
dt_cag[Predicted_Pop == "Ambiguous", Predicted_Pop := "Undefined"]

# Build palette: reuse REF continent colors; add grey for Undefined
ref_levels <- setdiff(names(continent_colors), "CaG")
palette_ext <- c(
  continent_colors[ref_levels],
  "Undefined" = "grey70"
)

# Factor levels (no Ambiguous), then drop unused to keep legend clean
dt_cag[, Predicted_Pop := factor(Predicted_Pop, levels = c(ref_levels, "Undefined"))]
dt_cag[, Predicted_Pop := droplevels(Predicted_Pop)]

# (Optional) pick up k / n_pcs / threshold for title if they exist
k            <- get0("k", ifnotfound = NA_integer_)
threshold_N  <- get0("threshold_N", ifnotfound = NA_integer_)
n_pcs_used   <- get0("max_pc_idx", ifnotfound = {
  all_pc_names <- grep("^PC\\d+$", names(dt), value = TRUE)
  min(get0("n_pcs", ifnotfound = 10L), length(all_pc_names))
})

title_txt <- if (is.na(k) || is.na(threshold_N)) {
  "CaG only — colored by inferred population"
} else {
  sprintf("CaG only — inferred population (k=%d, PCs=1-%d, thresh=%d)",
          k, n_pcs_used, threshold_N)
}

p_cag_only <- ggplot(dt_cag, aes(PC1, PC2, color = Predicted_Pop)) +
  geom_point(size = 2.4, alpha = 0.95) +
  scale_color_manual(values = palette_ext) +  # drop=TRUE by default
  theme_classic(base_size = 20) +
  theme(legend.position = "bottom") +
  labs(title = title_txt, x = "PC1", y = "PC2", color = "Inferred pop")

ggsave("plot_PCA/PCA_CaG_only_by_inferred_population.png",
       p_cag_only, width = 12, height = 10, dpi = 300)



# ---- CaG-only scatterplot matrix PC1..PC5 ----

# Color vector matching Predicted_Pop
cols <- setNames(palette_ext, names(palette_ext))
pt_col <- cols[as.character(dt_cag$Predicted_Pop)]

# pairs(): lower triangle only, blank diag/upper
png("plot_PCA/PCA_CaG_scattermatrix_PC1toPC5_base.png", width = 1200, height = 1200, res = 120)
pairs(
  dt_cag[, ..pc_cols],
  col = pt_col, pch = 16,
  upper.panel = NULL,
  diag.panel  = NULL,
  cex = 0.6
)
legend("topright", inset = 0.02,
       legend = levels(dt_cag$Predicted_Pop),
       col = cols[levels(dt_cag$Predicted_Pop)], pch = 16, cex = 0.9, bty = "n")
dev.off()

