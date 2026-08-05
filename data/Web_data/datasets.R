# datasets.R

library(readxl)
library(dplyr)

# setwd("Pezzini")
# list.files()

# Create dirs
dirs <- c("gProfiler", "STRING", "GSEA", "Reactome")
lapply(dirs, function(d) {
    if (!dir.exists(d)) dir.create(d)
})

# Read Pezzini data
dat <- readxl::read_xlsx("Pezzini2016_SHSY5Ycelldiff_DE_table_filtering.xlsx", sheet = "Full_DEResults_Table")
names(dat)[1:4] <- c("Ensemble_ID", "Gene_Symbol", "Undifferentiated", "logFC")

# Data sets with Ensemble IDs
ens1 <- dat %>%
    filter(abs(logFC) > log2(2) & FDR < 0.0001) %>%
    arrange(FDR) %>%
    pull(Ensemble_ID)
write.table(ens1, file = "gProfiler/DE_198set_Enemble_IDs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ens1, file = "STRING/DE_198set_Enemble_IDs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
ens2 <- dat %>%
    filter(abs(logFC) > log2(8) & FDR < 0.01) %>%
    arrange(FDR) %>%
    pull(Ensemble_ID)
write.table(ens2, file = "gProfiler/DE_378set_Enemble_IDs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ens2, file = "Reactome/DE_378set_Enemble_IDs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
ens3 <- dat %>%
    arrange(Ensemble_ID) %>%
    pull(Ensemble_ID)
write.table(ens3, file = "gProfiler/Bg_14420set_Enemble_IDs.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Data sets with Gene symbols
sym1 <- dat %>%
    filter(abs(logFC) > log2(2) & FDR < 0.0001) %>%
    arrange(FDR) %>%
    pull(Gene_Symbol)
write.table(sym1, file = "gProfiler/DE_198set_Gene_Symbols.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(sym1, file = "STRING/DE_198set_Gene_Symbols.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
sym2 <- dat %>%
    filter(abs(logFC) > log2(8) & FDR < 0.01) %>%
    arrange(FDR) %>%
    pull(Gene_Symbol)
write.table(sym2, file = "gProfiler/DE_378set_Gene_Symbols.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(sym2, file = "Reactome/DE_378set_Gene_Symbols.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
sym3 <- dat %>%
    arrange(Gene_Symbol) %>%
    pull(Gene_Symbol)
write.table(sym3, file = "gProfiler/Bg_14420set_Gene_Symbols.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# -------------------------------------------------------------------------------------------------------
# Generate TMM normalised data for Reactome analysis
library(limma)
library(edgeR)
library(jsonlite)
library(dplyr)
library(tibble)

# counts_file <- 'Pezzini2016_SHSY5Ycelldiff_NonStrandedCounts-withNames-proteinCoding.txt'
counts_file <- "SHSY5Ycelldiff_Pezzini2016.tsv"
output_dir <- "output_dir"
count_cols <- c("SRR3133925_nodiff_rep1", "SRR3133926_nodiff_rep2", "SRR3133927_nodiff_rep3", "SRR3133928_diff_rep1", "SRR3133929_diff_rep2", "SRR3133930_diff_rep3")
design <- matrix(c(c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 1)), ncol = 2, dimnames = list(c("SRR3133925_nodiff_rep1", "SRR3133926_nodiff_rep2", "SRR3133927_nodiff_rep3", "SRR3133928_diff_rep1", "SRR3133929_diff_rep2", "SRR3133930_diff_rep3"), c("Undifferentiated ", "Differentiated")))
cont.matrix <- matrix(c(c(-1, 1)), ncol = 1, dimnames = list(c("Undifferentiated ", "Differentiated"), c("Differentiated")))
export_cols <- c(c("Gene.ID", "Gene.Name", "SRR3133925_nodiff_rep1", "SRR3133926_nodiff_rep2", "SRR3133927_nodiff_rep3", "SRR3133928_diff_rep1", "SRR3133929_diff_rep2", "SRR3133930_diff_rep3"))

# Maybe filter out samples that are not used in the model
if (FALSE) {
    # Remove columns not used in the comparison
    use.samples <- rowSums((design %*% cont.matrix) != 0) > 0
    use.conditions <- colSums(design[use.samples, ] != 0) > 0

    count_cols <- count_cols[use.samples, drop = F]
    design <- design[use.samples, use.conditions, drop = F]
    cont.matrix <- cont.matrix[use.conditions, , drop = F]
}

# fileEncoding='UTF-8-BOM' should strip the BOM marker FEFF that some windows tools add
x <- read.delim(counts_file,
    sep = "\t",
    check.names = FALSE,
    colClasses = "character",
    na.strings = c(),
    skip = 0,
    fileEncoding = "UTF-8-BOM"
)

# Now re-read the first header line.  Workaround R problem that always has strip.white=T for read.table
colnames(x) <- scan(counts_file,
    what = "",
    sep = "	",
    nlines = 1,
    strip.white = F,
    quote = "\"",
    skip = "0",
    fileEncoding = "UTF-8-BOM"
)

# Check which rows have dodgy values in the counts columns
non_numeric_rows <- apply(x[, count_cols], 1, function(row) {
    any(is.na(suppressWarnings(as.numeric(row))))
})

if (any(non_numeric_rows)) {
    warning("Bad CSV, non-numeric counts in these rows : ", toString(paste(which(non_numeric_rows))))
}

x[, count_cols] <- apply(x[, count_cols], 2, function(v) as.numeric(v)) # Force numeric count columns
counts <- x[, count_cols]

# Keep rows based on string based filters of columns.  Rows must match all filters
filter_rows <- fromJSON("[]")
if (length(filter_rows) > 0) {
    keepRows <- apply(apply(filter_rows, 1, function(r) grepl(r["regexp"], x[, r["column"]], perl = T, ignore.case = T)), 1, all)
} else {
    keepRows <- rep(TRUE, nrow(x))
}

keepMin <- apply(counts, 1, max) >= 10.0
keepCpm <- rowSums(cpm(counts) > 0.0) >= 0 # Keep only genes with cpm above x in at least y samples

keep <- keepMin & keepCpm & keepRows
x <- x[keep, ]
counts <- counts[keep, ]



# nf <- calcNormFactors(counts)
# y <- voom(counts, design, plot = FALSE, lib.size = colSums(counts) * nf)
nf <- calcNormFactors(counts, method = "TMM")
y <- voom(counts, design, plot = FALSE, lib.size = colSums(counts) * nf)

y$genes <- data.frame(
    Gene.ID   = x$"Gene.ID",
    Gene.Name = x$"Gene.Name"
)

fit <- lmFit(y, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

out <- topTable(fit2, n = Inf, sort.by = "none")


# TMM Normalised -------------------------------------------------------------------------------------------------------
y$genes <- x$"Gene.ID"
normalized <- y$E %>% as.data.frame()
# grp1 <- grep("_diff_", colnames(normalized))
# grp2 <- grep("_nodiff_", colnames(normalized))
# normalized <- normalized[, c(grp1, grp2)] # Re-ordering exp data
row.names(normalized) <- y$genes
normalized <- normalized %>% rownames_to_column(var = "Ensemble_IDs")
write.table(normalized, file = "Reactome/TMM_normalized_data_Reactome.csv", row.names = FALSE, sep = ",", na = "", quote = FALSE)
# ----------------------------------------------------------------------------------------------------------------------


# Generate input data fro GSEA analysis
# Get the ranks --------------------------------------------------------------------------------------------------------
rnk1 <- out %>%
    arrange(-t) %>%
    dplyr::select(Gene.ID, t) %>%
    rename("Ensemble_ID" = "Gene.ID", "t-stat" = "t")
rnk2 <- out %>%
    arrange(-t) %>%
    dplyr::select(Gene.Name, t) %>%
    rename("Gene_Symbol" = "Gene.Name", "t-stat" = "t")
write.table(rnk1, file = "STRING/Pre_Ranked_List_t_stat_ENSEMBL", row.names = FALSE, sep = ",", na = "", quote = FALSE)
write.table(rnk2, file = "STRING/Pre_Ranked_List_t_stat_SYMBOL", row.names = FALSE, sep = ",", na = "", quote = FALSE)
rnk3 <- out %>%
    arrange(-logFC) %>%
    dplyr::select(Gene.ID, logFC) %>%
    rename("Ensemble_ID" = "Gene.ID", "logFC" = "logFC")
rnk4 <- out %>%
    arrange(-logFC) %>%
    dplyr::select(Gene.Name, logFC) %>%
    rename("Gene_Symbol" = "Gene.Name", "logFC" = "logFC")
write.table(rnk3, file = "STRING/Pre_Ranked_List_logFC_ENSEMBL", row.names = FALSE, sep = ",", na = "", quote = FALSE)
write.table(rnk4, file = "STRING/Pre_Ranked_List_logFC_SYMBOL", row.names = FALSE, sep = ",", na = "", quote = FALSE)
# ----------------------------------------------------------------------------------------------------------------------


# Create .gct file -----------------------------------------------------------------------------------------------------
nf <- calcNormFactors(counts, method = "TMM")
y <- voom(counts, design, plot = FALSE, lib.size = colSums(counts) * nf)

y$genes <- x$"Gene.ID"
normalized <- y$E
row.names(normalized) <- y$genes

normalizedDF <- cbind(description = "NA", normalized)
normalizedDF <- data.frame(normalizedDF) %>% rownames_to_column(var = "NAME")

# write.table(normalizedDF, file = "TMM_normalized.tsv", row.names = FALSE, sep = "\t", na = "", quote = FALSE)
# expressed_data <- readLines("TMM_normalized.tsv")
# OR
# Convert it to a character vector like a file's lines
tmp_text <- capture.output(
    write.table(normalizedDF, row.names = FALSE, sep = "\t", na = "", quote = FALSE)
)
conn <- textConnection(tmp_text)
expressed_data <- readLines(conn)
close(conn)


new_rows <- c("#1.2", paste(nrow(normalized), ncol(normalized), sep = "\t"))

GCT_FILE <- c(new_rows, expressed_data)
writeLines(GCT_FILE, "GSEA/Expression_Data.gct")


# Create .cls file
sampls <- ncol(normalized)
groups <- 2
labels <- c(rep("Nodiff", 3), rep("Diff", 3))
uniq <- unique(labels)
cls_content <- c(
    paste(sampls, groups, 1), # Number of samples, number of classes, and type
    paste("#", paste(uniq, collapse = " ")), # Class labels
    paste(labels, collapse = " ") # Phenotype labels for each sample
)
writeLines(cls_content, con = "GSEA/Phenotype_Labels.cls")
# ----------------------------------------------------------------------------------------------------------------------


# # No neeed to Run the following code for the current purpose, but kept for reference
# out2 <- cbind(
#     fit2$coef,
#     out[, c("P.Value", "adj.P.Val", "AveExpr")],
#     x[, export_cols]
# )

# if (FALSE) {
#     library(topconfects)
#     confect <- limma_confects(fit2, coef = 1, fdr = 0.0)
#     out2 <- cbind(out2, confect = confect$table$confect[order(confect$table$index)])
# }

# param_normalized <- ""
# normalized <- matrix(0)
# if (param_normalized == "") {
#     write.csv(out2, file = paste0(output_dir, "/output.txt"), row.names = FALSE, na = "")
# } else if (param_normalized == "backend") {
#     normalized <- y$E
# } else if (param_normalized == "remove-hidden") {
#     hidden_factors <- c()
#     if (length(hidden_factors) > 0) {
#         # Remove the batch effect (as done in removeBatchEffect)
#         beta <- fit$coefficients[, hidden_factors, drop = FALSE]
#         normalized <- as.matrix(y) - beta %*% t(design[, hidden_factors, drop = FALSE])
#     }
# }

# cat(
#     toJSON(list(
#         rank = fit2$rank, df_prior = fit2$df.prior,
#         design = data.frame(fit2$design), contrasts = data.frame(fit2$contrasts),
#         cov_coefficients = data.frame(fit2$cov.coefficients),
#         normalized = list(columns = colnames(normalized), values = normalized)
#     )),
#     file = paste0(output_dir, "/extra.json")
# )
