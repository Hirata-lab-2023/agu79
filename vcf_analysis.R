library(dplyr)
library(ggplot2)
library(openxlsx)

# Define output directory
family <- "agu79"

# Load sequencing depth information
depth_path <- sprintf("%s/depth/", family)
depth_files <- list.files(depth_path, pattern = "new_depth.txt", full.names = TRUE)

depth <- sapply(depth_files, function(file) {
  dat <- read.table(file, header = FALSE)
  dat <- dat[-26, ]
  sum(dat$V3) / sum(dat$V4)
})

depth_df <- data.frame(depth = depth)
rownames(depth_df) <- c("hetero", paste0("homo_", 1:5))
print(depth_df)

# Load and parse gene annotation from GTF file
gtf <- read.table("Danio_rerio.GRCz11.109.gff3", sep = "\t", header = FALSE)
gene <- gtf[gtf$V3 == "exon", ]
gene_all <- do.call(rbind, lapply(1:nrow(gene), function(i) {
  gene_1 <- unlist(strsplit(as.character(gene$V9[i]), "; "))
  gene_id <- unlist(strsplit(gene_1[6], " "))[2]
  transcript_id <- unlist(strsplit(gene_1[3], " "))[2]
  cbind(gene_id, transcript_id)
}))

gene_all <- as.data.frame(gene_all)
gene_all_new <- gene_all %>% distinct(V1, V2, .keep_all = TRUE)
colnames(gene_all_new) <- c("Gene_name", "Transcript")

# Count number of splicing variants per gene
Splicing_variant_name <- as.data.frame(table(gene_all_new$Gene_name))
colnames(Splicing_variant_name) <- c("Gene_name", "Splicing_variant")

# Load HIGH impact SNPs from heterozygous sample
HIGH_path <- sprintf("%s/vcf/HIGH", family)
HIGH_list <- list.files(HIGH_path, pattern = "ano_all_HIGH.vcf.gz$", full.names = TRUE)

hetero_high <- read.table(HIGH_list[1], sep = "\t", header = FALSE)
hetero_high <- hetero_high[grep("0/1", hetero_high$V10), ]

# Extract annotation fields
extract_ann <- function(vcf_row, impact_type) {
  p <- cbind(vcf_row$V1, vcf_row$V2, vcf_row$V4, vcf_row$V5)
  ann_field <- unlist(strsplit(vcf_row$V8, "\\|") )
  i <- grep(impact_type, ann_field)
  res <- NULL
  for (n in i) {
    values <- unlist(strsplit(vcf_row$V8, "\\|"))
    selected <- values[c(n-1, n+1, n+4, n+7, n+8, n+10)]
    res <- rbind(res, cbind(p, t(selected)))
  }
  return(res)
}

# Process HIGH variants in homozygous samples
homo_high_list <- list()
for (i in 2:6) {
  vcf <- read.table(HIGH_list[i], sep = "\t", header = FALSE)
  vcf <- vcf[grep("1/1", vcf$V10), ]
  ann <- do.call(rbind, lapply(1:nrow(vcf), function(j) extract_ann(vcf[j, ], "HIGH")))
  homo_high_list[[i - 1]] <- as.data.frame(unique(ann))
}

# Load MODERATE impact SNPs from heterozygous sample
MOD_path <- sprintf("%s/vcf/MOD", family)
MOD_list <- list.files(MOD_path, pattern = "ano_all_MOD.vcf.gz$", full.names = TRUE)

hetero_mod <- read.table(MOD_list[1], sep = "\t", header = FALSE)
hetero_mod <- hetero_mod[grep("0/1", hetero_mod$V10), ]

# Process MODERATE variants in homozygous samples
homo_mod_list <- list()
for (i in 2:6) {
  vcf <- read.table(MOD_list[i], sep = "\t", header = FALSE)
  vcf <- vcf[grep("1/1", vcf$V10), ]
  ann <- do.call(rbind, lapply(1:nrow(vcf), function(j) extract_ann(vcf[j, ], "MODERATE")))
  homo_mod_list[[i - 1]] <- as.data.frame(unique(ann))
}

# Save example output
write.xlsx(homo_high_list[[1]], file = sprintf("%s_HOMO_HIGH_sample.xlsx", family))
write.xlsx(homo_mod_list[[1]], file = sprintf("%s_HOMO_MOD_sample.xlsx", family))
