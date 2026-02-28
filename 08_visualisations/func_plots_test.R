func_table <- read.csv("07_functional_annotation/annotated_46_var.csv", header = TRUE, stringsAsFactors = FALSE)

signals <- c("chr7:129013434:C:CAAAA", "chr7:128984718:T:C", "chr7:128930474:CA:C")

coloc_eqtl <- func_table[func_table$category == "Coloc_eQTL", "variant_id"]

top_cs <- func_table[func_table$category == "Credible_set", ]
top_cs <- top_cs[order(-top_cs$PIP), ][1:10, "variant_id"]

priority_ids <- c(signals, coloc_eqtl, top_cs)
func_priority <- func_table[func_table$variant_id %in% priority_ids, ]

nrow(func_priority)
func_priority[, c("variant_id", "rsid", "category", "PIP")]

keep <- c(
  "chr7:129013434:C:CAAAA", "chr7:128984718:T:C", "chr7:128930474:CA:C",
  "chr7:128937860:C:CGCGGG", "chr7:128936032:T:C",
  "chr7:128932712:G:A", "chr7:128935498:C:G",
  "chr7:129034894:C:CAAA", "chr7:129046195:CA:C", "chr7:128958137:C:CTTTTT"
)
func_priority <- func_table[func_table$variant_id %in% keep, ]
nrow(func_priority)

func_priority$label <- ifelse(func_priority$rsid == "." | is.na(func_priority$rsid),
                              func_priority$variant_id,
                              func_priority$rsid)
func_priority$IRF5_eQTL <- ifelse(is.na(func_priority$IRF5_eQTL), FALSE, func_priority$IRF5_eQTL)

# fix empty strings to NA for cleaner display
func_priority[func_priority == ""] <- NA

# columns to plot
annot_cols <- c("cCRE_type", "ChromHMM_E029_Mono", "ChromHMM_E062_PBMC",
                "RegulomeDB_rank", "CADD_PHRED", "IRF5_eQTL")

# reshape to long format
library(tidyr)
func_priority[annot_cols] <- lapply(func_priority[annot_cols], as.character)
long <- pivot_longer(func_priority,
                     cols = all_of(annot_cols),
                     names_to = "annotation",
                     values_to = "value")
long$annotation <- factor(long$annotation, levels = annot_cols)
long$label <- factor(long$label,
                     levels = func_priority$label[order(func_priority$position)])

# separate continuous from categorical
long$is_cadd <- long$annotation == "CADD_PHRED"
long$cadd_val <- as.numeric(ifelse(long$is_cadd, long$value, NA))
long$cat_val <- ifelse(long$is_cadd, NA, long$value)

p_func <- ggplot(long, aes(x = annotation, y = label)) +
  geom_tile(data = subset(long, is_cadd),
            aes(fill = cadd_val), colour = "white") +
  scale_fill_gradient(low = "white", high = "#8B0000",
                      na.value = "grey90", name = "CADD\nPHRED") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(long, !is_cadd),
            aes(fill = cat_val), colour = "white") +
  scale_fill_discrete(na.value = "grey90", name = "Category") +
  geom_text(aes(label = ifelse(is_cadd, round(cadd_val, 1), cat_val)),
            size = 2.5, na.rm = TRUE) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

print(p_func)