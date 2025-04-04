#Aim1

#input the dataset
setwd("D:/analysis/project")
raw_data<-read.csv("D:/analysis/project/data_GSE255720.csv", header = TRUE, sep = ",")

library(tidyr) 
library(tidyverse)
library(dplyr)

#data wrangling
#change column names and remove additional rows and columns

colnames(raw_data)<- raw_data[1,]
clean_data<- raw_data[!(rownames(raw_data) == "Sample_type"), ]
colnames(clean_data)[2]<-"gene_symbol"
clean_data<- clean_data[,-1]
clean_data<- clean_data[,-2]
str(clean_data)

#remove not available in gene_symbol column
clean_data[clean_data==""]<-NA
clean_data <- drop_na(clean_data, gene_symbol)

#replace the gene_symbol with ENSEMBL_ID
ensembl_data<- read.csv("D:/analysis/project/ensembl_data.csv",header = TRUE, sep = "," )
clean_data_mix <- clean_data %>%
  left_join(ensembl_data, by = "gene_symbol")
any(is.na(clean_data_mix))
clean_data_mix <- drop_na(clean_data_mix, ENSEMBL.ID)
clean_data_mix <- clean_data_mix %>% select(-gene_symbol)

#determine ENSEMBL_ID as the row names
rownames(clean_data_mix)<- clean_data_mix$ENSEMBL.ID
data<-clean_data_mix %>% select(-ENSEMBL.ID)

# converted the count data values to integers and removed genes with a total count lower than 10
data <- data %>% mutate(across(everything(), as.numeric)) %>% round(0)
count<- data[which(rowSums(data)>10),]

``````````````````````````````````````````````````````````````
#B cell-related genes visualized in a heatmap 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
install.packages("limma")
library(limma)
install.packages("RColorBrewer")
library(RColorBrewer)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
install.packages("pheatmap")
library(pheatmap)

#x is my count data in which columns are samples(UC, non.IBD) and rows are genes 
x<- clean_data
logcounts <- log2(x + 1)
target_genes <- c( "CD19"  ,  "CD79B" ,"CD79A", "SDC1", "PTPRC") # select target genes
mypalette <- brewer.pal(11, "RdYlBu")
condition<- factor(c( "UC"   ,    "UC"  ,   "non.IBD"  , "non.IBD", "non.IBD" ,"UC"  ,   "UC"   ,  "non.IBD", "UC"  ,   "non.IBD")) #determine condition for columns
morecols <- colorRampPalette(mypalette)
sorted_indices <- order(condition)
logcounts_sorted <- logcounts[target_genes, sorted_indices]
condition_sorted <- factor(condition[sorted_indices], levels = c("UC", "non.IBD"))
annotation_col <- data.frame(Group = condition_sorted)
rownames(annotation_col) <- colnames(logcounts_sorted)

pheatmap(logcounts_sorted,
         color = rev(morecols(50)),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,  
         main = "Expression of B cell-related genes",
         scale = "row",
         fontsize = 12,
         fontface_row = "bold",
         border_color = "darkgray",
         legend = TRUE,
         legend_breaks = seq(-2, 2, 0.5),
         annotation_legend = TRUE,
         legend_position = "left",
         annotation_legend_side = "left",
         gaps_col = NULL,
         fontsize_col = 12,
         cellwidth = 25,
         cellheight = 25,
         margins = c(10, 10),
         show_colnames = TRUE
)
``````````````````````````````````````````````````````````````
#Gene Set Enrichment Analysis (GSEA) visualized in a bar chart 
#Differential gene expression
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2",force = TRUE )
library(DESeq2)
count<- data[which(rowSums(data)>10),]
condition<- factor(c( "UC"   ,    "UC"  ,   "non.IBD"  , "non.IBD", "non.IBD" ,"UC"  ,   "UC"   ,  "non.IBD", "UC"  ,   "non.IBD")) #determine condition for columns
coldata<- data.frame(row.names = colnames(count), condition)
dds<- DESeqDataSetFromMatrix (countData = count, colData = coldata, design = ~condition)
dds<-DESeq(dds)
# Compare UC vs non.IBD
res <- results(dds, contrast = c("condition", "UC", "non.IBD"))
sigs<-na.omit(res)
sigs<-sigs[sigs$padj<0.05,]

#GSEA
BiocManager:: install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
res1<- res[res$baseMean>10,]
res1<-res1[order(-res1$log2FoldChange),]
gene_list<- res1$log2FoldChange
names(gene_list) <- rownames(res1)               
gene_list <- sort(gene_list, decreasing = TRUE)  
head(gene_list)

library(clusterProfiler)
library(org.Hs.eg.db)
library(grid) 

gse<-gseGO(gene_list,
           ont = "BP",
           keyType = "ENSEMBL",
           OrgDb = "org.Hs.eg.db", 
           eps = 1e-300)
#GSEA bar chart
target_description <- c("antimicrobial humoral response",
                        "immunoglobulin production",
                        "humoral immune response", 
                        "B cell mediated immunity",
                        "T cell proliferation",
                        "T cell differentiation")

# Filter the rows based on the target descriptions
filtered_data <- gse_data[gse_data$Description %in% target_description, ]

# Create the bar chart using the filtered data
library(ggplot2)
##LINE bar chart
ggplot(data = filtered_data, aes(y = reorder(Description, enrichmentScore))) +
  geom_segment(aes(x = 0, xend = enrichmentScore, yend = Description), 
               color = "red", size = 1.5) +  
  geom_point(aes(x = enrichmentScore), 
             color = "red", size = 6) +  
  geom_text(aes(x = max(enrichmentScore) + 0.3, 
                label = format(pvalue, scientific = TRUE)), 
            hjust = 0, size = 5, color = "black") +  
  annotate("text", 
           x = max(filtered_data$enrichmentScore) + 0.4, 
           y = nrow(filtered_data) + 0.5, 
           label = "P-value", 
           hjust = 0, size = 4.5, fontface = "bold", color = "black") +  
  labs(
    title = " GSEA for UC vs non-IBD patients",
    x = "Normalised Enrichment Score",
    y = "Gene Set Description"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13,face = "bold", color = "black"), 
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 14, face = "bold"),     
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_x_continuous(expand = c(0.1, 0.1)) +  
  coord_cartesian(clip = "off")

``````````````````````````````````````````````````````````````
# Immunoglobulin gene expression visualised in a Volcano plot 
install.packages("ggplot2")
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
install.packages("textshaping")
library(textshaping)
library(org.Hs.eg.db)
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
sigs<-na.omit(res)
sigs.df<- as.data.frame(sigs)
sigs.df
sigs.df$symbol<- mapIds(org.Hs.eg.db, keys = row.names(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")
##indicate target genes:
target_genes <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4","IGHA1", "IGHA2", "IGHM", "IGHD", "IGHE" )
sigs.df$label <- ifelse(sigs.df$symbol %in% target_genes, sigs.df$symbol, NA)

EnhancedVolcano(
  sigs.df,
  x = "log2FoldChange",
  y = "pvalue",
  lab =NA,  
  pCutoff = 1e-4,       
  FCcutoff = 2,         
  title = "Volcano Plot",
  subtitle = "Immunoglobulin Differential Genes Expression",
  caption = "Condition: UC vs non-IBD",
  pointSize = 4.0,      
  labSize = 6.0        
) +
  # Adjust label positions for target genes
  geom_text_repel(
    data = subset(sigs.df, symbol %in% target_genes),
    aes(x = log2FoldChange, y = -log10(pvalue), label = symbol),
    size=6, 
    box.padding = 0.5,       
    point.padding = 0.3,     
    segment.color = "yellow",
    segment.size = 0.5,      
    max.overlaps = Inf       
  )+
  # Add halo around target genes
  geom_point(
    data = subset(sigs.df, symbol %in% target_genes),
    aes(x = log2FoldChange, y = -log10(pvalue)),
    size = 1.5,             
    shape = 21,             
    fill = "transparent",    
    color = "black",         
    stroke = 2               
  ) 


``````````````````````````````````````````````````````````````

#Aim2

#input the dataset
setwd("D:/analysis/project")
microarray_data<-read.csv("D:/analysis/project/GSE206171_Array_data_GEO.csv", header = TRUE, sep = ",")
#Violin plot for log2FC of FDCSP gene expression
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2",force = TRUE )
library(DESeq2)

#Data Wrangling
second_data<-2^(microarray_data)
rownames(second_data)<- second_data[,1]
count<- second_data[which(rowSums(data)>10),]
condition <- factor(c(rep("Control", 38), rep("UC", 114)))
coldata<- data.frame(row.names = colnames(count), condition)
dds<- DESeqDataSetFromMatrix (countData = count, colData = coldata, design = ~condition)
dds<-DESeq(dds)

``````````````````````````````````````````````````````````````
#Violin Plot of FDCSP Expression

library(ggplot2)
install.packages("ggstatsplot")
library(ggstatsplot)
install.packages("sysfonts")
install.packages("showtextdb")
library(sysfonts)
library(showtextdb)
library(ggpubr)
library(DESeq2)

plotCounts(dds, gene="ENSG00000181617", )
plot_data <- plotCounts(dds, gene = "ENSG00000181617", intgroup = "condition", returnData = TRUE)
plot_data$log_count <- log2(plot_data$count + 1)

plt <- ggbetweenstats(
  data = plot_data,
  x = condition,
  y = log_count,
  results.subtitle = FALSE,  
  ggtheme = theme_minimal()  
)

# Customize plot appearance
plt <- plt +
  labs(
    x = "Condition",
    y = "Log2 Normalized Count",
    title = "FDCSP Expression in UC vs. Control"
  ) +
  theme(
    text = element_text(family = "Roboto", size = 10, color = "black"),
    plot.title = element_text(
      family = "Lobster Two",
      size = 20,
      face = "bold",
      color = "#2a475e",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  )


print(plt)
grid.text(
  "P_adj=7.06e-08", 
  x = 0.5, y = 0.85, 
  gp = gpar(fontsize = 14, fontface = "bold")
) 

``````````````````````````````````````````````````````````````
#Regression Analysis of FDCSP with demographic information
# Set the row names properly
row.names(second_data)[1] <- "Mayo_score"

# Extract the Mayo score (first row) and FDCSP counts (specific row)
Mayo_score <- as.numeric(second_data["Mayo_score", ])
FDCSP <- as.numeric(second_data["ENSG00000181617", ])

threshold <- median(FDCSP, na.rm = TRUE)  # You can adjust this threshold as needed
FDCSP_binary <- ifelse(FDCSP > threshold, 1, 0)

df <- data.frame(Mayo_score = Mayo_score, FDCSP_binary = FDCSP_binary) # Combine them into a new data frame for logistic regression
logit_model <- glm(FDCSP_binary ~ Mayo_score, family = "binomial", data = df) # Fit logistic regression model

summary_stats <- summary(logit_model)# Summary of the logistic regression model
# Extract key statistics for annotations
coeff <- coef(summary_stats)  # Coefficients table
p_value <- round(coeff[2, "Pr(>|z|)"], 4)  # P-value of Mayo_score
odds_ratio <- exp(coeff[2, 1])  # Odds ratio for Mayo_score

ggplot(df, aes(x = Mayo_score, y = FDCSP_binary)) +
  geom_point(size = 2, alpha = 0.7, color = "darkblue") +
  geom_smooth(
    formula = y ~ x,
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE,
    color = "red"
  ) +
  labs(
    title = "Logistic Regression: Mayo Score vs FDCSP Expression",
    x = "Mayo Score",
    y = "FDCSP Binary"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey", linetype = "dotted"),
    panel.grid.minor = element_blank()
  ) +
  annotate(
    "text",
    x = max(Mayo_score) * 0.5,
    y = 0.8,
    label = paste(
      "P-value:", p_value, 
      "\nOdds Ratio:", round(odds_ratio, 2)
    ),
    size = 5,
    color = "darkred",
    hjust = 0
  )



``````````````````````````````````````````````````````````````
# correlation of B cell related genes and FDCSP&TNF genes
# Heatmap of correlation between two gene groups
library(Hmisc)
library(pheatmap)
B_cell <- c("CD19", "CD79A", "CD79B", "SDC1", "PTPRC")
TNF_family <- c("TNFSF13B", "TNF", "FDCSP")

# Extract expression data (Make sure x_symbole is a matrix or dataframe with rownames as gene names)
B_cell_gene <- x_symbole[rownames(x_symbole) %in% B_cell, ]
TNF_family_gene <- x_symbole[rownames(x_symbole) %in% TNF_family, ]

# Convert to numeric matrices to avoid errors
B_cell_gene_t <- as.matrix(t(B_cell_gene))  
TNF_family_gene_t <- as.matrix(t(TNF_family_gene))  
B_cell_gene_t <- apply(B_cell_gene_t, 2, as.numeric)
TNF_family_gene_t <- apply(TNF_family_gene_t, 2, as.numeric)

# Compute correlation only between B cell and TNF family genes
correlation_matrix <- cor(B_cell_gene_t, TNF_family_gene_t, use = "pairwise.complete.obs")  

# Generate heatmap
pheatmap(
  correlation_matrix,  # Use correlation values
  color = colorRampPalette(c("navy", "white", "firebrick"))(50), 
  cluster_rows = FALSE,              
  cluster_cols = FALSE,              
  main = " Correlation of B cell Genes vs FDCSP&TNF Genes",  
  display_numbers = TRUE,            
  fontsize = 12,                     
  fontsize_row = 14,                 
  fontsize_col = 13,                 
  border_color = "black",             
  legend = TRUE,                      
  legend_breaks = seq(-1, 1, 0.5),    
  annotation_legend = TRUE,           
  cellwidth = 30,                     
  cellheight = 30,                    
  angle_col = 45                      
)
