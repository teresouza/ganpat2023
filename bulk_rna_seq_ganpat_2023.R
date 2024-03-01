##Script to process RNA-seq RMS datasets
meta.rseq <- read.csv("/Users/terezinhadesouza/surfdrive/PrinsesMaxima/Other/maroussia/rna/metadata.csv", 
                      header = T, row.names = 1) %>%
 filter(condition %in% c("Fusion", "Fusion_PDX", "Luciferase", "TFE3"))
basedir <- "/Users/terezinhadesouza/surfdrive/PrinsesMaxima/Other/maroussia/rna/"
biomart <- readRDS("/Users/terezinhadesouza/surfdrive/PrinsesMaxima/files_analysis/biomart.RDS")
rfGenes <- read.csv("/Users/terezinhadesouza/surfdrive/PrinsesMaxima/Sarcoma/NRSTSSeq/rfGenes.csv", header = T)
kidneyData <- read.delim("../TitoRaw_annotated.txt", header = T) %>% 
  select(ensembl_gene_id,external_gene_name, all_of(metaKidney$ID)) %>% melt %>% 
  reshape2::dcast(.,  external_gene_name ~ variable, value.var = "value", fun.aggregate = sum)

##Function to transform read counts to transcripts per million 9TOM 
tpm <- function(counts, lengths) {
  rate <- counts /(lengths/1000)
  rate / sum(rate) * 1e6
}

##Process RNA-sequencing feature counts files
in.files <- list.files(basedir, pattern = "*_features.txt")
ids <- gsub("_features.txt", "", in.files)
rseq <- list()
for(id in ids){
  filename <- paste0(id, "_features.txt")
  all.file <- readLines(paste0(basedir, filename)) 
  dfs <- read.delim(paste0(basedir, filename), comment.char = "#", header = T)
  colnames(dfs)[7] <- "Counts"
  rseq[[id]] <-  dfs %>%
    dplyr::select("Geneid", "Length", "Counts") %>% drop_na() %>% 
    add_column("TPM" = tpm(.$Counts, .$Length)) %>% 
    add_column("log2TPM" = log(.$TPM+1, base = 2)) %>%
    add_column("EnsemblGene" =gsub("\\..*","", .$Geneid)) %>%
    add_column("GeneName" = biomart$hgnc[match(.$EnsemblGene, biomart$ensembl)])
}

#Flatten list to a data frame with all samples and save data file
rseq.df <- ldply(rseq, data.frame) %>% 
  filter(!GeneName %in% rfGenes$Var1)
rseq.df2 <- reshape2::dcast(rseq.df,  GeneName ~ .id, value.var = "Counts", 
                            fun.aggregate = sum) %>% 
  dplyr::filter(rowMeans(.[,-c(1,2)]) > 1) %>% na.omit() %>%
  filter(GeneName != "") %>%
  right_join(kidneyData, ., by = c("external_gene_name" = "GeneName")) %>%
  drop_na
#saveRDS(rseq.df2, "CR_RNA_Jun2022.RDS")

###Run DESEq2 to get DEGs
matrix.count <- as.matrix(rseq.df2)
matrix.count <- matrix.count[,colnames(matrix.count) %in% row.names(meta.rseq)]
matrix.count <- matrix.count[,row.names(meta.rseq)]
storage.mode(matrix.count) = "integer"
dds <- DESeqDataSetFromMatrix(countData = matrix.count,
                              colData = meta.rseq,
                              design = ~ condition)
keep <- rowMeans(counts(dds)) >= 5
dds <- dds[keep,]
dds3 <-DESeq(dds)
resultsNames(dds3)

vst<-vst(dds,blind=TRUE, fitType='local')
plotPCA(vst,
        intgroup = 'condition',
        returnData = FALSE)
##Write table with DE genes
degs <- list()
for(name1 in unique(meta.rseq$condition)){
  for(name2 in unique(meta.rseq$condition)){
    if(name1 != name2){
degs[[paste(name1, name2, sep = "_")]] <- data.frame(results(dds3, contrast = c("condition", name1, name2), 
                               pAdjustMethod = "BH", cooksCutoff=FALSE)) %>%
  filter(padj < 0.05, log2FoldChange > 0) %>% row.names
    }
    else { next;}
  }
}
degsAll <- ldply(degs, data.frame)


##Retrieve normalized data and add gene names
norm <- data.frame(log(counts(dds3, normalized=TRUE)+1, base = 2), 
                   check.names = F)
norm$rowNames <- row.names(norm)
norm <- norm[,-ncol(norm)]


tops <- list()
for(group in unique(meta.rseq$condition)){
  samples <- row.names(meta.rseq)[meta.rseq$condition == group]
  tops[[group]] <- norm %>% select(all_of(samples)) %>% 
    mutate(mean = rowMeans(.)) %>%
    top_n(., 500, mean) %>% row.names
}
topsAll <- ldply(tops, data.frame) %>% distinct(X..i..)


##log2tpm
rseq.df3 <- reshape2::dcast(rseq.df,  GeneName ~ .id, value.var = "log2TPM", 
                            fun.aggregate = sum) %>% 
  dplyr::filter(rowMeans(.[,-c(1,2)]) > 5) %>% na.omit() %>%
  filter(GeneName != "") %>%
  filter(!GeneName %in% rfGenes$Var1)
row.names(rseq.df3) <- rseq.df3$GeneName
rseq.df3 <- rseq.df3[,-1] %>% select(all_of(row.names(meta.rseq)))


##Top5000 variable genes
sel = order(apply(norm, 1, var), decreasing=TRUE)[1:500]
rv <- norm[sel,]


##Plots
#dataset <- rseq.df3 #norm[row.names(norm) %in% topsAll$X..i..,]
genesA <- data.frame(results(dds3, contrast = c("condition", "Fusion", "Luciferase"), 
                   pAdjustMethod = "BH",cooksCutoff=FALSE)) %>%
  dplyr::filter(padj < 0.05) %>% row.names
dataset <- norm %>% filter(row.names(.) %in%genesA) #%>%
  #dplyr::select(!all_of(c("O2978T2", "O2560M", "O45103T2")))
meta.rseq2 <- meta.rseq[,1,drop=F] %>% filter(condition != "MRTK")
meta.rseq2$condition <- factor(meta.rseq2$condition, 
                               levels = c("Fusion", "Fusion_PDX", "TFE3", "Luciferase", "RCC", "WT", "Healthy"))
mycolors <- c("#fcb969","#FE875E", "#5f71b6", "#d1d3d4", "#a34100", "#BEE5B0", "#d4f0f7")
names(mycolors) <- levels(meta.rseq2$condition)
mycolors <- list(condition = mycolors)
pheatmap(as.matrix(dataset), 
         show_rownames = F,
         annotation_col = meta.rseq2, scale = "row",
         cluster_rows = T, cluster_cols = T,
         show_colnames = T, clustering_method = "ward.D2", 
         annotation_colors = mycolors
         #cellwidth = 15, cellheight = 1.1,
         #width = 3, height=9.1
         )


##Analysis of enriched TFs
tfDf$.id <- factor(tfDf$.id, levels = c("Fusion", "RCC", "TFE3"))
tfrcc <- c(
           "USF1",
           "USF2",
           "TFE3",
           "BHLHE41",
           "MYOD1",
           "TFEC",
           "MLX",
           "SREBF2(var.2)"
           
)
tfall <- c("HEY2","ARNT2",
           "HES1", "HES2","TFEB","Arntl",
           "MAX::MYC","MAX","INSM1","HES5","BHLHE40",
           "HEY1","PAX5",
           "MYC"
)
tfDfint <- tfDf[tfDf$Transcription.Factor.Name %in% tfrcc,]
tfDfint$Transcription.Factor.Name <- factor(tfDfint$Transcription.Factor.Name, 
                                            levels = sort(unique(tfDfint$Transcription.Factor.Name), 
                                            decreasing = T))
#topNdf <- d %>% group_by(grp) %>% slice_max(order_by = Significance.Score, n = 10)
ggplot(tfDfint, 
       aes(y= Significance.Score, x = Transcription.Factor.Name, fill = .id)) + 
  geom_bar(stat = "identity", position =  position_dodge2(width = 0.9, preserve = "single")) + 
  coord_flip() + theme_bw() +
  scale_fill_manual(values=c("#fcb969","#5f71b6", "#a34100")) + xlab("Transcription factor") +
  ylab("Significance score") +
  theme(text = element_text(size = 20),legend.title = element_blank())   


##PCA all samples, raw data
pca <- prcomp(t(as.matrix(dataset)),  scale. = T, center = T)
pca2 <- merge(pca$x, meta.rseq, by = "row.names")
project.pca.proportionvariances <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100
ggplot(pca2[pca2$condition != "MRTK"], 
       aes(PC1, PC2, color = condition, size = 10)) + geom_point() +
  xlab(paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%")) +
  ylab(paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%")) +
  ggtitle(paste0("Normalized counts "))


##Generate pca from normalized data (see norm3 above for specific filterings)
pca <- prcomp(t(as.matrix(norm3)), scale. = T, center = T)
pca2 <- merge(pca$x, meta.rseq, by = "row.names")
project.pca.proportionvariances <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100
ggplot(pca2, aes(PC1, PC2, color = type, size = 10)) + geom_point() +
  xlab(paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%")) +
  ylab(paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%")) +
  ggtitle(paste0("Normalized read counts: ", nrow(norm3), " genes")) + 
  labs(fill = "Tumor type") + guides(size = "none") 

##Dot plot GO terms
top50 <- data.frame(pca$rotation) %>% arrange(desc(PC1)) %>% 
  dplyr::slice(1:50) %>% row.names(.)
ego2 <- enrichGO(gene         = gs$V1,
                 OrgDb         = org.Hs.eg.db,
                 universe = row.names(norm),
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.1)
dotplot(ego2, showCategory=15)

norm2 <- 
  melt(as.matrix(norm)) %>% merge(., meta.rseq, by.x ="Var2", 
                 by.y  = "row.names") %>% filter(condition != "MRTK")
norm2$condition <- factor(norm2$condition, levels =c("Healthy","WT", "RCC", "Fusion","Fusion_PDX", "TFE3", "Luciferase"))
mycolors2 <- c("#d4f0f7","#BEE5B0", "#a34100", "#fcb969","#FE875E", "#5f71b6", "#d1d3d4")
path.name <- "CST3"
ggplot(norm2[norm2$Var1 == path.name,], aes(y = value, x = condition, fill = condition)) + 
  geom_boxplot() + ggtitle(paste0(path.name, " normalized expression (log2)")) +
  scale_fill_manual(values=mycolors2) + theme_bw() + ylab("Normalized expression (log2)") +
  xlab(NULL) + guides(fill=guide_legend(title=NULL)) + theme(text = element_text(size = 16)) 

ggplot(norm2[norm2$Var1 %in% test$V1,], aes(condition, value)) + geom_boxplot() + facet_wrap(.~Var1)

rseq.dfmerge <- merge(rseq.df, meta.rseq, by.x = ".id", by.y = "row.names") %>%
  filter(GeneName != "")
rseq.dfmerge$condition <- factor(rseq.dfmerge$condition, levels =c("RCC", "Fusion", "Healthy", "Luciferase", "TFE3"))

ggplot(rseq.dfmerge[rseq.dfmerge$GeneName == path.name,], aes(y = log2TPM, x = condition, fill = condition)) + 
  geom_boxplot() + ggtitle(paste0(path.name, " normalized expression (log2)"))


genes.int <- setdiff(unique(degs$fusion_tfe3, degs$fusion_tfe3),
                     degs$fusion_RCC)

 ##Use a custom set of go terms
library(GSVA)
library(GSVAdata)
library(biomaRt)
##Using all GOS
gos <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")$GO

##Using a pre-dterminede list of GOs
gos <- c("GO:0016055", "GO:0072001", "GO:0030198", "GO:0030111", "GO:0198738")
names(gos) <- c("Wnt signaling pathway", "Renal system development", 
                "ECM organization", "Regulation of Wnt signaling pathway",
                "Cell-cell signaling by Wnt"
                )
customGOs <- list()
for(go in gos[!is.na(gos)]){
  allegs = get(go, org.Hs.egGO2ALLEGS)
  genes = unlist(mget(allegs,org.Hs.egSYMBOL))
  customGOs[[go]] <- genes
  #customGOs[[names(gos[gos == go])]] <- genes
}

gsva_results <- gsva(
  as.matrix(norm),
  customGOs,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Poisson",
  # Minimum gene set size
  min.sz = 5,
  # Maximum gene set size
  max.sz = 900,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)
pathwaysInt <- gsva_results[row.names(gsva_results) %in% c("GO:0001915",
                                                           "GO:0002246"),]
#res2 <- cor(t(pathwaysInt[,grepl("RT", colnames(pathwaysInt))]))
#pheatmap(res2, main = "RT")

pathwaysInt2 <- melt(pathwaysInt) %>% filter(!Var2 %in% c("O2978T2", "O2560M", "O45103T2"))
pathwaysInt2$condition <- meta.rseq$condition[match(pathwaysInt2$Var2, row.names(meta.rseq))]
pathwaysInt2$condition <- factor(pathwaysInt2$condition, 
                                 levels = c("Healthy","WT", "RCC", "Fusion","Fusion_PDX", "TFE3", "Luciferase"))
pathwaysInt2$Var1 <- factor(pathwaysInt2$Var1, levels = names(gos))

ggplot(pathwaysInt2, aes(condition, value)) + geom_boxplot() + 
  facet_wrap(.~Var1, labeller = label_wrap_gen(width=30))
ggboxplot(pathwaysInt2, x="type", y ="value") + 
  stat_compare_means() +
  facet_wrap(.~Var1, labeller = label_wrap_gen(width=30)) +
  ylab("GSVA score") +xlab(NULL)

##GSVA across samples, focus on fusion
gsva_all <- gsva_resultsAllMsig  %>% melt
gsva_all$condition <- meta.rseq$condition[match(gsva_all$Var2, row.names(meta.rseq))]
gsva_all2 <- gsva_all %>% group_by(Var1,condition) %>% dplyr::summarise(mean = mean(value))
ggplot(gsva_all[gsva_all$Var1 %in% c("GO:0043087", "GO:0001655", "GO:0072001", 
                                     "GO:0038127", "GO:0001822", "GO:0072006", 
                                     "GO:0016197"
                                     ),], aes(condition, value)) + geom_boxplot() + 
  facet_wrap(.~Var1)

##GO terms
##GSVA per sample
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "C2" # Only hallmark gene sets
)
hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)
hallmarks_list <- hallmarks_list[grep("REACTOME_*", names(hallmarks_list))]
gsva_resultsAllMsig <- gsva(
  as.matrix(norm),
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 5,
  # Maximum gene set size
  max.sz = 5000,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

gobp <- c("REACTOME_RAC1_GTPASE_CYCLE",
          "REACTOME_EGFR1",
          "REACTOME_SIGNAL_TRANSDUCTION",
          "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
          "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES",
          "REACTOME_DOWNSTREAM_SIGNAL_TRANSDUCTION",
          "REACTOME_SIGNALING_BY_VEGF",
          "REACTOME_SIGNALING_BY_MET"
          #"GOBP_MITOTIC_CELL_CYCLE" ,
          #"GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION",
          #"GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION"
)

pathwaysInt <- gsva_resultsAllMsig[row.names(gsva_resultsAllMsig) %in% gobp,]
#res2 <- cor(t(pathwaysInt[,grepl("WT", colnames(pathwaysInt))]))
#pheatmap(res2, main = "WT")


pathwaysInt2 <- melt(pathwaysInt)
pathwaysInt2$condition <- meta.rseq$condition[match(pathwaysInt2$Var2, row.names(meta.rseq))]
mycols <- c("#fcb969", "#5f71b6", "#d1d3d4")
pathwaysInt2$condition <- factor(pathwaysInt2$condition, 
                                 levels = c("Fusion", "Luciferase", "TFE3"))
ggplot(pathwaysInt2[pathwaysInt2$condition != "Fusion_PDX",], 
       aes(condition, value, fill = condition)) + geom_boxplot() + 
  scale_fill_manual(values=mycols) +
  facet_wrap(.~Var1, labeller = label_wrap_gen(width=5)) +
  ylab("GSVA score") + xlab(NULL) + theme_bw()

##Plots some genes
path.name ="C1GALT1"
ggboxplot(norm[norm$rowNames == path.name,], 
          x = "type.y", y = "value", width = 0.8, add= "boxplot") + 
  ylab('Normalized read counts (log2)') +
  ggtitle("IMDPH1 expression") + stat_compare_means()

genes2 <- c("ENSG00000103024",
            "ENSG00000025708",
            "ENSG00000198931",
            "ENSG00000112667",
            "ENSG00000178035")
genes3 <- c("TYMS", 
            "SHMT1", 
            "SHMT2",
            "MTHFD2",
            "MTHFD1", 'IMPDH1', 'IMPDH2')
ggboxplot(norm2[norm2$hgnc %in% genes3,], 
          x = "type.y", y = "value", width = 0.8, add= "boxplot") + 
  ylab('Normalized read counts (log2)') + stat_compare_means() +
  facet_wrap(.~hgnc)

ggboxplot(norm2[norm2$hgnc == "MTHFD1",], 
          x = "type.y", y = "value", width = 0.8, add= "boxplot") + 
  ylab('Normalized read counts (log2)') +
  facet_wrap(.~hgnc)
