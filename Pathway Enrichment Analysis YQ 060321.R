#Pathview/Gage Analysis for three different comparisons
#updated June 6th, 2021
########################First Comparision C-vs-P##########################
#set directory to the source file
setwd("/Volumes/GoogleDrive/Shared drives/Lab General Folder/Bioinformatics/LLC-PK1/YijunQi LLC-PK1 06032021 Pathvivew_GAGE/C-vs-P")

#check what file I have under the directory
dir()
#this is a test
#Using Significant-DEGs.csv
gene<-read.csv(file="Differential_expression_analysis_table.csv")
head(gene) #checking the gene list
dim(gene)
##loading necessary R packages (run this if it is your first time using it): 

# source("https://bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
# BiocManager::install(version = "3.12")
# BiocManager::install("gage")
# BiocManager::install("pathview")
# BiocManager::install("gageData")

library("pathview")
library("gage")
library("gageData")
# install.packages("dplyr")
library("dplyr")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(DOSE)
library(stringr)
BiocManager::install("org.Ss.eg.db")
BiocManager::install("org.Ss.eg.db")
library(org.Ss.eg.db) #Loading the library gene list for Sus scrofa 

kg.ss <- kegg.gsets(species="ssc")
kegg.sigmet.ss <- kg.ss$kg.sets[kg.ss$sigmet.idx]
kegg.dise.gs <- kg.ss$kg.sets[kg.ss$dise.idx]


kegg.sets.ssc <- kg.ss$kg.sets[kg.ss$sigmet.idx]#USE THIS ONE
head(kegg.sets.ssc,3)

#Using Significant-DEGs.csv (C vs P Comparison)
gene<-read.csv(file="Differential_expression_analysis_table.csv")
head(gene) #checking the gene list
gene$ID # the GAGE function does not support this ID as a fromType, I need to change it to ENSEMBL

#Mapping for the Ensembl database online
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl")
listDatasets(ensembl)# try to find the ss gene dataset
#sscrofa_gene_ensembl # find it
mart <- useDataset("sscrofa_gene_ensembl", useMart("ensembl"))
listAttributes(mart) 
G_list <- getBM(filters= "ensembl_gene_id" , attributes= c("ensembl_gene_id","hgnc_symbol", "description"),values=gene$ID,mart= mart)
geneid=as.data.frame(gene$ID)
gene=as.data.frame(gene)
G_list=as.data.frame(G_list)
CvsP_DEG=merge(gene,G_list,by.x="ID",by.y="ensembl_gene_id")
write.table(CvsP_DEG, file="Differential_expression_analysis_table_annotated.csv",sep=",")

colnames(CvsP_DEG)
gene.df<-bitr(geneID=gene$Gene.name, fromType = "SYMBOL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Ss.eg.db)

head(gene.df)
foldchanges = gene$log2FoldChange
names(foldchanges)= gene.df$ENTREZID

head(foldchanges)
gene.DF=as.data.frame(gene.df)
colnames(gene.DF)
colnames(Gene)
Gene=as.data.frame(gene)

M_gene=merge(gene.DF,Gene,by.x="SYMBOL",by.y= "Gene.name",row.names=NULL,na.rm=TRUE)
colnames(M_gene)
head(M_gene)
#rownames(M_gene)=as.vector(M_gene[,1])
exprs=M_gene[ ,c(2,7:25)]
head(exprs)
E=data.matrix(exprs, rownames.force = T)
#rownames(E)=E[ ,1]
exprs=E
colnames(E)
dim(E)
library(org.Ss.eg.db)
org.Ss.eg.db
#keggres = gage(exprs, ref=c(2:6,12:16),samp=c(7:11,17:20) ,compare = "unpaired",gsets = kegg.sets.ssc)
#keggres = gage(exprs, ref=c(2:6,12:16),gsets = kegg.sets.hs) #used the human one, ssc doesnt work
#fc.kegg.p <- gage(exprs, gsets =kegg.ssc, ref=c(2:6,12:16))
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]


kegg.sets.hs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]#USE THIS ONE, HAVE Hippo pathway in it. Looks like the one I needed.

##Here we look at the expression changes towards one direction (either up- or down- regulation) 
##in the gene sets. The results of such 1-directional analysis is a list of two matrices, corresponding to the two directions. 
##Each result matrix contains multiple columns of test statistics and p-/q-values for all tested gene sets. Among them, p.val is the global 
##p-value and q.val is the corresponding FDR q-value. Gene sets are ranked by significance. Type ?gage for more information.

####Corrected######
keggres = gage(exprs, ref=c(2:6,12:16),samp=c(7:11,17:20),compare = "unpaired",gsets = kegg.sets.hs) 
head(keggres,3)

# col <- colorRampPalette(c("red","yellow","darkgreen"))(30) 

# png("Heatmap_pathwayGreater.png")
# #keggres$greater[ ,c(6:14)]
# heatmap.2(as.matrix(a), Rowv = T, Colv = FALSE, 
#           dendrogram = "row", scale = "row", col = col, density.info = "none", trace = "none", margins = c(7, 15),cluster_rows=FALSE)
# dev.off()

#heatmap(data.matrix(keggres$greater[ c(5:10), ]))
#dev.off()
#kegg.sets.ssc
#"ssc"  
#org.Ss.eg.db
keggres=keggres[!is.na(keggres)]
as.matrix(keggres)

write.table(keggres, file="Pathway Analysis CvsP.csv",sep=",")
summary(keggres)
# #kegg.gsets works with 3000 KEGG species,for examples:
# data(korg)
# head(korg[,1:3])
# 
# 
# head(gene)
# dim(gene)
# control_count=gene[,c(5,6,7,8,9,15,16,17,18,19)]
# sample_count=gene[,-c(1:4,5,6,7,8,9,15,16,17,18,19,24)]
# 
# gage(gse16873.sym, gsets = kegg.gs.sym,
#      ref = hn, samp = dcis)

####################################################################
#?gage
#BMC Bioinformatics 2009, 
##10:161 doi:10.1186/1471-2105-10-161BMC Bioinformatics 2009, 10:161 doi:10.1186/1471-2105-10-161
# gage(exprs, gsets, ref = NULL, samp = NULL, set.size = c(10, 500),
#      same.dir = TRUE, compare = "paired", rank.test = FALSE, use.fold = TRUE,
#      FDR.adj = TRUE, weights = NULL, full.table = FALSE, saaPrep = gagePrep,
#      saaTest = gs.tTest, saaSum = gageSum, use.stouffer=TRUE, ...)
#keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE,na.rm=TRUE)
# gres = gage(foldchanges, gsets =go.cc.gs, same.dir = TRUE,na.rm=TRUE)
# keggres=keggres[!is.na(keggres)]
# keggres

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

#pick the top and bottom five pathways to draw. you can change the number depends on your interests. 
keggrespathwaysGreater = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tibble::as_tibble() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()

keggrespathwaysLess = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tibble::as_tibble() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()

keggrespathwaysGreater
keggrespathwaysLess

# Get the IDs.
keggresidsGreater = substr(keggrespathwaysGreater, start=1, stop=8)
summary(keggresidsGreater)

keggresidsLess = substr(keggrespathwaysLess, start=1, stop=8)
summary(keggresidsLess)

keggresids=c(keggresidsGreater,keggresidsLess) #prepare to output all of the figures at the sametime

data(gene.idtype.list); gene.idtype.list
#set up matrix

M_gene_path=M_gene[ ,c(2,4)]
M_gene_path=data.matrix(M_gene_path)
rownames(M_gene_path)=M_gene_path[ ,1]
# define the function
plot_pathway = function(pid) 
  pathview(
    gene.data=M_gene_path[ ,2], pathway.id=keggresids, species="hsa", new.signature=FALSE)

# output mutiple pathways and output to the working directory automatically
tmp = sapply(keggresids, plot_pathway ) # takes a bit to finish running, be patient;). 
# note: the human pathway package is a bit outdated, so 
#if it doesnt have the pathways you wanted, you can go to https://pathview.uncc.edu/analysis, somehow this website's library is a bit dated.

#making a heat map for C vs P with gene counts and symbol

#
#H=M_gene[ ,c(1,7:ncol(M_gene))]
head(M_gene)
H=M_gene[ ,c(7:ncol(M_gene))]
head(H)
H=as.matrix(H)
colnames(H)=c("CTL-Set1-33", "CTL-Set1-52","CTL-Set1-86","CTL-Set1-92","CTL-Set1-112",
                       "PKD2-Set1-47","PKD2-Set1-54","PKD2-Set1-64","PKD2-Set1-96","PKD2-Set1-103",
                       "CTL-Set2-33","CTL-Set2-52", "CTL-Set2-86","CTL-Set2-92","CTL-Set2-112",
                       "PKD2-Set2-47","PKD2-Set2-54","PKD2-Set2-96","PKD2-Set2-103")
#filiter genes for heatmap plot with the most significant genes
M_gene_f=which(M_gene$padj<0.01) #select rows'number which has a padj value smaller than 0.01
head(M_gene)
M_gene_f=M_gene[c(M_gene_f), ]
head(M_gene_f)
O=order(M_gene_f$padj) #given the row number, ranked by padj from smallest to largest (0.01)
H=M_gene_f[c(O), ]
head(H)
H_count=H[ ,-c(2:6)]
dim(H_count)
head(H_count)
H_count=as.data.frame(H_count)
RN=H_count[ ,1] #save row names for H_count
H_count=data.matrix(H_count)
H_count=H_count[ ,c(2:ncol(H_count))]
head(H_count)
rownames(H_count)=RN #adding gene names to the selected count data
head(H_count)
H_count=H_count[ ,c(1:5,11:15,6:10,16:ncol(H_count))]
head(H_count)
summary(H_count)
colnames(H_count)=c("CTL-Set1-33", "CTL-Set1-52","CTL-Set1-86","CTL-Set1-92","CTL-Set1-112",
                    "CTL-Set2-33","CTL-Set2-52",
                    "CTL-Set2-86","CTL-Set2-92","CTL-Set2-112",
                    "PKD2-Set1-47","PKD2-Set1-54","PKD2-Set1-64",
                    "PKD2-Set1-96","PKD2-Set1-103",
                    "PKD2-Set2-47","PKD2-Set2-54","PKD2-Set2-96","PKD2-Set2-103")
#filiter genes for heatmap
#H=as.matrix(H)
head(H_count)
dim(H_count)
H_count=as.matrix(H_count)
H_count=data.matrix(H_count)
colnames(H_count)
png("Heatmap_counts(adjP-value ranked).png",width = 4, height = 15, units = 'in', res = 600)
#heatmap(H_count,cluster_rows=F, cluster_cols=F,dendrogram='none',Rowv=FALSE, Colv=FALSE,trace="none") ##
heatmap(H_count) ##you could take out "PKD2-Set1-64" (-13)
colnames(H_count)
dev.off()

#try using heatmap.2 function in gplots package
# # load package
# library(gplots)
# 
# my.image=structure(H_count[1:50, ],.Dim=c(50,19), .dimnames=list(c(RN[1:50]),colnames(H_count))) 
# heatmap.2(my.image, density.info="none", trace="none", dendrogram='none', 
#           Rowv=FALSE, Colv=FALSE)

#filiter genes for heatmap plot with the most significant genes, then low expressed
M_gene_f=which(M_gene$padj<0.01) #select rows'number which has a padj value smaller than 0.01
head(M_gene)
M_gene_f=M_gene[c(M_gene_f), ]
head(M_gene_f)
O=order(M_gene_f$log2FoldChange) 
H=M_gene_f[c(O), ]
head(H) #the most negative
tail(H) #the most positive
H_count=H[ ,-c(2:6)]
dim(H_count)
head(H_count)
H_count=as.data.frame(H_count)
RN=H_count[ ,1] #save row names for H_count
H_count=data.matrix(H_count)
H_count=H_count[ ,c(2:ncol(H_count))]
head(H_count)
rownames(H_count)=RN #adding gene names to the selected count data
head(H_count)
H_count=H_count[ ,c(1:5,11:15,6:10,16:ncol(H_count))]
head(H_count)
summary(H_count)
colnames(H_count)=c("CTL-Set1-33", "CTL-Set1-52","CTL-Set1-86","CTL-Set1-92","CTL-Set1-112",
                    "CTL-Set2-33","CTL-Set2-52",
                    "CTL-Set2-86","CTL-Set2-92","CTL-Set2-112",
                    "PKD2-Set1-47","PKD2-Set1-54","PKD2-Set1-64",
                    "PKD2-Set1-96","PKD2-Set1-103",
                    "PKD2-Set2-47","PKD2-Set2-54","PKD2-Set2-96","PKD2-Set2-103")
#filiter genes for heatmap
#H=as.matrix(H)
head(H_count)
dim(H_count)
H_count=as.matrix(H_count)
H_count=data.matrix(H_count)
colnames(H_count)
png("Heatmap_counts(ranked via fold2change).png",width = 4, height = 15, units = 'in', res = 600)
#heatmap(H_count,cluster_rows=F, cluster_cols=F,dendrogram='none',Rowv=FALSE, Colv=FALSE,trace="none") ##
heatmap(H_count) ##you could take out "PKD2-Set1-64" (-13)
colnames(H_count)
dev.off()
col <- colorRampPalette(c("red","yellow","darkgreen"))(30) 

png("Heatmap_most positive top 10.png")
heatmap.2(H_count[c(c(nrow(H_count)-10):nrow(H_count)), ], Rowv = T, Colv = FALSE, dendrogram = "row", scale = "row", col = col, density.info = "none", trace = "none", margins = c(7, 15) )
          dev.off()
png("Heatmap_most negative top 10.png")
heatmap.2(H_count[1:10, ], Rowv = T, Colv = FALSE, dendrogram = "row", scale = "row", col = col, density.info = "none", trace = "none", margins = c(7, 15) )
dev.off()
          
          