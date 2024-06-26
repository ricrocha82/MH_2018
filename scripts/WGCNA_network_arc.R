# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# https://ramellose.github.io/networktutorials/wgcna.html
# https://deneflab.github.io/HNA_LNA_productivity/WGCNA_analysis.html#
# https://github.com/horto2dj/GLCW/blob/master/WGCNA_glcw.R
# https://figshare.com/articles/dataset/Untitled_ItemSupplementary_Files_for_Henson_et_al_2018_L_O/12055947
# https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/

#############################################
# -----------   Network analysis -----------#
#############################################

# Notes: 
# node = each OTU
# edges = links
# hub = highly connected nodes = high degree centrality 
# degree = distribution of the number of links 
# power law means there is a degree of order in the universe. It is related to
    # free-scale networks (more similar to nature) = many nodes with only a few links, 
        # and a few hubs with large number of links
# modules = subset of nodes or clusters of nodes
# An eigengene is 1st principal component of a module expression matrix and represents a suitably defined average OTU community.

#################
# --- WGCNA ----#
#################

# load packages
library(WGCNA)
library(tidyverse)
library(phyloseq)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#-----------------------------
# 1 - Data input and cleaning 
#------------------------------

# Pull out the information and clean up the data
pseq <- pseq.clr.list[[2]]
# Clust contains the samples we want to keep
paste("There are", ntaxa(pseq), "OTUs in the dataset.")
paste("There are", nsamples(pseq), "samples in the dataset.")

tax.net <- pseq %>% tax_table() %>% data.frame() %>% rownames_to_column("OTU")
kingdom <- pseq %>% tax_table() %>% data.frame() %>% pull(Kingdom) %>% unique()

# Keep in mind that WGCNA needs to have the taxa as columns.
clr_otu <- pseq %>% otu_table() %>% t() %>% as.data.frame()
# select only the numeric variables
meta.pseq <- pseq %>% sample_data() %>% 
  as.tibble() %>% 
  unite("ID",Sample_ID, Layer, sep = "_") %>%
  column_to_rownames("ID") %>%
  select(where(is.numeric)) %>%
  select(where(~!any(is.na(.)))) %>%
  select(Temperature:NH4)

# data cleaning
# Below is an analysis to find out whether there are outliers.
# Check for OTUs and samples with too many missing values
good_otus <- goodSamplesGenes(clr_otu, verbose = 3)
good_otus$allOK # no outliers!!
# if outliers, check https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

# via clustering
sampleTree <- hclust(dist(clr_otu), method = "ward.D2")
sampleTree$labels
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(meta.pseq, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(meta.pseq),
                    main = paste0("Sample dendrogram based on microbial composition (",kingdom,")\n and environmental variables heatmap"))
# the plot shows how environmental parameters relate to the samples dendrogram based on microbial composition

#-----------------------------
# 2 - Automatic network construction and module detection 
#------------------------------

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=50, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(clr_otu, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))

# Set some parameters
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

paste0("The soft threshold is ", sft$powerEstimate,"!")

# Create adjacency matrix by raising OTU matrix by beta and identify subnetworks (modules)
otu_WGCNA2 <- as.matrix(clr_otu)
mode(otu_WGCNA2)
# class(otu_WGCNA2) <- "numeric"

##################33
# The shoftpower is NA - so there is no way to go further since it is impossible to make the 
# 
##################333

# Check that the network ensures scale-free topology at that power
# R should be close to 1 (R > 0.85, I believe), should see a straight line.
##### scaleFreePlot #####
# here we define the adjacency matrix using soft thresholding with beta=14
ADJ1=abs(cor(otu_WGCNA2,use="p"))
# When you have relatively few genes (<5000) use the following code
k = as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
#k=softConnectivity(datE=otu_WGCNA2,power=14)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
scaleFreeFitIndex(k)

#R^2 of 0.87, this suggests we meet the assumption of scale-free topol.

# One-step network construction and module detection
# block-wise network construction and module detection
# power of 4 chosen based on powerEstimate from 'pst'
net = blockwiseModules(clr_otu, power=14, minModuleSize=30,
                       corType = "pearson", saveTOMs = TRUE, 
                       saveTOMFileBase = "blockwiseTOM", pamStage=FALSE, verbose=5)
# Plot the dendrogram
moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    main = paste0("Cluster Dendrogram - ",kingdom))

# Multi-dimensional scaling plots
# Calculate topological overlap. 
# Calculate during module detection, but calculating again here:
dissTOM <- 1-TOMsimilarityFromExpr(clr_otu, power = 14) 

cmd1 = cmdscale(as.dist(dissTOM),2)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(moduleColors), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

#-----------------------------
# 3 - Relating modules to environmental traits
#------------------------------

# quantify module-trait associations
# Identify Eigenvalues for subnetworks by sample
# define numbers of Otus and samples
nPop = ncol(clr_otu)
nSamples = nrow(clr_otu)

# Recalculate MEs with color labels
MEsO = moduleEigengenes(clr_otu, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
names(MEs) <- str_replace_all(names(MEs), "ME","")
save(MEs, moduleLabels, moduleColors, geneTree,file = paste0(folder_path,"data/all/network/Module_eigen_values_",kingdom,".Rdata"))

# boxplots to see if Module eigengene values can separate samples based on Hsbitat and water layer
box.MEs <- MEs %>% rownames_to_column("ID") %>% 
  left_join(meta %>% rownames_to_column("ID") %>% select(ID, Habitat_Layer)) %>%
  select(-ID) %>%
  mutate(Habitat_Layer = fct_relevel(Habitat_Layer, c("ocean_S", "MH_S","MH_P","GR_S"))) %>%
  pivot_longer(!Habitat_Layer, names_to = "Module", values_to = "Eigengene") %>%
  ggplot(aes(x =Habitat_Layer, y = Eigengene, color = Habitat_Layer)) +
    geom_boxplot() +
    scale_color_manual(values = cols_hab_layer) +
    geom_jitter(alpha=0.2, position = position_jitterdodge())+
    geom_point(size=1,alpha=0.2)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~Module, scales = "free_y") +
  labs(title = paste0("Module Eigengene and Sample Locations - ",kingdom), x = "Sites", y = "Eigengene values") +
  theme_Publication_3() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
box.MEs 
# https://www.researchgate.net/publication/40906845_Gene_expression_profiling_in_C57BL6J_and_AJ_mouse_inbred_strains_reveals_gene_networks_specific_for_brain_regions_independent_of_genetic_background/figures?lo=1

# Save data
write.csv(MEs,file=paste0(folder_path,"data/all/network/Module_eigen_values_",kingdom,".csv"))
write.csv(net$colors,file=paste0(folder_path,"data/all/network/Module_composition_",kingdom,".csv"))

# Correlate Eigenvalues to metadata and create heatmap
moduleTraitCor <- cor(MEs, meta.pseq, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# PLOT
sizeGrWindow(10,6)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(meta.pseq),
               yLabels = names(MEs), 
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait Relationships -",kingdom,sep = " "))

# set default plot margin area
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


# OTU significance and module membership
#----------------------------
# Now make a plot for specific module <-> trait (metadata component) pairings
# This allows us to explore the structure of submodule OTU correlations with a given metadata component
# Here we will use "env variable" as trait and "color" as module
# First get the links between all modules and this trait
vecOTUnames <- list()
for (n in names(meta.pseq)) {
  env.par <- n

# Define variable weight containing the weight column of datTrait

env.par.df <- meta.pseq %>% select(env.par)

# Calculate the correlations between modules
geneModuleMembership <- as.data.frame(WGCNA::cor(clr_otu, MEs, use = "p"))

# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# What's the correlation for the trait: bacterial production?
geneTraitSignificance <- as.data.frame(cor(clr_otu, env.par.df, use = "p"))

# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))

# Define OTU significance (GS, originally called gene significance) as the absolute value of the correlation 
# between the gene and the trait.
names(geneTraitSignificance) <- paste("GS.", names(env.par.df), sep = "")
names(GSPvalue) <- paste("p.GS.", names(env.par.df), sep = "")

# Summary output of network analysis 
# Prepare pvalue df
GSpval <- GSPvalue %>% rownames_to_column(var = "OTU") 

gMM_df <- geneModuleMembership %>%
  tibble::rownames_to_column(var = "OTU") %>%
  pivot_longer(!OTU,names_to = "moduleColor", values_to = "moduleMemberCorr") 

# Prepare gene significance df
GS_bacprod_df <- geneTraitSignificance %>%
  data.frame() %>% rownames_to_column(var = "OTU")

# Put everything together 
allData_df <- gMM_df %>%
  left_join(GS_bacprod_df, by = "OTU") %>%
  left_join(GSpval, by = "OTU") %>%
  left_join(tax.net, by = "OTU")

# Write a file 
write.csv(allData_df, file = paste0(folder_path,"data/all/network/corr_",kingdom,"_",env.par,"_WGCNA.csv"))


# Fix the names so that they match the actual color
#names(geneModuleMembership) <- substring(names(geneModuleMembership), 3) # Remove the first two letters in string


par(mfrow = c(2,3))  

# Initialize for loop to plot each module vs the env parameter
# NEED: modNames
vecOTUnames[[n]] <- list()
for (i in names(geneModuleMembership)) {
  
  # Pull out the module we're working on
  module <- i
  print(paste(module,env.par,sep = "-"))   

  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors %>% magrittr::extract(. == module) %>% as.data.frame() %>% 
    rownames_to_column("OTU") %>% rename(module = ".")
  print(paste("There are ", nrow(moduleGenes), " OTUs in the ", module, " module.", sep = ""))
  moduleGenes$OTU %>% print()
  # NOTE: This makes hidden variables with the OTU names: brown_OTUs, yellow_OTUs
  vecOTUnames[[n]][[i]] = moduleGenes$OTU 
  #moduleGenes %>% unite("OTU_module", module, OTU, sep = "_") %>% pull() %>% str_remove('_[[:digit:]]+') 
  # prepare the values to plot
  member <- geneModuleMembership %>% rownames_to_column("OTU") %>% 
    filter(OTU %in% moduleGenes$OTU) %>% pull(module)
  sig <- geneTraitSignificance %>% rownames_to_column("OTU") %>% 
    filter(OTU %in% moduleGenes$OTU) %>% select(-OTU) %>% pull()
  # Make the plot
  verboseScatterplot(abs(member), 
                     abs(sig),
                     xlab = paste0("Module Membership in ", module, "module"),
                     ylab = paste("Gene significance for",env.par,sep = " "),
                     main = paste0(kingdom,"- Module (",module,") membership vs. OTU significance \n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
}    

# Visualizing the network of eigengenes
# Get the similarity between eigenvalues and weight
MET = orderMEs(cbind(MEs, env.par.df))
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with ",env.par), 
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = TRUE,
                      xLabelsAngle = 90)
}

# 5 - OTUs - Module Membership

mmdf <- geneModuleMembership %>%
  tibble::rownames_to_column(var = "OTU") %>%
  mutate(blue = round(abs(blue), digits = 3), 
         grey = round(abs(grey), digits = 3),  
         brown = round(abs(brown), digits = 3), 
         turquoise = round(abs(turquoise), digits = 3),  
         yellow = round(abs(yellow), digits = 3)) %>%
  left_join(tax.net, by = "OTU")
otu.member <-list()
for (i in names(geneModuleMembership)) {
otu.member[[i]] <- mmdf %>% filter(OTU %in% vecOTUnames$Temperature[[i]])
}

otu.member %>% 
  names(.) %>%
  walk(~ write_csv(otu.member[[.]], paste0(folder_path,"data/all/network/otu_module_membership_",kingdom,"_", ., ".csv")))

#############################################################
#------ Moving to machine learning PLS and VIP scores ------
#############################################################

# metadata component (e.g. Nox, designated "weight/env.par" here as above) is the same as before, we just replicate the row names for pls
# the below r2 threshold changes with user input, good correlation above "0.x"
# this is important, because VIP scores don't really mean anything without a good correlation

# set/choose thresholds, module and trait and OTU co-correlation threshold (thr_con)
th_r2 <- 0.5

# 1 - bacteria
#---------------
# nox = blue
env.par <- "Nox"
module <- "blue"
thr_con <- 0.5 # bacteria - Nox/blue

# NH4 = turquoise
env.par <- "NH4"
module <- "turquoise"
thr_con <- 0.2 # bacteria - NH4/turquoise

# 2 - archaea
#--------------

# 3- euk
#--------------

# make a string to be use latter
GS.envpar <- paste0("GS.",env.par)

# get the otus related to determined module (subnetwork)
moduleGenes = moduleColors %>% magrittr::extract(. == module) %>% as.data.frame() %>% 
  rownames_to_column("OTU") %>% rename(module = ".") 
subnetwork <- clr_otu %>% select(any_of(moduleGenes %>% pull(OTU)))

# the env parameter (trait)
weight <- meta.pseq %>% select(env.par)
weight <- as.matrix(weight)
subnetwork <- as.matrix(subnetwork)
class(weight) <- "numeric"
class(subnetwork) <- "numeric"

# run the pls
pls_result <- plsr(weight ~ subnetwork, validation="LOO",method="oscorespls")

# chceck if the R^2 if above the threshold
r2_vector <- R2(pls_result)
max <- 0
max_comp <- -1
for (j in 1:length(r2_vector$val)){
  if(r2_vector$val[j] > th_r2){         # We will only look at the PLS if the correlation is better than th_r2
    if(r2_vector$val[j] > max){
      max <- r2_vector$val[j]
      max_comp <- r2_vector$comp[j]
    }
  }
}
print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep="",".pdf"))

if(max==0){
  print ("No good correlation, we stop here")
} else{
  print("Good correlation, we check the VIP!")
}

# Checking the VIP
vip_result <- VIP (pls_result)
dim(subnetwork)
pls.component <- paste("Comp",max_comp,sep = " ")
vip_components <- vip_result %>% as.data.frame() %>% rownames_to_column("comp") %>% 
  filter(comp == pls.component)  %>% 
  pivot_longer(!comp, names_to = "OTU", values_to = "VIP") %>% arrange(desc(VIP)) %>%
  mutate(to_print = paste0("Rank ",row.names(.)," we have ",OTU," with a VIP of ",VIP))
write.csv(vip_components, paste0(folder_path,"data/all/network/",kingdom,"/Vip_scores_",env.par,"_",module,"_",kingdom,".csv"))

pls.pred <- pls_result[["validation"]][["pred"]] %>% as.data.frame() %>%
  select(all_of(max_comp)) %>% 
  rownames_to_column("ID") 

weight_2 <- weight %>% as.data.frame() %>% rename(measured = env.par) %>% rownames_to_column("ID")
df <- weight_2 %>% left_join(pls.pred) %>% column_to_rownames("ID")
colnames(df)<-c("x","y")

df %>% ggplot() + 
  geom_point(aes(x=x , y=y)) + 
  geom_smooth(aes(x=x,y=y),method=lm) + 
  xlab("Measured") + ylab("Predicted") + 
  ggtitle(paste0("Comparison of ",env.par," measured vs predicted for module ",module)) + 
  theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black")) +
  ggpubr::stat_cor(aes(x=x , y=y), method="pearson", p.accuracy = 0.001)

# Establish the correlation between predicted and modeled
# This is the data to report with the figure (R2, CI, signif, etc.)
# cor.test(df$x,df$y)


# Identify node centrality based on co-occurence data for each OTU in the module

TOM = TOMsimilarityFromExpr(clr_otu, power = 14, corType = "pearson")
# Select submodule of interest based on high correlation and signficance
#module <-"blue" # <- this changes with module color being currently explored
# Select module probes
moduleGenes = moduleColors %>% magrittr::extract(. == module) %>% as.data.frame() %>% 
  rownames_to_column("OTU") %>% rename(module = ".") 
probes = names(clr_otu)
inModule = (moduleColors==module)
modProbes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
write.csv(modTOM, paste0(folder_path,"data/all/network/",kingdom,"/nodeconnections_",env.par,"_",module,"_",kingdom,".csv"))

# Number of cells above X threshhold <-- this number is flexible and should change with your data
# i.e significant Spearman correlations above the threshold as edges
# this is done for visualization purposes since WGCNA subnetworks (based on the Topology Overlap Measure (TOM) between nodes) are hyper-connected.

x <- as.data.frame(rowSums(modTOM > thr_con)) %>% rename(connectivity = 'rowSums(modTOM > thr_con)')
write.csv(x, paste0(folder_path,"data/all/network/",kingdom,"/number_nodes_above_threshold_",env.par,"_",module,"_",kingdom,".csv"))
x %>% arrange(desc(connectivity)) %>% head()

# make scatter hive plots
# You will need to make the Nodeworksheet by combining OTUinfo table for submodule of interest, 
# See workflow for more details. (link above - Henson et al 2018)
# OTUinfo
otuinfo <- read.csv(paste0(folder_path,"data/all/network/",kingdom,"/corr_",kingdom,"_",env.par,"_WGCNA.csv"), header = T, row.names = 1)
otuinfo <- otuinfo %>% filter(moduleColor == module)  %>% filter(OTU %in% moduleGenes$OTU)

# VIP scores
vip.scores <- read.csv(paste0(folder_path,"data/all/network/",kingdom,"/Vip_scores_",env.par,"_",module,"_",kingdom,".csv"))
vip.scores <- vip.scores %>% select(-X, -comp, -to_print) 
n_top <- 10 # to plot the names 
top <- vip.scores %>% slice_max(n = n_top, VIP) %>% left_join(tax.net) %>%
  mutate(names = paste0(OTU,"(",Class,", ",Family,")")) %>% select(OTU, names)
# Nodeconnections. 
node_cent <- read.csv(paste0(folder_path,"data/all/network/",kingdom,"/number_nodes_above_threshold_",env.par,"_",module,"_",kingdom,".csv")) 
node_cent <- node_cent %>% rename(OTU = X)

# join the datasets 
hive.df <- vip.scores %>% left_join(otuinfo) %>% 
  left_join(node_cent) %>%
  left_join(top)
write.csv(hive.df,paste0(folder_path,"data/all/network/",kingdom,"/hive_plot_",env.par,"_",module,"_",kingdom,".csv")) 

# prepare the dataframe to plot
# OTU, connectivity, GS.env.par, VIP, taxa, names
# filter for better visualization
hive <- read.csv(paste0(folder_path,"data/all/network/",kingdom,"/hive_plot_",env.par,"_",module,"_",kingdom,".csv"), header = T, row.names = 1) 
thr_cor_env <- 0.5
thr_VIP <- 0.5
hive <- hive %>% filter(.data[[GS.envpar]] > thr_cor_env & VIP > thr_VIP) 
hive <- hive %>% mutate(OTU = fct_relevel(OTU, unique(OTU)))

# set colors
library(RColorBrewer)
colourCount <- hive %>% pull(Phylum) %>% n_distinct()
# colourCount <-length(unique(ps.glom[,length(ps.glom)]))
getPalette <- colorRampPalette(brewer.pal(colourCount, "Paired"))

# Network visualization and results of PLS analysis on the subnetwork most correlated with env.par.
p <- hive %>% 
  ggplot(aes(x= connectivity, y= .data[[GS.envpar]])) +
  geom_point(aes(size = VIP, colour = Phylum, alpha = 0.5)) +
  scale_size_area(max_size= 10) +
  scale_color_manual(values = getPalette(colourCount)) +
  scale_alpha(guide=FALSE) +
  ggrepel::geom_text_repel(aes(label = names), force = 1) +
  labs(x="Node centrality", y=paste0("Correlation to ",env.par ), color="Phylum",
       size = "VIP") +
  theme_Publication_3() + 
  theme(legend.margin=margin(t = 0, unit='cm'),
        legend.position = "bottom",
        legend.key.size = unit(-0.5,"cm")) 

p



################################################
### Export to Cytoscape for Network Building ###
################################################
# Select modules
module = "blue"
# Select module probes
moduleGenes = moduleColors %>% magrittr::extract(. == module) %>% as.data.frame() %>% 
  rownames_to_column("OTU") %>% rename(module = ".") 
probes = names(clr_otu)
inModule = is.finite(match(moduleColors, module))
modProbes = probes[inModule]
modGenes = tax.net %>% filter(OTU %in% modProbes) %>% unite("name", Phylum, Family, sep = "_")
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste0(folder_path,"data/all/network/",kingdom,"/CytoscapeInput_edges_",module,".txt"),
                               nodeFile = paste0(folder_path,"data/all/network/",kingdom,"/CytoscapeInput_nodes_",module,".txt"),
                               weighted = TRUE,
                               threshold = 0.5,
                               nodeNames = modProbes,
                               altNodeNames = modGenes$name,
                               nodeAttr = moduleColors[inModule])
?exportNetworkToCytoscape
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(modGenes$name, modGenes$name))
