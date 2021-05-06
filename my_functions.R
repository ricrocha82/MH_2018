# tips
# {{Rank}} for function - to wrap the function arguments 
# .data[[Rank]] for string inputs
# https://www.tidyverse.org/blog/2019/06/rlang-0-4-0/
# https://www.tidyverse.org/blog/2020/02/glue-strings-and-tidy-eval/


# some functions were based on microbiomeSeq package 
# more at: <http://userweb.eng.gla.ac.uk/umer.ijaz/projects/microbiomeSeq_Tutorial.html>
# Alfred Ssekagiri, William T. Sloan, * Umer Zeeshan Ijaz (* Correspondence: Umer.Ijaz@glasgow.ac.uk)


#'ANOVA of environmental variables
#'
#'This functions applies analysis of variance on the measured environmental variables
#'in specified groups and constructs a plot showing how the environmental variables vary in the groups.
#'The plot is annotated with significance level as obtained from \code{anova}.
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Name of a categorical variable that is preffered for grouping the.
#'        information, this should be one of the components of grouping vector.
#' @param pValueCutoff. p-value threshold for significance of anova results, default set to 0.05.
#' @param select.variables. A vector of character strings(s) specifying environmental variable(s) to be analysed. If
#' not supplied, all numeric variables are analysed.
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#' @examples
#' data(pitlatrine)
#' physeq<-pitlatrine
#' p1<-plot_anova_env(physeq,grouping_column =  "Country",select.variables=c("Temp","pH"))
#' print(p1)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export plot_anova_env
#'
plot_anova_env <- function(physeq, grouping_column, pValueCutoff=0.05,select.variables=NULL, print.lines = TRUE){
  
  #get meta data from phyloseq object
  meta_table <- as.data.frame(sample_data(physeq))
  #pick numerical variables of environmental data
  env_table <- meta_table[,sapply(meta_table,is.numeric)]
  df<- reshape2::melt(as.matrix(env_table))
  names(df)<-c("sample","measure","value")
  #Incorporate categorical data in df
  df<-data.frame(df,(meta_table[, grouping_column])[as.character(df$sample),])
  
  #do anova of environmental variables between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  df <- anova_res$df #get updated environmental measure information
  
  #pick selected variables
  if(!is.null(select.variables)){
    df <- df[which(df$measure%in%select.variables),]
    df_pw<-df_pw[which(df_pw$measure%in%select.variables),]
  }
  #Draw the boxplots
  p<-ggplot2::ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+ggplot2::geom_boxplot()+geom_jitter(position = position_jitterdodge(), alpha = 0.3)
  p<-p+ggplot2::theme_bw()+geom_point(size=1,alpha=0.2)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+ggplot2::facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Groups")
  p<-p+ggplot2::theme(strip.background = element_rect(fill = "white"))
  #This loop will generate the lines and signficances
  if(print.lines == FALSE){
    return(p)
  } else {
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
     p<-p+ ggplot2::geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
     p<-p+ ggplot2::geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
  }
  }
  return(p)
}
}

##################################3
##################################

#'pairwise ANOVA
#'
#'This function performs ANOVA of a given measure in specified groups, in addition it
#'also performs pairwise ANOVA of the measure between possible pairs of levels in the grouping variable. It returns
#'p-values obtained.
#'
#' @param df (Required). A \code{data.frame} containg measures a measure to be analysed.
#' @param meta_table (Required). A data frame containing atleast one variable (grouping variable).
#' @param grouping_column (Required). A character string specifying name of the grouping variable in the supplied meta table.
#' @return It returns a list of two data frames: One being p-values obtained from the pairwise ANOVA of
#' the measure and levels of grouping variable , the other is containing updated measure.
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @import data.table
#'
#' @export perform_anova
#'

perform_anova <- function(df,meta_table,grouping_column,pValueCutoff){
  
  dt<-data.table::data.table(data.frame(df,.group.=meta_table[,grouping_column]))
  #specifying a p-value cutoff for the ggplot2 strips
  pval<-dt[, list(pvalue = sprintf("%.2g",
                                   tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
           by=list(measure)]
  #Filter out pvals that we are not significant
  pval<-pval[!pval$pvalue=="",]
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]
  
  #using sapply to generate significances for pval$pvalue using the cut function.
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, pValueCutoff, Inf),label=c("***", "**", "*", "")))})
  
  #Update df$measure to change the measure names if the grouping_column has more than three classes
  if(length(unique(as.character(meta_table[,grouping_column])))>2){
    df$measure<-as.character(df$measure)
    if(dim(pval)[1]>0){
      for(i in seq(1:dim(pval)[1])){
        df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
      }
    }
    df$measure<-as.factor(df$measure)
  }
  #Get all possible pairwise combination of values in the grouping_column
  s<-combn(unique(as.character(df[,grouping_column])),2)
  
  #df_pw will store the pair-wise p-values
  df_pw<-NULL
  for(k in unique(as.character(df$measure))){
    #We need to calculate the coordinate to draw pair-wise significance lines
    #for this we calculate bas as the maximum value
    bas<-max(df[(df$measure==k),"value"])
    #Calculate increments as 10% of the maximum values
    inc<-0.1*bas
    #Give an initial increment
    bas<-bas+inc
    for(l in 1:dim(s)[2]){
      #Do a pair-wise anova
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      #Ignore if anova fails
      if(!is.na(as.numeric(tmp[length(tmp)]))){
        #Only retain those pairs where the p-values are significant
        if(as.numeric(tmp[length(tmp)])<pValueCutoff){
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
          #Generate the next position
          bas<-bas+inc
        }
      }
    }
  }
  if(!is.null(df_pw)){
    df_pw<-data.frame(row.names=NULL,df_pw)
    names(df_pw)<-c("measure","from","to","y","p")
  }
  out <- list("df_pw"=df_pw, "df"=df)
  return(out)
}

###############################################3
#################################################
################################################3


#'Alpha diversity measure
#'
#'This function calculates alpha diversity of provided community data using
#'selected indices/method.
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "richness", "fisher", "simpson", "shannon" and "evenness".
#' @return It returns a data frame of diversity measure and corresponding indices/methods
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export alpha_div
#'
alpha_div <- function(physeq,method){
  #==check for validity of selected methods
  method<- match.arg(method,c("richness", "fisher", "simpson", "shannon", "evenness","pd", "invsimpson"), several.ok = TRUE)
  
  abund_table <- otu_table(physeq)
  df <- NULL
  if("richness"%in%method){
    R<- vegan::rarefy(abund_table,min(rowSums(abund_table)))
    df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))
    if(is.null(df)){
      df<-df_R}
    else {
      df<-rbind(df,df_R)}
  }
  if("fisher"%in%method){
    alpha <- vegan::fisher.alpha(abund_table)
    df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))
    if(is.null(df)){
      df<-df_alpha}
    else {
      df<-rbind(df,df_alpha)}
  }
  if("simpson"%in%method){
    simp <- vegan::diversity(abund_table, "simpson")
    df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))
    if(is.null(df)){
      df<-df_simp}
    else {
      df<-rbind(df,df_simp)}
  }
  if("invsimpson"%in%method){
    invsimp<- vegan::diversity(abund_table, "invsimpson")
    df_invsimp<-data.frame(sample=names(invsimp),value=invsimp,measure=rep("InvSimpson",length(invsimp)))
    if(is.null(df)){
      df<-df_H}
    else {
      df<-rbind(df,df_invsimp)}
  }
  if("shannon"%in%method){
    H<- vegan::diversity(abund_table)
    df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
    if(is.null(df)){
      df<-df_H}
    else {
      df<-rbind(df,df_H)}
  }
  if("evenness"%in%method){
    H<-vegan::diversity(abund_table)
    S <- specnumber(abund_table)
    J <- H/log(S)
    df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))
    if(is.null(df)){
      df<-df_J}
    else {
      df<-rbind(df,df_J)}
  }
  if("pd"%in%method){
    otu_tree <- phyloseq::phy_tree(physeq)
    PD <- pd(abund_table, otu_tree ,include.root = TRUE)
    df_PD<-data.frame(sample=names(PD),value=PD,measure=rep("PD",length(PD)))
    if(is.null(df)){
      df<-df_PD}
    else {
      df<-rbind(df,df_PD)}
  }
  return(df)
}


###############################################3
##############################################3
##################################################

#' Co-occurence network
#'
#'This function uses co-occurence pattern analysis to identify co-occuring features/taxa in community
#'data under specified environmental conditions. Co-occurence is measured as positive correlation whose threshold(s)
#'can be specified as indicated in arguments section. Amongst these features, pairwise co-occurences which are outstanding within sub communities
#'are detected. p-values generated during pairwise correlations are adjusted for multiple comparisons
#'by false discovery rate. The network statistics used to assign importance of taxa/features include betweenness, closeness,
#'and eigenvector centrality.
#'
#' @param physeq(Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method. A character string for the correlation method to be used; options include: "cor" and "bicor"
#' @param rhos (required). A vector specifying thresholds for correlation between co-occuring pairs of features.
#' @param select.condition (optional). A character string speifying name of a desired condition/group. If not supplied, co-occuence analysis
#'       is performed amongst all conditions present in the grouping variable.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information.
#' @param qval_threshold Cut off for "fdr" adjusted p-values.
#'
#' @return Files of visual representation of the network showing subcommunities (identified by colors),
#' network statistics and file containing pairwise corrrelations of taxa in under different conditions.
#' @examples
#' data(pitlatrine)
#' physeq <- taxa_level(pitlatrine, "Phylum")
#' co_occr <- co_occurence_network(physeq, grouping_column = "Country", rhos = 0.35, select.condition = "V", scale.vertex.size=3, scale.edge.width=15)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references [Ryan J. Williams, Adina Howe and Kirsten S. Hofmockel.
#'            Demonstrating microbial co-occurrence pattern analyses within and between ecosystems,
#'            Frontiers in Microbial Ecology, 5:358, 2014].
#'
#' @author  Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export co_occurence_network
#'

co_occurence_network <- function(physeq ,qval_threshold=0.05, grouping_column,rhos=c(-0.75,-0.5,0.5,0.75),select.condition=NULL,method="cor", filename=NULL, ...){
  
  corr_results <- network_correlation(physeq, grouping_column, select.condition, method, filename)
  
  comm_results<-comm_stats(corr_results$results, rhos)
  
  edge_lists<-co_occur_pairs(corr_results$results, rhos)
  
  edge_lists<-edge_lists[edge_lists$qval<=qval_threshold,]
  
  net <- construct_network(edge_lists, comm_results, ...)
  
  out <- list(net=net, comm_data=corr_results$comm.data, meta_table=corr_results$meta_table)
  
  return(out)
}

#computing correlation of taxa/features
network_correlation <- function(physeq, grouping_column, select.condition, method, filename){
  abund_table <- as.data.frame(otu_table(physeq))
  meta_table <- data.frame(sample_data(physeq))
  meta_table$Groups <- meta_table[,grouping_column]
  
  if(!is.null(select.condition)){meta_table <- subset(meta_table, Groups==select.condition)}
  
  # === rarefy
  comm.data<-vegan::rrarefy(abund_table,min(rowSums(abund_table)))
  trts<-unique(meta_table$Groups)
  
  results<-matrix(nrow=0,ncol=7)
  for(a in 1:length(trts)){
    trt.temp<-trts[a]
    #subset the dataset for those treatments
    temp<-subset(comm.data, meta_table$Groups==trt.temp)
    #remove empty species
    temp<-temp[,colSums(temp)>0]
    
    new.row<-NULL
    if(method=="bicor"){
      #Biweight midcorrelation
      res<- WGCNA::bicorAndPvalue(temp,use="p")
      
      res$bicor[upper.tri(res$bicor)] <- NA
      diag(res$bicor)<-NA
      res$p[upper.tri(res$p)]<-NA
      diag(res$p)<-NA
      
      new.row<-data.frame(rep(trts[a],dim(temp)[2]),reshape2::melt(as.matrix(res$bicor)),reshape2::melt(as.matrix(res$p))[,3])
      
    } else if (method=="cor"){
      #Pearson correlation
      res<- WGCNA::corAndPvalue(temp,use="p")
      
      res$bicor[upper.tri(res$cor)] <- NA
      diag(res$cor)<-NA
      res$p[upper.tri(res$p)]<-NA
      diag(res$p)<-NA
      
      new.row<-data.frame(rep(trts[a],dim(temp)[2]), reshape2::melt(as.matrix(res$cor)),reshape2::melt(as.matrix(res$p))[,3])
    }
    colnames(new.row)<-c("trt","taxa1","taxa2","rho","p.value")
    new.row<-data.frame(new.row,ab1=as.vector(sapply(as.character(new.row$taxa1),function(x){colSums(temp)[x]})),ab2=as.vector(sapply(as.character(new.row$taxa2),function(x){colSums(temp)[x]})))
    results<-rbind(results,new.row)
  }
  
  results$taxa1<-as.character(results$taxa1)
  results$taxa2<-as.character(results$taxa2)
  #We remove the upper triangle and diagonal from the correlation just calculated using NA assignments done before
  results<-results[complete.cases(results),]
  results[(results$ab1 <= 1) | (results$ab2 <= 1),"rho"]<-0
  results[(results$ab1 <= 1) | (results$ab2 <= 1),"p.value"]<-1
  
  if(!is.null(filename)){write.csv(results,paste(filename,"_co-occurence","_",paste(unique(meta_table$Groups),collapse="_vs_"),"_correlations.csv",sep=""))}
  
  out <- list(comm.data=comm.data, results=results, meta_table=meta_table)
  return(out)
}


#######################################3
#########################################
###########################################


#'Taxa-Environment correlation
#'
#'This function finds the relationship between most abundant taxa
#'and numerical environmental variables based on correlation. The abundance of each feature/taxa is
#'correlated with each of the environmental variables. A correlation test is performed and associated
#'p-values are adjusted for multiple testing. The scheme of adjustment is elaborated in the arguments section. The function returns a data.frame with raw p-values, corrected p-values,
#'and correlation results.
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information, this should be one of the components of grouping vector.
#' @param method A character string indicating which correlation coefficient is to be computed, available options are "pearson" which is also the default, "kendall" and "spearman".
#' @param pvalue.threshold. Cut off p-value for significance of correlation between taxa abundance and environmental variables, default is 0.05.
#' @param padjust.method. Method for adjusting p-values.
#' @param adjustment (optional). An integer to specify how the p-values should be adjusted;
#'        \itemize{
#'        \item 1 - donot adjust
#'        \item 2 - adjust environmental variables + Groups (column on the correlation plot)
#'        \item 3 - adjust Taxa + Groups (row on the correlation plot for each Groups)
#'        \item 4 - adjust Taxa (row on the correlation plot)
#'        \item 5 - adjust environmental variables (panel on the correlation plot)
#'        }
#' @param select.variables (optional). Character string for the names(s) of environmental variables to be considered in correlation compution.
#'        if not given, all numerical variables are considered.
#' @return A \code{data.frame} of correlation results.
#'
#' @examples
#' data(pitlatrine)
#' physeq <- pitlatrine
#' physeq <- taxa_level(physeq, "Phylum")
#' tax_env_cor <- taxa.env.correlation(physeq,"Country")
#' plot_taxa_env(tax_env_cor)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export taxa.env.correlation
#' @export plot_taxa_env

taxa.env.correlation <- function(physeq, grouping_column, method="pearson", pvalue.threshold=0.05,
                                 padjust.method="BH", adjustment=1, num.taxa=50, select.variables=NULL){
  
  method<- match.arg(method,c("pearson", "kendall", "spearman"),several.ok = F)
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- as.data.frame(otu_table(physeq))
  #abund_table <- t(otu_table(physeq))
  meta_table <- data.frame(sample_data(physeq))
  #get grouping information
  groups<-meta_table[,grouping_column]
  #select variables to show (exclude) in (from) correlation plot.
  if(!is.null(select.variables)){
    meta_table <- subset(meta_table,select=select.variables)
  }
  #pick the numerical environmental variables since correlation function only accepts numerical variable
  mt_env <- meta_table[,sapply(meta_table,is.numeric)]
  
  #Now get a filtered abundance table based on selected variables
  abund_table_filt<-abund_table[rownames(mt_env),]
  
  # #pick top most  num.taxa taxa
  # select.top.taxa <- top.taxa(abund_table, num.taxa)
  # abund_table_filt <- select.top.taxa$abund_table
  abund_table_filt<-abund_table_filt[,order(colSums(abund_table_filt),decreasing=TRUE)]
  #Extract list of top N Taxa
  taxa_list<-colnames(abund_table_filt)[1:num.taxa]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
  abund_table_filt<-data.frame(abund_table_filt[,colnames(abund_table_filt) %in% taxa_list])
  
  #Now calculate the correlation between individual Taxa and the environmental data
  df <- tables.correlate(abund_table_filt, mt_env, groups, method)
  colnames(df)<-c("Taxa","Env","Correlation","Pvalue","Type")
  df$Pvalue<-as.numeric(as.character(df$Pvalue))
  df$Correlation<-as.numeric(as.character(df$Correlation))
  # add column for adjusted p-values
  df$AdjPvalue<-rep(0,dim(df)[1])
  
  # correct pvalues for multiple testing
  df <- p.adjust.cor(df, adjustment, padjust.method)
  #Now we generate the labels for signifant values
  df$Significance<-cut(df$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  #We ignore NAs
  df<-df[complete.cases(df),]
  return(df)
}

tables.correlate<-function(table1, table2, groups=NULL, method){
  df<-NULL
  for(i in colnames(table1)){
    for(j in colnames(table2)){
      
      if(!is.null(groups)){
        for(k in unique(groups)){
          a<-table1[groups==k,i,drop=F]
          b<-table2[groups==k,j,drop=F]
          tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
          
          if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}
        }
      }
      else{
        
        a<-table1[,i,drop=F]
        b<-table2[,j,drop=F]
        tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value)
        
        if(is.null(df)){df<-tmp} else{df<-rbind(df,tmp)}
        
      }
      
    }
  }
  
  df<-data.frame(row.names=NULL,df)
  return(df)
}


plot_taxa_env <- function(df){
  p <-ggplot2::ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=df)
  p <- p + ggplot2::geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C")
  p<-p+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  p<-p+ ggplot2::geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)
  p<-p+ ggplot2::facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")
  p<-p+ ggplot2::xlab("Groups")
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))
  return(p)
}

#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
#adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")

# df is a data frame
p.adjust.cor <- function(df,adjustment=1,padjust.method="BH"){
  if(adjustment==1){
    df$AdjPvalue<-df$Pvalue
  } else if (adjustment==2){
    for(i in unique(df$Env)){
      for(j in unique(df$Type)){
        sel<-df$Env==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      }
    }
  } else if (adjustment==3){
    for(i in unique(df$Taxa)){
      for(j in unique(df$Type)){
        sel<-df$Taxa==i & df$Type==j
        df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
      }
    }
  } else if (adjustment==4){
    for(i in unique(df$Taxa)){
      sel<-df$Taxa==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
    }
  } else if (adjustment==5){
    for(i in unique(df$Env)){
      sel<-df$Env==i
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method=padjust.method)
    }
  }
  return(df)
}



#############################################333
################################################33
##################################################


#'Taxa level collating
#'
#'This function takes a physeq object and returns a physeq object with taxa at a specified taxonomic level.
#'It extends the \link[phyloseq]{tax_glom}  to include names of corresponding taxa
#'at a specified taxonomy level and  updating tree tip labels accordingly.
#'
#' @param physeq (Required). A \link[phyloseq]{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param which_level (Required). Character string specifying taxonomic level.
#' @return physeq object at specified taxonomic level.
#'
#' @examples
#'
#' data(pitlatrine)
#' physeq<-data(physeq)
#' physeq <- taxa_level(physeq = physeq,which_level = "Family")
#'
#'@export taxa_level
#'
taxa_level <- function(physeq,which_level){
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}


#'Alpha diversity with ANOVA
#'
#'This function calculates alpha diversity of provided community data using
#'selected indices/method(s). It performs pair-wise ANOVA of diversity measures between groups
#'and outputs a plot for each of the selected methods(indices) annotated with significance levels.
#' @importFrom phyloseq t
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Required). A list of character strings specifying \code{method} to be used to calculate for alpha diversity
#'        in the data. Available methods are: "richness", "fisher", "simpson", "shannon" and "evenness".
#' @param grouping_column (Required). A character string specifying the name of a categorical variable containing  grouping information.
#' @return Returns a ggplot object which can further be manipulated further.
#'
#' @examples
#' data(pitlatrine)
#' physeq <- pitlatrine
#' p<-plot_anova_diversity(physeq, method = c("richness","simpson"),grouping_column =  "Country",pValueCutoff=0.05)
#' plot_anova_diversity(physeq, method = c("richness","shannon"), grouping_column = "Depth")
#' print(p)
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com},Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export plot_anova_diversity
#'

plot_anova_diversity <- function(physeq, method, grouping_column,pValueCutoff=0.05, outfile="anova_diversity.csv", print.lines = TRUE)
{
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  
  #get diversity measure using selected methods
  div.df <- alpha_div(physeq,method)
  
  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  
  #perform anova of diversity measure between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  write.csv(df_pw, file = outfile)
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("Samples")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("Groups")
  
  #This loop will generate the lines and signficances
  if(print.lines == FALSE){
    return(p)
  } else {
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
    }
  }
  return(p)
}
}

#'Data normalisation
#'
#'This function normalises taxa abundance
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method (Required). A \code{method} to used to normalise the data. Available methods include:
#'        "edgeRnorm", "varstab", "randomsubsample", "proportion" .
#' @return Returns a phyloseq object whose abundance data is normalised.
#' @examples
#' data(pitlatrine)
#' physeq<-data(pitlatrine)
#' physeq <- normalise_data(physeq,method = "randomsubsample")
#' physeq <- normalise_data(physeq,method = "edgeRnorm")
#' physeq <- normalise_data(physeq,method = "proportion")
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export normalise_data
#'
normalise_data <- function(physeq, norm.method, ...){
  norm.method = match.arg(norm.method,c("edgernorm","varstab","randomsubsample",
                                        "proportion","relative","log-relative","scale"))
  switch(norm.method,
         "randomsubsample"=randomsubsample(physeq),
         "proportion"=proportion(physeq),
         "varstab"=deseq_varstab(physeq, ...),
         "edgernorm"=edgeRnorm(physeq, ...),
         "log-relative"=log_relative(physeq, ...),
         "relative"=relative(physeq, ...),
         "scale"=scale.meta(physeq, ...)
  )
}

#===== normalisation functions =========================#
#Each of the normalisation functions takes up a phloseq object and returns
#returns a physeq object whose otu-table is transformed.
#1. edgernorm
edgeRnorm = function(physeq, ...){
  
  abund_table <- otu_table(physeq)
  
  # Enforce orientation.
  if(!taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }
  x = as(abund_table, "matrix")
  # See if adding a single observation, 1, 
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify 
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems. 
  # Can the 1 be reduced to something smaller and still work?
  x = x + 1
  # Now turn into a DGEList
  y = edgeR::DGEList(counts=x, remove.zeros=TRUE)
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, ...)
  # A check that we didn't divide by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors")
  }
  otu_table(physeq) <- otu_table(z$counts, taxa_are_rows=TRUE)
  return(physeq)
}
#2. variance stabilisation
deseq_varstab = function(physeq, sampleConditions=rep("A", nsamples(physeq)), ...){
  
  abund_table <- otu_table(physeq)
  
  # Enforce orientation.
  if(!taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }
  
  x = as(abund_table, "matrix")
  # The same tweak as for edgeR to avoid NaN problems
  # that cause the workflow to stall/crash.
  x = x + 1
  
  cds = DESeqDataSetFromMatrix(x, DataFrame(sampleConditions), ~ 1)
  # First estimate library size factors
  cds = estimateSizeFactors(cds)
  # Variance estimation, passing along additional options
  cds = estimateDispersions(cds, ...)
  #vsmat = varianceStabilizingTransformation(cds)
  vsmat <- varianceStabilizingTransformation(cds)
  otu_table(physeq) <- otu_table(assay(vsmat), taxa_are_rows=TRUE)
  return(physeq)
}
#3.proportion
# Normalize total sequences represented 

# Scale by dividing each variable by its standard deviation.
#physeq = transform_sample_counts(physeq, function(x) x/sd(x))
# Center by subtracting the median
#physeq = transform_sample_counts(physeq, function(x) (x-median(x)))
proportion = function(physeq){
  normf = function(x, tot=max(sample_sums(physeq))){ tot*x/sum(x) }
  physeq = transform_sample_counts(physeq, normf)
  return(physeq)
}

#4.random sampling 
randomsubsample = function(physeq, smalltrim=0.15, replace=TRUE,meta=F){
  # Set the minimum value as the smallest library quantile, n`smalltrim` 
  samplemin = sort(sample_sums(physeq))[-(1:floor(smalltrim*nsamples(physeq)))][1]
  physeqr = rarefy_even_depth(physeq, samplemin, rngseed=TRUE,replace=replace, trimOTUs=TRUE)
  return(physeqr)
}

#5.relative transformation
relative <- function(physeq,norm.meta=F,select.variables=NULL){
  if(norm.meta){
    get.vars <- get.num.variables(physeq)
    
    norm.variables <- get.vars$num.variables/rowSums(get.vars$num.variables)
    
    meta_table <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
    
    sample_data(physeq) <- meta_table
  }
  else if(!norm.meta){
    abund_table <- otu_table(physeq)
    otu_table(physeq) <- abund_table/rowSums(abund_table)
  }
  return(physeq)
}
#6. log relative transformation
log_relative <- function(physeq, norm.meta=F, select.variables=NULL){
  if(norm.meta){
    get.vars <- get.num.variables(physeq)
    norm.variables <- log(get.vars$num.variables/rowSums(get.vars$num.variables))
    meta_table <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
    sample_data(physeq) <- meta_table
  }
  else if(!norm.meta){
    abund_table <- otu_table(physeq)
    otu_table(physeq) <- log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
  }
  return(physeq)
}

#7.
scale.meta <- function(physeq,type="scale",select.variables=NULL){
  
  get.vars <- get.num.variables(physeq)
  num.variables <- get.vars$num.variables
  
  type <- match.arg(type, c("scale","log","sqrt"))
  
  if(type=="scale"){
    norm.variables <- data.frame(scale(num.variables))
  }
  else if(type=="log"){
    norm.variables <- log2(num.variables)
  }
  else if(type=="sqrt"){
    norm.variables <- sqrt(num.variables)
  }
  
  sample_data(physeq) <- select.vars(norm.variables, get.vars$notnum.variables, select.variables)
  
  return(physeq)
}


# This function is not part of microbiomeSeq
abund.table <- function(pseq, taxo, stats = FALSE) {
  otu_df <- pseq %>% otu_table() %>% as.data.frame() %>% rownames_to_column(var = "otu_id")
  taxo_df <- pseq %>% tax_table() %>% as.data.frame() %>% rownames_to_column(var = "otu_id") %>%
    mutate(Family.Genus = paste(Family, Genus, sep = "_")) 
  otu_df <- left_join(otu_df, taxo_df)
  x <- colnames(otu_table(pseq))
  otu_df <- gather(otu_df, "Sample", "Abundance", x)
  df <- otu_df %>% 
    group_by_at(vars(all_of(taxo), Sample)) %>%
    summarise(Abundance = sum(Abundance)) %>% 
    ungroup() %>% 
    spread(Sample, Abundance) %>%
    as.data.frame() %>%
    'rownames<-'(.[,1]) %>%
    select(-1)
  if(stats == TRUE){ # return a table with statistics summary
    df.stats <- df  %>%  
      mutate(sum = rowSums(.[1:length(.)]),
             mean = rowMeans(.[1:length(.)]),
             stdev= matrixStats::rowSds(as.matrix(.[1:length(.)])),
             freq = round(100*rowSums(.[1:length(.)])/sum(rowSums(.[1:length(.)])),2)) %>%
      select(mean, stdev, freq) 
    row.names(df.stats) <- row.names(df) 
    df.stats <- df.stats  %>% rownames_to_column(var = "taxa") %>% arrange(desc(freq))
    return(df.stats)
  }else{
    return(df)
  }
  rm(otu_df,taxo_df,x)
}

###########################################
# ------- FASTER PSMELT - using dplyr ----
###########################################

psmelt.dplyr = function(physeq) {
  sd = data.frame(sample_data(physeq)) %>% rownames_to_column("Sample")
  TT = data.frame(tax_table(physeq)) %>% rownames_to_column("OTU")
  otu.table = data.frame(otu_table(physeq), check.names = FALSE) %>% rownames_to_column("OTU")
  otu.table %>%
    pivot_longer(!OTU, names_to = "Sample", values_to = "Abundance") %>%
    left_join(sd) %>%
    left_join(TT) 
}

#'Local Contribution to Beta Diversity (LCBD)
#'
#'This function finds the most abundant taxa  for each sample. It also calculates local
#'contribution to beta  diversity using a selected dissimilarity coefficient.
#'It returns a ggplot object which a visual representation of the most abundant
#'taxa for each of the samples. The number of top taxa can be suggested as an argument. The visual
#'representation is limited to 21 top taxa, more information can be availed by supplying a file name to which
#'details of LCBD can be written.
#' @importFrom adespatial beta.div
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param method. A \code{method} to used to calculate dissimilarity coefficents.
#'        Available methods are: "hellinger", "chord", "chisquare", "profiles",
#'        "percentdiff", "ruzicka", "divergence", "canberra",
#'        "whittaker", "wishart", "kulczynski", "jaccard", "sorensen", "ochiai",
#'        "ab.jaccard", "ab.sorensen", "ab.ochiai", "ab.simpson", "euclidean". the function uses
#'        "hellinger" as a default.
#' @param grouping_column(Required). Name of a categorical variable that is preffered for grouping the.
#'        information.
#' @return Returns a ggplot object. This can further be manipulated as preferred by user.
#'
#' @examples plot_taxa(physeq, "Country")
#'
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#'
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#'
#' @export plot_taxa
#'

plot_taxa <- function(physeq,taxo,grouping_column,method="hellinger",filename=NULL){
  
  #==extract components of the phyloseq object
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))
  #Enforce orientation of the phyloseq object
  if(taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }
  
  #===Calculate beta diversity and extract measure for local contribution to beta diversity
  beta_div<-adespatial::beta.div(abund_table,method=method,sqrt.D=F,samp=T,nperm=9999)
  df_LCBD<-data.frame(Sample=names(beta_div$LCBD),LCBD=beta_div$LCBD,p.LCBD=beta_div$p.LCBD)
  
  #=== add grouping information to the LCBD results
  df_LCBD<-data.frame(df_LCBD,Groups=meta_table[rownames(df_LCBD),grouping_column])
  
  if(!is.null(filename)){
    write.csv(df_LCBD,paste(filename,"_LCBD",".csv",sep=""))
  }
  
  #arrange data for plotting in a format compatible to ggplot
  sample_data(physeq)$LCBD <- df_LCBD$LCBD
  sample_data(physeq)$p.LCBD <- df_LCBD$p.LCBD
  sample_data(physeq)$Groups <- df_LCBD$Groups
  sample_data(physeq)$Groups <- factor(sample_data(physeq)$Groups,levels=unique(sample_data(physeq)$Groups))
  
  #arrange data for plotting in a format compatible to ggplot
  ps.glom <- physeq %>%
    #speedyseq::tax_glom(taxrank = taxo) %>%                     # agglomerate at phylum level
    psmelt.dplyr() %>%                                         # Melt to long format
    arrange(taxo, Abundance)                              # Sort data frame alphabetically by phylum
  
  # Creating the  <1% Abundance
 # ps.glom[,length(ps.glom)] <- as.character(ps.glom[,length(ps.glom)]) #convert to character
  
  #simple way to rename phyla with < 1% abundance
 # ps.glom[,length(ps.glom)][ps.glom$Abundance < 0.01] <- " < 1% abund."
  
 # Use !! with the definition operator := as mentioned here, 
  # to set a variable name as the column name. 
  ps.glom <- ps.glom %>% 
    mutate(!!taxo := case_when(Abundance < 0.01 ~ "< 1% abund.",
                            TRUE ~ as.character(.[,taxo])))
  
  # https://dplyr.tidyverse.org/articles/programming.html
  
  # plot the data
  # set color
  library(RColorBrewer)
  colourCount <- ps.glom %>% select(!!taxo) %>% n_distinct()
  # colourCount <-length(unique(ps.glom[,length(ps.glom)]))
  getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
  
  ord_plot = c('ocean_S','C_S','B_S','A_S','C_P','B_P','A_P','GR_S')
  ord_plot2 = c('ocean_S','MH_S','MH_P')
 
  # Plotting
  p <- ps.glom %>% mutate(Groups = factor(Site_Layer, levels = ord_plot)) %>%
  # p <- ps.glom %>% mutate(Groups = factor(Habitat_Layer, levels = ord_plot2)) %>%
  #p <- ps.glom %>%
    ggplot(aes(x = Sample, y = Abundance, fill = ps.glom[,length(ps.glom)])) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = getPalette(colourCount)) +
    facet_grid(. ~ Groups, drop=TRUE,scale="free",space="free_x") +
    guides(keywidth = 1, keyheight = 1) + 
    ylab("Relative Abundance (%)")
  p <- p + guides(fill=guide_legend(ncol=1)) + xlab("Samples")
  p <- p + scale_y_continuous(expand = c(0.02,0)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  p <- p + geom_point(aes(Sample,-0.02, size=LCBD)) + theme(strip.background = element_rect(fill = "white"))
  p <- p + guides(fill=guide_legend(title="Taxa"))
  return(p)
  
} 

##################################################
#------------ Frequency table -------------------
##################################################
# Summarize taxa (https://github.com/joey711/phyloseq/issues/418)

# to summarize a group and get frequency table based on all community
summarize_taxa = function(physeq, Rank, GroupBy = NULL, arrange = FALSE, top_rank = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with psmelt
  ps.melt <- physeq %>%
    #speedyseq::tax_glom(taxrank = Rank) %>%                     # agglomerate at 'Rank' level
    psmelt.dplyr() %>%                                         # Melt to long format
    arrange(Rank)                                     # arrange by 'Rank'
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    summ.df <- ps.melt %>% 
      # group by Rank and Group
      group_by(.data[[GroupBy]], .data[[Rank]]) %>% 
      # summarise and create a column with the relative abundance
      dplyr::summarise(Abundance = sum(Abundance)) %>%
      ungroup() %>%
      mutate(freq = paste0(round(100 * Abundance / sum(Abundance), 2), "%"))
  }else{# only with 'Rank'
    summ.df <- ps.melt %>%
      group_by(.data[[Rank]]) %>%
      dplyr::summarise(Abundance = sum(Abundance)) %>%
      mutate(freq = paste0(round(100 * Abundance / sum(Abundance), 2), '%'))
  }
  if(arrange == FALSE){
    return(summ.df)
    }else{# sort data with the most abundant in the first line
    summ.df.ar <- summ.df %>% arrange(desc(Abundance), .data[[GroupBy]]) %>%
    slice_head(n = 5)
    return(summ.df.ar)
    }

}

# to summarize a group and get frequency table based on within specific taxa group
summarize_by_subtaxa = function(physeq, taxa, Rank, GroupBy){

total <- physeq %>%  
  psmelt.dplyr() %>%
  summarise(Abundance = sum(Abundance)) %>% pull(Abundance)
df <- physeq %>%  
  psmelt.dplyr() %>% 
  group_by(Phylum,.data[[GroupBy]], .data[[Rank]]) %>%
  filter_all(any_vars(. == taxa)) %>%
  summarise(Abundance = sum(Abundance)) %>%
  #  ungroup() %>%
  mutate(freq.total = paste0(round(100 * Abundance / total, 2), '%')) %>% 
  mutate(freq.within.group = paste0(round(100 * Abundance / sum(Abundance), 2), '%')) %>% 
  arrange(desc(Abundance), .data[[GroupBy]])

return(df)
}


##############################
# ------ Box plot -------
##############################

boxplot.pseq <- function(pseq, Rank, GroupBy){
  # Start with psmelt
  ps.melt <- pseq %>%
    #speedyseq::tax_glom(taxrank = Rank) %>%                     # agglomerate at 'Rank' level
    psmelt.dplyr() %>%                                         # Melt to long format
    arrange(Rank)                                     # arrange by 'Rank'
  # x <- levels(ps.melt[,GroupBy])
  Abund <- "Abundance"   # transform as a string factor
  # Let’s generate box plots according to group and 
  # facet them by phylum using the raw counts. 
  box.p <- ps.melt %>% 
    ggplot(aes_(x = as.name(GroupBy), y = as.name(Abund))) +
    geom_boxplot(outlier.shape  = NA) +
    geom_jitter(aes(color = OTU), height = 0, width = .2) +
    labs(x = "", y = "Relative Abundance (%)\n") +
    scale_y_continuous(labels=scales::percent)+
    # scale_x_discrete(labels=x) +
    facet_wrap(as.formula(paste("~", Rank)), scales = "free") +
    theme(legend.position = 'none')
  return(box.p)  
}

########################
# ===== Ordination =====
########################

#############################
# Ordination using factoextra
# -----------------------------

my_ord <- function(pseq, GroupBy, leg){
  # get otu table
  clr_otu <- as.data.frame(otu_table(pseq))
  # performe principal components Analysis
  pcx.abund <- prcomp(t(clr_otu))
  # get the group
  colour <- sample_data(pseq)[,GroupBy]
  colour <- colour[[GroupBy]]
  # set color
  library(RColorBrewer)
  colourCount <- length(unique(colour))
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  # plot
  p1 <- factoextra::fviz_pca_ind(pcx.abund, 
                     # Individuals (samples)
                     geom.ind = "point",
                     habillage=colour, col.ind = "black",
                     pointshape = 19, pointsize = 2,
                     palette = getPalette(colourCount),
                     # addEllipses = TRUE,
                     #  ellipse.level=0.75,   # data ellipses encompass 75% of the points in a group.
                     #ellipse.type ="confidence",
                     geom.var = F,
                     legend.title = leg)
  return(p1)
}

#####################################
# ==== Kmean and RF clustering  ====
#####################################
# return a list with three plots (one MDS based on RF proximity matrix and 2 dendrograms: RF and Euclidean)
# also return the names of the samples and the related cluster based on #clusters (k) selected

plot.cluster.pseq <- function(pseq, k = 3){
  
  clr_otu <- as.data.frame(otu_table(pseq))
  pcx.abund <- prcomp(t(clr_otu))
  King <- unique(data.frame(tax_table(pseq)[,1]))  
  # run the unsupervised RF
  rf2 <- randomForest::randomForest(x = t(clr_otu), ntree = 501, proximity = TRUE)
  # get the distance matrix
  distance.matrix <- as.dist(1-rf2$proximity)
  # create a dataframe to plot the MDS plot
  mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
  ## calculate the percentage of variation that each MDS axis accounts for...
  mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
  
  ## now make a fancy looking plot that shows the MDS axes and the variation:
  mds.values <- mds.stuff$points
  mds.data <- data.frame(Sample=rownames(mds.values),
                         X=mds.values[,1],
                         Y=mds.values[,2],
                         Status=sample_data(pseq)$region)
  p1 <- ggplot(data=mds.data, aes(x=X, y=Y, label=Status)) + 
    geom_text(aes(color=Status)) +
    xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
    ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
    theme_Publication_2()+
    ggtitle(paste("MDS plot using (1 - Random Forest Proximities): ", King))
  
  ## Ward Hierarchical Clustering of the RF output
  hclust.rf <- hclust(as.dist(1-rf2$proximity), method = "ward.D2")
  names.rf.clust = as.factor(cutree(hclust.rf, k=k))
  p2 <- factoextra::fviz_dend(hclust.rf, cex = 0.6,  
                              rect = TRUE, rect_fill = TRUE, 
                              main = paste("Dendrogram using (1 - Random Forest Proximities) - Ward.D2: ", King),
                              xlab = "Samples", ylab = "Distance", sub = "")
  # Using Euclidean Distance and Wrad.D2 Hierarchical Clustering
  hc <- hclust(dist(t(clr_otu), method = 'euclidean'), method = "ward.D2")
  names.hc.cluster = as.factor(cutree(hc, k=k))
  # plot 
  p3 <- factoextra::fviz_dend(hc, cex = 0.6,
                              rect = TRUE, rect_fill = TRUE, 
                              main = paste("Dendrogram using Euclidean Distance - Ward.D2: ", King),
                              xlab = "Samples", ylab = "Distance", sub = "")
  
  # return a list with plots and a vector with the names of clusters
  return(list(p1,p2,p3,names.rf.clust,names.hc.cluster))
}


# performing PerMANOVA to test the hypothesis that bacterial communities within each group
# are more similar to each other than those under other groups. 
# perform PerMANOVA using Euclidean distances (for RDA)
#Generate distance matrix

permanova.pseq <- function(pseq, GroupBy, plot.taxa.coeff = FALSE){
  set.seed(seed)
  clr_dist_matrix <- phyloseq::distance(pseq, method = "euclidean")
  group <- data.frame(phyloseq::sample_data(pseq))[,GroupBy]
  
  #ADONIS test (tests homogeneity of dispersion among groups)
  perm.df <- vegan::adonis(t(abundances(pseq)) ~ group, permutations = 9999, pairwise = TRUE, method = 'euclidean')
  adonis.tb <- perm.df$aov.tab %>% rename(p.value = `Pr(>F)`) %>% as.data.frame() %>% drop_na()
  row.names(adonis.tb)[1] <-'PERMANOVA'
  adonis.tb <- adonis.tb %>% rownames_to_column(var = 'pairs') %>% select(-MeanSqs) 
  adonis.tb$p.adjusted <- NA
  # Pairwise Adonis test 
  pair.df <- pairwiseAdonis::pairwise.adonis(clr_dist_matrix, group, perm = 9999,p.adjust.m = "bonferroni")
  pair.tb <- pair.df %>% select(-sig) 
  #pair.df[[1]] <- NULL
  #pair.tb <- data.table::rbindlist(pair.df, use.names=TRUE, idcol=TRUE) %>% 
  # drop_na() %>% rename(Comparision = .id , p_value = `Pr(>F)`)
  perm.tb <- rbind(adonis.tb, pair.tb) %>%
    mutate(across(where(is.numeric),round, 2)) %>%
    # mutate_if(is.numeric, ~if_else(. < 0.01, "< 0.01", as.numeric(.x)))
    mutate(across(where(is.numeric), ~if_else(. < 0.01, recode(. , .default ="< 0.01"), as.character(.))))
  
  #Dispersion test and plot (it tests whether composition among groups is similar or not)
  dispr <- vegan::betadisper(clr_dist_matrix, group)
  anova(dispr)
  #plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
  #boxplot(dispr, main = "", xlab = "")
  betad <- permutest(dispr,permutations = 9999, pairwise = TRUE)
  betad.tb <- betad$tab %>% rename(p_value_betadisper = `Pr(>F)`) %>% select(p_value_betadisper) %>% drop_na() 
  row.names(betad.tb)[1] <-'PERM_BETAD'
  pair.betad <- betad$pairwise$observed %>% as.data.frame()
  row.names(betad.tb)[1] <-'PERM_BETAD'
  colnames(pair.betad) <- 'p_value_betadisper'
  pair.betad <- pair.betad %>% rownames_to_column('id') %>% mutate(id = str_replace(id, "-"," vs ")) %>% column_to_rownames('id')
  perm.betad.tb <- rbind(betad.tb, pair.betad) %>% 
    rownames_to_column('pairs') %>% mutate(p_value_betadisper = round(p_value_betadisper, 2)) %>%
    # mutate_if(is.numeric, ~if_else(. < 0.01, "< 0.01", as.numeric(.x)))
    mutate(across(where(is.numeric), ~if_else(. < 0.01, recode(. , .default ="< 0.01"), as.character(.))))
  
  # Tukey's HSD
  beta.disper.tukey <- dispr %>% TukeyHSD() %>% tidy() %>% select(-term) %>%
                      mutate(contrast = str_replace(contrast, "-"," vs ")) %>% 
                      mutate(across(where(is.numeric), round, 2)) %>% add_row() 
  # create a list for the output
  out <- list()
  out$permutest <- perm.betad.tb 
  out$adonis <- perm.tb
  out$tukey <- beta.disper.tukey
  if(plot.taxa.coeff == TRUE){
    # check which taxa contribute most to the community differences
    King <- unique(data.frame(tax_table(pseq)[,1]))
    df.coef <- perm.df$coefficients
    for (j in row.names(df.coef)) {
      coef <- perm.df$coefficients[j,]
      tax.tb <- data.frame(tax_table(pseq)) %>% rownames_to_column(var = 'OTU')
      coef.tb <- as.data.frame(coef) %>% rownames_to_column(var = 'OTU') %>% left_join(tax.tb)
      coef.tb$coef <- coef
      top.coef <- coef[rev(order(abs(coef)))[1:20]]
      x <- data.frame(tax_table(pseq)[names(top.coef), c("Phylum","Class","Family", "Genus")]) %>% 
        rownames_to_column(var = 'OTU') %>%
        mutate(id = paste(OTU,Phylum,Class,Family, sep = '_'))
      x$coef <- top.coef
      # plot the top 20 taxa which contribute most to the community differences.
      p <- x %>% ggplot(aes(x = reorder(id, coef), y = coef, fill = coef > 0)) +
        geom_bar(stat="identity") +
        coord_flip() +
        labs(x = 'Taxa id', y = 'Coefficients',
             title = paste0('Top taxa separating the groups - ',King,' ',j)) +
        #subtitle = 'Top taxa separating the groups') 
        theme_Publication_2() + guides(fill = F)
      print(p)
    }
    return(out)
  }else{
    return(out) 
  }
}


#=====================================
# ----- Constrained ordination -----
#=====================================

# dbRDA
contr.ord.pseq <- function(pseq, Rank, GroupBy){
  clr_otu <- otu_table(pseq)
  ord_meta <- data.frame(sample_data(pseq))
  King <- unique(data.frame(tax_table(pseq)[,1]))
  cond <- apply(ord_meta, 2, function(x) any(is.na(x))) # check for NAs
  X <- ord_meta %>% select(where(~is.numeric(.) && all(!is.na(.))))
  X <- vegan::decostand(X, method='standardize')
  Y <- as.data.frame(t(clr_otu)) 
  
  # RDA
  # res <- rda(Y ~ ., data=X)
  # plot(res, xlab=QsRutils::ord_labels(res)[1], ylab=QsRutils::ord_labels(res)[2])
  
  # dbRDA 
  # note the result is the same as for RDA when using euclidean distances
  res <- capscale(Y ~ ., data=X, distance='euclidean')
  # plot(res, xlab=QsRutils::ord_labels(res)[1], ylab=QsRutils::ord_labels(res)[2])
  
  # exploring statistically
  # summary of the results
  # summary(res, display=NULL)
  
  ## select only the variables that are explaining variation efficiently.
  #------------------------------------------------------------------------
  # calculate variance inflation factors (VIF)
  # variables with scores >10 are redundant
  sort(vif.cca(res))
  
  # Fitting environmental vectors/factors onto an ordination 
  # calculate fit
  bio.fit <- envfit(Y ~ ., data=X)
  
  # Stepwise selection (ordistep)
  # set up full and null models for 'ordistep'
  # full model
  rda1 <- rda(Y ~ ., data=X)
  # intercept-only (null) model
  rda0 <- rda(Y ~ 1, data=X)
  
  # perform forward and backward selection of explanatory variables
  # output not shown
  step.env <- ordistep(rda0, Pin = 0.01, scope=formula(rda1), direction='both')
  
  # look at the significant variables 
  #step.env$anova
  # For the sake of argument, let’s include the variables
  # with low VIF scores along with those identified by ordistep.
  # code to get variable names from 'ordistep' and 'envfit' results
  vars.ordistep <- step.env$anova %>% row.names() %>% str_replace('^. ', '')
  vars.biofit <- bio.fit$vectors$pvals %>% as.data.frame() %>% filter(. <=0.01) %>% row.names()
  vars.envfit <- vif.cca(res) %>% as.data.frame() %>% filter(. <=10) %>% row.names() 
  vars <- unique(c(vars.ordistep, vars.envfit,vars.biofit))
  
  # make a dataframe with pvalues and envpars
  vars.df <- data.frame(step.env$anova) %>% 
    rownames_to_column("env_par") %>%
    mutate(env_par = str_trim(str_replace(env_par, "[+]", ""))) %>% 
    # mutate(env_par = str_replace(env_par, "\"", "")) %>% 
    rename(p_value.ordistep = Pr..F.) %>% 
    select(env_par,p_value.ordistep) %>% 
    full_join(data.frame(bio.fit$vectors$pvals) %>% filter(. <=0.01) %>%
                rename(p_value.env.fit = bio.fit.vectors.pvals) %>%
                rownames_to_column("env_par")
    ) %>%
    full_join(as.data.frame(vif.cca(res)) %>% filter(. <=10) %>%
                rename(vif.score = 'vif.cca(res)') %>%
                rownames_to_column("env_par")
    ) %>%
    full_join(bio.fit[["vectors"]][["r"]] %>% data.frame() %>% 
                rename(r_envfit = '.') %>%
                rownames_to_column("env_par") %>%
                mutate(r_envfit = round(r_envfit,2))
    )
  
  #vars 
  # select variables to keep from table 'Y'
  X1 <- X %>% select(matches(vars))
  #str(X1)
  
  # RDA
  res <- rda(Y ~ ., data=X1)
  
  # summary of the results
  # summary(res, display=NULL)
  # permutation test of statistical significance
  # anova(res)
  
  # plotting
  # load library
  
  # set up dataframes for plotting the results
  sit <- ord_meta %>% 
    left_join(vegan::scores(res, display='sites', choices=c(1:2)) %>% 
                                  as.data.frame() %>% 
                                  rownames_to_column("Sample_ID")) %>%
    column_to_rownames("Sample_ID")
  spp <- data.frame(tax_table(pseq)) %>% rownames_to_column('otu') %>% 
    left_join(vegan::scores(res, display='species') %>% 
                data.frame() %>% 
                rownames_to_column('otu')) %>%
                column_to_rownames("otu")
  vec <- vegan::scores(res, display='bp') %>% data.frame()
  
  # use these to adjust length of arrows and position of arrow labels
  adj.vec <- 2
  adj.txt <- 2.5
  arrow_style <- arrow(length = unit(.05, "inches"),
                       type = "closed")
  
  # 'site' ordination
  p1 <- sit %>%
    ggplot(aes(x=RDA1, y=RDA2)) +
    geom_point(aes_string(color= GroupBy, shape = 'Layer') , size=5, alpha = 0.5) +
    geom_segment(data = vec,
                 inherit.aes=F,
                 mapping=aes(x=0, y=0, xend=adj.vec*RDA1, yend=adj.vec*RDA2),
                 arrow = arrow_style) +
    geom_text(data=vec, inherit.aes=F, 
              mapping=aes(x=adj.txt*RDA1, y=adj.txt*RDA2, 
                          label=c(rownames(vec))), 
              hjust = 0, 
              vjust = 1,
              size = 5, 
              color = '#0A537D') +
    ggtitle(paste0('dbRDA: ',King)) +
    # scale_color_manual(values = cols_site_layer) + 
    geom_hline(yintercept=0, linetype="dashed", alpha = 0.3) +
    geom_vline(xintercept=0, linetype="dashed", alpha = 0.3) +
    #labs(colour = GroupBy)+ #another way to set the labels, in this case, for the colour legend
    theme_Publication_3() +
    theme(legend.position = 'bottom') 
  p1$labels$x <- QsRutils::ord_labels(res)[1]
  p1$labels$y <- QsRutils::ord_labels(res)[2]
  
  # 'species' ordination
  # set color
  color.plot <- randomcoloR::distinctColorPalette(k = length(unique(spp[,Rank])))
  
  p2 <- spp %>%
    ggplot(aes(x=RDA1, y=RDA2)) + 
    geom_point(aes_string(color = Rank),size=4, alpha=0.7) + 
    scale_color_manual(values = color.plot) +
    ggtitle(paste0('dbRDA: Taxa - ',King,' - ', Rank)) +
    theme_Publication_3() +
    theme(legend.position = 'bottom') 
  
  p2$labels$x <- QsRutils::ord_labels(res)[1]
  p2$labels$y <- QsRutils::ord_labels(res)[2]
  
  # p3 <- egg::ggarrange(p1,p2, labels = c('A', 'B'), nrow = 2, top = paste0("dbRDA: ",King," - ",Rank))
  
  #--------------------------
  # Variation Partitioning
  #--------------------------
  # extract the sample data from the 'phyloseq' object
  # then remove the 'SampleID' column
  # for only categorical predictors - standardize
  dat <- ord_meta[, c(GroupBy,vars)] 
  
  # use standardised categorical and continuous predictors in the constraint (X)
  X <- ade4::dudi.mix(dat, scannf=F, nf=2)$tab
  # define the partitions
  # MH variables from ordistep and envfit analyses
  X1 <- X %>% select(matches(vars))
  # Group variables
  X2 <- X %>% select(-vars)
  
  # partition variation and plot the result
  vp <- varpart(Y, X1, X2)
  part.plot <- plot(vp, Xnames=c('Env', GroupBy), bg = c("hotpink","skyblue"), main=paste0('Variation partitioning - ',King,' - ', Rank))
  
  # return sit to calculate and plot cluster ellipses
  return(list(p1,p2,vars.df))
  
  
}

###########################
#-----Mantel Test ---------
###########################
# Does Environment Select?
mantel.env.pseq <- function(pseq, scatter.plot = FALSE) {
  abund = data.frame(t(microbiome::abundances(pseq)))
  #abundance data frame - euclidean distance
  dist.abund = dist(abund, method = "euclidean")
  # get the metadata
  meta_mantel <- sample_data(pseq) %>% data.frame()
  # get only numerical variables
  cond <- apply(meta_mantel, 2, function(x) any(is.na(x))) # check for NAs
  env <- meta_mantel %>% select_if(function(col) any(!is.na(col)) & is.numeric(col)) 
  #scale data 
  scale.env = scale(env, center = TRUE, scale = TRUE)
  #create distance matrix of scaled data
  dist.env = dist(scale.env, method = "euclidean") 
  # run mantel with cumulative env 
  abund_env = mantel(dist.abund, dist.env, method = "spearman", permutations = 9999, na.rm = TRUE)
  # run mantel for each env parameter
  mantel.results <- list()
  for (i in unique(names(env))) {
    env.par = meta_mantel[,i]
    dist.par = dist(env.par, method = "euclidean")
    abund_env[[i]] = mantel(dist.abund, dist.par, method = "spearman", permutations = 9999, na.rm = TRUE)
    mantel.results[[i]]$Mantel_r <- abund_env[[i]]$statistic
    mantel.results[[i]]$p.value = abund_env[[i]]$signif
  }
  
  # include the stats and pvalue from the cumulative Mantel
  mantel.results$cumulative$Mantel_r <- abund_env$statistic
  mantel.results$cumulative$p.value <- abund_env$signif
  # extract the r and p.values for each variable from the list
  df <- mantel.results %>% map_df(as.data.frame)
  row.names(df) <- names(mantel.results)
  
  if(scatter.plot == TRUE){
    aa = as.vector(dist.abund)
    tt = as.vector(dist.env)
    
    #new data frame with vectorized distance matrices
    mat = data.frame(aa,tt)
    
    #abundance vs env paramters
    lm = lm(aa~tt, data = mat)
    summary(lm)
    my.formula <- aa~tt
    mm = ggplot(mat, aes(y = aa, x = tt)) + 
      geom_point(size = 3, alpha = 0.5) + 
      geom_smooth(method = "lm", colour = "blue", alpha = 0.2) +
      ggpmisc::stat_poly_eq(formula = my.formula, 
                            aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                            parse = TRUE,
                            label.x.npc = "left", label.y.npc = 75) +   
      ggpmisc::stat_fit_glance(method = 'lm',
                               geom = 'text',
                               aes(label = paste("P-value = ", signif(..p.value.., digits = 4),
                                                 sep = ""))) +
      theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
             axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
             axis.title= element_text(face = "bold", size = 14, colour = "black"), 
             panel.background = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black"),
             legend.position = "top") + 
      labs(x = "Difference in Env Factors", y = "Euclidean Distance")
    print(mm)
    output <- list()
    output$mantel <- df
    output$plot <- mm
    return(output)
  }else{
    return(df)
  }
}

# Is there dispersal limitation in the community?
# Latitude and Longitude need to be named as Lat and Long 
mantel.geo.pseq <- function(pseq, scatter.plot = FALSE) {
  abund = data.frame(t(microbiome::abundances(pseq)))
  #abundance data frame - euclidean distance
  dist.abund = dist(abund, method = "euclidean")
  # #longitude and latitude
  geo <- sample_data(pseq) %>% data.frame() %>% select(Long, Lat)
  #geographic data frame - haversine distance 
  d.geo = geosphere::distm(geo, fun = geosphere::distHaversine)
  dist.geo = as.dist(d.geo)
  # run mantel: abundance vs geographic 
  abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
  abund_geo 
  # include the stats and pvalue from the geo Mantel
  df <- NULL
  df$Mantel_r <- abund_geo$statistic
  df$p.value <- abund_geo$signif
  df <- do.call(rbind, df) %>% t()
  row.names(df) <- 'abund_vs_geo'
  if(scatter.plot == TRUE){
    aa = as.vector(dist.abund)
    tt = as.vector(dist.geo)
    
    #new data frame with vectorized distance matrices
    mat = data.frame(aa,tt)
    
    #abundance vs geographic distance
    lm = lm(aa~tt, data = mat)
    summary(lm)
    my.formula <- aa~tt
    mm = ggplot(mat, aes(y = aa, x = tt)) + 
      geom_point(size = 3, alpha = 0.5) + 
      geom_smooth(method = "lm", colour = "blue", alpha = 0.2) +
      ggpmisc::stat_poly_eq(formula = my.formula, 
                            aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                            parse = TRUE,
                            label.x.npc = "left", label.y.npc = 75) +   
      ggpmisc::stat_fit_glance(method = 'lm',
                               geom = 'text',
                               aes(label = paste("P-value = ", signif(..p.value.., digits = 4),
                                                 sep = ""))) +
      theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
             axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
             axis.title= element_text(face = "bold", size = 14, colour = "black"), 
             panel.background = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black"),
             legend.position = "top") + 
      labs(x = "Physical Separation (Km)", y = "Euclidean Distance")
    print(mm)
    output <- list()
    output$mantel <- df
    output$plot <- mm
    return(output)
  }else{
    return(df)
  }
}

###################################
##### Correlation Matrix ##########
###################################
# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
# the results will be displayed in html or latex format
# source: http://www.sthda.com/english/wiki/elegant-correlation-table-using-xtable-r-package
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 

###########################################
##### Flatten Correlation Matrix ##########
###########################################
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


#### The VIP function bellow was based on Horton et al 2019.
#### https://github.com/horto2dj/GLCW
### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
### `pls' package.
### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $

### Copyright © 2006,2007 Bjørn-Helge Mevik
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License version 2 as
### published by the Free Software Foundation.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.

### A copy of the GPL text is available here:
### http://www.gnu.org/licenses/gpl-2.0.txt

### Contact info:
### Bjørn-Helge Mevik
### bhx6@mevik.net
### Rødtvetvien 20
### N-0955 Oslo
### Norway

### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
### some variable selection methods when multicollinearity is present,
### Chemometrics and Intelligent Laboratory Systems 78, 103--112.

## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}
