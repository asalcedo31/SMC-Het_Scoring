library(BoutrosLab.plotting.general)
library(plyr)
library(reshape2)
#tsv.dir = "scoring_metric_data/text_files/"
#plot.dir = "/u/asalcedo/SMC-HET-Paper1/plots/"
tsv.dir = "~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/"
#tsv.dir = "/u/asalcedo/SMC-Het-Challenge-master/smc_het_eval/scoring_metric_data/text_files/"
plot.dir = "~/Shad_scoring/smc_het_eval/plots/"

metric_behaviour_multi_2A <- function(infile=NULL, outfile=NULL, prop_small=2.5, prop_big=12, ordering=NULL, penalty="spearman"){
  d = read.csv(file=infile, sep="\t",header=F,stringsAsFactors=FALSE)

  # make the header more readable
  header <- sapply(d[1,], as.character)
  colnames(d) <- header
  d <- d[-1,]
	d$iter <- unlist(alply(c(1:100),1, function(x) rep(x,6)))

  CaseName <- c("SplitClusterBot", "MergeClusterBot", "OneCluster", "NCluster", "SmallExtra", "BigExtra");
  CaseOrd<- CaseName[c(4,3,6,2,1,5)]
  if (is.null(ordering)){ordering <- match(CaseOrd,d$Case)}

  weight.diff <- function(d_sub,penalty="spearman"){
  sapply(header[-1], function(col){
    d.wght <- d_sub[,c('Case',col)]
    colnames(d.wght) <- c('Case', 'Metric')

    d.wght$Metric <- as.numeric(d.wght$Metric) # change the metric column from a factor to a numeric value
    d.wght$Metric[d.wght$Metric=='NaN'] <- 0
    ideal.order <- ordering
 #   d.wght <- merge(rank, d.wght, by="Case")
    actual.order <- order(d.wght[,"Metric"])
 #   ideal.order <- order(order(d.wght[ordering,]))# order the columns based on the given value

    if(penalty=="abs"){
      diff <- sum(abs(actual.order - ideal.order))
    } else if(penalty == 'sq'){
      diff <- sqrt(sum((actual.order - ideal.order)^2))
    } else if(penalty == 'spearman'){
      n <- length(actual.order)
      diff <- 1-(6 * (sum((actual.order - ideal.order)^2) / (n * (n^2 - 1))))
    }
    return(diff)
  }) }

  prop_match <- function(d_sub){
    sapply(header[-1], function(col){
    d.wght <- d_sub[,c('Case',col)]
    colnames(d.wght) <- c('Case', 'Metric')
    d.wght$Metric <- as.numeric(d.wght$Metric) # change the metric column from a factor to a numeric value
    d.wght$Metric[d.wght$Metric=='NaN'] <- 0
    ideal.order <- ordering
    actual.order <- order(d.wght[,"Metric"])
    ifelse(identical(ideal.order,actual.order), return(1),return(0))
  })
	}

  corr.rep <- daply(d, .(iter), function(x) weight.diff(x))
  prop.rep  <- daply(d,.(iter), function(x) prop_match(x))
  colnames(corr.rep) <- header[-1]
  colnames(prop.rep) <- header[-1]
  print(colSums(prop.rep))
  js.cast <- dcast(iter~Case,data=d[,c(1,14,21)],value.var="js_divergence")
  js.cast <- data.frame(apply(js.cast,2,function(x) as.numeric(x)))
  js.counts <- data.frame(apply(js.cast[-1],2, function(x) return(table(cut(x,breaks=seq(0,1,0.02),include.lowest=TRUE)))))
  js.counts$bin <- seq(0.01,0.99,0.02)
  js.counts.long <- melt(js.counts,id.vars="bin")
  js.counts.long$value[js.counts.long$value==0] <- NA

  pearson.cast <- dcast(iter~Case,data=d[,c(1,20,21)],value.var="pearson")
  pearson.cast <- data.frame(apply(pearson.cast,2,function(x) as.numeric(x)))
  pearson.counts <- data.frame(apply(pearson.cast[-1],2, function(x) return(table(cut(x,breaks=seq(0,1,0.02),include.lowest=TRUE)))))
  pearson.counts$bin <- seq(0.01,0.99,0.02)
  pearson.counts.long <- melt(pearson.counts,id.vars="bin")
  pearson.counts.long$value[pearson.counts.long$value==0] <- NA


  comp.cast <- dcast(iter~Case,data=d[,c(1,16,21)],value.var="js_divergence+mcc+pearson")
  comp.cast <- data.frame(apply(comp.cast,2,function(x) as.numeric(x)))
  comp.counts <- data.frame(apply(comp.cast[-1],2, function(x) return(table(cut(x,breaks=seq(0,1,0.02),include.lowest=TRUE)))))
  comp.counts$bin <- seq(0.01,0.99,0.02)
  comp.counts.long <- melt(comp.counts,id.vars="bin")
  comp.counts.long$value[comp.counts.long$value==0] <- NA

  avg.cast <- dcast(iter~Case,data=d[,c(1,17,21)],value.var="js_divergence+pearson")
  avg.cast <- data.frame(apply(avg.cast,2,function(x) as.numeric(x)))
  avg.counts <- data.frame(apply(avg.cast[-1],2, function(x) return(table(cut(x,breaks=seq(0,1,0.02),include.lowest=TRUE)))))
  avg.counts$bin <- seq(0.01,0.99,0.02)
  avg.counts.long <- melt(avg.counts,id.vars="bin")
  avg.counts.long$value[avg.counts.long$value==0] <- NA

  corr.long <- melt(corr.rep[,c(13,19,14,15)])

  metric.cov.data <- data.frame(js = c('black','white','black','black'),pcc=c('white','black','black','black'), mcc= c('white','white','white','black'),stringsAsFactors=FALSE)

  cov.heatmap <- create.heatmap(x=metric.cov.data,
                  #              filename="~/test_col.png",
                                clustering.method='none',
                                input.colours=TRUE,
                                print.colour.key=FALSE,
                                yaxis.lab= c('JS','PCC','MCC'))

  compare_rho <- create.boxplot(
    value ~ Var2,
    data=corr.long,
    add.stripplot=TRUE,
    xaxis.lab=c('','','',''),
    xat= c(),
  #  filename="~/metric_boxplot.png",
    xaxis.rot = 40,
    xaxis.tck=0,
    xaxis.cex = 0.8,
    ylab.label = expression(bold(paste("Spearman's ",rho))),
    xlab.label='',
    ylab.cex=2,
    yaxis.cex=1.2,
    yaxis.fontface='plain',
    resolution=200)

  compare_rho_mpp <- create.multipanelplot(
                      plot.objects = list(compare_rho,cov.heatmap),
                      plot.objects.heights=c(1, 0.12),
                      resolution=200,
                   #   filename="~/rho_metric_multiplot.png"
                      )

  xat=seq(0.5,50.5,10)

  col.df <- data.frame(case = levels(js.counts.long$variable), col=default.colours(6),stringsAsFactors=FALSE)
  rownames(col.df) <- col.df$case
  col.df <- col.df[CaseName[ordering],]

  nice_names <- c("Split Cluster", "Merge Cluster", "One Cluster", "N Clusters", paste0("Extra Cluster: ", prop_small*100,"% of ssms"), paste0("Extra Cluster: ", prop_big*100,"% of ssms"));
  js_hist <- create.barplot( value~ bin,
    data = js.counts.long,
    groups=js.counts.long$variable,
    border.col = default.colours(6),
    col= 'transparent',
    xlab.label = 'JS Divergence',
    ylab.label = 'Count',
    xaxis.lab = c(unique(js.counts.long$bin),1.01)[xat+0.5]-0.01, 
  #  filename="~/js_hist.png",
    xat=xat,
    box.ratio=10,
    
    resolution=200)

  pearson_hist <- create.barplot( value~ bin,
    data = pearson.counts.long,
    groups=pearson.counts.long$variable,
    border.col = default.colours(6),
    col= 'transparent',
    xlab.label = 'PCC/MCC',
    ylab.label = 'Count',
    xaxis.lab = c(unique(pearson.counts.long$bin),1.01)[xat+0.5]-0.01, 
   # filename="~/pearson_hist.png",
    xat=xat,
    box.ratio=10,   
    resolution=200)

  avg_hist <- create.barplot( value~ bin,
    data = avg.counts.long,
    groups=avg.counts.long$variable,
    border.col = default.colours(6),
    col= 'transparent',
    xlab.label = 'JS+PCC',
    ylab.label = 'Count',
    xaxis.lab = c(unique(avg.counts.long$bin),1.01)[xat+0.5]-0.01, 
   # filename="~/pearson_hist.png",
    xat=xat,
    box.ratio=10,   
    resolution=200)

  comp_hist <- create.barplot( value~ bin,
    data = comp.counts.long,
    groups=comp.counts.long$variable,
    border.col = default.colours(6),
    col= 'transparent',
    xlab.label = 'JS+PCC+MCC',
    ylab.label = 'Count',
    xaxis.lab = c(unique(comp.counts.long$bin),1.01)[xat+0.5]-0.01, 
  #  filename="~/comp_hist.png",
    xat=xat,
    box.ratio=10,
    
    resolution=200)

legendG <- legend.grob(
  list( 
    legend = list(   colours = col.df$col,
                      labels = nice_names[ordering],
                      title = 'Case',
                      title.cex =1,
                      label.cex = 1,
                      border='transparent',
                      size =3 )
    ,
    legend = list(title='Metric',
                  colours=c('black','white'),
                  labels = c('included','not included')
                  )
  ))
  create.multipanelplot(
    plot.objects = list(js_hist,pearson_hist,comp_hist, avg_hist),
    layout.width=2,
    layout.height=3,
    x.spacing=1,
    y.spacing = c(0,-1.2),

    legend = list(right=list(fun=legendG)),
    layout.skip=c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE),
    plot.objects.heights = c(1,1,0.2),
    filename=outfile,
    width=12,
    resolution=200
    )

}

infile=paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000.tsv")
outfile="~/2A_1k_6clust.png"
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.1_sep0.02.tsv"), "~/2A_1k_6clust_bep0.1_sep0.02.png",prop_big=0.1, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.15_sep0.02.tsv"), "~/2A_1k_6clust_bep0.15_sep0.02.png",prop_big=0.15, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.1_sep0.01.tsv"), "~/2A_1k_6clust_bep0.1_sep0.01.png",prop_big=0.1, prop_small=0.01)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.1_sep0.05.tsv"), "~/2A_1k_6clust_bep0.1_sep0.05.png",prop_big=0.1, prop_small=0.05)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc3_s1000_bep0.1_sep0.02.tsv"), "~/2A_1k_3clust_bep0.1_sep0.02.png",prop_big=0.1, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.6_sep0.02.tsv"), "~/2A_1k_6clust_bep0.6_sep0.02.png",prop_big=0.6, prop_small=0.02)

metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc3_s1000_bep0.15_sep0.02.tsv"), "~/2A_1k_3clust_bep0.15_sep0.02.png",prop_big=0.15, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc3_s1000_bep0.1_sep0.01.tsv"), "~/2A_1k_3clust_bep0.1_sep0.01.png",prop_big=0.1, prop_small=0.01)


metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.1_sep0.02.tsv"), "~/2A_1k_6clust_bep0.1_sep0.02_avg.png",prop_big=0.1, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.15_sep0.02.tsv"), "~/2A_1k_6clust_bep0.15_sep0.02_avg.png",prop_big=0.15, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.1_sep0.01.tsv"), "~/2A_1k_6clust_bep0.1_sep0.01_avg.png",prop_big=0.1, prop_small=0.01)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.1_sep0.05.tsv"), "~/2A_1k_6clust_bep0.1_sep0.05_avg.png",prop_big=0.1, prop_small=0.05)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc3_s1000_bep0.1_sep0.02.tsv"), "~/2A_1k_3clust_bep0.1_sep0.02_avg.png",prop_big=0.1, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc6_s1000_bep0.6_sep0.02.tsv"), "~/2A_1k_6clust_bep0.6_sep0.02_avg.png",prop_big=0.6, prop_small=0.02)

metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc3_s1000_bep0.15_sep0.02.tsv"), "~/2A_1k_3clust_bep0.15_sep0.02_avg.png",prop_big=0.15, prop_small=0.02)
metric_behaviour_multi_2A(paste0(tsv_dir,"weights2A_all_cases_js_divergence_mcc_pearson_nc3_s1000_bep0.1_sep0.01.tsv"), "~/2A_1k_3clust_bep0.1_sep0.01_avg.png",prop_big=0.1, prop_small=0.01)
