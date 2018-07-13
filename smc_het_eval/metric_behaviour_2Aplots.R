library(BoutrosLab.plotting.general)
library(plyr)
library(reshape2)
library(png)

std <- function(x) sd(x)/sqrt(length(x))

#tsv.dir = "scoring_metric_data/text_files/"
#plot.dir = "/u/asalcedo/SMC-HET-Paper1/plots/"
tsv.dir = "~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/"
#tsv.dir = "/u/asalcedo/SMC-Het-Challenge-master/smc_het_eval/scoring_metric_data/text_files/"
plot.dir = "~/Shad_scoring/smc_het_eval/plots/"

get.ranking <- function(scenarios.i=NA, aggr.meth='Copeland', filename="aggregate_scenario_rankings.csv", plot=TRUE){
  rank <- read.csv(file=paste("/u/asalcedo/SMC-HET-Paper1/smc_het_eval/scoring_metric_data/text_files/", filename, sep=""), sep=",", header=TRUE)
  if(is.na(scenarios.i)){ # if scenarios.i is NA then include all rows/mistake scenarios
    scenarios.i <- 1:dim(rank)[1]
  }
  colnames(rank)[1] <- 'Case'
  rownames(rank) <- rank$Case
  rank <- rank[scenarios.i,!(names(rank) %in% c('Copeland', 'Borda', 'Std.Dev', 'SD'))] # consider only the mistake scenarios specified
  
   # TODO: remove this and name the cross validation tsv files correctly
  
  num.sc <- dim(rank)[1] # number of scenarios
  
  aggr.rank <- rep(0, num.sc) # aggregate score
  for(col in 2:dim(rank)[2]){
    for(i in 1:(num.sc-1)){
      for(j in (i+1):num.sc){
        if(rank[i, col] > rank[j, col]){ # if i is better than j
          aggr.rank[i] <- aggr.rank[i] + 1
          if(aggr.meth == 'Copeland'){
            aggr.rank[j] <- aggr.rank[j] - 1
          }
        } else if(rank[i, col] < rank[j, col]){ # if j is better than i
          aggr.rank[j] <- aggr.rank[j] + 1
          if(aggr.meth == 'Copeland'){
            aggr.rank[i] <- aggr.rank[i] - 1
          }
        } else if(aggr.meth == 'Borda'){ # if i and j are equally good
          aggr.rank[i] <- aggr.rank[i] + 0.5
          aggr.rank[j] <- aggr.rank[j] + 0.5
        }
      }
    }
  }
  ranking <- data.frame(rank$Case, aggr.rank, stringsAsFactors=FALSE)
  colnames(ranking) <- c('Case', aggr.meth)
  return(ranking[order(ranking[,2]),])
}


  weight.diff <- function(d_sub, ordering,penalty="spearman"){
  sapply(colnames(d_sub)[3:ncol(d_sub)], function(col){
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

  prop_match <- function(d_sub, ordering){
    sapply(colnames(d_sub)[3:ncol(d_sub)], function(col){
    d.wght <- d_sub[,c('Case',col)]
    colnames(d.wght) <- c('Case', 'Metric')
    d.wght$Metric <- as.numeric(d.wght$Metric) # change the metric column from a factor to a numeric value
    d.wght$Metric[d.wght$Metric=='NaN'] <- 0
    ideal.order <- ordering
    actual.order <- order(d.wght[,"Metric"])
    ifelse(identical(ideal.order,actual.order), return(1),return(0))
  })
  }
  

  contrast_stats <- function(infile,sc,plot=TRUE){
    d = read.csv(file=infile, sep="\t",header=F,stringsAsFactors=FALSE)
    CaseName <- c("SplitClusterBotSame", "MergeClusterMid&BotMultiChild", "OneCluster", "NClusterOneLineage", "SmallExtraMid", "BigExtraMid");
    d <- d[d$V1 %in% CaseName,]
    if(sc == '2B' | sc =='3B'){
      colnames(d) <- c("Case","iter", "js_divergence","mcc", "aupr", "pearson")
      d$js_mcc_aupr_pearson <- rowMeans(d[,3:6])
      d$js_mcc_pearson <- rowMeans(d[,c(3:4,6)])
      d$mcc_pearson <- rowMeans(d[,c(4,6)])
    } else{
      colnames(d) <- c("Case","iter", "js_divergence","mcc", "aupr")    
    }
    ord <- get.ranking(CaseName)
    ordering <- match(ord$Case, d$Case)
    print(ord$Case[ordering])
    d$js_mcc_aupr <- rowMeans(d[,3:5])
    d$js_mcc <- rowMeans(d[,3:4])
    d$js_mcc_mcc <- rowMeans(d[,3:4])
    print(head(d))

    d$Case <- as.factor(d$Case)
    worst <-c(1/4,1/4,-0.5,-0.5,1/4,1/4)
    big <- c(-1,0.5,0,0,0,0.5)
    small <- c(0,-0.5,0,0,1,-0.5)
    other1 <- c(0,0,0,-1,1)
    other2 <- c(0,1,0,0,0,-1)
    
    # contrasts.m <- rbind(rep(1/6,6), worst,big,small,other1,other2)
    contrasts.m <- cbind(worst, big,small)
    # contrasts.m <- rbind(rep(1/6,6), worst,big,small,o1,02)
    # c.inv <- solve(contrasts.m)
    contrasts(d$Case) <- contrasts.m
    mod <- alply(c(3:ncol(d)), 1, function(x){f <- as.formula(paste(colnames(d)[x],"~Case")); lm(f, data=d, contrasts=list(Case=contrasts.m))}) 
    names(mod) <- colnames(d)[3:ncol(d)]
    return(mod)
  }

  correlation_stats <- function(infile,sc,plot=TRUE){
    d = read.csv(file=infile, sep="\t",header=F,stringsAsFactors=FALSE)
    if (sc=='2A' | sc == '2B'){
       CaseName <- c("SplitClusterBotSame", "MergeClusterMid&BotMultiChild", "OneCluster", "NClusterOneLineage", "SmallExtraMid", "BigExtraMid");
       d <- d[d$V1 %in% CaseName,]
    } else {
      CaseName <- NA
    }
    if(sc == '2B' | sc =='3B'){
      colnames(d) <- c("Case","iter", "js_divergence","mcc", "aupr", "pearson")
      d$js_mcc_aupr_pearson <- rowMeans(d[,3:6])
      d$js_mcc_pearson <- rowMeans(d[,c(3:4,6)])
      d$mcc_pearson <- rowMeans(d[,c(4,6)])
    } else{
      colnames(d) <- c("Case","iter", "js_divergence","mcc", "aupr")    
    }
    print(head(d))
    ord <- get.ranking(CaseName)
    ordering <- match(ord$Case, d$Case)
    print(ord$Case[ordering])
    d$js_mcc_aupr <- rowMeans(d[,3:5])
    d$js_mcc <- rowMeans(d[,3:4])
    d$js_mcc_mcc <- rowMeans(d[,c(3,4,4)])
    
    corr.rep <- daply(d, .(iter), function(x) weight.diff(x,ordering))
    prop.rep  <- daply(d,.(iter), function(x) prop_match(x,ordering))
    
    colnames(corr.rep) <- colnames(d)[3:ncol(d)]
    corr.long <- melt(corr.rep)
    anova.corr <- anova(lm(value~Var2, corr.long))

  if (plot==TRUE){
     compare_rho <- create.boxplot(
      value ~ Var2,
      data=corr.long,
      add.stripplot=TRUE,
      xaxis.lab= TRUE,
      xat= c(),
      filename=paste("~/metric_boxplot",sc,".png",sep=""),
      xaxis.rot = 40,
      xaxis.tck=0,
      xaxis.cex = 0.8,
      ylab.label = expression(bold(paste("Spearman's ",rho))),
      xlab.label='',
      ylab.cex=2,
      yaxis.cex=1.2,
      yaxis.fontface='plain',
      resolution=200)

      # metric.cov.data <- data.frame(js = c('black','white','black','black'),pcc=c('white','black','black','black'), mcc= c('white','white','white','black'),stringsAsFactors=FALSE)

      # cov.heatmap <- create.heatmap(x=metric.cov.data,
      #                 #              filename="~/test_col.png",
      #                               clustering.method='none',
      #                               input.colours=TRUE,
      #                               print.colour.key=FALSE,
      #                               yaxis.lab= c('JS','PCC','MCC'))

      # compare_rho_mpp <- create.multipanelplot(
      #                     plot.objects = list(compare_rho,cov.heatmap),
      #                     plot.objects.heights=c(1, 0.12),
      #                     resolution=200,
      #                  #   filename="~/rho_metric_multiplot.png"
      #                     )
  }
    return(list(anova.corr,corr.rep))
 }


test_contrasts <- contrast_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='2A')
contrasts_2B <- contrast_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='2B')

test3B <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='3B', plot=FALSE)
test2B <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='2B', plot=FALSE)
test3A <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='3A', plot=FALSE)
test2A <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='2A', plot=FALSE)


test3B <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='3B')
test2B <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='2B')
test3A <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='3A')
test2A <- correlation_stats(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='2A')

rho_table <- ldply(list(test2A[[2]],test2B[[2]],test3A[[2]],test3B[[2]]), function(x) {means <- colMeans(x); se <- apply(x, 2, function(y) std(y)); return(rbind(means,se))})
rho_table <- apply(rho_table,2, function(x) round(x,3))
rho_table <- rho_table[,c()]
write.table(rho_table,"~/spearmans_rho_table.txt",sep="\t", quote=FALSE, row.names=FALSE)
metric_behaviour_sp <- function(infile=NULL, outfile=NULL, sc='2A', xlab=NULL, print.legend=FALSE){
 
  d = read.csv(file=infile, sep="\t",header=F,stringsAsFactors=FALSE)
  names(d) <- c("x","nssm","js2A","mcc2A","aupr2A","js3A","mcc3A","aupr3A")
  if(sc=='2A'){
    w <- d[,c(1,3,4,5)]
    } else if (sc=='3A'){
    w <- d[,c(1,6,7,8)]
  }
  l <- melt(w,id.vars="x")
  l$variable <- gsub(sc,"",l$variable)
 # print(head(l))
  l_means <- ddply(l, .(x,variable), summarize, mean=mean(value), se=sd(value)/sqrt(length(value)))
  l_means$variable <- as.factor(l_means$variable)
  print(l_means)
 metric.key = list(
    text = list(
      lab = c("AUPR","JS Divergence","Pearson/MCC"),
      cex = 1,
      col = 'black'
      ),
    points = list(
      pch = 19,
      col =  default.colours(3),
      cex = 0.75
      ),
    x = 0.04,
    y = 0.95,
    padding.text = 1
    )
legend.print=NULL
  if(print.legend == TRUE){
    legend.print = list(inside = list(fun = draw.key, args = list(key=metric.key), x=0.65, y = 0.95))
  }


  create.scatterplot(mean ~ x,
  data = l_means,
  groups = l_means$variable,
  col = default.colours(3),
  y.error.up = l_means$se,
  y.error.bar.col = default.colours(3),
  ylab.label = "Score",
  xlab.label = xlab,
  file=outfile,
  legend = legend.print,
  xaxis.cex = 1.2,
  yaxis.cex =1.2,
  ylab.axis.padding=2,

  resolution=200)
}

metric_behaviour_sp(infile="out_prop2A3A.txt", outfile="prop_extra_test.png", sc='2A', xlab="Proportion of mutations", print.legend=TRUE)
metric_behaviour_sp(infile="out_prop2A3A.txt", outfile="prop_extra_test_3A.png", sc='3A', xlab="Proportion of mutations", print.legend=FALSE)

metric_behaviour_sp(infile="out_prop_merged2A3A.txt", outfile="prop_merged_test.png", sc='2A', xlab="Proportion of mutations", print.legend=FALSE)
metric_behaviour_sp(infile="out_prop_split2A3A.txt", outfile="prop_split_test.png", sc='2A', xlab="Proportion of mutations", print.legend=FALSE)

metric_behaviour_sp(infile="out_nssm_small_extra_mid_2A3A.txt", outfile="prop_nssm_test.png", sc='2A', xlab="Number of mutations", print.legend=FALSE)
metric_behaviour_sp(infile="out_nssm_small_extra_mid_2A3A.txt", outfile="prop_nssm_test3A.png", sc='3A', xlab="Number of mutations", print.legend=FALSE)



metric_behaviour_sp_f2 <-  function(infile=NULL, outfile=NULL, sc='2A', xlab=NULL, print.legend=FALSE){
  d = read.csv(file=infile, sep="\t",header=F,stringsAsFactors=FALSE)
  d = d[1:50,]
  if (sc=='2A' | sc == '2B'){
       CaseName <- c("SplitClusterBotSame", "MergeClusterMid&BotMultiChild", "OneCluster", "NClusterOneLineage", "SmallExtraCurBot", "BigExtraCurBot");

      xlim=c(0,7)
      xlab=c(1,2,4,5,6,8)

    } else {
        CaseName <-  c("NClusterOneLineage", "OneCluster", "ParentIsSiblingWithChildren", "SmallExtraCurBot", "BigExtraCurBot", "MergeClusterTop&Mid", "MergeClusterMid&BotMultiChild", "SplitClusterBotSame")
        xlim=c(0,9)
        xlab=seq(1:8)
      }

    d <- d[d$V1 %in% CaseName,]
    if(sc == '2B' | sc =='3B'){
      colnames(d) <- c("CaseName","iter", "js_divergence","mcc", "aupr", "pearson")
    } else{
      colnames(d) <- c("CaseName","iter", "js_divergence","mcc", "aupr")    
    }
    # print(head(d))
    ord <- get.ranking(CaseName)
    ordering <- match(ord$Case, d$CaseName)
   CaseOrd <- ord$Case   
  #CaseOrd<- CaseName[c(4,3,6,2,5,1)]
  d_long <- melt(d,id.vars=c('CaseName','iter'))
  d_means <- ddply(d_long, .(variable, CaseName), summarise, mean=mean(value), se=sd(value)/sqrt(length(value)))
  d_means$ord <- match(d_means$CaseName, CaseOrd)
  if (sc=='2A' | sc == '3A'){

   metric.key = list(
    points = list(
      pch = 19,
      col =  rev(default.colours(3)),
      cex = 0.3
      ),
    text = list(
      lab = c("JS Divergence","MCC/PCC","AUPR"),
      cex = 0.5,
      col = 'black'
      ),
    x = 0,
    y = 0.97,
    padding.text = 1
    )
  d_means$variable <- factor(d_means$variable, levels=c("aupr","mcc","js_divergence"), labels = c(1,2,3))
   n_col=3
  }  else{
   cols = default.colours(4)
   metric.key = list(
    points = list(
      pch = 19,
      col =  cols[c(3,2,4,1)],
      cex = 0.3
      ),
    text = list(
      lab = c("JS Divergence", "MCC","PCC","AUPR"),
      cex = 0.55,
      col = 'black'
      ),
    x = 95,
    y = 0.3,
    padding.text = 1
    )
  d_means$variable <- factor(d_means$variable, levels=c("aupr","mcc","js_divergence","pearson"), labels = c(1,2,3,4))
   n_col=4
 }
legend.print=NULL
  if(print.legend == TRUE){
    legend.print = list(inside = list(fun = draw.key, args = list(key=metric.key), x=0.48, y = 0.3))
  }

  d_means <- d_means[order(d_means$ord, d_means$variable),]
  sp <- create.scatterplot(mean ~ ord,
  data = d_means,
  groups = d_means$variable,
  col = default.colours(n_col),
  y.error.up = d_means$se,
  y.error.bar.col = default.colours(n_col),
  ylab.label = '', #"Score",
  xlab.label = '', #xlab,
  file=outfile,
  xlim=xlim,
  xaxis.lab=xlab,
  cex = 0.3,
  xaxis.tck=c(0.3,0),
  yaxis.tck=c(0.3,0),
  xat = seq(1,(xlim[2]-1)),
  type=c('p','l'),
  legend = legend.print,
  xaxis.cex = 0.5,
  yaxis.cex =0.5,
  error.bar.length=0.05,
  ylab.axis.padding=2,
  xaxis.fontface = 'plain',
  yaxis.fontface = 'plain',
  #filename="~/test_2A.png",
  resolution=200)

  return(sp)
}


sp3B <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='3B',print.lege)
#sp2B <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='2B', print.legend=FALSE)
#sp3A <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='3A')
#sp2A <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='2A',print.legend=FALSE)

cols=default.colours(4)
mpp <- create.multipanelplot(
  plot.objects = list(sp2A, sp3A,sp2B, sp3B),
  layout.width = 2, 
  layout.height=2,
  #plot.objects.heights=c(1, 0.12),
  y.spacing=0.7,
  x.spacing=0.7,
  left.padding=0,
  right.padding=0,
  ylab.axis.padding=0,
  xlab.axis.padding=0,
  xlab.cex=0.7,
  ylab.cex=0.7,
  resolution=200,
  # ylab.label='Score',
  # xlab.label='Case',
  height=5,
  width=5,
  #filename="~/F2_multiplot.png"
  )
metric.leg  = list(
    legend=list(points = list(
      pch = 19,
      col =  cols[c(3,2,4,1)],
      cex = 0.3
      ),
    text = list(
      lab = c("JS Divergence", "MCC","PCC","AUPR"),
      cex = 0.5,
      col = 'black'
      )
    ))

legendG <- legend.grob(
  list( 
    metric.leg
  ))


mp.image <- readPNG("~/cluster/Shad_scoring/smc_het_eval/plots/figure.2.5.png",native=TRUE)
#mp.Raster <- as.raster(mp.image)
png("~/cluster/Shad_scoring/smc_het_eval/plots/figure2_multi.png", height=7, width=7, unit="in",res=400)
grid.newpage()
lay0 <- grid.layout(2,1,heights=unit(c(2.5,4.5),c("in","in")),widths=unit(c(7),c("in","in")))
vp0 <- viewport(layout=lay0,name="top.vp",height=unit(7,"in"),width=unit(7,"in"))
pushViewport(vp0)

pushViewport(viewport(layout.pos.row=1,layout.pos.col=1,name="scenarios", height=2.5, width=7))

grid.raster(mp.image,interpolate=TRUE,x=unit(0.5,"npc"),y=unit(0,"npc"),just=c("center","bottom"))
grid.text("a",x=unit(0.15,"npc"), y=unit(0.9,"npc"), gp=gpar(font.family='Arial', fontface='bold',fontsize=12))

upViewport()
lay1 <- grid.layout(2,2,heights=unit(c(2.25,2.25),c("in","in")),widths=unit(c(2.25,2.25),c("in","in")))

pushViewport(viewport(layout.pos.row=2,layout.pos.col=1,layout=lay1,name="scatterplots",height=unit(4.5,"in"), width=unit(7,"in"), clip=FALSE))
#grid.text("Score", y=unit(0.55,"npc"),x=unit(0.12, "npc"),rot=90,gp=gpar(font.family='Arial', fontface='bold', fontsize=18))
grid.text("Case", y=unit(0.05,"npc"),x=unit(0.55, "npc"),gp=gpar(font.family='Arial', fontface='bold', fontsize=18))
#grid.rect(gp=gpar(col='red'))
pushViewport(viewport(layout.pos.row=2, layout.pos.col=2,name="sp3B",height=unit(2.25,"in"), width=unit(3.5,"in")))
grid.text("e",x=unit(0.15,"npc"), y=unit(1.04,"npc"), gp=gpar(font.family='Arial', fontface='bold',fontsize=12))
sp3B <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='3B', print.legend=TRUE)
print(sp3B,newpage=FALSE,panel.width=list(1.75,"in"),panel.height=list(1.75,"in"), position=c(0.2,0.1,1.2,1.1))
upViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=1,name="sp2B",height=unit(2.25,"in"), width=unit(3.5,"in")))
sp2B <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2B_cases__nc6_rep50_s1200_bep0.15_sep0.05.txt",sc='2B', print.legend=FALSE)
print(sp2B,newpage=FALSE,panel.width=list(1.75,"in"),panel.height=list(1.75,"in"), just=c("bottom","right"),position=c(-0.1,0.1,1.1,1.1) )
grid.text("d",x=unit(0,"npc"), y=unit(1.04,"npc"), gp=gpar(font.family='Arial', fontface='bold',fontsize=12))
grid.text("Score",x=unit(0.0,"npc"), y=unit(0.6,"npc"), rot=90,gp=gpar(font.family='Arial', fontface='bold',fontsize=10))

upViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=2,name="sp3A",height=unit(2.25,"in"), width=unit(3.5,"in")))
grid.text("c",x=unit(0.2,"npc"), y=unit(1.05,"npc"), gp=gpar(font.family='Arial', fontface='bold',fontsize=12))
sp3A <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring3A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='3A')
print(sp3A,newpage=FALSE, panel.width=list(1.75,"in"),panel.height=list(1.75,"in"), just=c("bottom","left"),position=c(0.2,0.07,1.2,1.07))
upViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=1,name="sp2A",height=unit(2.25,"in"), width=unit(3.5,"in")))
sp2A <- metric_behaviour_sp_f2(infile="~/Shad_scoring/smc_het_eval/scoring_metric_data/text_files/scoring2A_cases__nc6_s1200_bep0.15_sep0.05.txt",sc='2A',print.legend=FALSE)
print(sp2A,newpage=FALSE, panel.width=list(1.75,"in"),panel.height=list(1.75,"in"), just=c("bottom","left"),position=c(-0.1,0.07,1.1,1.07) )
grid.text("b",x=unit(0,"npc"), y=unit(1.05,"npc"), gp=gpar(font.family='Arial', fontface='bold',fontsize=12))
grid.text("Score",x=unit(0.01,"npc"), y=unit(0.6,"npc"), rot=90,gp=gpar(font.family='Arial', fontface='bold',fontsize=10))
dev.off()



metric_behaviour_multi_2A <- function(infile=NULL, outfile=NULL, prop_small=2.5, prop_big=12, ordering=NULL, penalty="spearman"){
  d = read.csv(file=infile, sep="\t",header=F,stringsAsFactors=FALSE)

  # make the header more readable
  header <- sapply(d[1,], as.character)
  colnames(d) <- c("CaseName", "js_divergence","mcc")
  d <- d[-1,]
	d$iter <- unlist(alply(c(1:100),1, function(x) rep(x,6)))

  CaseName <- c("SplitClusterBot", "MergeClusterBot", "OneCluster", "NCluster", "SmallExtra", "BigExtra");
  CaseOrd<- CaseName[c(4,3,6,2,1,5)]
  if (is.null(ordering)){ordering <- match(CaseOrd,d$Case)}
  d_sub <- d[d$CaseName %in% CaseOrd,]

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
