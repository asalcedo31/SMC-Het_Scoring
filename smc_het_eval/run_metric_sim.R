source('~/Shad_scoring/smc_het_eval/metric_test_sim.R')

n.iter=500

sims <- run.sim(n.iter =100,
		C.values = 2:6,
		K.u.values = 2:10,
		K.n.values = 1:2,
		e1.values=seq(0,0.1,by=0.01),#c(0.02,0.06,0.08,0.1),
		e2.values= seq(0,0.2,0.02)#c(0.066,0.133,0.2))
)

sims <- run.sim(n.iter =100,
        C.values = 5,
        K.u.values = 2:5,
        K.n.values = 0:2,
        e1.values=seq(0,0.1,by=0.01),#c(0.02,0.06,0.08,0.1),
        e2.values= seq(0,0.2,0.02)#c(0.066,0.133,0.2))
)


sims <- run.sim(n.iter =100,
        C.values = 5,
        K.u.values = 6:10,
        K.n.values = 0:2,
        e1.values=seq(0,0.1,by=0.01),#c(0.02,0.06,0.08,0.1),
        e2.values= seq(0,0.2,0.02)#c(0.066,0.133,0.2))
)

sims <- read.sim(n.iter =100,
		C.values = 4,
		K.u.values = 2:9,
		K.n.values = 1:3,
		e1.values=seq(0,0.4,by=0.05),#c(0.02,0.06,0.08,0.1),
		e2.values= seq(0,0.4,0.05)#c(0.066,0.133,0.2))
)


run.sim(n.iter =1,
		C.values = 2,
		K.u.values = 3,
		K.n.values = 1,
		e1.values=0.1,#c(0.02,0.06,0.08,0.1),
		e2.values= 0.2)

plot.sim.batch(sims)

p1 <- test.p1( C=3:6,K.u.values=2:9,K.n.values=3,e1.values=0.4,e2.values=0.4,metrics=c("Pseudo.V","MCC","AUPR","X2A_metric"))
p2 <- test.p2( C=2:6,K.u.values=2:9,K.n.values=3,e1.values=0.4,e2.values=0.4,metrics=c("Pseudo.V","MCC","AUPR","X2A_metric"))
p3 <- test.p3( C=4,K.u.values=2:9,K.n.values=3,e1.values=seq(0,0.4,by=0.05),e2.values=0.4,metrics=c("Pseudo.V","MCC","AUPR","X2A_metric"))
p4 <- test.p4( C=4,K.u.values=6,K.n.values=c(1,2,3),e1.values=0.4,e2.values=seq(0,0.4,0.05),metrics=c("Pseudo.V","MCC","AUPR","X2A_metric"))

prop.heatmap <- plot.properties.heatmap(p1,p2,p3,p4,"scoring_metric_data/sim_plots/metric_properties_heatmap2.png")

useful <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pearson","Pseudo.V","MCC","AUPR","X2A_metric"), filename="scoring_metric_data/sim_plots/useful2.png",param.interest="K.u",print.legend=FALSE)

useful.prop <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pearson","Pseudo.V","MCC","AUPR","X2A_metric"), filename="scoring_metric_data/sim_plots/useful_prop2.png",param.interest="epsilon1")
useful.noise <- plot.sim(sims,c(6,NA,0.2,NA),metrics =c("Pearson","Pseudo.V","MCC","AUPR","X2A_metric"), filename="scoring_metric_data/sim_plots/useful_noise2.png",param.interest='epsilon2',print.legend=FALSE)


panel.size = 2.0
font.y = 0.95
font.cex = 0.8
pdf("scoring_metric_data/sim_plots/all_sim3.pdf",height=7,width=7.1)
grid.newpage()
lay1 <- grid.layout(2,2,heights=unit(c(3.5,3.5),c("in","in")),widths=unit(c(3.5,3.5),c("in","in")))
vp1 <- viewport(layout=lay1,name="top.vp",height=unit(7,"in"),width=unit(7,"in"))
pushViewport(vp1)
pushViewport(viewport(layout.pos.row=1,layout.pos.col=1,name="useful_prop"))
useful.prop <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric"), filename="scoring_metric_data/sim_plots/useful_prop.png",param.interest="epsilon1")
print(useful.prop,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.1,0,1,0.9))
grid.text('a',x=unit(0.02,'npc'),y=unit(font.y-0.1,'npc'),gp=gpar(cex=font.cex,fontface="bold"))
upViewport()
pushViewport(viewport(layout.pos.row=1,layout.pos.col=2,name="useful_noise"))
useful.noise <- plot.sim(sims,c(6,NA,0.2,NA),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric"), filename="scoring_metric_data/sim_plots/useful_noise.png",param.interest='epsilon2',print.legend=FALSE)
print(useful.noise,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.05,0,0.95,0.9))
grid.text('b',x=unit(0.02,'npc'),y=unit(font.y-0.1,'npc'),gp=gpar(cex=font.cex,fontface="bold"))
upViewport()
pushViewport(viewport(layout.pos.row=2,layout.pos.col=2,name="useful"))
prop.heatmap <- plot.properties.heatmap(p1,p2,p3,p4,"scoring_metric_data/sim_plots/metric_properties_heatmap.png")
print(prop.heatmap,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.01,0.1,0.89,1))
grid.text('d',x=unit(0.02,'npc'),y=unit(font.y,'npc'),gp=gpar(cex=font.cex,fontface="bold"))
upViewport()
pushViewport(viewport(layout.pos.row=2,layout.pos.col=1,name="heatmap"))
useful <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric"), filename="scoring_metric_data/sim_plots/useful.png",param.interest="K.u",print.legend=FALSE)
print(useful,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.1,0.1,1,1))
grid.text('c',x=unit(0.02,'npc'),y=unit(font.y,'npc'),gp=gpar(cex=font.cex,fontface="bold"))
dev.off()



p1 <- test.p1( C=3:6,K.u.values=2:9,K.n.values=3,e1.values=0.4,e2.values=0.4,metrics=c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"))
p2 <- test.p2( C=2:6,K.u.values=2:9,K.n.values=3,e1.values=0.4,e2.values=0.4,metrics=c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"))
p3 <- test.p3( C=4,K.u.values=2:9,K.n.values=3,e1.values=seq(0,0.4,by=0.05),e2.values=0.4,metrics=c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"))
p4 <- test.p4( C=4,K.u.values=6,K.n.values=c(1,2,3),e1.values=0.4,e2.values=seq(0,0.4,0.05),metrics=c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"))

useful.prop <- plot.sim(sims,c(6,NA,0.2,NA),limits=TRUE,metrics =c("Pearson","Pseudo.V","MCC","AUPR","X2A_metric","SCC","SSE"), lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"), filename="scoring_metric_data/sim_plots/useful_prop_all.pdf",param.interest="epsilon2",key.x=5,key.y=-0.5)

prop.heatmap <- plot.properties.heatmap(p1,p2,p3,p4,col='skyblue1',textcol='grey39',yaxis.lab= c("JS divergence","Pearson/MCC","AUPR","2A score", "SCC","SSE"),"scoring_metric_data/sim_plots/all_metric_properties_heatmap.png")


cairo_pdf("scoring_metric_data/sim_plots/all_metrics_all_sim3.pdf",height=6.3,width=7)
# CairoPDF("scoring_metric_data/sim_plots/all_metrics_all_sim3.pdf",height=7,width=7.1)
grid.newpage()
lay1 <- grid.layout(2,2,heights=unit(c(3.5,3.5),c("in","in")),widths=unit(c(3.5,3.5),c("in","in")))
vp1 <- viewport(layout=lay1,name="top.vp",height=unit(7,"in"),width=unit(7,"in"))
pushViewport(vp1)
pushViewport(viewport(layout.pos.row=1,layout.pos.col=1,name="useful_prop"))
# dev.off()
useful.prop <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"), limits=TRUE, lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"),filename=NULL,param.interest="epsilon1",print.legend=FALSE)
print(useful.prop,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.05,0,0.95,0.9))
grid.text("\u03b5",x=unit(0.535,'npc'),y=unit(0.1,'npc'),gp=gpar(cex=font.cex,fontface="bold", fontfamily="DejaVu Sans"))
grid.text('a',x=unit(0.02,'npc'),y=unit(font.y-0.1,'npc'),gp=gpar(cex=font.cex,fontface="bold", fontfamily="Arial"))
upViewport()
pushViewport(viewport(layout.pos.row=1,layout.pos.col=2,name="useful_noise"))
useful.noise <- plot.sim(sims,c(6,NA,0.2,NA),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"),limits=TRUE, lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"), filename=NULL,param.interest='epsilon2',print.legend=FALSE)
print(useful.noise,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.05,0,0.95,0.9))
grid.text("\u03b5",x=unit(0.535,'npc'),y=unit(0.1,'npc'),gp=gpar(cex=font.cex,fontface="bold", fontfamily="DejaVu Sans"))
grid.text('b',x=unit(0.02,'npc'),y=unit(font.y-0.1,'npc'),gp=gpar(cex=font.cex,fontface="bold", fontfamily="Arial"))
upViewport()
pushViewport(viewport(layout.pos.row=2,layout.pos.col=2,name="useful"))
prop.heatmap <- plot.properties.heatmap(p1,p2,p3,p4,col='skyblue1',textcol='grey39',yaxis.lab= c("JS divergence","Pearson/MCC","AUPR","2A score", "SCC","SSE"),filename=NULL)
print(prop.heatmap,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.01,0.1,0.89,1))
grid.text('d',x=unit(0.02,'npc'),y=unit(font.y,'npc'),gp=gpar(cex=font.cex,fontface="bold", fontfamily="Arial"))
upViewport()
pushViewport(viewport(layout.pos.row=2,layout.pos.col=1,name="heatmap"))
useful <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"), limits=TRUE, lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"), filename=NULL, param.interest="K.u",print.legend=TRUE, key.y=0.4, key.x=0.5)
print(useful,newpage=FALSE,panel.height=list(c(panel.size),c("inches")),position=c(0.05,0.1,0.95,1))
grid.text('c',x=unit(0.02,'npc'),y=unit(font.y,'npc'),gp=gpar(cex=font.cex,fontface="bold", fontfamily="Arial"))
dev.off()


useful.prop <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"), limits=TRUE, lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"),filename="~/SMCHET_R1/SuppFig1_a.tiff",param.interest="epsilon1",print.legend=FALSE)
useful.noise <- plot.sim(sims,c(6,NA,0.2,NA),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"),limits=TRUE, lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"), filename="~/SMCHET_R1/SuppFig1_b.tiff",param.interest='epsilon2',print.legend=FALSE)
prop.heatmap <- plot.properties.heatmap(p1,p2,p3,p4,col='skyblue1',textcol='grey39',yaxis.lab= c("JS divergence","Pearson/MCC","AUPR","2A score", "SCC","SSE"),filename="~/SMCHET_R1/SuppFig1_c.tiff")
useful <- plot.sim(sims,c(NA,0,NA,0),metrics =c("Pseudo.V","MCC","AUPR","X2A_metric", "SCC","SSE"), limits=TRUE, lab=c("AUPR","Pearson/MCC","JS divergence","SCC", "SSE","2A score"), filename="~/SMCHET_R1/SuppFig1_d.tiff", param.interest="K.u",print.legend=TRUE, key.y=0.4, key.x=0.5)
