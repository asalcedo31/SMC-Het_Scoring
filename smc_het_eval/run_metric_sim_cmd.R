source('metric_test_sim.R')
library(getopt)

cli.params <- c(
               'C', 'c', 0, 'numeric',
			   'Ku', 'u', 0, 'numeric',
			   'Kn', 'n', 0, 'numeric',
			   'e1', 'e', 0, 'numeric'
                );

params = matrix(cli.params, ncol = 4, byrow = TRUE);
opt <- getopt(params)
sims <- run.sim(n.iter =25,
        C.values = opt$C,
        K.u.values = opt$Ku,
        K.n.values = opt$Kn,
        e1.values=opt$e1,#c(0.02,0.06,0.08,0.1),
        e2.values= seq(0,0.4,0.05)#c(0.066,0.133,0.2))
)
