source('~/Shad_scoring/smc_het_eval/metric_test_sim.R')
library(pryr)
library(getopt)
library(doParallel)
registerDoParallel(cores=8)


cli.params <- c(
               'C', 'c', 0, 'numeric',
			   'Ku', 'u', 0, 'numeric',
			   'Kn', 'n', 0, 'numeric',
			   'e1', 'e', 0, 'numeric',
               'e2', 'f', 0, 'numeric'
                );

params = matrix(cli.params, ncol = 4, byrow = TRUE);
opt <- getopt(params)
# sims <- run.sim.par(n.iter =10,
#         C.values = opt$C,
#         K.u.values = opt$Ku,
#         K.n.values = opt$Kn,
#         e1.values=opt$e1,#c(0.02,0.06,0.08,0.1),
#         e2.values= c(0.066,0.133,0.2))
# )
# print("mem used before calling")
# print(mem_used())
sims <- run.sim.single(n.iter =100,
        C = opt$C,
        K.u = opt$Ku,
        K.n = opt$Kn,
        e1=opt$e1,#c(0.02,0.06,0.08,0.1),
        e2= opt$e2#c(0.066,0.133,0.2))
)
