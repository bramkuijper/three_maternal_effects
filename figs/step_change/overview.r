# plot how an environmental perturbation works for
# different maternal effects

the.data <- read.table("summary_three.csv"
                        ,sep=";"
                        ,header=T)


pdf("overview_perturb.pdf")
print(levelplot(generation ~ rate_t0 * intptb |
                    mu_m_g * mu_m_e * mu_m_m * mu_b,
                    data=the.data,
                    strip=function(strip.levels,...) {
                        strip.default(strip.levels=T,...)
                    },
                    col.regions=matlab.like))
dev.off()
