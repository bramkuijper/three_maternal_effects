# read in the data
the.data <- read.table("summary_three_types_new_mge.csv",sep=";",header=T)

# plot the m_ijs
pdf(file="overview_three_m.pdf")
print(xyplot(
                meanm_m + meanm_g + meanm_e ~
                    rate_t0 | omega2 * mu_m_g * mu_m_e * mu_m_m * mu_b * tau
                ,data=the.data,
                ,strip=function(strip.levels,...) 
                        {
                            strip.default(strip.levels=T,...)
                        }
                ,auto.key=T
                )
        )
dev.off()

