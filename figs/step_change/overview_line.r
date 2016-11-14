the.data <- read.table("summary_three_x.csv",
                        sep=";",
                        header=T)

the.data <- the.data[the.data$rate_t0 < pi,]

freq.u <- sort(unique(the.data$rate_t0))



if (!exists("comb.table"))
{
    # now make table of combinations
    comb.table <- as.data.frame(
            table(
                    the.data[,c("rate_t0","mu_b","mu_m_m","mu_m_e","mu_m_g")]
                    )
            )

    # remove factor thingies
    for (name_i in names(comb.table)) {
    comb.table[,name_i] <- as.numeric(
            as.character(comb.table[,name_i]))
    }

    for (i in 1:nrow(comb.table))
    {
        print(i)

        rate_t0 <- comb.table[i,"rate_t0"]
        mu_b <- comb.table[i,"mu_b"]
        mu_m_m <- comb.table[i,"mu_m_m"]
        mu_m_e <- comb.table[i,"mu_m_e"]
        mu_m_g <- comb.table[i,"mu_m_g"]


        max_delta <- max(
                the.data[
                    the.data$rate_t0 == rate_t0
                    & the.data$mu_b == mu_b 
                    & the.data$mu_m_m == mu_m_m 
                    & the.data$mu_m_e == mu_m_e 
                    & the.data$mu_m_g == mu_m_g
                    & the.data$generation == 50000
                    ,"intptb"])

        comb.table[i,"max_delta"] <- max_delta
    }
}



pdf("max_delta_line.pdf")
print(xyplot(max_delta ~ rate_t0 |
                    mu_m_g * mu_m_e * mu_m_m * mu_b,
                data=comb.table,
                strip=function(strip.levels,...) {
                    strip.default(strip.levels=T,...)
                },
                type="l"
                )
        )
dev.off()

