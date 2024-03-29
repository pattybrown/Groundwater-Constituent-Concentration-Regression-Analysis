### Build predictive model to project when constituents of concern will fall below groundwater protection
### standards using groundwater data collected since April 2021. 

require(dplyr)
require(readxl)
require(moderndive)
require(ggplot2)
require(chemCal)
require(mbdm)
require(RobustLinearReg)

# Import and view data:
data <- structure(list(monitoring_well = c("MW", "MW", "MW", 
                                                "MW", "MW", "MW", "MW", "MW", "MW", "MW", 
                                                "MW", "MW"), date = structure(c(1617807900, 1622638800, 
                                                                                        1628069100, 1632839400, 1652187300, 1654530900, 1659519000, 1664985300, 
                                                                                        1681563300, 1686230700, 1690973700, 1696265100), class = c("POSIXct", 
                                                                                                                                                   "POSIXt"), tzone = "UTC"), antimony = c(0.0116, 0.0131, 0.0136, 
                                                                                                                                                                                           0.0124, 0.0107, 0.0112, 0.011, 0.0107, 0.0115, 0.00971, 0.00938, 
                                                                                                                                                                                           0.00882), arsenic = c(0.125, 0.142, 0.161, 0.172, 0.193, 0.186, 
                                                                                                                                                                                                                 0.183, 0.173, 0.17, 0.165, 0.139, 0.132), fluoride = c(7.94, 
                                                                                                                                                                                                                                                                        8.75, 8.99, 9.13, 7.28, 5.62, 7.22, 7.31, 5.58, 7.86, 6.94, 6.96
                                                                                                                                                                                                                 ), pH = c(11.03, 10.7, 10.77, 11.1, 10.77, 10.7, 10.95, 10.74, 
                                                                                                                                                                                                                           10.9, 10.6, 10.46, 11.7), sample = c(1, 2, 3, 4, 5, 6, 7, 8, 
                                                                                                                                                                                                                                                                9, 10, 11, 12)), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                         -12L))
View(data)

# Build linear model for antimony: 
ols_antimony <- lm(antimony ~ sample, data)
coef(ols_antimony)
summary(ols_antimony)
olsregtable <- get_regression_table(ols_antimony)
get_regression_points(ols_antimony)

# Plot antimony model residuals:
antimony_residuals <- resid(ols_antimony)
antimony_residuals <- resid(ols_antimony)
qqnorm(antimony_residuals)
qqline(antimony_residuals)
plot(density(antimony_residuals))

# Shaprio-Wilk test for antimony model residual normality:
shapiro.test(rstandard(ols_antimony))
    # P-value = 0.9086: null hypothesis of normality accepted.

# Return prediction for y = 0.006 (antimony GWPS):
inverse.predict(ols_antimony, 0.006)
    # 21.82: Antimony will fall below the GWPS on this sample event. 


# Build linear model for fluoride:
ols_fluoride <- lm(fluoride ~ sample, data)
coef(ols_fluoride)
summary(ols_fluoride)
olsregtable <- get_regression_table(ols_fluoride)
get_regression_points(ols_fluoride)

# Plot fluoride model residuals:
fluoride_residuals <- resid(ols_fluoride)
fluoride_residuals <- resid(ols_fluoride)
qqnorm(fluoride_residuals)
qqline(fluoride_residuals)
plot(density(fluoride_residuals))

# Shaprio-Wilk test for fluoride model residual normality:
shapiro.test(rstandard(ols_fluoride))
    # P-value = 0.4548: null hypothesis of normality accepted.

# Return prediction for y = 4 (fluoride GWPS):
inverse.predict(ols_fluoride, 4)
    # 25.92: Fluoride will fall below the GWPS on this sample event.

# Plot OLS regression line for antimony: 
ggplot(data, aes(sample, antimony)) + geom_point() +
  geom_abline(aes(slope = -0.0003356294, intercept = 0.0133240909, color = "OLS Regression Line")) +
  ylim(0.005, 0.014) +
  scale_x_continuous(breaks=seq(0,24,by=2), limits = c(0, 24)) +
  geom_hline(aes(yintercept=0.006, color = "GWPS (0.006 mg/L)")) +
  ggtitle("Antimony OLS Trend at MW") + 
  ylab("Antimony (mg/L)") +
  xlab("Sample Event") +
  scale_color_manual(name = element_blank(), values=c("red", "blue")) +
  theme(legend.position = "bottom", plot.title = element_text(size = 13, hjust = 0.5)) +
  geom_point(aes(21.82, 0.006)) +
  geom_text(aes(21.82, 0.006), label = "21.82", size = 3, hjust = -0.2, vjust = -1) 

# Plot OLS regression line for fluoride:
ggplot(data, aes(sample, fluoride)) + geom_point() +
  geom_abline(aes(slope = -0.1784615, intercept = 8.6250000, color = "OLS Regression Line")) +
  ylim(3.5, 10) +
  scale_x_continuous(breaks=seq(0,27,by=3), limits = c(0, 27)) +
  geom_hline(aes(yintercept=4, color = "GWPS (4 mg/L)")) +
  ggtitle("Fluoride OLS Trend at MW") + 
  ylab("Fluoride (mg/L)") +
  xlab("Sample Event") +
  scale_color_manual(name = element_blank(), values=c("red", "blue")) +
  theme(legend.position = "bottom", plot.title = element_text(size = 13, hjust = 0.5)) +
  geom_point(aes(25.92, 4)) +
  geom_text(aes(25.92, 4), label = "25.92", size = 3, hjust = -0.2, vjust = -1)

