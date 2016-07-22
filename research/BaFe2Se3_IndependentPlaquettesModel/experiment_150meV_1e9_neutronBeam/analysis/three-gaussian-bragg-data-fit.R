# R script for fitting three Gaussian curves to a set of data (in this case, three Bragg peaks)
# inside the R scripting environment, use the following command:  source("three-gaussian-bragg-data-fit.R")
# data is read from "bragg_firstExperiment_revised.csv"

# load library for plotting
library(ggplot2)

# load csv file containing data points of this order:  "Grazing Angle","scaled_pk-pk","scaled_sigma_pk-pk"
#bragg <- read.csv(file="bragg_firstExperiment_revised.csv",sep=",")
bragg <- read.csv(file="integrated_intensity.csv",sep=",")
# NOTE: R will change the columns to "Grazing.Angle", "scaled_pk.pk", "scaled_sigma_pk.pk" (apparently hyphens not allowed)

#debugger code
#print(bragg)
#readline("Press <Enter> to continue")

#create R lists (or vectors, apparently same thing) from dataframe
energy <- bragg[['energy']]
I <- bragg[['intensity']]
#errors <- bragg[['scaled_sigma_pk.pk']]

# create new dataframe with new labels
df <- data.frame(energy, I)

# use nls function to create fit
# nls() can be fussy about choosing "reasonable" starting guesses for the parameters,
# producing convergence errors.  The default algorithm Gauss-Newton can be even fussier
#trigauss <- nls( signal ~ (C1 * exp(-(angle-mean1)**2/(2 * sigma1**2)) +
#                           C2 * exp(-(angle-mean2)**2/(2 * sigma2**2)) +
#                           C3 * exp(-(angle-mean3)**2/(2 * sigma3**2)) ),
#                           data=df,
#                           start=list(C1=36000, mean1=0, sigma1=4,
#                                      C2=3700, mean2=24, sigma2=2,
#                                      C3=5010, mean3=29, sigma3=2), algorithm="port",
#                           weights=errors
#)

trigauss <- nls( I ~ (C1 * exp(-(energy-mean1)**2/(2 * sigma1**2)) +
                           C2 * exp(-(energy-mean2)**2/(2 * sigma2**2)) +
                           C3 * exp(-(energy-mean3)**2/(2 * sigma3**2)) ),
                           data=df,
                           start=list(C1=0.06, mean1=0, sigma1=1,
                                      C2=0.078, mean2=88, sigma2=1,
                                      C3=0.071, mean3=108, sigma3=1), algorithm="port"#,
                           #weights=errors
)

# Predict the fitted model to a denser grid of angle (x coordinate) values
#df_fit <- data.frame(angle=seq(0,54,0.1)) # for 100 planes (first set)
#df_fit <- data.frame(energy=seq(18,38,0.1)) # for 110 planes (second set)
#df_fit$I <- predict(trigauss, newdata=df_fit)
df_fit <- data.frame(energy=seq(-20,200,0.1))
df_fit$I <- predict(trigauss, newdata=df_fit)
#df <- predict(trigauss, newdata=df)

# Plot the data with the model superimposed
print(ggplot(df, aes(x=energy, y=I)) + geom_point() + geom_smooth(data=df_fit, stat="identity", color="red", size=1.5) +
      xlab("Energy (meV)") + ylab("Intensity (arb)") + ggtitle("Scattered Intensity vs Energy"))

# save plot to a file
#ggsave("100-planes-triple-gauss-fit-R.pdf", plot=last_plot())
ggsave("150meV-flatDispersions-triple-gauss-fit-R.pdf", plot=last_plot())

# disable scientific notation
options(scipen=999)

# print model coefficients that were calculated for the trigauss fit object
print(coef(trigauss))
