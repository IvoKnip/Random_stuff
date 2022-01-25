######################################################
####          Flood frequency analysis            ####
######################################################

# Flood frequency analysis from extreme value distribution
# with L-Moments parameter fitting.

# Helpful literature: 
# Parkes, 2015 Chapter 6 (https://kclpure.kcl.ac.uk/portal/files/51218432/2015_Parkes_Brandon_0838465_ethesis.pdf) 
# Hoskin & Wallis, 1997 (https://www.cambridge.org/core/books/regional-frequency-analysis/8C59835F9361705DAAE1ADFDEA7ECD30)

# Part 1: Selection of the best distribution

# The distribution that describes the relationship between
# flood magnitude and recurrence interval best is unknown 
# prior to the modelling. The most commonly used distributions
# are fitted and the performance is evaluated

library(reshape2)
library(extremeStat)
library(goftest)
library(lmomco)

# function for fitting distributions and calculating their goodness of fit
frequency_curves <- function(amax) {
  
  # Candidate distributions 
  distributions <- c("wak", "wei", "pe3", "gum", "kap", "gpa", "gev", "glo", "exp") 
  # For available distributions see:
  #dist.list
  #prettydist("gev")
  
  # fit the functions
  dlf <- distLfit(amax, ks=TRUE, selection = distributions)
  
  # Conduct Anderson-Darling test 
  kap <- goftest::ad.test(amax, "cdfkap", parkap(lmoms(amax))) # Kappa
  wak <- goftest::ad.test(amax, "cdfwak", parwak(lmoms(amax))) # Wakeby
  wei <- goftest::ad.test(amax, "cdfwei", parwei(lmoms(amax))) # Weibull
  pe3 <- goftest::ad.test(amax, "cdfpe3", parpe3(lmoms(amax))) # Pearson Type III
  gum <- goftest::ad.test(amax, "cdfgum", pargum(lmoms(amax))) # Gumbel
  gpa <- goftest::ad.test(amax, "cdfgpa", pargpa(lmoms(amax))) # Generalized pareto
  gev <- goftest::ad.test(amax, "cdfgev", pargev(lmoms(amax))) # GEV
  glo <- goftest::ad.test(amax, "cdfglo", parglo(lmoms(amax))) # GLO
  exp <- goftest::ad.test(amax, "cdfexp", parexp(lmoms(amax))) # Exponential
  
  ## store goodness of fit statistics
  
  # Kolmogorov-Smirnov
  df_stats <- dlf$gof
  df_stats$distribution <- rownames(df_stats)
  
  # Anderson-Darling
  df_stats$an <- NA
  df_stats$an_p <- NA
  df_stats$an[df_stats$distribution == "pe3"] <- pe3$statistic
  df_stats$an_p[df_stats$distribution == "pe3"] <- pe3$p.value
  df_stats$an[df_stats$distribution == "wak"] <- wak$statistic
  df_stats$an_p[df_stats$distribution == "wak"] <- wak$p.value
  df_stats$an[df_stats$distribution == "kap"] <- kap$statistic
  df_stats$an_p[df_stats$distribution == "kap"] <- kap$p.value
  df_stats$an[df_stats$distribution == "wei"] <- wei$statistic
  df_stats$an_p[df_stats$distribution == "wei"] <- wei$p.value
  df_stats$an[df_stats$distribution == "gev"] <- gev$statistic
  df_stats$an_p[df_stats$distribution == "gev"] <- gev$p.value
  df_stats$an[df_stats$distribution == "gpa"] <- gpa$statistic
  df_stats$an_p[df_stats$distribution == "gpa"] <- gpa$p.value
  df_stats$an[df_stats$distribution == "gum"] <- gum$statistic
  df_stats$an_p[df_stats$distribution == "gum"] <- gum$p.value
  df_stats$an[df_stats$distribution == "glo"] <- glo$statistic
  df_stats$an_p[df_stats$distribution == "glo"] <- glo$p.value
  df_stats$an[df_stats$distribution == "exp"] <- exp$statistic
  df_stats$an_p[df_stats$distribution == "exp"] <- exp$p.value
  
  out <- list("dlf" = dlf, "stats" = df_stats)
  
  return(out)
}

# Prepare dummy data
amax <- runif(40, 7000, 15000)
years <- seq(2020, 2059, 1)
plot(years, amax, type = "b")

# fit curve and calculate goodness of fit
out <- frequency_curves(amax = amax)
stats <- out$stats
dlf <- out$dlf

# Plot the curves
plotLfit(distLfit(amax), cdf=TRUE, nbest = nrow(stats))


# Part 2: Calculation of the return periods for the best distribution

# Define function for the calculation of return periods
# In this example the Wakeby (wak) distribution is fitted.
# The distribution can be adjusted in the function.
# For this, the functions parwak and cdfwak need to be replaced 
# by their, e.g. pargev and cdfgev for generalised extreme value distribution
fun_return_periods <- function(amax, return_period) {
  
  # parameter of the CDF (best performing CDF - Wakeby)
  param <- parwak(lmoms(amax))
  # for other distributions see:
  #dist.list
  #prettydist("wak") # -> parwak()
  
  # generate synthetic flow period from amax
  low_bnd <- min(amax) / 2
  high_bnd <- max(amax) * 2
  flows <- seq(low_bnd, high_bnd, (high_bnd - low_bnd) / 1999)
  
  # compute non-exceedance probability for synthetic flows
  p_non <- cdfwak(flows, param)
  
  # calculate exceedance probability
  df <- data.frame("flow" = flows, "p_non" = p_non)
  df$p_exceed <- 1 - df$p_non
  df$rp <- 1/df$p_exceed
  
  # calculate flow for return periods
  df_rp <- data.frame("return_period" = return_period, "flow" = NA)
  for (i in 1:length(return_period)) {
    df_rp$flow[i] <- round(df$flow[which(abs(df$rp-return_period[i])==min(abs(df$rp-return_period[i])))], 0)
  }
  
  return(df_rp)
}

# Define retunr periods
return_periods <- c(1, 2, 5, 10, 20, 50, 100)
df_rp <- fun_return_periods(amax = amax, return_period = return_periods)

# Extrapolation beyond length of record (in this case 40 years) is uncertain and needs to be done carefully.
# Some publications recommend a max return period of half of the length of the data!
