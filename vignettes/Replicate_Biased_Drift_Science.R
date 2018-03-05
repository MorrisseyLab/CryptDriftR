## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, message=FALSE-------------------------------------------
library(CryptDriftR)

## ------------------------------------------------------------------------
data(PulseChaseData)
ls()

## ------------------------------------------------------------------------
# Get time points from colnames
time_interval = data_SI_ahCreER_WT %>% colnames() %>% as.numeric()
# Fit model
fit_out_WT   = fitNeutralDrift(data_SI_ahCreER_WT, time_interval)

## ---- fig.width = 7------------------------------------------------------
pp_convergence = plotsConvergence_Neutral(fit_out_WT)
plot(pp_convergence)

## ---- fig.width = 8.5, fig.height= 7-------------------------------------
pp_posterior = plotPosterior_Neutral(fit_out_WT) 
plot(pp_posterior)

## ------------------------------------------------------------------------
print(fit_out_WT$MAP_est)

## ---- fig.width = 10, fig.height= 4--------------------------------------
plotsNeutralDrift_Fit(fit_out_WT)

## ------------------------------------------------------------------------
time_interval = data_SI_ahCreER_KRAS %>% colnames() %>% as.numeric()
# We use the WT parameters inferred from the WT data
fit_out_mut   = fitBiasedDrift(data_SI_ahCreER_KRAS, time_interval, WT_params = fit_out_WT$MAP_est)

## ---- fig.width = 7------------------------------------------------------
pp_convergence =  plotsConvergence_Biased(fit_out_mut)
plot(pp_convergence)

## ---- fig.width = 6------------------------------------------------------
pp_posterior = plotPosterior_Pr(fit_out_mut) 
plot(pp_posterior)

## ------------------------------------------------------------------------
print(fit_out_mut$Pr_quant)

## ---- fig.width = 10, fig.height= 5--------------------------------------
pp_KRAS = plotsBiasDrift_Fit(fit_out_mut)
plot(pp_KRAS)

## ---- fig.width = 6, fig.height= 6---------------------------------------
time_interval = data_SI_ahCreER_APChet %>% colnames() %>% as.numeric()
fit_out_het   = fitBiasedDrift(data_SI_ahCreER_APChet, time_interval, WT_params = fit_out_WT$MAP_est)
time_interval = data_SI_ahCreER_APChom %>% colnames() %>% as.numeric()
fit_out_hom   = fitBiasedDrift(data_SI_ahCreER_APChom, time_interval, WT_params = fit_out_WT$MAP_est)
pp_het = plotPosterior_Pr(fit_out_het) + ggtitle("Apc-Het")
pp_hom = plotPosterior_Pr(fit_out_hom) + ggtitle("Apc-Hom")
grid.arrange(pp_het, pp_hom)

## ------------------------------------------------------------------------
time_interval    = data_Colon_LGR5CreER_WT %>% colnames() %>% as.numeric() 
fit_out_colon_WT = fitNeutralDrift(data_Colon_LGR5CreER_WT, time_interval)
print(fit_out_colon_WT$MAP_est)

## ---- fig.width = 6------------------------------------------------------
time_interval   = data_Colon_LGR5CreER_p53 %>% colnames() %>% as.numeric()
fit_out   = fitBiasedDrift(data_Colon_LGR5CreER_p53, time_interval, WT_params = fit_out_colon_WT$MAP_est) 
pp_Pr_p53 = plotPosterior_Pr(fit_out) 
plot(pp_Pr_p53)

## ---- fig.width = 8.5, fig.height= 7-------------------------------------
time_interval   = data_DSSColon_LGR5CreER_WT %>% colnames() %>% as.numeric()
fit_out   = fitNeutralDrift(data_DSSColon_LGR5CreER_WT, time_interval) 
pp_post = plotPosterior_Neutral(fit_out)
plot(pp_post)

## ---- fig.width = 6------------------------------------------------------
time_interval = data_DSSColon_LGR5CreER_p53 %>% colnames() %>% as.numeric()
fit_out_mut   = fitBiasedDrift(data_DSSColon_LGR5CreER_p53, time_interval, WT_params = fit_out$MAP_est) 
pp_Pr_p53     = plotPosterior_Pr(fit_out_mut) 
plot(pp_Pr_p53)

