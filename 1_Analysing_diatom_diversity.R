#########################################################################
## R code for data analyses in: Long-term trends in diatom diversity and palaeoproductivity: a 16 000-year multidecadal study from Lake Baikal, southern Siberia (2022). Clim. Past, 18, 1-18, https:/doi.org/10.5194/cp-18-1-2022 
## Authors: Anson W. Mackay, Vivian A. Felde, David W. Morley, Natalia Piotrowska, Patrick Rioual, Alistair W. R. Seddon, and George E. A. Swann
#########################################################################

#########################################################################
## 1. Prepare R workspace and data
#########################################################################

# 1.1 Load libraries 
library(tidyverse)
library(vegan)
library(ggvegan)
library(mgcv)
library(mvpart)
library(GGally)
library(ggrepel)

# 1.2 Source R script with functions not in libraries
source("functions/hill_diversity_estimate.R")
source("functions/renumgr.R")

# 1.3 Import diatom and palaoproductivity data (DOI:)
plankton <- read.csv("data/raw_plankton.csv")
productivity <- read.csv("data/productivity.csv")

# 1.4 Create regional climatic periods using the chronological ages
period <- cut(plankton$AgeBac20, breaks = c(-64, 4200, 8200, 11700, 12900, 14700, 16000), labels = c("Late Holocene", "Middle Holocene", "Early Holocene", "Younger Dryas", "B-A", "pre B-A"))

# 1.5 Extract taxa only
plankton_taxa <- plankton %>% 
  dplyr::select(Aulacoseira.baicalensis:Diatoma.vulgaris)

#########################################################################
## 2. Diatom analysis: Estimate planktonic composition using DCA and CCA 
#########################################################################

# 2.1 Data transformation
plankton.percent <- plankton_taxa/rowSums(plankton_taxa)*100
plankton.hel <- decostand(plankton_taxa, method = "hellinger") 
plankton.norm <- decostand(plankton_taxa, method = "normalize")


# 2.2 Get DCA scores
plankton.dca <- decorana(plankton.hel)
plankton.dca

plankton.dca.sc <- scores(plankton.dca, choices = 1:2)

# 2.3 CCA with time as constraining variable
plankton.cca <- cca(plankton.hel~plankton$AgeBac20)
plankton.cca


plankton.cca.fort <- fortify(plankton.cca)

autoplot(plankton.cca) + 
  theme_bw() +
  geom_text_repel(data = plankton.cca.fort %>% 
                    filter(Score %in% "species"), 
                  aes(x = CCA1, y = CA1, label = Label), size = 3) +
  geom_text_repel(data = plankton.cca.fort %>% 
                    filter(Score %in% "sites"), 
                  aes(x = CCA1, y = CA1, label = Label), size = 3)


plankton.cca.sc <- scores(plankton.cca, display = "sites", scaling = "sites", choices = 1:2)

#########################################################################
## 3. Run Multivariate classification trees (MCT) 
#########################################################################

# 3.1 MCT constrained by the regional climatic periods
set.seed(100)

mvpart.period <- mvpart(data.matrix(plankton.norm)~period, 
                        xv = "1se", 
                        xval = nrow(plankton.norm),
                        xvmult = 100, 
                        data = plankton.norm)

summary(mvpart.period)
mvpart.gr.period <- renumgr(mvpart.period$where)

# 3.2 Save table summary results
write.csv(mvpart.period$cptable, "cptable.periods.csv")
write.csv(mvpart.period$splits, "mct.splits.period.csv")

# 3.3 Extract taxa with the highest contribution for each split/group
taxa_gr4 <- apply(plankton.norm[mvpart.gr.period ==4, colSums(plankton.norm[mvpart.gr.period==4, ])>0 ],2, sum)

sort(taxa_gr4[taxa_gr4 > 4], decreasing = TRUE)

taxa_gr3 <- apply(plankton.norm[mvpart.gr.period==3, colSums(plankton.norm[mvpart.gr.period==3,])>0 ],2, sum)

sort(taxa_gr3[taxa_gr3 > 8], decreasing = TRUE)

taxa_gr2 <- apply(plankton.norm[mvpart.gr.period==2, colSums(plankton.norm[mvpart.gr.period==2,])>0 ],2, sum)

sort(taxa_gr2[taxa_gr2 > 10], decreasing = TRUE)

taxa_gr1 <- apply(plankton.norm[mvpart.gr.period==1, colSums(plankton.norm[mvpart.gr.period==1,])>0 ],2, sum)

sort(taxa_gr1[taxa_gr1 > 5], decreasing = TRUE)


#########################################################################
## 4. Estimate diatom richness, diversity, and evenness
#########################################################################

min(rowSums(plankton_taxa))
plankton_rarediv_estimate <- hill_diversity_estimate(plankton_taxa) # smallest sample size = 150

#########################################################################
## 5. Combine data
#########################################################################

# 5.1 Combine plankton diversity and palaeoproductivity, and recalculate diversity divided on sedimentation rate
plankton_diversity <- 
  data.frame(depth = plankton$depth..cm, 
                                 age = plankton$AgeBac20, 
                                 period = period,
                                 accumRate = plankton$accumRate, 
                                 counts = rowSums(plankton_taxa),
                                 mct.period = mvpart.gr.period,
                                 plankton_rarediv_estimate,
                                 plankton.dca.sc) %>%
  left_join(productivity, 
            by = c("age" = "AgeBac20")) %>%
  mutate(N0 = N0/accumRate,
         N1 = N1/accumRate,
         N2 = N2/accumRate,
         N2.N0 = N2/N0,
         N2.N1 = N2/N1,
         TotValveBVAR_log = log(TotValveBVAR))

# 5.2 Save results
write.csv(plankton_diversity, "plankton_diversity.csv")

# Note: Stratigraphical figures are produced using software C2 Data Analysis Version 1.7.7 (Juggins, 2014).

#########################################################################
## 6. Palaeoproductivity-diversity relationships
#########################################################################

# 6.1 Inspect variables
plankton_diversity %>%
  dplyr::select(age, period, N0, N1, N2, TotValveBVAR, TotValveBVAR_log) %>%
  ggpairs()

# 6.2 Filter out outliers in palaeoproductivity
plankton_diversity2 <- plankton_diversity %>% 
  filter(TotValveBVAR_log > 6) 

# 6.3 Explore GAM models of the variables N2 and TotValveBVAR_log using factor smooth interaction
model1 <- gam(N2 ~ period + s(TotValveBVAR_log, by = period, bs = "tp"), method = "REML", family = Gamma(link = "log"), data = plankton_diversity2)

summary(model1)
gam.check(model1)
plot(model1, pages = 1, shade = TRUE, scale = 0)
acf(model1$residuals)
pacf(model1$residuals)


# simple autocorrelation structure
# valRho <- acf(resid(model1), plot=FALSE)$acf[2]
# 
# plankton_diversity2 <- plankton_diversity2 %>% mutate(start.AR = 1:nrow(plankton_diversity2) == 1)
# 
# model2 <- gam(N2 ~ period + s(TotValveBVAR_log, by = period, bs = "tp"), 
#           rho = valRho, AR.start = start.AR, 
#           method = "REML", family = Gamma(link = "log"), data = plankton_diversity2)
# 
# summary(model2)
# plot(model2, pages = 1, shade = TRUE, scale = 0)
# acf(model2$residuals)
# pacf(model2$residuals)
# gam.check(model2)

# try different temporal autocorrelation structures
# model3 <- gamm(N2 ~ period + s(TotValveBVAR_log, by = period, bs = "tp"),
#           correlation = corAR1(form = ~ age),
#           method = "REML", family = Gamma(link = "log"),  data = plankton_diversity2)
# 
# summary(model3$gam)
# plot(model3$gam, pages = 1, shade = TRUE, scale = 0)
# gam.check(model3$gam)
# acf(model3$gam$residuals)
# pacf(model3$gam$residuals)
# acf(resid(model3$lme, type = "normalized"))
# pacf(resid(model3$lme, type = "normalized"))
# 
# model4 <- gamm(N2 ~ period + s(TotValveBVAR_log, by = period, bs = "tp"),
#            correlation = corARMA(form = ~ age, p = 2, q = 0),
#            method = "REML", family = Gamma(link = "log"),  data = plankton_diversity2)
# 
# acf(model4$gam$residuals)
# pacf(model4$gam$residuals)
# acf(resid(model4$lme, type = "normalized"))
# pacf(resid(model4$lme, type = "normalized"))
# 
# 
# model5 <- gamm(N2 ~ period + s(TotValveBVAR_log, by = period, bs = "tp"),
#            correlation = corARMA(form = ~ age, p = 2, q = 1),
#            method = "REML", family = Gamma(link = "log"),  data = plankton_diversity2)
# 
# acf(model5$gam$residuals)
# pacf(model5$gam$residuals)
# acf(resid(model5$lme, type = "normalized"))
# pacf(resid(model5$lme, type = "normalized"))
# 
# model6 <- gamm(N2 ~ period + s(TotValveBVAR_log, by = period, bs = "tp"),
#            correlation = corARMA(form = ~ age, p = 2, q = 2),
#            method = "REML", family = Gamma(link = "log"),  data = plankton_diversity2)
# 
# acf(model6$gam$residuals)
# pacf(model6$gam$residuals)
# acf(resid(model6$lme, type = "normalized"))
# pacf(resid(model6$lme, type = "normalized"))


# 6.4 Get fitted values from the model and produce figure
new_data <- tidyr::expand(plankton_diversity2  %>%
                          rowid_to_column("sample.id"), 
                          nesting(period, TotValveBVAR_log), 
                          sample.id = unique(sample.id))

model1_pred <- bind_cols(new_data, 
                         as.data.frame(
                           predict(model1, 
                                   newdata = new_data, 
                                   type = "response", 
                                   se.fit = TRUE)))

palaeoprod.diversity.relationships.fig <- 
  ggplot(model1_pred %>%
         mutate(lower = fit - 1.96 * se.fit, 
                upper = fit + 1.96 * se.fit),
         aes(x = TotValveBVAR_log, y = fit, color = period)) + 
  geom_point(data = plankton_diversity2, 
             aes(x = TotValveBVAR_log, y = N2), alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = period), alpha = 0.1) +
  geom_line(size = 1) +
  facet_wrap(~period) +
  theme_bw() +
  labs(y = expression("Fitted N2 ("~cm^-2~yr^-1~")"), 
       x = expression("Total Valve Biovolume (log BVAR" ~"\u00b5m\u00b3"~cm^-2~yr^-1~")"))

palaeoprod.diversity.relationships.fig

ggsave("palaeoprod.diversity.relationships.fig.pdf")

