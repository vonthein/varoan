#
# function for multimodel analysis with
#
# numerical scores
# at baseline and as endpoint 
# for sum score estimands and ICC
#
sumscorefit <- function(daten = daten,  # preprocessed to different formats
                        items = 35:39, # vector of integer referring to daten$rs
                        treat = "arm", # character for allocation to two groups
                        endpoint = "GOALS", 
                        id = "Proband", # cluster variable
                        covar, # vector of character of covariable names for adjustment
                        svg = TRUE,    # same plots as vector graphics?
                        Deutsch = FALSE # language switch
    ){
# for tests
 # items = 35:39 # vector of integer referring to daten$ss
 # treat = "arm"; endpoint = "GOALS"; id = "Proband"
 # covar = c("Dauer der Erkrankung (Tage)", "Alter.Patient", "BMI");
 # svg = TRUE; Deutsch = FALSE

    # Transformation of bounded data
  logit <- function(x, inf = 0, sup = 1){
    log((x-inf)/(sup-x))
  }
  # Backtransformation
  expit <- function(x, inf = 0, sup = 1){
    inf + (sup - inf) * exp(x) / (1 + exp( x))
  }
  
  times <- length(unique(daten$ssFAS$time)) 
  raters <- length(unique(daten$ssFAS$rater)) 
  item_n <- length(items)
  
  datenFAS   <- daten$FAS
  datenPP    <- daten$PP
  datensFAS  <- daten$sFAS
  datensPP   <- daten$sPP
  datenssFAS <- daten$ssFAS
  datenssPP  <- daten$ssPP
  datensrFAS <- daten$srFAS
  datensrPP  <- daten$srPP
# unify variable names
  datenFAS$GOALS1 <- datenFAS[, paste0(endpoint, 1)]
  datenFAS$GOALS2 <- datenFAS[, paste0(endpoint, times)]
  datenFAS$GOALSchange <- datenFAS[, paste0(endpoint, "change")]
  datenFAS$logit.GOALS1 <- datenFAS[, paste0("logit.", endpoint, 1)]
  datenFAS$logit.GOALS2 <- datenFAS[, paste0("logit.", endpoint, times)]
  datenFAS$logit.GOALSchange <- datenFAS[, paste0("logit.", endpoint, "change")]
  datenFAS$arm <- datenFAS[, treat]
  datenPP$GOALS1    <- daten$PP[, paste0(endpoint, 1)]
  datenPP$GOALS2    <- daten$PP[, paste0(endpoint, 1)]
  datenPP$GOALSchange    <- daten$PP[, paste0(endpoint, 1)]
  datenPP$logit.GOALS1    <- daten$PP[, paste0("logit.", endpoint, 1)]
  datenPP$logit.GOALS2    <- daten$PP[, paste0("logit.", endpoint, 1)]
  datenPP$logit.GOALSchange    <- daten$PP[, paste0("logit.", endpoint, 1)]
  datenPP$arm <- datenPP[, treat]
  datensFAS$GOALS  <- daten$sFAS[, endpoint]
  datensFAS$logit.GOALS  <- daten$sFAS[, paste0("logit.", endpoint)]
  datensPP$logit.GOALS   <- daten$sPP[, paste0("logit.", endpoint)]
  datensPP$GOALS   <- daten$sPP[, endpoint]
  datenssFAS$GOALS <- daten$ssFAS[, endpoint]
  datenssPP$GOALS  <- daten$ssPP[, endpoint]
  datensrFAS$GOALS <- daten$srFAS[, endpoint]
  datensrPP$GOALS  <- daten$srPP[, endpoint]
  datenssFAS$logit.GOALS <- daten$ssFAS[, paste0("logit.", endpoint)]
  datenssPP$logit.GOALS  <- daten$ssPP[, paste0("logit.", endpoint)]
  # datensrFAS$logit.GOALS <- daten$srFAS[, paste0("logit.", endpoint)]
  # datensrPP$logit.GOALS  <- daten$srPP[, paste0("logit.", endpoint)]
  datenssFAS$Gesamt1 <- daten$ssFAS[, paste0(endpoint, 1)]
  datenssPP$Gesamt1  <- daten$ssPP[, paste0(endpoint, 1)]
  datensrFAS$Gesamt1 <- daten$srFAS[, paste0(endpoint, 1)]
  datensrPP$Gesamt1  <- daten$srPP[, paste0(endpoint, 1)]
  # datenssFAS$GesamtE1 <- daten$ssFAS[, paste0(endpoint, "E1")]
  # datenssPP$GesamtE1  <- daten$ssPP[, paste0(endpoint, "E1")]
  datensrFAS$GesamtE1 <- daten$srFAS[, paste0(endpoint, "E1")]
  datensrPP$GesamtE1  <- daten$srPP[, paste0(endpoint, "E1")]
  # datenssFAS$logit.GesamtE1 <- daten$ssFAS[, paste0("logit.", endpoint, "E1")]
  # datenssPP$logit.GesamtE1  <- daten$ssPP[, paste0("logit.", endpoint, "E1")]
  # datensrFAS$logit.GesamtE1 <- daten$srFAS[, paste0("logit.", endpoint, "E1")]
  # datensrPP$logit.GesamtE1  <- daten$srPP[, paste0("logit.", endpoint, "E1")]
  datenssFAS$GesamtD1 <- daten$ssFAS[, paste0(endpoint, "D1")]
  datenssPP$GesamtD1  <- daten$ssPP[, paste0(endpoint, "D1")]
  datensrFAS$GesamtD1 <- daten$srFAS[, paste0(endpoint, "D1")]
  datensrPP$GesamtD1  <- daten$srPP[, paste0(endpoint, "D1")]
  
  Toolbox <- labels(datenFAS[,treat])[1] 
  Waiting <- labels(datenFAS[,treat])[2]
   
  # FAS
# Raw endpoint
t2 <- t.test(GOALS2 ~ arm, data = datenFAS)
est2 <- c(t2$estimate, t2$estimate %*%c(-1, 1), -t2$conf.int[2:1], t2$p.value)
# Raw change
t2d <- t.test(GOALSchange ~ arm, data = datenFAS)
est2d <- c(t2d$estimate, t2d$estimate %*%c(-1, 1), -t2d$conf.int[2:1], t2d$p.value)
# Baseline adjusted
fit <- lm(GOALS2 ~ GOALS1 + arm, data = datenFAS)
sumfit <- summary(fit)
confit <- confint(fit)
prefit <- predict(fit, 
                  newdata = data.frame(arm = c(Toolbox, Waiting),
                                       GOALS1 = rep(mean(datenFAS$GOALS1, 
                                                         na.rm = TRUE), 
                                                    2)
                  )
)
t2adj <- sumfit$coefficients["arm" == substr(rownames(sumfit$coefficients),1,3),
                             c("Estimate", "Pr(>|t|)")]
est2adj <- c(prefit, t2adj[1], confit["arm" == substr(rownames(confit),1,3),], t2adj[2])

# mixed model
library(lme4)
mfit <- lmer(GOALS ~ time * arm + (1 | Proband), data = datensFAS)
mfit0 <- lmer(GOALS ~ time + arm + (1 | Proband), data = datensFAS)
summfit <- summary(mfit)
conffit <- confint(mfit)
newdata <- data.frame(time = c("Baseline", "Final", NA, "Baseline", "Final"),
                      arm = c(Waiting, Waiting, NA, Toolbox, Toolbox))
# pre <- data.frame(newdata, 
#                   mean = predict(mfit, re.form = NA, newdata = newdata))
pre <- unique(data.frame(rbind(datensFAS[,c("time", "arm")], NA),
                         mean = c(predict(mfit, re.form = NA), NA)))[c(1,3,5,2,4),]
est2m <- c(pre$mean[5], pre$mean[2], 
           summfit$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"],
           conffit[paste0("timeFinal:arm", Waiting),],
           anova(mfit, mfit0)["mfit", "Pr(>Chisq)"])

# adjusted
mfita <- lmer(formula(paste0("GOALS ~ time * arm + (1 | Proband)", 
                paste0(" + `", covar, collapse = "`"), "`")), 
              data = datensFAS)
mfit0a <- lmer(formula(paste0("GOALS ~ time + arm + (1 | Proband)", 
                              paste0(" + `", covar, collapse = "`"), "`")), 
               data = datensFAS)
summfita <- summary(mfita)
conffita <- confint(mfita)
means <- apply(datensFAS[, covar], 2, mean, na.rm = TRUE)
mean.data <- data.frame(time = c("Baseline", "Final", "Baseline", "Final"),
                        arm = c(Waiting, Waiting, Toolbox, Toolbox),
                        rbind(means, means, means, means))
names(mean.data) <- c("time", "arm", names(means))
prea <- data.frame(mean.data, 
                   mean = predict(mfita, newdata = mean.data, re.form = NA))
est2ma <- c(prea$mean[4], prea$mean[2], 
            summfita$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"],
            conffita[paste0("timeFinal:arm", Waiting),],
            anova(mfita, mfit0a)["lmfit", "Pr(>Chisq)"])
aest2ma <- cbind(summfita$coefficients[4:(3 + length(covar)), "Estimate"],
                              -conffita[6:(5 + length(covar)), ])


# mixed model with assessor agreement
mfitsr <- lmer(GOALS ~ arm + rater + (1 | Proband), data = datensrFAS)
mfit0sr <- lmer(GesamtE1 ~ arm + rater + (1 | Proband), data = datensrFAS)
mfit0srRA <- lmer(GesamtE1 ~ arm + (1 | rater) + (1 | Proband), data = datensrFAS)
mfit1sr <- lmer(GesamtE1 ~ Gesamt1 + arm + rater + (1 | Proband), data = datensrFAS)
mfit2srRA <- lmer(GesamtE1 ~ Gesamt1 + arm + (1 | rater) + (1 | Proband), data = datensrFAS)

mfit0ss <- lmer(GOALS ~ time + arm + rater + (1 | Proband), data = datenssFAS)
mfitss <- lmer(GOALS ~ time * arm + rater + (1 | Proband), data = datenssFAS)
mfit0ssi <- lmer(GOALS ~ time + arm + rater + (1+time | Proband), data = datenssFAS)
mfitssi <- lmer(GOALS ~ time * arm + rater + (1+time | Proband), data = datenssFAS)
mfit0ssR <- lmer(GOALS ~ time + arm + Proband + (1 | rater), data = datenssFAS)
mfitssR <- lmer(GOALS ~ time * arm + Proband + (1 | rater), data = datenssFAS)
# ICC <- mfitssR@theta
mfit0ssRA <- lmer(GOALS ~ time + arm + (1 | Proband) + (1 | rater), data = datenssFAS)
mfitssRA <- lmer(GOALS ~ time * arm + (1 | Proband) + (1 | rater), data = datenssFAS)
mfit0ssRAi <- lmer(GOALS ~ time + arm + (1+time | Proband) + (1 | rater), data = datenssFAS)
# icc <- mfitssRAi@theta[4]
mfitssRAi <- lmer(GOALS ~ time * arm + (1+time | Proband) + (1 | rater), data = datenssFAS)
mfit0ssF <- lm(GOALS ~ time + arm + rater + Proband, data = datenssFAS)
mfitssF <- lm(GOALS ~ time * arm + rater + Proband, data = datenssFAS)
anova(mfitss, mfit0ssRA, mfitssRA, mfit0ss)
anova(mfit0ss, mfitss, 
      mfit0ssR, mfitssR, 
      #     mfit0ssRA, mfitssRA, 
      mfit0ssRAi, mfitssRAi, 
      mfit0ssi, mfitssi, 
      mfit0ssF, mfitssF)
summfitmm1 <- summary(mfitss)
conffitmm1 <- confint(mfitss)
summfitmm2 <- summary(mfitssi)
conffitmm2 <- confint(mfitssi)
summfitmm3 <- summary(mfitssF)
conffitmm3 <- confint(mfitssF)
summfitmm4 <- summary(mfitssRAi)
conffitmm4 <- confint(mfitssRAi)
newdatar <- data.frame(time = c("Baseline", "Final", NA,"Baseline", "Final"),
                       arm = c(Waiting, Waiting, NA, Toolbox, Toolbox), 
                       rater = as.factor(c(1, 1, NA, 1, 1)))
newdataF <- data.frame(time = c("Baseline", "Final", NA,"Baseline", "Final"),
                       arm = c(Waiting, Waiting, NA, Toolbox, Toolbox), 
                       rater = as.factor(c(1, 1, NA, 1, 1)), 
                       Proband = as.factor(c(2, 2, NA, 2, 2)))
# premm1 <- data.frame(newdatar, 
#                   mean = predict(mfitss, re.form = NA, newdata = newdatar))
premm1 <- unique(data.frame(rbind(datenssFAS[,c("time", "arm", "rater")], NA),
                            mean = c(predict(mfitss, re.form = NA), NA)))[c(1,3,13,2,4,13,
                                                                            5,7,13,6,8,13,
                                                                            9,11,13,10,12),]
# premm2 <- data.frame(newdatar, 
#                   mean = predict(mfitssi, re.form = NA, newdata = newdatar))
premm2 <- unique(data.frame(rbind(datenssFAS[,c("time", "arm", "rater")], NA),
                            mean = c(predict(mfitssi, re.form = NA), NA)))[c(1,3,13,2,4,13,
                                                                             5,7,13,6,8,13,
                                                                             9,11,13,10,12),]
# premm3 <- data.frame(newdataF, 
#                   mean = predict(mfitssF, re.form = NA, newdata = newdataF))
premm3 <- unique(data.frame(rbind(datenssFAS[,c("time", "arm", "rater")], NA),
                            mean = c(predict(mfitssF, re.form = NA), NA)))[c(1,3,13,2,4,13,
                                                                             5,7,13,6,8,13,
                                                                             9,11,13,10,12),]
# premm4 <- data.frame(newdatar, 
#                   mean = predict(mfitssRAi, re.form = NA, newdata = newdatar))
premm4 <- unique(data.frame(rbind(datenssFAS[,c("time", "arm", "rater")], NA),
                            mean = c(predict(mfitssRAi, re.form = NA), NA)))[c(1,3,13,2,4,13,
                                                                               5,7,13,6,8,13,
                                                                               9,11,13,10,12),]
est2mm1 <- c(premm1$mean[5], premm1$mean[2], 
             summfitmm1$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"],
             conffitmm1[paste0("timeFinal:arm", Waiting),],
             anova(mfitss, mfit0ss)["mfitss", "Pr(>Chisq)"])
est2mm2 <- c(premm2$mean[5], premm2$mean[2], 
             summfitmm2$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"],
             conffitmm2[paste0("timeFinal:arm", Waiting),],
             anova(mfitssi, mfit0ssi)["mfitssi", "Pr(>Chisq)"])
est2mm3 <- c(premm3$mean[5], premm3$mean[2], 
             summfitmm3$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"],
             conffitmm3[paste0("timeFinal:arm", Waiting),],
             anova(mfitssF, mfit0ssF)[2, "Pr(>F)"])
est2mm4 <- c(premm4$mean[5], premm4$mean[2], 
             summfitmm4$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"],
             conffitmm4[paste0("timeFinal:arm", Waiting),],
             anova(mfitssRAi, mfit0ssRAi)["mfitssRAi", "Pr(>Chisq)"])
# cat(pander(cbind(est2mm1,est2mm2,est2mm3,est2mm4)))
# ICC another way
library(irr)
ic <- function(vari = "GOALS", daten = datenssFAS, unit = "single"){
  data <- data.frame(a = datenssFAS[1:36, vari],
                     b = datenssFAS[37:72, vari], 
                     c = datenssFAS[73:108, vari])
  ic <- icc(data, model = "twoway", type = "agreement", unit = unit)
  out <- round(data.frame(ICC = ic$value, LCL = ic$lbound, UCL = ic$ubound), 2)
  return(out)
}
vars <- c("Tiefenwahrnehmung1", "Bim.Geschicklichkeit1", "Effizient1", "Umgang mit Gewebe1", "GOALS")
ICC <- t(sapply(vars, ic))
# rownames(ICC) <- vars[-5] 
ICCm <- t(sapply(vars, ic, unit = "average"))
# rownames(ICC) <- dimensions[-5] 
library(pander)
if(Deutsch) {
  colnames(ICC) <- c("ICC", "UKG", "OKG")
  colnames(ICCm) <- c("ICCm", "UKG", "OKG")
  cat(pander(ICC, caption = "Intraklasskorrelation in der Zweifachklassifikation (ICC) mit 3 Befundern und 95% Konfidenzintervall. UKG = untere Konfidenzgrenze, OKG = obere Konfidenzgrenze"))
  cat(pander(ICCm, caption = "Intraklasskorrelation in der Zweifachklassifikation (ICC) des Mittels von 3 Befundern mit 95% Konfidenzintervall. UKG = untere Konfidenzgrenze, OKG = obere Konfidenzgrenze"))
  } else {
  cat(pander(ICC, caption = "Two-way agreement intraclass correlation (ICC) of 3 raters with 95% confidence interval. LCL = lower confidence limit, UCL = upper confidence limit"))
  cat(pander(ICCm, caption = "Two-way agreement intraclass correlation (ICC) of 3 raters when the average of several measurements is regarded with 95% confidence interval. LCL = lower confidence limit, UCL = upper confidence limit"))
  }

# Logit transformed endpoint
# Raw endpoint
lt2 <- t.test(logit.GOALS2 ~ arm, data = datenFAS)
lest2 <- c(expit(lt2$estimate, 5, 25), 
           exp(lt2$estimate %*%c(1, -1)),
           exp(lt2$conf.int), 
           lt2$p.value)
# Raw change
lt2d <- t.test(logit.GOALSchange ~ arm, data = datenFAS)
lest2d <- c(expit(lt2d$estimate, 5, 25), 
            exp(lt2d$estimate %*%c(1, -1)),
            exp(lt2d$conf.int), 
            lt2d$p.value)
# Baseline adjusted
lfit <- lm(logit.GOALS2 ~ logit.GOALS1 + arm, data = datenFAS)
lsumfit <- summary(lfit)
lconfit <- confint(lfit)
lprefit <- predict(lfit, 
                   newdata = data.frame(arm = c(Toolbox, Waiting),
                                        logit.GOALS1 = rep(mean(datenFAS$logit.GOALS1, 
                                                                na.rm = TRUE), 
                                                           2)
                   )
)
lt2adj <- lsumfit$coefficients["arm" == substr(rownames(lsumfit$coefficients),1,3),
                               c("Estimate", "Pr(>|t|)")]
lest2adj <- c(expit(lprefit, 5, 25), exp(-lt2adj[1]),
              exp(-lconfit["arm" == substr(rownames(lconfit),1,3),2:1]),
              lt2adj[2])

# mixed model
lmfit <- lmer(logit.GOALS ~ time * arm + (1 | Proband), data = datensFAS)
lmfit0 <- lmer(logit.GOALS ~ time + arm + (1 | Proband), data = datensFAS)
lsummfit <- summary(lmfit)
lconffit <- confint(lmfit)
# lpre <- data.frame(newdata, 
#                   mean = predict(lmfit, re.form = NA, newdata = newdata))
lpre <- unique(data.frame(rbind(datensFAS[,c("time", "arm")], NA),
                          mean = c(predict(lmfit, re.form = NA), NA)))[c(1,3,5,2,4),]
lest2m <- c(expit(lpre$mean[5], 5, 25), expit(lpre$mean[2], 5, 25), 
            exp(lsummfit$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"]),
            exp(lconffit[paste0("timeFinal:arm", Waiting),]),
            anova(lmfit, lmfit0)["lmfit", "Pr(>Chisq)"])
# adjusted
lmfita <- lmer(logit.GOALS ~ time * arm + (1 | Proband) + `Dauer der Erkrankung (Tage)` + Alter.Patient + BMI, data = datensFAS)
lmfit0a <- lmer(logit.GOALS ~ time + arm + (1 | Proband) + `Dauer der Erkrankung (Tage)` + Alter.Patient + BMI, data = datensFAS)
lsummfita <- summary(lmfita)
lconffita <- confint(lmfita)
medians <- apply(datensFAS[, covar], 2, median, na.rm = TRUE)
med.data <- data.frame(time = c("Baseline", "Final", "Baseline", "Final"),
                       arm = c(Toolbox, Toolbox, Waiting, Waiting),
                       rbind(medians, medians, medians, medians))
names(med.data) <- c("time", "arm", covar)

lprea <- data.frame(med.data, 
                    mean = predict(lmfita, newdata = med.data, re.form = NA))
lest2ma <- c(expit(lprea$mean[4], 5, 25), expit(lprea$mean[2], 5, 25), 
             exp(lsummfita$coefficients[paste0("timeFinal:arm", Waiting), "Estimate"]),
             exp(lconffita[paste0("timeFinal:arm", Waiting),]),
             anova(lmfita, lmfit0a)["lmfit", "Pr(>Chisq)"])
alest2ma <- cbind(exp(lsummfita$coefficients[4:(3+length(covar)), "Estimate"]),
              exp(lconffita[6:(5+length(covar)), ]))

# table
estimands <- rbind(#est2fas, 
                    est2,  est2d,  est2adj,  est2m,  est2ma,
                   lest2, lest2d, lest2adj, lest2m, lest2ma) 
if(Deutsch) {
  colnames(estimands) <- c("Intervention", "Kontrolle", "Differenz oder OR", "95%-confidence", "interval", "P")
  rownames(estimands) <- c(#"Mediane Differenz",
                           "Finale Erwartung", 
                           "Erwartete Änderung", 
                           "Erwartete Änderung adj.f. Baseline", 
                           "Behandlung-Zeit-Wechselwirkung an mittlerer Baseline", 
                           "Adjustierte Behandlung-Zeit-Wechselwirkung an mittlerer Baseline", 
                           "Finaler Median (OR)", 
                           "Mediane Änderung (OR)", 
                           "Mediane Änderung adj.f. Baseline (OR)", 
                           "Behandlung-Zeit-Wechselwirkung an medianer Baseline (OR)", 
                           "Adjustierte Behandlung-Zeit-Wechselwirkung an medianer Baseline (OR)")
  cat(pander(estimands, caption = "Tabelle 3: Estimanden des Behandlungseffekts"))
  } else {
    colnames(estimands) <- c("Treatment", "Control", "Difference or OR", "95%-confidence", "interval", "P")
rownames(estimands) <- c(#"Median difference",
                         "Final mean", 
                         "Mean change", 
                         "Mean change adj.f. baseline", 
                         "Treatment-time interaction at mean baseline", 
                         "Adjusted treatment-time interaction at mean baseline", 
                         "Final median  (OR)", 
                         "Median change (OR)", 
                         "Median change adj.f. baseline (OR)", 
                         "Treatment-time interaction at median baseline (OR)", 
                         "Adjusted treatment-time interaction at median baseline (OR)")
cat(pander(estimands, caption = "Table 3: Estimands for the Toolbox effect"))
} # if not Deutsch

# Plot the data
if(svg) {
    svg(filename = paste0(endpoint, "_FAS.svg"),
        width = 3, height = 3)
}
cols <- c("grey70", "pink")
ot <- par(mar = c(3,2,0,0) + 0.1)
plot(GOALS ~ when, col = cols[datensFAS$Randomisierung], 
     data = datensFAS, 
     xlim = c(0.7, 2.6), xaxp=c(1, 2, 1),
     ylim = c(5, 25), xaxt = "n")
lapply(split(datensFAS,datensFAS$Proband), function(x){
  lines(x$when, x$GOALS, col = cols[x$Randomisierung[1]], lwd = 1)
})
lines(pre$mean[1:2], col = "black", lwd = 3)
lines(pre$mean[4:5], col = "red", lwd = 3)
axis(1, c("Baseline", "Final"), at = 1:2)
text(x=1.5, y=25, endpoint)
text(x=2.3, y=pre$mean[2], Waiting, col = "black")
text(x=2.3, y=pre$mean[5], Toolbox, col = "red")
if(svg) {dev.off()}
# # PP
# if(svg) {
#   svg(filename = paste0(endpoint, "PP.svg"),
#       width = 3, height = 3)
# }
# plot(GOALS ~ when, col = cols[datensPP$Randomisierung], 
#      data = datensPP, 
#      xlim = c(0.7, 2.6), xaxp=c(1, 2, 1),
#      ylim = c(5, 25), xaxt = "n")
# lapply(split(datensPP,datensPP$Proband), function(x){
#   lines(x$when, x$GOALS, col = cols[x$Randomisierung[1]], lwd = 1)
# })
# lines(prePP$mean[1:2], col = "black", lwd = 3)
# lines(prePP$mean[4:5], col = "red", lwd = 3)
# axis(1, c("Baseline", "Final"), at = 1:2)
# text(x=1.5, y=25, endpoint)
# text(x=2.3, y=prePP$mean[2], Waiting, col = "black")
# text(x=2.3, y=prePP$mean[5], Toolbox, col = "red")
par(ot)
# if(svg) {dev.off()}

library(forestplot)
extra <- FALSE
ot <- par(mar=c(4, 1, 0, 0) + 0.1)
# difference
breit <- 6#; hoch <- 3.5
myticks                 <- -2:15
attr(myticks, "labels") <- c("-2", "", "0", "", "2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14", "")
horm <- "ITT"
datentabelle <- estimands[1:5, 3:5]
colnames(datentabelle) <- c("mean", "lower", "upper")
hoch <- 2 + nrow(datentabelle) / 6
if(svg) {
  svg(filename = paste0(endpoint, horm, "forrestplot_Diff.svg"),
      width = breit, height = hoch)
}
# cat(pander(datentabelle),"\n")
forestplot(datentabelle,
           new_page = TRUE,
           is.summary = rep(FALSE,  nrow(datentabelle)),
           xlab = "Difference",
           clip = c(floor(min(datentabelle[,"lower"], na.rm = TRUE)), 
                     ceiling(max(datentabelle[,"upper"], na.rm = TRUE))), 
           boxsize = 0.25, 
           zero = 0,
           xlog = FALSE, 
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
           vertices = TRUE,
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 0.8), 
                            xlab  = gpar(cex = 0.8), 
                            label = gpar(cex = 0.8)),
           xticks = myticks, 
           grid = TRUE)
if(svg) {dev.off()}
# if(extra) {system(paste0("open ", endpoint, horm, "forrestplot_Diff.svg")) }
# 
# horm <- "PP"
# datentabelle <- estimandsPP[1:5, 3:5]
# colnames(datentabelle) <- c("mean", "lower", "upper")
# hoch <- 2+nrow(datentabelle)/6
# if(svg) {
#   svg(filename = paste0(horm, "forrestplot_GoalsDiff.svg"),
#       width = breit, height = hoch)
# }
# forestplot(datentabelle,
#            new_page = TRUE,
#            is.summary = rep(FALSE,  nrow(datentabelle)),
#            xlab = "Difference",
#            clip = c(floor(min(datentabelle[,"lower"], na.rm = TRUE)), 
#                   ceiling(max(datentabelle[,"upper"], na.rm = TRUE))), 
#            boxsize = 0.25, 
#            zero = 0,
#            xlog = FALSE, 
#            col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
#            vertices = TRUE,
#            txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 0.8), 
#                             xlab  = gpar(cex = 0.8), 
#                             label = gpar(cex = 0.8)),
#            xticks = myticks, 
#            grid = TRUE)
# if(svg) {dev.off()}
# if(extra) {system(paste0("open ", horm, "forrestplot_GoalsDiff.svg")) }

# OR
myORticks                 <- log(c(1, 1.5, 2, 3, 5, 10, 15))
attr(myORticks, "labels") <- c("1", "1.5", "2", "3", "5", "10", "15")
horm <- "ITT"
datentabelle <- estimands[6:10, 3:5]
colnames(datentabelle) <- c("mean", "lower", "upper")
hoch <- 2 + nrow(datentabelle) / 6
if(svg) {
  svg(filename = paste0(endpoint, horm, "forrestplot_OR.svg"),
      width = breit, height = hoch)
}
forestplot(datentabelle,
           new_page = TRUE,
           is.summary = rep(FALSE,  nrow(datentabelle)),
           xlab = "Odds Ratio",
           #           clip = c(0.5, 12.5), 
           boxsize = 0.25, 
           zero = 1,
           xlog = TRUE, 
           col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
           vertices = TRUE,
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 0.8), 
                            xlab  = gpar(cex = 0.8), 
                            label = gpar(cex = 0.8)),
           xticks = myORticks, 
           grid = exp(myORticks))
if(svg) {dev.off()}
if(extra) {system(paste0("open ", endpoint, horm, "forrestplot_OR.svg")) }
#
# horm <- "PP"
# datentabelle <- estimandsPP[7:11, 3:5]
# colnames(datentabelle) <- c("mean", "lower", "upper")
# if(svg) {
#   svg(filename = paste0(endpoint, horm, "forrestplot_OR.svg"),
#       width = breit, height = hoch)
# }
# forestplot(datentabelle,
#            new_page = TRUE,
#            is.summary = rep(FALSE,  nrow(datentabelle)),
#            xlab = "Difference",
#            clip = c(floor(min(datentabelle[,"lower"], na.rm = TRUE)), 
#                   ceiling(max(datentabelle[,"upper"], na.rm = TRUE))), 
#            boxsize = 0.25, 
#            zero = 1,
#            xlog = FALSE, 
#            col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
#            vertices = TRUE,
#            txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 0.8), 
#                             xlab  = gpar(cex = 0.8), 
#                             label = gpar(cex = 0.8)),
#            xticks = myticks, 
#            grid = TRUE)
# if(svg) {dev.off()}
# if(extra) {system(paste0("open ", endpoint, horm, "forrestplot_Diff.svg")) }
par(ot)


# Output ICC
if(Deutsch) {
  ICCde <- ICC
  colnames(ICCde) <- c("ICC", "UKI", "OKI")
  cat("Der ICC für den Gesamtscore war ",
      unlist(ICC["GOALS", 1]), " (95%KI ",
      unlist(ICC["GOALS", 2]), " bis ",
      unlist(ICC["GOALS", 3]), "). ", 
      "Der ICC für einzelne Aspekte reichte von  ",
      min(unlist(ICC[1:4, 1])), " bis ", 
      max(unlist(ICC[1:4, 1])),".")
  cat(pander(ICCde, caption = paste0("Der Intraklasskorrelationskoeffizient (ICC) von ",
                                   raters, " in der zweifachen Klassifikation mit 95% Konfidenzintervall. UKL = untere Konfidenzgrenze, OKL = obere Konfidenzgrenze")))
} else {
  cat("The ICC for the GOALS total score was ",
      unlist(ICC["GOALS", 1]), " (95%CI ",
      unlist(ICC["GOALS", 2]), " to ",
      unlist(ICC["GOALS", 3]), "). ", 
      "The ICC for the single aspects ranged from ",
      min(unlist(ICC[1:4, 1])), " to ", 
      max(unlist(ICC[1:4, 1])),".")
  cat(pander(ICC, caption = paste0("Two-way agreement intraclass correlation (ICC) of ",
                                   raters, " raters with 95% confidence interval. LCL = lower confidence limit, UCL = upper confidence limit")))
}
} # sumscorefit
# sumscorefit(daten, svg = FALSE, covar = c("Dauer der Erkrankung (Tage)", "Alter.Patient", "BMI"))

