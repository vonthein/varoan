#
# function for multimodel analysis with
#
# GLM for positive variable (time)
# and percentage 
# measured at baseline and endpoint
# for mean and median change
#
glmfit <- function(daten = daten, # preprocessed in different formats
                   plot_boxcox = TRUE, # plot all ?
                   treat = "arm",      # character name of dichotomous variable
                   Y1 = c("Operationszeit.E", "`OP-Zeit Proband:in.E`"),      # endpoints, vector of character
                   Y0 = c("Operationszeit", "`OP-Zeit Proband:in`"),      # baselines in order of endpoints, vector of character
                   X = c("`Dauer der Erkrankung (Tage)`", "Alter.Patient", "BMI"),  # covariates besides baseline, vector of character
                   conf.int = 0.95 # degree between 0 and 1
                   ){
  
  library(lme4)
  library(MASS)
  library(gamlss)

  datenFAS <- daten$FAS
  datensFAS <- daten$sFAS
  
  # datenFAS$id <- 1:nrow(datenFAS)
  # datenFAS$base <- rep(0, nrow(datenFAS)) 
  # datenFAS$ende <- rep(1, nrow(datenFAS))
  # 
  # if(X=="") {
  #   datensFAS <- data.frame(time = c(datenFAS$base, datenFAS$ende),
  #                             id = c(datenFAS$id,   datenFAS$id),
  #                          treat = c(datenFAS$treat,datenFAS$treat)
  #   )
  #   datensFAS <- data.frame(datensFAS, rbind(datenFAS[,Y0], datenFAS[,Y1]))
  #     
  #     # rbind(datenFAS[,c("id", Y0, treat, "base")],
  #     #                  datenFAS[,c("id", Y1, treat, "ende")])
  # } else {
  #   datensFAS <- rbind(datenFAS[,c("id", Y0, treat, "base", X)],
  #                      datenFAS[,c("id", Y1, treat, "ende", X)])
  # } 
  

  for(y in 1:length(Y1)) {
  lhs <- paste(Y1[y], "~", treat, "+", Y0[y])
  Lhs <- paste(Y1[y], "~", treat, "* time + ", Y0[y])
  rhs <- paste(X, collapse = "+")
  formulae <- vector()
  formulae[1] <- lhs
  formulae[2] <- paste(lhs, "+", rhs)
  formulae[3] <- paste(lhs, "+ (1|id)")
  formulae[4] <- paste(lhs, "+", rhs, "+ (1|id)")
  boxcox(data = datenFAS, formula = formulae[1])
  boxcox(data = datenFAS, formula = formulae[2])

  fit.norm    <- glm(data = datenFAS, formula = formulae[1], gaussian("identity"))
  fit.gamma   <- glm(data = datenFAS, formula = formulae[1], Gamma("identity"))
  fit.invnorm <- glm(data = datenFAS, formula = formulae[1], inverse.gaussian("identity"))
  fit.norma   <- glm(data = datenFAS, formula = formulae[2], gaussian("identity"))
  fit.gammaa  <- glm(data = datenFAS, formula = formulae[2], Gamma("identity"))
  fit.invnorma<- glm(data = datenFAS, formula = formulae[2], inverse.gaussian("identity"))
  fit.lognorm <- gamlss(formulae[1], family = LOGNO2(mu.link = "identity"))

  }
  
  fit.mix <- lmer(`OP-Zeit Proband:in` ~ time * arm + (1 | Proband), data = datensFAS)
ci.mix <- confint(fit.mix, "timeFinal:armWaiting")
fit.imix <- glmer(`OP-Zeit Proband:in` ~ time * arm + (1 | Proband),inverse.gaussian("identity"), data = datensFAS)
ci.imix <- confint(fit.imix, "timeFinal:armWaiting", method = "Wald")
# fit.imixa <- glmer(`OP-Zeit Proband:in` ~ time * arm + (1 | Proband) + `Dauer der Erkrankung (Tage)` + Alter.Patient + BMI,
#                    inverse.gaussian("identity"), data = datensFAS)
# ci.imixa <- confint(fit.imixa, "timeFinal:armWaiting", method = "Wald")
fit.mixa <- lmer(`OP-Zeit Proband:in` ~ time * arm + (1 | Proband) + `Dauer der Erkrankung (Tage)` + Alter.Patient + BMI, data = datensFAS)
ci.mixa <- confint(fit.mixa, "timeFinal:armWaiting")
fit.ia <- glm(`OP-Zeit Proband:in` ~ time * arm + `Dauer der Erkrankung (Tage)` + Alter.Patient + BMI, data = datensFAS)
ci.ia <- confint(fit.ia, "timeFinal:armWaiting")
ci.norm <- confint(fit.norm,"datenFAS$Randomisierung")
ci.invnorm <- confint(fit.invnorm,"datenFAS$Randomisierung")
ci.norma <- confint(fit.norma,"datenFAS$Randomisierung")
ci.invnorma <- confint(fit.invnorma,"datenFAS$Randomisierung")
ci.lognorm <- confint(fit.lognorm)["mu.datenFAS$Randomisierung",]
ci.gamma <- confint(fit.gamma,"datenFAS$Randomisierung")
fitnorm <- glm(datenFAS$`Operationszeit.E`~datenFAS$Operationszeit+datenFAS$Randomisierung, gaussian("identity"))
fitgamma <- glm(datenFAS$`Operationszeit.E`~datenFAS$Operationszeit+datenFAS$Randomisierung, Gamma("identity"))
fitinvnorm <- glm(datenFAS$`Operationszeit.E`~datenFAS$Operationszeit+datenFAS$Randomisierung, inverse.gaussian("identity"))
fitlognorm <- gamlss(datenFAS$`Operationszeit.E`~datenFAS$Operationszeit+datenFAS$Randomisierung, family = LOGNO2(mu.link = "identity"))
# mixed model
# library(lme4)
fitmix <- lmer(Operationszeit ~ time * arm + (1 | Proband), data = datensFAS)
cimix <- confint(fitmix, "timeFinal:armWaiting")
cinorm <- confint(fitnorm,"datenFAS$Randomisierung")
ciinvnorm <- confint(fitinvnorm,"datenFAS$Randomisierung")
cilognorm <- confint(fitlognorm)["mu.datenFAS$Randomisierung",]
cigamma <- confint(fitgamma,"datenFAS$Randomisierung")
# mixed adjusted
fitmixa <- lmer(Operationszeit ~ time * arm + (1 | Proband) + `Dauer der Erkrankung (Tage)` + Alter.Patient + BMI, 
                data = datensFAS)
cimixa <- confint(fitmixa, "timeFinal:armWaiting")
# fitimixa <- glmer(Operationszeit ~ time * arm + (1 | Proband) + 
#                     I(`Dauer der Erkrankung (Tage)` - mean(`Dauer der Erkrankung (Tage)`)) + 
#                     I(Alter.Patient - mean(Alter.Patient)) + 
#                     I(BMI - mean(BMI)), 
#                 inverse.gaussian("identity"), data = datensFAS)
fitia <- glm(Operationszeit ~ time * arm +
               I(`Dauer der Erkrankung (Tage)` - mean(`Dauer der Erkrankung (Tage)`)) + 
               I(Alter.Patient - mean(Alter.Patient)) + 
               I(BMI - mean(BMI)), 
             inverse.gaussian("identity"), data = datensFAS)
ciia <- confint(fitia, "timeFinal:armWaiting")

CIs <- round(rbind(-cimix[2:1], -cimixa[2:1], cinorm, cilognorm, cigamma, ciinvnorm, -ciia[2:1], 
                   -ci.mix[2:1], -ci.mixa[2:1], ci.norm, ci.lognorm, ci.gamma, ci.invnorm, -ci.ia[2:1]), 1)
cis <- apply(CIs, 1, paste0, collapse = " to ")
if(Deutsch) {
  worte <- c("Finale Mediane",
             "Mediane Änderung",
             "Wechselwirkung bei Normalverteilung",
             "Adjustierte Wechselwirkung bei Normalverteilung",
             "Erwartungswerte von Normalverteilungen*",
             "Mediane von Lognormalverteilungen*",
             "Erwartungswerte von Gamma-Verteilungen*",
             "Erwartungswerte von inversen Gauß-Verteilungen*",
             "Adjustierte Wechselwirkung bei inverser Gauß-Verteilung"
  )
} else {
  worte <- c("Final medians",
             "Medians of changes",
             "Interaction of normals",
             "Adjusted interaction of normals",
             "Means of normals*",
             "Medians of lognormals*",
             "Means of gammas*",
             "Means of inverse gaussians*",
             "Adjusted interaction of inverse gaussians"
  )
}
time.tab <- data.frame(Time = c(rep("Operation", 9), rep("Resident surgeon", 9)),
                       Difference = rep(worte, 2),
                       Minutes = round(c(erg2[1,3],
                                         erg2[3,3],
                                         -fixed.effects(fitmix)["timeFinal:armWaiting"],
                                         -fixed.effects(fitmixa)["timeFinal:armWaiting"],
                                         fitnorm$coefficients["datenFAS$Randomisierung"],
                                         fitlognorm$mu.coefficients["datenFAS$Randomisierung"],
                                         fitgamma$coefficients["datenFAS$Randomisierung"],
                                         fitinvnorm$coefficients["datenFAS$Randomisierung"],
                                         -fitia$coefficients["timeFinal:armWaiting"],
                                         erg2[2,3],
                                         erg2[4,3],
                                         -fixed.effects(fit.mix)["timeFinal:armWaiting"],
                                         -fixed.effects(fit.mixa)["timeFinal:armWaiting"],
                                         fit.norm$coefficients["datenFAS$Randomisierung"],
                                         fit.lognorm$mu.coefficients["datenFAS$Randomisierung"],
                                         fit.gamma$coefficients["datenFAS$Randomisierung"],
                                         fit.invnorm$coefficients["datenFAS$Randomisierung"],
                                         -fit.ia$coefficients["timeFinal:armWaiting"]
                       ), 1),
                       ` 95% conf. int.` = c(d2[1,2], d2[3,2], cis[1:7],
                                             d2[2,2], d2[4,2], cis[8:14]))
row.names(time.tab) <- NULL
if(Deutsch) {
  cat(pander(time.tab), "Zeitdifferenz adjustiert nach der Baseline im FAS. * adjustiert nach der Baseline.")
} else {
cat(pander(time.tab), "Time differences adjusted for baseline in the full analysis set. * adjusted for baseline.")               
} # if
} # glmfit

glmfit(datenFAS = data$FAS, # data.frame, wide format
                   plot_boxcox = TRUE, # plot all ?
                   treat = "arm",   # character name of dichotomous variable
                   Y1 = c("`Operationszeit.E","`OP-Zeit Proband:inE`"),      # endpoints, vector of character
                   Y0 = c("Operationszeit","`OP-Zeit Proband:in`"),      # baselines in order of endpoints, vector of character
                   X = c("BMI", "Alter.Patient", "Dauer der Erkrankung (Tage)"),  # covariates besides baseline, vector of character
                   conf.int = 0.95 # degree between 0 and 1
)