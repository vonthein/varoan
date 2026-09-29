# function for data rearrangement
# aggregations and logits

# take declaration of roles to make different stacked data.frames 
sumscore_preprocess <- function(daten = daten, # data.frame
                            groups = 2,
                            times = 2, 
                            items = 5,
                            raters = 3, 
                            case_vars = 1:7,
                            case_sums = 8:25,
                            rep_vars = c(26:33, 52:59),
                            rep_sums = c(34:36, 60:62),
                            item_vars = c(seq(37, 36 + items * raters, 1),
                                          seq(63, 62 + items * raters, 1)),
                            PP = "Ausschlussgrund" # character variable name that is 0 for PP 
) {
# for tests
# data = daten; groups = 2; times = 2; items = 5; raters = 3;case_vars = 1:7; case_sums = 8:25
# rep_vars = c(26:33, 52:59); rep_sums = c(34:36, 60:62); item_vars = c(seq(37, 36 + items * raters, 1), seq(63, 62 + items * raters, 1))
# PP = "Ausschlussgrund"
  
  # Transformation of bounded data
  logit <- function(x, inf = 0, sup = 1){
    log((x-inf)/(sup-x))
  }
  # Backtransformation
  expit <- function(x, inf = 0, sup = 1){
    inf + (sup - inf) * exp(x) / (1 + exp( x))
  }
  
  N <- nrow(daten)
  K <- ncol(daten)
  interim <- list()
  rep_var_n <- length(rep_vars) / times
  rep_sum_n <- length(rep_sums) / times
  items_n  <- length(item_vars) / times # items * raters
  rep_n <- rep_var_n + rep_sum_n + items_n
  
  # Sum scores
  # wide
  rep_sum_start <- length(c(case_vars, case_sums)) + length(rep_vars) / times
  rep_sum_startE <- K - items_n - raters # final time is primary
  for(a in 1:raters) {
    daten[, paste0("logit.Gesamt", a)] <- logit(daten[,rep_sum_start  + a], 5 - 1/3, 25 + 1/3)
  } # for raters
  for(a in 1:raters) {
    daten[, paste0("logit.GesamtE",a)] <- logit(daten[,rep_sum_startE + a], 5 - 1/3, 25 + 1/3)
  } # for raters
  K_new <- ncol(daten)
  # aggregation to minimize technical variability
  # only baseline and primary
  daten$GOALS1 <- apply(daten[,rep_sum_start  + (1:raters)], 1, median, na.rm = TRUE)
  daten$GOALS2 <- apply(daten[,rep_sum_startE + (1:raters)], 1, median, na.rm = TRUE)
  # for(z in 1:2) { # only baseline and primary
  for(z in 1:times) {
    daten[, paste0("logit.GOALS", z)] <- apply(daten[,(K_new - z * raters + 1):(K_new - (z - 1) * raters)], 
                                               1, mean, na.rm = TRUE)
    }
  daten$GOALSchange <- daten[, paste0("GOALS", times)] - daten$GOALS1
  daten$logit.GOALSchange <- daten[, paste0("logit.GOALS", times)] - daten$logit.GOALS1
  
  
  # stack data for points in time
    for(r in 1:times) {
    interim[[r]] <- daten[, c(case_vars, case_sums, 
                              rep_vars[((r-1)*rep_var_n+1):(r*rep_var_n)],
                              rep_sums[((r-1)*rep_sum_n+1):(r*rep_sum_n)],
                              item_vars[((r-1)*items_n+1):(r*items_n)])] }
  datens <- interim[[1]]
  for(r in 2:times) {
    datens2 <- interim[[r]]
    colnames(datens2) <- colnames(datens)
    datens <- rbind(datens, datens2)
  }
  datens$Zeit <- rep(1:times, each = nrow(datens2))

  # single stacked
  kt <- ncol(datens) # teporary number of columns
  repsum_start <- length(c(case_vars, case_sums)) + rep_var_n 
    for(a in 1:raters) {
      datens[, paste0("logit.Gesamt", a)] <- logit(datens[,repsum_start + a], 5 - 1/3, 25 + 1/3)
    }
  datens$GOALS <- apply(datens[,repsum_start + (1:raters)], 1, 
                        median, na.rm = TRUE)
  datens$logit.GOALS <- apply(datens[,(kt + 1):(kt + raters)], 1, 
                            mean, na.rm = TRUE)
# first case variable is caseID per default
# second case variable is assigned group per default
  daten$arm <- as.factor(daten[,2])
  datens$arm <- as.factor(datens[,2])
# last time is primary per default
datens$time <- as.factor(ifelse(datens$Zeit == 1, "Baseline", 
                                ifelse(datens$Zeit == max(datens$Zeit), "Final", 
                                       datens$Zeit))) 
# jitter in plot
n2 <- dim(datens)[1]
datens$when <- datens$Zeit + 
  0.1 * as.numeric(datens$arm)  - rep(0.15, n2) +
  rnorm(n2, 0, 0.03)

# doubly stacked data for assessor agreement among three assessors
# case and observation covariates
general <- c(1:repsum_start, which(names(datens)== "arm"),
             which(names(datens)== "time"), which(names(datens)== "Zeit")) 
assessor <- c(seq(repsum_start + 1, repsum_start + (items + 1) * raters, raters)) # items

datenss <- datens[, c(general, assessor)] 
for(i in 1:(raters - 1)) {
  ass <- datens[, c(general, assessor + i)]
  colnames(ass) <- colnames(datenss)
  datenss <- rbind(datenss, ass)
}
datenss$GOALS <- datenss$GOALS1
datenss$logit.GOALS <- logit(datenss$GOALS, 5, 25)
datenss$rater <- as.factor(rep(1:3, each = n2))
# table(datenss$time, datenss$arm, datenss$rater)


# intermediate data set
# assessors stacked, times not
general <- c(1:rep_sum_start, 
             (rep_sum_startE - rep_var_n + 1):rep_sum_startE,
             which(names(daten)== "arm"),
             which(names(daten)== "time"), 
             which(names(daten)== "Zeit"))
assessor <- c(seq(rep_sum_start + 1, rep_sum_start + (rep_sum_n + items) * times, raters),
              seq(rep_sum_startE + 1, rep_sum_startE + (rep_sum_n + items) * raters, raters))
#general <- c(1:16, 33:39, 78:83)
#assessor <- c(17, 20, 23, 26, 30, 40, 43, 46, 49, 53, 56, 59, 62, 65, 68, 69, 72, 75)

datensr <- daten[, c(general, assessor)]
for(a in 1:(raters - 1)) {
  ass <- daten[, c(general, assessor + a)]
  colnames(ass) <- colnames(datensr)
  datensr <- rbind(datensr, ass)
}
datensr$rater <- as.factor(rep(1:3, each = N))
datensr$GOALS <- datensr$GOALSE1 - datensr$GOALS1

# definition of analysis data sets
nFAS <- is.na(daten$GOALS2)
PP  <- daten[, PP] == 0
data <- list()
data$FAS <- daten[!nFAS,]
data$PP  <- daten[PP,]
data$sFAS <- datens[!rep(nFAS, times),]
data$sPP  <- datens[rep(PP, times),]
data$ssFAS <- datenss[!rep(nFAS, times * raters),]
data$ssPP  <- datenss[rep(PP, times * raters),]
data$srFAS <- datensr[!rep(nFAS, raters),]
data$srPP  <- datensr[rep(PP, raters),]

return(data)
} # sumscore_preprocess

daten <- sumscore_preprocess(daten = daten)
