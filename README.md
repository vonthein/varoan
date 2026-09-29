Analysis of variance  
R package by Reinhard Vonthein to produce the analyses in DOI 10.1097/JS9.0000000000002304  
under development  
After downloading that article's data in DOI 10.5281/zenodo.10428829,  
prepare this specific file for further processing. 
daten <- file_preprocess(datafile = "../Data/NOVICE_rawDataAnonymized.xlsx")
Preprocess it further to have all formats, transformations, and summaries needed. 
daten <- sumscore_preprocess(daten = daten)
Then produce the first tables of estimates meaning nearly the same   
sumscorefit(daten = daten, svg = FALSE, covar = c("Dauer der Erkrankung (Tage)", "Alter.Patient", "BMI"))
Produce another table of estimates with glmfit(). 
glmfit(datenFAS = data$FAS, 
       plot_boxcox = TRUE, 
       treat = "arm",  
       Y1 = c("`Operationszeit.E","`OP-Zeit Proband:inE`"),      # endpoints
       Y0 = c("Operationszeit","`OP-Zeit Proband:in`"),      # baselines 
       X = c("BMI", "Alter.Patient", "Dauer der Erkrankung (Tage)")) # covariates 
