# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# computes artifact multiplier based on eye tracker accuracy, saves
# back the data with the multiplier, also figures with correction factors
# and a table with eye tracker specifiation for the manuscript


# ----------------------------------------------------------------------
# Loading data
# ----------------------------------------------------------------------

# specify your path here if you want to use this script interactively, 
# and uncomment the line:
# setwd("/home/hstojic/Research/project/attention_meta/scripts")

# housekeeping
rm(list = ls())

# import packages and functions
source("utils.R")

# loading data
data = as.data.table(read_csv2(
	file.path(dataDir, "EMMA_ES_data.csv")
))
ET_specs = as.data.table(read_csv2(
	file.path(dataDir, "EMMA_ET_specs_data.csv")
))


# ----------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------

# -----
# Prepare data
# ----

# impute eye tracker "Unknown" when eye tracker name not identified
data$Eye.tracker[is.na(data$Eye.tracker)] = "Unknown"

# merge with data file for precision values of ET
data = merge(data, ET_specs, by = c("Eye.tracker"), all.x = TRUE)

# long format for analysis of accuracy and precision on effect size attenuation
data_long = as.data.table(
	gather(data, DV, R, fix.count:dwell.count, factor_key = TRUE)
)
data_long = data_long[data_long$R != "NA"]

# the models are computed on the absolute effect sizes
data_long$R_abs = abs(data_long$R)  

# averaging studies with multiple effect sizes
data_long = data_long[, 
	list(
		R_abs = mean(R_abs), 
		Precision = unique(Precision),  
	    Accuracy = unique(Accuracy), 
	    Eye.tracker = unique(Eye.tracker), 
	    IV = unique(IV)
	), 
	by = c("Study", "IV")
]

# chi square test
table = table(data_long$IV, data_long$Eye.tracker)
chisq.test(table)


# -----
# Test whether eye tracker accuracy or precision attenuate effect sizes 
# -----

# mixed effects linear models
model1 <- lmer(R_abs~ + (1|Eye.tracker), data_long)
model2 <- lmer(R_abs~ + (1|Study), data_long)
model3 <- lmer(R_abs~ + (1|Eye.tracker) + (1|Study), data_long)
model21 <- lmer(R_abs~ Accuracy + (1|Study), data_long)
model212 <- lmer(R_abs~ Accuracy + IV + (1|Study), data_long)
model22 <- lmer(R_abs~ Precision + (1|Study), data_long)
model222 <- lmer(R_abs~ Precision + IV + (1|Study), data_long)
model23 <- lmer(R_abs~ Precision + Accuracy + (1|Study), data_long)

# winning model21 has lowest BIC
# includes accuracy, but not precision
anova(model1,model2,model3,model21,model22,model23,model212,model222) 

# compute artefact multiplier based on accuracy for psychometric meta-analysis
b0 = coef(summary(model21))[1,1] 
b1 = coef(summary(model21))[2,1]
data$a_acc = (b0 + b1*data$Accuracy)/b0 # see eq. 1 in manuscript


# -----
# Save corrected data
# ----

write_csv(data, file.path(dataDir, "EMMA_ES_data_corrected.csv"))


# -----
# Plot of accuracy on effect size
# ----

# strangely, has two columns with same name??
data_long[,7] <- NULL

# perhaps add slope coeff?
figure <- 
	ggplot(data = data_long, aes(Accuracy, R_abs)) +
	geom_point(alpha = .5, size = pointSize) +
	# geom_smooth(method = "lm", color = "black", size = lineSize*2) +
	geom_abline(intercept = b0, slope = b1, size = lineSize*1) +
	scale_x_continuous("Eye tracker accuracy", breaks = seq(0.4, 1, by=.1)) +
	ylab("Observed effect size") +
	mytheme
filename <- file.path(figsDir, "ET_accuracy_effectsize.pdf")
savePlots(figure, filename, fd_SI_1x1.5)
  

# -----
# Make accuracy results table for manuscript
# -----

acc_table = data.table(summary(model21)$coef) 
obs = data.table("Number of observations =",NROW(data_long),NA,NA,NA)
Ll = data.table("Log likelihood =",as.numeric(summary(model21)$logLik),NA,NA,NA)
AIC = data.table("AIC =",AIC(model21),NA,NA,NA)
BIC = data.table("BIC =",BIC(model21),NA,NA,NA)
Raneff = data.table("SD(study) =",data.frame(summary(model21)$varcor)[1,5],NA,NA,NA)
l = list(acc_table,obs,Ll,AIC,BIC,Raneff)
acc_table = rbindlist(l)
acc_table = acc_table %>% mutate_at(vars(2:5), round, 3)
setnames(acc_table, c(2:5), c("SE", "df","t","p"))
write_csv(acc_table, file.path(tablesDir, "accuracy_winner_model.csv"))

# latex version
print(
	xtable(acc_table), 
	include.rownames = FALSE,
	file = file.path(tablesDir, "accuracy_winner_model.tex")
)


# -----
# Make eye tracker specifications table for manuscript
# -----

# make eye tracker specifications table for manuscript appendix
ET_specs_final = data[, list(A = unique(a_acc), accuracy = unique(Accuracy), precision = unique(Precision)), by= c("Eye.tracker")]
ET_specs_final$A <- round(ET_specs_final$A, 4)
setnames(ET_specs_final, c(1:2), c("Eye tracker model", "$a_a$"))
write_csv(ET_specs_final, file.path(tablesDir, "eyetracker_specs.csv"))

# latex version
tab_caption <- "Eye tracker specifications table"
tab_label <- "tab:eyetracker_specifications"
tab_note <- paste0("\\hline \n \\multicolumn{4}{l}",
           "{\\scriptsize{\\textit{Note.} $a_a$ = artifact multiplier.}} \n")
print(
	xtable(ET_specs_final, caption=tab_caption, label=tab_label), 
	include.rownames = FALSE,
	caption.placement = "top", 
	hline.after = c(-1, 0),
	add.to.row = list(
		pos = list(nrow(ET_specs_final)),
        command = tab_note
    ),
    sanitize.text.function = function(x){x},
    file = file.path(tablesDir, "eyetracker_specs.tex")
)
