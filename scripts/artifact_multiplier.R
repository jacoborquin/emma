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
# setwd("/Users/au161118/Dropbox/ASB/Admin stuff/Posters & Papers/PAPERS/EMMA/scripts/emma/scripts")

# housekeeping
rm(list = ls())

# import packages and functions
source("./scripts/utils.R")

# loading data
data = as.data.table(read_excel(
	file.path(dataDir, "EMMA_ES_data.xlsx")
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

# impute values when ET "unknown". we use the average accuracy and precision
ET_specs$Accuracy[is.na(ET_specs$Accuracy)] = mean(ET_specs$Accuracy, na.rm = T)
ET_specs$Precision[is.na(ET_specs$Precision)] = mean(ET_specs$Precision, na.rm = T)

# merge with data file for precision values of ET
data = merge(data, ET_specs, by = c("Eye.tracker"), all.x = TRUE)

# long format for analysis of accuracy and precision on effect size attenuation
data_long = as.data.table(
	gather(data, DV, R, fix.count:dwell.count, factor_key = TRUE)
)
data_long = data_long[data_long$R != "NA"]

# averaging studies with multiple effect sizes
data_long = data_long[, 
	list(
		R = FisherZInv(mean(FisherZ(R))),
		Rz = mean(FisherZ(R)),
		Precision = unique(Precision),  
	    Accuracy = unique(Accuracy), 
	    Eye.tracker = unique(Eye.tracker)
	), 
	by = c("Study", "IV")
]

# verify
any(table(data_long$Study, data_long$IV) > 1) == FALSE

# -----
# Test whether eye tracker accuracy and precision attenuate effect sizes 
# While precision is a better predictor we will use accuracy for artifact multiplier
# the reason is that we have more exact numbers on accuracy than precision
# -----

# mixed effects linear models
model1 <- lmer(Rz ~ + (1|Study), data_long)
model2 <- lmer(Rz ~ Accuracy + (1|Study), data_long)
model3 <- lmer(Rz ~ Accuracy + IV + (1|Study), data_long)
model4 <- lmer(Rz ~ Precision + (1|Study), data_long)
model5 <- lmer(Rz ~ Precision + IV + (1|Study), data_long)
model6 <- lmer(Rz ~ Accuracy + Precision + IV + (1|Study), data_long)

# includes accuracy, but not precision
#anova(model1,model2,model3,model4,model5,model6) 
#summary(model3)

# compute artefact multiplier based on accuracy for psychometric meta-analysis
b0 = coef(summary(model3))[1,1] 
b1 = coef(summary(model3))[2,1]
data$a_acc = (b0 + b1*data$Accuracy)/b0 

# paste result for manuscript
result <- paste0(
  "$\\beta_0=", 
  round(coef(summary(model3))[1,1], 3), 
  "$, $\\SE=", 
  round(coef(summary(model3))[1,2], 3), 
  "$, $t=", 
  round(coef(summary(model3))[1,4], 3), 
  "$, $", 
  ifelse(round(coef(summary(model3))[1,5], 3) == 0, "p<0.001", paste0("p=", round(coef(summary(model3))[1,5], 3))), 
  "$, $\\beta_{\\textrm{accuracy}} =",
  round(coef(summary(model3))[2,1], 3), 
  "$, $\\SE=", 
  round(coef(summary(model3))[2,2], 3), 
  "$, $t=", 
  round(coef(summary(model3))[2,4], 3), 
  "$, $", 
  ifelse(round(coef(summary(model3))[2,5], 3) == 0, "p<0.001", paste0("p=", round(coef(summary(model3))[2,5], 3))), 
  "$;"
)
cat(result, file = file.path(tablesDir, "artifactregresult.tex"))

# -----
# Save corrected data
# ----

write_csv(data, file.path(dataDir, "EMMA_ES_data_corrected.csv"))

# -----
# Plot of accuracy on effect size
# ----

# create ES-accuracy  scatter plot
figure <- 
	ggplot(data = data_long, aes(Accuracy, Rz)) +
	geom_point(alpha = .5, size = pointSize) +
	# geom_smooth(method = "lm", color = "black", size = lineSize*2) +
	geom_abline(intercept = b0, slope = b1, size = lineSize*1) +
	scale_x_continuous("Eye tracker accuracy (visual angle)", breaks = seq(0, 1.5, by=.1)) +
	ylab("Effect size (z)") +
	mytheme
filename <- file.path(figsDir, "ET_accuracy_effectsize.pdf")
savePlots(figure, filename, fd_SI_1x1.5)
  
# -----
# Make eye tracker specifications table for manuscript
# -----

# make eye tracker specifications table for manuscript appendix
ET_specs_final = data[, list(A = unique(a_acc), accuracy = unique(Accuracy), precision = unique(Precision)), by= c("Eye.tracker")]
ET_specs_final$A <- round(ET_specs_final$A, 4)
setnames(ET_specs_final, c(1:4), c("Eye tracker model", "$a_a$", "Accuracy", "Precision"))
write_csv(ET_specs_final, file.path(tablesDir, "eyetracker_specs.csv"))

# latex version
tab_caption <- "Eye tracker specifications table, with accuracy and precision for each eye tracker as extracted from the manufacturer website, and computed artifact multiplier used for correcting for a bias in effect size estimates."
tab_label <- "tab:eyetracker_specifications"
tab_note <- paste0("\\hline \n \\multicolumn{4}{l}",
           "{\\scriptsize{\\textit{Note.} $a_a$ = artifact multiplier.}} \n")
print(
	xtable(
		ET_specs_final, 
		caption = tab_caption, 
		label = tab_label, 
		align = "llccc",
		digits = c(0,0,4,2,2)
	), 
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
