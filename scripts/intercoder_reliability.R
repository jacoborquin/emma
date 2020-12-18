# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# extracts stats about intercoder reliability


# ----------------------------------------------------------------------
# Loading data
# ----------------------------------------------------------------------

# specify your path here if you want to use this script interactively, 
# and uncomment the line:
# setwd("/home/hstojic/research/project/attention_meta/scripts")
# setwd("/Users/au161118/Dropbox/ASB/Admin stuff/Posters & Papers/PAPERS/EMMA/scripts/emma/scripts")

# housekeeping
rm(list = ls())

# import packages and functions
source("utils.R")

# loading data
coder1 = read_csv(
	file.path(dataDir, "data_intercoder_reliability_1.csv")
)
coder2 = read_csv(
	file.path(dataDir, "data_intercoder_reliability_2.csv") 
)
coder3 = read_csv(
  file.path(dataDir, "data_intercoder_reliability_3.csv") 
)
coder4 = read_csv(
  file.path(dataDir, "data_effect_sizes.csv") 
)


# ----------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------

# join coder realiability data for old data
ICCdata = left_join(coder1, coder2, by = "Paper")

# kappa for categorical variables
DV_kappa = kappa2(ICCdata[,c("DV.x","DV.y")], "unweighted")
IV_kappa = kappa2(ICCdata[,c("IV.x","IV.y")], "unweighted")
ET_kappa = kappa2(ICCdata[,c("EyeTracker.x","EyeTracker.y")], "unweighted")
domain_kappa = kappa2(ICCdata[,c("Domain.x","Domain.y")], "unweighted")

# ICC for continuous variables
ES_ICC = icc(ICCdata[,c("ES.x","ES.y")], model="oneway", type="agreement")
N_ICC = icc(ICCdata[,c("N.x","N.y")], model="oneway", type="agreement")

# paste summary of inter coder reliability for manuscript
result <- paste0(
	"effect size, $\\textrm{ICC} = ", round(ES_ICC$value, 3),"$,",
	" sample size, $\\textrm{ICC} = ", round(N_ICC$value, 3),"$,",
	" research domain, $\\kappa = ", round(domain_kappa$value, 3),"$,",
	" eye tracker model, $\\kappa = ", round(ET_kappa$value, 3),"$,",
	" dependent variable, $\\kappa = ", round(DV_kappa$value, 3),"$,", 
	" independent variable, $\\kappa = ", round(IV_kappa$value, 3),"$"
)
cat(result, file = file.path(tablesDir, "intercoder_reliability.tex"))

# join coder reliability for recoded ES data (revision round 2)
coder3 = melt(as.data.table(coder3), id.vars = c(2,3,10,11), measure.vars = c(4:7))
coder4 = melt(as.data.table(coder4), id.vars = c(2,3,10,11), measure.vars = c(4:7))
ESdata = data.table(merge(coder3, coder4, by = c("Study", "IV", "Research.Domain", "Alt.att", "variable"), all.x = TRUE))
ESdata = ESdata[is.na(ESdata$value.x) == FALSE & is.na(ESdata$value.y) == FALSE]
ES_ICC2 = icc(ESdata[,c("value.x","value.y")], model="oneway", type="agreement")

result <- paste0("$\\textrm{ICC} = ", round(ES_ICC2$value, 3),"$")
cat(result, file = file.path(tablesDir, "intercoder_reliability2.tex"))

