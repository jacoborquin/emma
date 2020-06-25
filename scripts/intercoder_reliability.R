# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# extracts stats about intercoder reliability


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
source("utils.R")

# loading data
coder1 = read_csv(
	file.path(dataDir, "EMMA_intercoder_reliability_data_1.csv")
)
coder2 = read_csv(
	file.path(dataDir, "EMMA_intercoder_reliability_data_2.csv") 
)


# ----------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------

# join coder realiability data
ICCdata = left_join(coder1, coder2, by = "Paper")
# print.data.frame(ICCdata)

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
	" research domain, $\\textrm{ICC} = ", round(domain_kappa$value, 3),"$,",
	" eye tracker model, $\\textrm{ICC} = ", round(ET_kappa$value, 3),"$,",
	" dependent variable, $\\kappa = ", round(DV_kappa$value, 3),"$,", 
	" independent variable, $\\kappa = ", round(IV_kappa$value, 3),"$"
)
cat(result, file = file.path(tablesDir, "intercoder_reliability.tex"))
