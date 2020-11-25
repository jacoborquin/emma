# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# executes whole analysis in a correct sequence
# figures are placed in ../figs folder
# tables and textual output in ../tables folder


# ----------------------------------------------------------------------
# Loading data
# ----------------------------------------------------------------------

# specify your path here if you want to use this script interactively, 
# and uncomment the line:
# setwd("/home/hstojic/Research/project/attention_meta/scripts")
 setwd("/Users/au161118/Dropbox/ASB/Admin stuff/Posters & Papers/PAPERS/EMMA/scripts/emma/scripts/")

# adds eye tracker artifact multiplier to the data
# produces associated figures and tables
source("./scripts/artifact_multiplier.R")

# adds metric correction factors to the data
# produces associated figures and tables
source("./scripts/metric_correction.R")

# psychometric meta-analysis
# produces tables with main results, forest and funnel plots
source("meta_analysis.R")

# import packages and functions
source("intercoder_reliability.R")
