# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# compute correction factor am for transforming effect sizes in different 
# metrics
# correction factor ratio = M_to / M_from
# M_. = FisherZ(sum(Fix*N)/N) see eq. 2 in manuscript


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
data = as.data.table(read_csv(
	file.path(dataDir, "EMMA_ES_data_corrected.csv")
))

# ----------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------

# -----
# Prepare data
# ----

# correction factor for fixation count to fixation likelihood and vice versa
TEMP = data[, fix.count:N]
TEMP = TEMP[TEMP$fix.count != "NA" & TEMP$fix.like != "NA"]
FL_to_FC = FisherZ(sum(TEMP$fix.count*TEMP$N)/sum(TEMP$N))/
           FisherZ(sum(TEMP$fix.like*TEMP$N)/sum(TEMP$N))
FC_to_FL = FisherZ(sum(TEMP$fix.like*TEMP$N)/sum(TEMP$N))/
           FisherZ(sum(TEMP$fix.count*TEMP$N)/sum(TEMP$N))

# correction factor for total fixation duration to fixation count
TEMP = data[, fix.count:N]
TEMP = TEMP[TEMP$fix.count != "NA" & TEMP$TFD != "NA"]
TFD_to_FC = FisherZ(sum(TEMP$fix.count*TEMP$N)/sum(TEMP$N))/
            FisherZ(sum(TEMP$TFD*TEMP$N)/sum(TEMP$N))

# correction factor for total fixation duration to fixation likelihood
TEMP = data[, fix.count:N]
TEMP = TEMP[TEMP$fix.like != "NA" & TEMP$TFD != "NA"]
TFD_to_FL = FisherZ(sum(TEMP$fix.like*TEMP$N)/sum(TEMP$N))/
            FisherZ(sum(TEMP$TFD*TEMP$N)/sum(TEMP$N))

# correction factor for total fixation duration to fixation likelihood and fixation count
TEMP = data[, fix.count:N]
TEMP = TEMP[TEMP$dwell.count != "NA" & TEMP$TFD != "NA"]
DC_to_TFD = FisherZ(sum(TEMP$TFD*TEMP$N)/sum(TEMP$N))/
            FisherZ(sum(TEMP$dwell.count*TEMP$N)/sum(TEMP$N))

# no cases for computing DC to FC or FL hence the current approach
DC_to_FC = DC_to_TFD*TFD_to_FC 
DC_to_FL = DC_to_TFD*TFD_to_FL
  
# correcting effect sizes to fixation count and fixation likelihood based on observed metric type
data$fix.count.m = FisherZInv(
  ifelse(is.na(data$fix.count) & is.na(data$fix.like) == F, FisherZ(data$fix.like)*FL_to_FC, 
  ifelse(is.na(data$fix.count) & is.na(data$TFD) == F, FisherZ(data$TFD)*TFD_to_FC, 
  ifelse(is.na(data$fix.count) & is.na(data$dwell.count) == F, FisherZ(data$dwell.count)*DC_to_FC, FisherZ(data$fix.count)))))

data$fix.like.m = FisherZInv(
  ifelse(is.na(data$fix.like) & is.na(data$fix.count) == F, FisherZ(data$fix.count)*FC_to_FL, 
  ifelse(is.na(data$fix.like) & is.na(data$TFD) == F, FisherZ(data$TFD)*TFD_to_FL, 
  ifelse(is.na(data$fix.like) & is.na(data$dwell.count) == F, FisherZ(data$dwell.count)*DC_to_FL, FisherZ(data$fix.like)))))

# -----
# Save corrected data
# ----

write_csv(data, file.path(dataDir, "EMMA_ES_data_corrected.csv"))

# -----
# Make correction factor table for manuscript
# -----

metric_correction_factors <- data.frame(
	c("Fixation count", "Fixation likelihood", "Total fixation duration", 
		"Total fixation duration", "Dwell count", "Dwell count"),
	c("Fixation likelihood", "Fixation count", "Fixation likelihood", 
		"Fixation count", "Fixation likelihood", "Fixation count"),
	round(c(FC_to_FL, FL_to_FC, TFD_to_FL, TFD_to_FC, DC_to_FL, DC_to_FC), 3)
)
colnames(metric_correction_factors)[1:3] <- c(
	"Correcting from", "Correcting to", "$a_m$"
)
write_csv(metric_correction_factors, file.path(tablesDir, "metric_correction.csv"))

# latex version
tab_caption <- "Metric correction factor $a_m$ when correcting to either fixation count or fixation likelihood. These correction factors were used to make sure all dependent variables are comparable."
tab_label <- "tab:metric_correction"
print(
	xtable(
		metric_correction_factors, 
		caption = tab_caption, 
		label = tab_label, 
		align = "lllc",
		digits = 3
	), 
	include.rownames = FALSE,
	caption.placement = "top", 
    sanitize.text.function = function(x){x},
    file = file.path(tablesDir, "metric_correction.tex")
)

# -----
# Scatter plots
# ----

# showing the relationship between effect sizes expressed in 
# fixation likelihood and fixation count

# fix like to fix count plot
FLFCfig = data[, fix.count:N]
FLFCfig = FLFCfig[FLFCfig$fix.count != "NA" & FLFCfig$fix.like != "NA"]
FLFC_plot = genScatter(
	FLFCfig, fix.count, fix.like,
	"Fixation count", "Fixation likelihood"
)

# fix like to TFD plot
FLTFDfig = data[, fix.count:N]
FLTFDfig = FLTFDfig[FLTFDfig$TFD != "NA" & FLTFDfig$fix.like != "NA"]
FLTFD_plot = genScatter(
	FLTFDfig, TFD, fix.like,
	"Total fixation duration", "Fixation likelihood"
)

# fix count to TFD plot
FCTFDfig = data[, fix.count:N]
FCTFDfig = FCTFDfig[FCTFDfig$TFD != "NA" & FCTFDfig$fix.count != "NA"]
FCTFD_plot = genScatter(
	FCTFDfig, TFD, fix.count,
	"Total fixation duration", "Fixation count"
)

# dwell count to TFD plot
DCTFDfig = data[, fix.count:N]
DCTFDfig = DCTFDfig[DCTFDfig$TFD != "NA" & DCTFDfig$dwell.count != "NA"]
DCTFD_plot = genScatter(
	DCTFDfig, TFD, dwell.count,
	"Total fixation duration", "Dwell count"
)

# arrange plots in panel for manuscript
figure <- plot_grid(
	FLFC_plot, FLTFD_plot, FCTFD_plot, DCTFD_plot, 
	labels = LETTERS[1:4], 
	ncol = 2
)
filename <- file.path(figsDir, "metric_correction.pdf")
savePlots(figure, filename, fd_SI_2x2)
