# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Loading data
# ----------------------------------------------------------------------

# housekeeping
rm(list = ls())

# import packages and functions
source("./scripts/utils.R")

# loading data
rawdata = as.data.table(read_excel(
  file.path(dataDir, "EMdata.xlsx")
))

# ----------------------------------------------------------------------
# corrections to raw data
# ----------------------------------------------------------------------

# selected wrong author name
rawdata$author[rawdata$author_confirm == "bagger 3"] = "Bagger 2016 Study 3"

# coding sheet to be deleted
rawdata = rawdata[rawdata$author_confirm != "glaholt 2009b 1"] # initial coding missed the research domain 

# correct coded value
rawdata$EM1_2[rawdata$author == "Gidloef et al. 2017"] = 0.8845253
rawdata$EM1_2[rawdata$author == "Hwang & Lee 2017"] = 798.7778
rawdata$EM1_7[rawdata$author == "Hwang & Lee 2017"] = 652.3333
rawdata$EM3_2[rawdata$author == "Peschel et al. 2019"] = 247.2184
rawdata$EM3_7[rawdata$author == "Peschel et al. 2019"] = 237.4305

# ----------------------------------------------------------------------
# Processing of conditional AOIs
# ----------------------------------------------------------------------

# select relevant columns (conditional AOIs)
data = rawdata[-1, c(20,27,29,31,36:45,
                     47,49,54:63,
                     65,67,72:81,
                     83,85,90:99,
                     101,103,108:117,
                     119,121,126:135,
                     137,139,144:153,
                     155,157,162:171)]

# melt one IV/DV at a time and append
iv1 = melt(data, id.vars = c(1,2,3,4), measure.vars = c(5:14))
setnames(iv1, c(3,4),c("IV","DV"))
iv2 = melt(data, id.vars = c(1,2,15,16), measure.vars = c(17:26))
setnames(iv2, c(3,4),c("IV","DV"))
iv3 = melt(data, id.vars = c(1,2,27,28), measure.vars = c(29:38))
setnames(iv3, c(3,4),c("IV","DV"))
iv4 = melt(data, id.vars = c(1,2,39,40), measure.vars = c(41:50))
setnames(iv4, c(3,4),c("IV","DV"))
iv5 = melt(data, id.vars = c(1,2,51,52), measure.vars = c(53:62))
setnames(iv5, c(3,4),c("IV","DV"))
iv6 = melt(data, id.vars = c(1,2,63,64), measure.vars = c(65:74))
setnames(iv6, c(3,4),c("IV","DV"))
iv7 = melt(data, id.vars = c(1,2,75,76), measure.vars = c(77:86))
setnames(iv7, c(3,4),c("IV","DV"))
iv8 = melt(data, id.vars = c(1,2,87,88), measure.vars = c(89:98))
setnames(iv8, c(3,4),c("IV","DV"))
EMlong = rbind(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8)

# create condition and AOI variables
EMlong$aoi = sapply(strsplit(as.character(EMlong$variable), "_"), `[`, 2)
EMlong$condition = ifelse(EMlong$aoi == 2 | EMlong$aoi == 3 | EMlong$aoi == 4 | EMlong$aoi == 5 | EMlong$aoi == 6, 0, 1)
EMlong$AOI = ifelse(EMlong$aoi == 2, 1,
                  ifelse(EMlong$aoi == 3, 2,
                         ifelse(EMlong$aoi == 4, 3,
                                ifelse(EMlong$aoi == 5, 4,
                                       ifelse(EMlong$aoi == 6, 5,
                                              ifelse(EMlong$aoi == 7, 1,
                                                     ifelse(EMlong$aoi == 8, 2,
                                                            ifelse(EMlong$aoi == 9, 3,
                                                                   ifelse(EMlong$aoi == 10, 4, 5)))))))))

# remove NA cells and temp variables
EMlong = EMlong[is.na(EMlong$value) == F, -c(5,7)]

# get IVs where missing
ESdata = as.data.table(read_excel(
  file.path(dataDir, "EMMA_ES_data.xlsx")
))

# rename variables and IV levels in ESdata
ESdata$IV = ifelse(ESdata$IV == "Pref.view", "Preferential viewing",
                   ifelse(ESdata$IV == "Center.position", "Central position",
                          ifelse(ESdata$IV == "Choice.bias", "Choice bias",
                                 ifelse(ESdata$IV == "LR.position", "Left vs right position",
                                        ifelse(ESdata$IV == "Task", "Task instruction",
                                               ifelse(ESdata$IV == "Size", "Surface size",
                                                      ifelse(ESdata$IV == "Setsize", "Set size", ESdata$IV)))))))

setnames(ESdata, old=c("Study", "Research.Domain"), new=c("author", "research_domain"))
ESdata$research_domain = ifelse(ESdata$research_domain == "pref. Consumer choice", "preferential consumer choice",
                                ifelse(ESdata$research_domain == "inf. Consumer choice", "inferential consumer choice",
                                       ifelse(ESdata$research_domain == "pref. Non-consumer choice", "preferential Non-consumer choice",
                                              ifelse(ESdata$research_domain == "inf. Non-consumer choice", "inferential Non-consumer choice", "Risky Gamble"))))

# combine multiple AOIs by taking average
EMlong$value = as.numeric(EMlong$value)
EMlong = EMlong[, list(value = mean(value)), by = c("author", "research_domain", "IV", "DV", "condition")]

# ------------------------------
# merge EM and ES data
# the complicated merging is due to papers being coded on different features and therefore
# they have to be merged on those features as well
# ------------------------------

# make subsets depending on the coding elements
EMlongNoIV_RD = EMlong[is.na(EMlong$IV) == T & is.na(EMlong$research_domain) == T]
EMlongIV = EMlong[is.na(EMlong$IV) == F & is.na(EMlong$research_domain) == T]
EMlongIV_RD = EMlong[is.na(EMlong$IV) == F & is.na(EMlong$research_domain) == F]

# merge subsets with ES data
EMlongIV_RD = merge(EMlongIV_RD, ESdata, by = c("author", "research_domain", "IV"))
EMlongIV = merge(EMlongIV, ESdata, by = c("author", "IV"))
EMlongNoIV_RD = merge(EMlongNoIV_RD, ESdata, by = c("author"), all.x = T)

# align subset col names etc
EMlongIV[, research_domain.x := NULL]
setnames(EMlongIV, "research_domain.y", "research_domain")
setcolorder(EMlongIV, c("Authors", "author", "IV", "DV", "research_domain", 
                        "Eye.tracker", "Alt.att", "N", "condition", "value", 
                        "fix.like", "fix.count", "TFD", "dwell.count"))
EMlongNoIV_RD[, c("research_domain.x", "IV.x") := NULL]
setnames(EMlongNoIV_RD, c("IV.y", "research_domain.y"), c("IV", "research_domain"))
setcolorder(EMlongNoIV_RD, c("Authors", "author", "IV", "DV", "research_domain", 
                        "Eye.tracker", "Alt.att", "N", "condition", "value", 
                        "fix.like", "fix.count", "TFD", "dwell.count"))
setcolorder(EMlongIV_RD, c("Authors", "author", "IV", "DV", "research_domain", 
                        "Eye.tracker", "Alt.att", "N", "condition", "value", 
                        "fix.like", "fix.count", "TFD", "dwell.count"))

EMlong = rbind(EMlongNoIV_RD, EMlongIV_RD, EMlongIV)

# rename TDT level
EMlong$DV = ifelse(EMlong$DV == "Total dwell time (total fixation duration)", "Total dwell time", EMlong$DV)

# ----------------------------------------------------------------------
# scatter plots
# ----------------------------------------------------------------------

# arranging data with condition 1 and 2 in separate columns
value.x = EMlong$value[EMlong$condition == 0]
value.y = EMlong$value[EMlong$condition == 1]
temp3 = unique(EMlong[, c("value", "condition") := NULL])
EMwide = cbind(temp3, cbind(value.x,value.y))

# aligning x,y data points so that the effect lands on the lower half of the scatter plot
EMwide$maxval = ifelse(EMwide$value.x > EMwide$value.y, EMwide$value.x, EMwide$value.y)
EMwide$minval = ifelse(EMwide$value.x > EMwide$value.y, EMwide$value.y, EMwide$value.x)
EMwide = EMwide[, - c("value.x", "value.y")]

# plot function
# EMplot = function(data, DV){
#   if (DV == "Fixation likelihood") {
#     xymax = 1 
#   } else {
#     xymax = max(data$maxval)}
#   ggplot(data, aes(maxval,minval))+
#     geom_point()+
#     geom_segment(aes(x=0,y=0,xend=xymax,yend=xymax), linetype = "dashed")+
#     xlab(paste(DV, "Condition 1"))+
#     ylab(paste(DV, "Condition 2"))+
#     facet_grid(~IV)+
#     mytheme
# }
# 
# flplot = EMplot(EMwide[EMwide$DV == "Fixation likelihood"], DV = "Fixation likelihood")
# fcplot = EMplot(EMwide[EMwide$DV == "Fixation count"], DV = "Fixation count")
# TDTplot = EMplot(EMwide[EMwide$DV == "Total dwell time"], DV = "Total dwell time")
# dcplot = EMplot(EMwide[EMwide$DV == "Dwell count"], DV = "Dwell count")
# 
# plot_grid(flplot, fcplot, TDTplot, dcplot, nrow = 4)

# ----------------------------------------------------------------------
# Examine relation between EM and ES
# ----------------------------------------------------------------------

# compute transformations of EM differences across conditions
# fixation likelihood is transformed into logits, fix count and TDT are log transformed
EMwide$diff = EMwide$maxval - EMwide$minval
EMwide$diffPercent = (EMwide$maxval - EMwide$minval)/EMwide$minval
EMwide$diffLog = log(EMwide$maxval) - log(EMwide$minval)
EMwide$diffLogit = log(EMwide$maxval/(1 - EMwide$maxval)) - log(EMwide$minval/(1 - EMwide$minval)) 
EMwide$diffLogit = ifelse(EMwide$diffLogit == Inf, max(EMwide$diffLogit[is.na(EMwide$diffLogit) == F & EMwide$diffLogit != Inf]), EMwide$diffLogit)

# make vectors for correlations between ES and EM data
fldiff = EMwide$diffLogit[is.na(EMwide$fix.like) == F & EMwide$DV == "Fixation likelihood"]
flES = EMwide$fix.like[is.na(EMwide$fix.like) == F & EMwide$DV == "Fixation likelihood"]
fcdiff = EMwide$diffLog[is.na(EMwide$fix.count) == F & EMwide$DV == "Fixation count"]
fcES = EMwide$fix.count[is.na(EMwide$fix.count) == F & EMwide$DV == "Fixation count"]
tdtdiff = EMwide$diffLog[is.na(EMwide$TFD) == F & EMwide$DV == "Total dwell time" & EMwide$author != "Orquin et al. 2019a Study 1"]
tdtES = EMwide$TFD[is.na(EMwide$TFD) == F & EMwide$DV == "Total dwell time" & EMwide$author != "Orquin et al. 2019a Study 1"]

# correlations between EM and ES data for FL, FC, and TDT
cor(fldiff,flES)
cor(fcdiff,fcES)
cor(tdtdiff,tdtES)

# plot with transformations
fldata = EMwide[is.na(EMwide$fix.like) == F & EMwide$DV == "Fixation likelihood"]
flfig = ggplot(fldata, aes(fix.like, diffLogit))+
  geom_point()+
  #geom_text(aes(label=author),hjust=-.01, vjust=-.2, size = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Effect size (r)")+
  ylab("Logit difference in fixation likelihood")+
  mytheme

fcdata = EMwide[is.na(EMwide$fix.count) == F & EMwide$DV == "Fixation count"]
fcfig = ggplot(fcdata, aes(fix.count, diffLog))+
  geom_point()+
  #geom_text(aes(label=author),hjust=-.01, vjust=-.2, size = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Effect size (r)")+
  ylab("Log difference in fixation count")+
  mytheme

tdtdata = EMwide[is.na(EMwide$TFD) == F & EMwide$DV == "Total dwell time"]
tdtfig = ggplot(tdtdata, aes(TFD, diffLog))+
  geom_point()+
  #geom_text(aes(label=author),hjust=-.01, vjust=-.2, size = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Effect size (r)")+
  ylab("Log difference in total dwell time")+
  mytheme

# ----------------------------------------------------------------------
# Processing of overall AOIs
# ----------------------------------------------------------------------

# select relevant columns (overall AOIs)
overalldata = rawdata[-1, c(20,29,31,35,
                            47,49,53,
                            65,67,71,
                            83,85,89,
                            101,103,107,
                            119,121,125,
                            137,139,143,
                            155,157,161)]

# melt one IV/DV at a time and append
iv1 = melt(overalldata, id.vars = c(1,2,3), measure.vars = c(4))
setnames(iv1, c(2,3),c("IV","DV"))
iv2 = melt(overalldata, id.vars = c(1,5,6), measure.vars = c(7))
setnames(iv2, c(2,3),c("IV","DV"))
iv3 = melt(overalldata, id.vars = c(1,8,9), measure.vars = c(10))
setnames(iv3, c(2,3),c("IV","DV"))
iv4 = melt(overalldata, id.vars = c(1,11,12), measure.vars = c(13))
setnames(iv4, c(2,3),c("IV","DV"))
iv5 = melt(overalldata, id.vars = c(1,14,15), measure.vars = c(16))
setnames(iv5, c(2,3),c("IV","DV"))
iv6 = melt(overalldata, id.vars = c(1,17,18), measure.vars = c(19))
setnames(iv6, c(2,3),c("IV","DV"))
iv7 = melt(overalldata, id.vars = c(1,20,21), measure.vars = c(22))
setnames(iv7, c(2,3),c("IV","DV"))
iv8 = melt(overalldata, id.vars = c(1,23,24), measure.vars = c(25))
setnames(iv8, c(2,3),c("IV","DV"))
EMoverall = rbind(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8)

# remove NA cells and temp variables
EMoverall = EMoverall[is.na(EMoverall$value) == F, -c(4)]
EMoverall$value = as.numeric(EMoverall$value)
EMoverall = EMoverall[, list(value = unique(value)), by = c("author", "DV")]

# compute overall EMs from conditional AOIs
EMwide$overall = (EMwide$maxval + EMwide$minval)/2
EMoverallComp = EMwide[, list(value = mean(overall)), by = c("author", "DV")]
EMoverall = merge(EMoverall, EMoverallComp, by = c("author", "DV"), all = T)
EMoverall$average = ifelse(is.na(EMoverall$value.x), EMoverall$value.y, EMoverall$value.x)
EMoverall[, c(3:4) := NULL]

# figures with overall distributions
flavg = EMoverall[DV == "Fixation likelihood"]
flavgfig = ggplot(flavg, aes(DV, average))+
  geom_violin()+
  geom_boxplot(width = .3)+
  xlab("")+
  ylab("Average fixation likelihood")+
  mytheme+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

fcavg = EMoverall[DV == "Fixation count"]
fcavgfig = ggplot(fcavg, aes(DV, average))+
  geom_violin()+
  geom_boxplot(width = .3)+
  xlab("")+
  ylab("Average fixation count")+
  mytheme+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

tdtavg = EMoverall[DV == "Total dwell time"]
tdtavgfig = ggplot(tdtavg, aes(DV, average))+
  geom_violin()+
  geom_boxplot(width = .3)+
  xlab("")+
  ylab("Average total dwell time")+
  mytheme+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

EMtoES = plot_grid(flavgfig, fcavgfig, tdtavgfig,
                   flfig, fcfig, tdtfig, labels = c("A", "B", "C", "D", "E", "F"), ncol = 3)

filename <- file.path(figsDir, "EMtoES.pdf")
savePlots(EMtoES, filename, c(1.5*height, 2*width))

# ----------------------------------------------------------------------
# ES to EM conversion table
# ----------------------------------------------------------------------

# overall EM means
flmean = mean(EMoverall$average[EMoverall$DV == "Fixation likelihood"])
fcmean = mean(EMoverall$average[EMoverall$DV == "Fixation count"])
tdtmean = mean(EMoverall$average[EMoverall$DV == "Total dwell time"])

# overall EM logits and logs
flmean = mean(EMoverall$average[EMoverall$DV == "Fixation likelihood"])
fcmean = mean(EMoverall$average[EMoverall$DV == "Fixation count"])
tdtmean = mean(EMoverall$average[EMoverall$DV == "Total dwell time"])
flmeanLogit = log(flmean/(1 - flmean))
fcmeanLog = log(fcmean)
tdtmeanLog = log(tdtmean)

# ES to EM conversion models
flmodel = lmer(diffLogit ~ fix.like + (1|author), data = EMwide[EMwide$DV == "Fixation likelihood"])
fcmodel = lmer(diffLog ~ fix.count + (1|author), data = EMwide[EMwide$DV == "Fixation count"])
tdtmodel = lmer(diffLog ~ TFD + (1|author), data = EMwide[EMwide$DV == "Total dwell time"])

# print .tex model output
twocoefTex(flmodel, "FLtoLogitModel.tex")
twocoefTex(fcmodel, "FCtoLogModel.tex")
twocoefTex(tdtmodel, "TDTtoLogModel.tex")

# get main results expressed in correlations
EMresults = as.data.table(read.csv(
  file.path(tablesDir, "main_results.csv")
))

EMresults = EMresults[is.na(EMresults$Group) == F & is.na(EMresults$X..rho.) == F, c(1,4)]
EMresults$rho = as.numeric(as.character(EMresults$X..rho.))
EMresults[, "X..rho." := NULL]

# transformation functions
logitTransform = function(IV, model){
  ES = EMresults$rho[EMresults$Group == IV]
  ESasLogit = summary(model)$coef[1,1] + summary(model)$coef[2,1] * ES
  upperLogit = ESasLogit + flmeanLogit
  EMupper = exp(upperLogit)/(exp(upperLogit) + 1)
  EMupper
}

logTransform = function(IV, model, mean){
  ES = EMresults$rho[EMresults$Group == IV]
  ESasLog = summary(model)$coef[1,1] + summary(model)$coef[2,1] * ES
  upperLog = ESasLog + mean
  EMupper = exp(upperLog)
  EMupper
}

FL = c(logitTransform("Salience", flmodel),
logitTransform("Surface size", flmodel),
logitTransform("Left vs right position", flmodel),
logitTransform("Center position", flmodel),
logitTransform("Set size", flmodel),
logitTransform("Task instructions", flmodel),
logitTransform("Preferential viewing", flmodel),
logitTransform("Choice bias", flmodel))

FC = c(logTransform("Salience", fcmodel, fcmeanLog),
logTransform("Surface size", fcmodel, fcmeanLog),
logTransform("Left vs right position", fcmodel, fcmeanLog),
logTransform("Center position", fcmodel, fcmeanLog),
logTransform("Set size", fcmodel, fcmeanLog),
logTransform("Task instructions", fcmodel, fcmeanLog),
logTransform("Preferential viewing", fcmodel, fcmeanLog),
logTransform("Choice bias", fcmodel, fcmeanLog))

TDT = c(logTransform("Salience", tdtmodel, tdtmeanLog),
logTransform("Surface size", tdtmodel, tdtmeanLog),
logTransform("Left vs right position", tdtmodel, tdtmeanLog),
logTransform("Center position", tdtmodel, tdtmeanLog),
logTransform("Set size", tdtmodel, tdtmeanLog),
logTransform("Task instructions", tdtmodel, tdtmeanLog),
logTransform("Preferential viewing", tdtmodel, tdtmeanLog),
logTransform("Choice bias", tdtmodel, tdtmeanLog))

EMresults$'FL mean' = round(flmean, digits = 3)
EMresults$'FL increase' = round(FL, digits = 3)
EMresults$'FC mean' = round(fcmean, digits = 3)
EMresults$'FC increase' = round(FC, digits = 3)
EMresults$'TDT mean' = round(tdtmean, digits = 3)
EMresults$'TDT increase' = round(TDT, digits = 3)

# latex version
tab_caption <- "Main effects expressed as absolute changes in the fixation likelihood, fixation count, and total dwell time."
tab_label <- "tab:em_results"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.95\\textwidth}}",
                   "{\\scriptsize{\\textit{Note.} FL = fixation likelihood, FC = fixation count, TDT = total dwell time}} \n")
print(
  xtable(
    EMresults, 
    caption = tab_caption, 
    label = tab_label,
    # align = "llp{0.03\\linewidth}p{0.05\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}",
    align = "lp{2cm}ccccccc"
    # digits = c(0,0,0,0,3,3,3,3,3,3,3)
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0),
  add.to.row = list(
    pos = list(nrow(EMresults)),
    command = tab_note
  ),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "em_results.tex")
)

