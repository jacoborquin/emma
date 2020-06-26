# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# main analysis + trim and fill analysis
# all analyses are performed on fixation count


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
data = as.data.table(read_csv(
	file.path(dataDir, "EMMA_ES_data_corrected.csv")
))


# ----------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------

# compute corrected (unattenuated) effect sizes 
# (performed on z transformed values)
# for psychometric meta-analysis
data$yi.c.FL = FisherZInv(FisherZ(data$fix.like.m)/data$a_acc)
data$yi.c.FC = FisherZInv(FisherZ(data$fix.count.m)/data$a_acc)


# -----
# Meta analyses - visual factors
# -----

# psychometric meta-analysis and trim and fill analysis of 
# visual salience
saldata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Salience"], vtype="AV")[,c(1:19,20)])
saldata$vi.c = saldata$vi/saldata$a_acc^2 # compute corrected variances based on artefact multiplier
salres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=saldata, method="HS")
salTrim = trimfill(salres)
# salresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=saldata)
# salTrimZ = trimfill(salresZ)
# salTrimZ$b = FisherZInv(salTrimZ$b)
# salTrimZ$ci.lb = FisherZInv(salTrimZ$ci.lb)
# salTrimZ$ci.ub = FisherZInv(salTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# surface size
sizedata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Size"], vtype="AV")[,c(1:19,20)])
sizedata$vi.c = sizedata$vi/sizedata$a_acc^2
sizeres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=sizedata, method="HS")
sizeTrim = trimfill(sizeres)
# sizeresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=sizedata)
# sizeTrimZ = trimfill(sizeresZ)
# sizeTrimZ$b = FisherZInv(sizeTrimZ$b)
# sizeTrimZ$ci.lb = FisherZInv(sizeTrimZ$ci.lb)
# sizeTrimZ$ci.ub = FisherZInv(sizeTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# left v right position
LRdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "LR.position"], vtype="AV")[,c(1:19,20)])
LRdata$vi.c = LRdata$vi/LRdata$a_acc^2
LRres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=LRdata, method="HS")
LRTrim = trimfill(LRres)
# LRresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=LRdata)
# LRTrimZ = trimfill(LRresZ)
# LRTrimZ$b = FisherZInv(LRTrimZ$b)
# LRTrimZ$ci.lb = FisherZInv(LRTrimZ$ci.lb)
# LRTrimZ$ci.ub = FisherZInv(LRTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# centrality position
centerdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Center.position"], vtype="AV")[,c(1:19,20)])
centerdata$vi.c = centerdata$vi/centerdata$a_acc^2
centerres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=centerdata, method="HS")
centerTrim = trimfill(centerres)
# centerresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=centerdata)
# centerTrimZ = trimfill(centerresZ)
# centerTrimZ$b = FisherZInv(centerTrimZ$b)
# centerTrimZ$ci.lb = FisherZInv(centerTrimZ$ci.lb)
# centerTrimZ$ci.ub = FisherZInv(centerTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# set size
setdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize"], vtype="AV")[,c(1:19,20)])
setdata$vi.c = setdata$vi/setdata$a_acc^2
setres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata, method="HS")
setTrim = trimfill(setres)
# setresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=setdata)
# setTrimZ = trimfill(setresZ)
# setTrimZ$b = FisherZInv(setTrimZ$b)
# setTrimZ$ci.lb = FisherZInv(setTrimZ$ci.lb)
# setTrimZ$ci.ub = FisherZInv(setTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# set size moderator analysis - effect of attribute vs alternatives
setmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=setdata, method="HS") 
setmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata[setdata$Alt.att == "attribute",], method="HS") 
setmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata[setdata$Alt.att == "alternative",], method="HS") 
setmod_attTrim = trimfill(setmod_att)
setmod_altTrim = trimfill(setmod_alt)

# setmod_attZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=setdata[setdata$Alt.att == "attribute",])
# setmod_attTrimZ = trimfill(setmod_attZ)
# setmod_attTrimZ$b = FisherZInv(setmod_attTrimZ$b)
# setmod_attTrimZ$ci.lb = FisherZInv(setmod_attTrimZ$ci.lb)
# setmod_attTrimZ$ci.ub = FisherZInv(setmod_attTrimZ$ci.ub)

# setmod_altZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=setdata[setdata$Alt.att == "alternative",])
# setmod_altTrimZ = trimfill(setmod_altZ)
# setmod_altTrimZ$b = FisherZInv(setmod_altTrimZ$b)
# setmod_altTrimZ$ci.lb = FisherZInv(setmod_altTrimZ$ci.lb)
# setmod_altTrimZ$ci.ub = FisherZInv(setmod_altTrimZ$ci.ub)


# -----
# Meta analyses - cognitive factors
# -----

# psychometric meta-analysis and trim and fill analysis of 
# task instruction 
taskdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task"], vtype="AV")[,c(1:19,20)])
taskdata$vi.c = taskdata$vi/taskdata$a_acc^2
taskres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata, method="HS")
taskTrim = trimfill(taskres)
# taskresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=taskdata)
# taskTrimZ = trimfill(taskresZ)
# taskTrimZ$b = FisherZInv(taskTrimZ$b)
# taskTrimZ$ci.lb = FisherZInv(taskTrimZ$ci.lb)
# taskTrimZ$ci.ub = FisherZInv(taskTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# task instruction moderator analysis - effect of attribute vs alternative 
taskmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=taskdata, method="HS")
taskmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata[taskdata$Alt.att == "alternative",], method="HS")
taskmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata[taskdata$Alt.att == "attribute",], method="HS")
taskmod_altTrim = trimfill(taskmod_alt)
taskmod_attTrim = trimfill(taskmod_att)

# taskmod_attZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=taskdata[taskdata$Alt.att == "attribute",])
# taskmod_attTrimZ = trimfill(taskmod_attZ)
# taskmod_attTrimZ$b = FisherZInv(taskmod_attTrimZ$b)
# taskmod_attTrimZ$ci.lb = FisherZInv(taskmod_attTrimZ$ci.lb)
# taskmod_attTrimZ$ci.ub = FisherZInv(taskmod_attTrimZ$ci.ub)

# taskmod_altZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=taskdata[taskdata$Alt.att == "alternative",])
# taskmod_altTrimZ = trimfill(taskmod_altZ)
# taskmod_altTrimZ$b = FisherZInv(taskmod_altTrimZ$b)
# taskmod_altTrimZ$ci.lb = FisherZInv(taskmod_altTrimZ$ci.lb)
# taskmod_altTrimZ$ci.ub = FisherZInv(taskmod_altTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing
prefdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view"], vtype="AV")[,c(1:19,20)])
prefdata$vi.c = prefdata$vi/prefdata$a_acc^2
prefres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata, method="HS")
prefTrim = trimfill(prefres)
# prefresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=prefdata)
# prefTrimZ = trimfill(prefresZ)
# prefTrimZ$b = FisherZInv(prefTrimZ$b)
# prefTrimZ$ci.lb = FisherZInv(prefTrimZ$ci.lb)
# prefTrimZ$ci.ub = FisherZInv(prefTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing moderator analysis - effect of alternative vs attribute 
prefmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=prefdata, method="HS")
prefmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata[prefdata$Alt.att == "alternative",], method="HS")
prefmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata[prefdata$Alt.att == "attribute",], method="HS")
prefmod_altTrim = trimfill(prefmod_alt)
prefmod_attTrim = trimfill(prefmod_att)

# prefmod_attZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=prefdata[prefdata$Alt.att == "attribute",])
# prefmod_attTrimZ = trimfill(prefmod_attZ)
# prefmod_attTrimZ$b = FisherZInv(prefmod_attTrimZ$b)
# prefmod_attTrimZ$ci.lb = FisherZInv(prefmod_attTrimZ$ci.lb)
# prefmod_attTrimZ$ci.ub = FisherZInv(prefmod_attTrimZ$ci.ub)

# prefmod_altZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=prefdata[prefdata$Alt.att == "alternative",])
# prefmod_altTrimZ = trimfill(prefmod_altZ)
# prefmod_altTrimZ$b = FisherZInv(prefmod_altTrimZ$b)
# prefmod_altTrimZ$ci.lb = FisherZInv(prefmod_altTrimZ$ci.lb)
# prefmod_altTrimZ$ci.ub = FisherZInv(prefmod_altTrimZ$ci.ub)

# psychometric meta-analysis and trim and fill analysis of 
# choice bias moderator analysis - effect of inferential vs preferential choice
choicedata_mod = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Choice.bias"], vtype="AV")[,c(1:19,20)])
choicedata_mod$vi.c = choicedata_mod$vi/choicedata_mod$a_acc^2
choicedata_mod$inf_prefmod = ifelse(choicedata_mod$Research.Domain == "Risky Gamble", "preferential", # create inferential vs preferential task moderator variable
                         ifelse(choicedata_mod$Research.Domain == "pref. Consumer choice", "preferential",
                         ifelse(choicedata_mod$Research.Domain == "pref. Non-consumer choice", "preferential",
                         ifelse(choicedata_mod$Research.Domain == "inf. Consumer choice", "inferential",
                         ifelse(choicedata_mod$Research.Domain == "inf. Non-consumer choice", "inferential", NA)))))
choicemod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ inf_prefmod, data=choicedata_mod, method="HS")

# psychometric meta-analysis and trim and fill analysis of 
# choice bias main analysis (averaging studies with more effect sizes)
choicedata = choicedata_mod[, list(yi.c.FC = mean(yi.c.FC), vi.c = mean(vi.c), IV = unique(IV), N = unique(N)), by = Study]
choiceres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=choicedata, method="HS")
choiceTrim = trimfill(choiceres)
# choiceresZ = rma(FisherZ(yi.c.FC), 1/sqrt(N-3), data=choicedata)
# choiceTrimZ = trimfill(choiceresZ)
# choiceTrimZ$b = FisherZInv(choiceTrimZ$b)
# choiceTrimZ$ci.lb = FisherZInv(choiceTrimZ$ci.lb)
# choiceTrimZ$ci.ub = FisherZInv(choiceTrimZ$ci.ub)

# -----
# Moderator results text for manuscript
# -----

result <- paste0(
	"$Q_M(1)=", round(setmod$QM, 3),"$, $p=", round(setmod$QMp, 3), "$"
)
cat(result, file = file.path(tablesDir, "moderator_setsize.tex"))

result <- paste0(
	"$Q_M(1)=", round(taskmod$QM, 3),"$, $p=", round(taskmod$QMp, 3), "$"
)
cat(result, file = file.path(tablesDir, "moderator_task.tex"))

result <- paste0(
	"$Q_M(1)=", round(prefmod$QM, 3),"$, $p=", round(prefmod$QMp, 3), "$"
)
cat(result, file = file.path(tablesDir, "moderator_pref.tex"))

result <- paste0(
	"$Q_M(1)=", round(choicemod$QM, 3),"$, $p=", round(choicemod$QMp, 3), "$"
)
cat(result, file = file.path(tablesDir, "moderator_choicebias.tex"))


# -----
# Table with main results for manuscript
# -----

extractMain <- function(name, mainres=centerres, data) {
    res <- c(
        name,
        paste0(mainres$k),
        paste0(sum(data$N)),
        paste0(round(coef(summary(mainres))$estimate, 2)),
        paste0(round(coef(summary(mainres))$se, 2)),
        paste0(round(coef(summary(mainres))$zval, 2)),
        ifelse(coef(summary(mainres))$pval < 0.001, "<0.001",
            paste0(round(coef(summary(mainres))$pval, 3))),
        paste0(round(coef(summary(mainres))$ci.lb, 2)),
        paste0(round(coef(summary(mainres))$ci.ub, 2)),
        round(mainres$I2, 2)
    )
    return(res)
}

extractTrim <- function(trimres) {
    res <- c(
        NA,
        paste0("(", trimres$k0, ")"),
        NA,
        paste0("(", round(coef(summary(trimres))$estimate, 2), ")"),
        paste0("(", round(coef(summary(trimres))$se, 2), ")"),
        paste0("(", round(coef(summary(trimres))$zval, 2), ")"),
        ifelse(coef(summary(trimres))$pval < 0.001, "(<0.001)",
            paste0("(", round(coef(summary(trimres))$pval, 3), ")")),
        paste0("(", round(coef(summary(trimres))$ci.lb, 2), ")"),
        paste0("(", round(coef(summary(trimres))$ci.ub, 2), ")"),
        NA
    )
    return(res)
}

mainresults = data.frame(rbind(
  c("\\textbf{Visual factors}", rep(NA, 9)),
  extractMain("Salience",salres,saldata),
  extractTrim(salTrim),
  extractMain("Surface size",sizeres,sizedata),
  extractTrim(sizeTrim),
  extractMain("Left vs right position",LRres,LRdata),
  extractTrim(LRTrim),
  extractMain("Center position",centerres,centerdata),
  extractTrim(centerTrim),
  extractMain("Set size",setres,setdata),
  extractTrim(setTrim),
  c("\\textbf{Cognitive factors}", rep(NA, 9)),
  extractMain("Task instructions",taskres,taskdata),
  extractTrim(taskTrim),
  extractMain("Preferential viewing",prefres,prefdata),
  extractTrim(prefTrim),
  extractMain("Choice bias",choiceres,choicedata),
  extractTrim(choiceTrim)
), stringsAsFactors = FALSE)

# format main results table, e.g. rounding, variable naming etc.
setnames(mainresults, c(1:10), c("Group","$k$","$N$","$\\rho$","SE","$Z$","$p$","$\\textrm{CI}_{95}$ LL","$\\textrm{CI}_{95}$ UL","$I^2$"))
write_csv(mainresults, file.path(tablesDir, "main_results.csv"))

# latex version
tab_caption <- "Main results of the meta-analysis, divided into visual and cognitive factor groups, and individual factors within them. The most important values are the corrected effect size estimate, $\\rho$, and the associated heterogeneity, $I^2$. Results of trim and fill analysis are in the parentesis."
tab_label <- "tab:main_results"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.95\\textwidth}}",
           "{\\scriptsize{\\textit{Note.} $k$ = number of studies (for trim and fill analysis number of imputed studies); $N$ = number of participants; $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}_{95}$ LL = lower limit of the 95\\% confidence interval; $\\textrm{CI}_{95}$ UL = upper limit of the 95\\% confidence interval, $I^2$ = within-group heterogeneity.}} \n")
print(
	xtable(
		mainresults, 
		caption = tab_caption, 
    label = tab_label,
    # align = "llp{0.03\\linewidth}p{0.05\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}",
    align = "llccccccccc"
    # digits = c(0,0,0,0,3,3,3,3,3,3,3)
  ), 
  size = "\\small",
	include.rownames = FALSE,
	caption.placement = "top", 
	hline.after = c(-1, 0),
	add.to.row = list(
		pos = list(nrow(mainresults)),
        command = tab_note
    ),
    sanitize.text.function = function(x){x},
    file = file.path(tablesDir, "main_results.tex")
)


# -----
# Table with moderator results for manuscript
# -----

modresults = data.frame(rbind(
    c("\\textbf{Set size}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        setmod_alt,setdata[setdata$Alt.att == "alternative",]),
    extractTrim(setmod_altTrim),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        setmod_att,setdata[setdata$Alt.att == "attribute",]),
    extractTrim(setmod_attTrim),
    
    c("\\textbf{Task instruction}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        taskmod_alt,taskdata[taskdata$Alt.att == "alternative",]),
    extractTrim(taskmod_altTrim),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        taskmod_att,taskdata[taskdata$Alt.att == "attribute",]),
    extractTrim(taskmod_attTrim),

    c("\\textbf{Preferential viewing}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        prefmod_alt,prefdata[prefdata$Alt.att == "alternative",]),
    extractTrim(prefmod_altTrim),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        prefmod_att,prefdata[prefdata$Alt.att == "attribute",]),
    extractTrim(prefmod_attTrim)

), stringsAsFactors = FALSE)

# format main results table, e.g. rounding, variable naming etc.
setnames(modresults, c(1:10), c("Group","$k$","$N$","$\\rho$","SE","$Z$","$p$","$\\textrm{CI}_{95}$ LL","$\\textrm{CI}_{95}$ UL","$I^2$"))
write_csv(modresults, file.path(tablesDir, "mod_results.csv"))

# latex version
tab_caption <- "Moderator analysis results. The most important values are the corrected effect size estimate, $\\rho$, and the associated heterogeneity, $I^2$. Results of trim and fill analysis are in the parentesis."
tab_label <- "tab:mod_results"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.9\\textwidth}}",
           "{\\scriptsize{\\textit{Note.} $k$ = number of studies (for trim and fill analysis number of imputed studies); $N$ = number of participants; $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}_{95}$ LL = lower limit of the 95\\% confidence interval; $\\textrm{CI}_{95}$ UL = upper limit of the 95\\% confidence interval, $I^2$ = within-group heterogeneity.}} \n")
print(
  xtable(
    modresults, 
    caption = tab_caption, 
    label = tab_label,
    align = "llccccccccc"
    # align = "llp{0.03\\linewidth}p{0.05\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}",
    # digits = c(0,0,0,0,3,3,3,3,3,3,3)
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0),
  add.to.row = list(
    pos = list(nrow(modresults)),
        command = tab_note
    ),
    sanitize.text.function = function(x){x},
    file = file.path(tablesDir, "mod_results.tex")
)


# -----
# Test of task instruction vs preferential viewing vs choice bias 
# -----

# perform wald test to compute z score difference between effect sizes
# then compute the significance level of the z score
# http://www.metafor-project.org/doku.php/tips:comp_two_independent_estimates
task_vs_pref_z = as.numeric(taskres$b - prefres$b)/sqrt(taskres$se^2 + prefres$se^2)
task_vs_pref_p = dnorm(task_vs_pref_z)
task_vs_choice_z = as.numeric(taskres$b - choiceres$b)/sqrt(taskres$se^2 + choiceres$se^2)
task_vs_choice_p = dnorm(task_vs_choice_z)
pref_vs_choice_z = as.numeric(prefres$b - choiceres$b)/sqrt(prefres$se^2 + choiceres$se^2)
pref_vs_choice_p = dnorm(pref_vs_choice_z)

# paste results for manuscript
result <- paste0(
	"$z=", round(task_vs_pref_z, 3), "$, $p=", round(task_vs_pref_p, 3), "$"
)
cat(result, file = file.path(tablesDir, "difftest_task_pref.tex"))

result <- paste0(
	"$z=", round(task_vs_choice_z, 3), "$, $p=", round(task_vs_choice_p, 3), "$"
)
cat(result, file = file.path(tablesDir, "difftest_task_choice.tex"))

result <- paste0(
	"$z=", round(pref_vs_choice_z, 3), "$, $p=", round(pref_vs_choice_p, 3), "$"
)
cat(result, file = file.path(tablesDir, "difftest_pref_choice.tex"))

# -----
# Table with raw data for appendix
# -----

overviewtabel = data[, c(2,4,9,14,15,1)]
tab_caption <- "Overview of studies"
tab_label <- "tab:overview"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.95\\textwidth}}",
                   "{\\scriptsize{\\textit{Note.} $k$ = number of studies (for trim and fill analysis number of imputed studies); $N$ = number of participants; $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}_{95}$ LL = lower limit of the 95\\% confidence interval; $\\textrm{CI}_{95}$ UL = upper limit of the 95\\% confidence interval, $I^2$ = within-group heterogeneity.}} \n")
print(
  xtable(
    overviewtabel, 
    caption = tab_caption, 
    label = tab_label,
    # align = "llp{0.03\\linewidth}p{0.05\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}p{0.07\\linewidth}",
    #align = "llccccccccc"
    # digits = c(0,0,0,0,3,3,3,3,3,3,3)
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0),
  add.to.row = list(
    pos = list(nrow(overviewtabel)),
    command = tab_note
  ),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "overviewtabel.tex")
)

# -----
# Forrest plots 
# -----

# generating plots for each subgroup
salplot = genForest(salres,saldata,"Salience","Salience", "yi.c.FC")
sizeplot = genForest(sizeres,sizedata,"Size","Surface size", "yi.c.FC")
LRplot = genForest(LRres,LRdata,"LR.position","Left vs right position", "yi.c.FC")
centerplot = genForest(centerres,centerdata,"Center.position","Center position", "yi.c.FC")
setplot = genForest(setres,setdata,"Setsize", "Set size", "yi.c.FC")
taskplot = genForest(taskres,taskdata,"Task", "Task instructions", "yi.c.FC")
prefplot = genForest(prefres,prefdata,"Pref.view", "Preferential viewing", "yi.c.FC")
choiceplot = genForest(choiceres,choicedata,"Choice.bias","Choice bias", "yi.c.FC")

# att vs alt
setplot_alt = genForest(setmod_alt,setdata[Alt.att == "alternative"],"Setsize", "Set size - alternative", "yi.c.FC")
setplot_att = genForest(setmod_att,setdata[Alt.att == "attribute"],"Setsize", "Set size - attribute", "yi.c.FC")
taskplot_alt = genForest(taskmod_alt,taskdata[Alt.att == "alternative"],"Task", "Task instructions - alternative", "yi.c.FC")
taskplot_att = genForest(taskmod_att,taskdata[Alt.att == "attribute"],"Task", "Task instructions - attribute", "yi.c.FC")
prefplot_alt = genForest(prefmod_alt,prefdata[Alt.att == "alternative"],"Pref.view", "Preferential viewing - alternative", "yi.c.FC")
prefplot_att = genForest(prefmod_att,prefdata[Alt.att == "attribute"],"Pref.view", "Preferential viewing - attribute", "yi.c.FC") 

# arange all forest plots in panel plot for manuscript - main ones
# BUleftside = plot_grid(salplot,LRplot,setplot_alt,labels = c("A","C","E"), ncol = 1, rel_heights = c(9,5,7))
# BUrightside = plot_grid(sizeplot,centerplot,setplot_att,labels = c("B","D","F"), ncol = 1, rel_heights = c(7,11,6))
# BUplot = plot_grid(BUleftside, BUrightside, ncol = 2)
# TDleftside = plot_grid(taskplot_alt,prefplot_alt,choiceplot,labels = c("G","I","K"), ncol = 1, rel_heights = c(8,6,11))
# TDrightside = plot_grid(taskplot_att,prefplot_att,labels = c("H","J","F"), nrow = 3, rel_heights = c(11,12,9))
# TDplot = plot_grid(TDleftside,TDrightside, ncol = 2)
# forestpanel = plot_grid(BUplot,TDplot, nrow = 2, rel_heights = c(2,3))
# filename <- file.path(figsDir, "forest_plots.pdf")
# savePlots(forestpanel, filename, fd_1c_3x2)

# arange all forest plots in panel plot for manuscript - visual factors
forestpanel = plot_grid(
	salplot, centerplot,
	LRplot, sizeplot,
	setplot, 
    labels = LETTERS[1:5], 
    ncol = 2
)
filename <- file.path(figsDir, "forest_plots_visual.pdf")
savePlots(forestpanel, filename, fd_1c_3x2)

# arange all forest plots in panel plot for manuscript - visual factors
forestpanel = plot_grid(
	taskplot, prefplot, choiceplot, 
    labels = LETTERS[1:3], 
    ncol = 1
)
filename <- file.path(figsDir, "forest_plots_cognitive.pdf")
savePlots(forestpanel, filename, fd_1c_3x1)

# arange all forest plots in panel plot for manuscript - alt vs att ones
forestpanel = plot_grid(
	setplot_alt, setplot_att,
	taskplot_alt, taskplot_att,
	prefplot_alt, prefplot_att, 
    labels = LETTERS[1:6], 
    ncol = 2
)
filename <- file.path(figsDir, "forest_plots_altatt.pdf")
savePlots(forestpanel, filename, fd_SI_3x2)


# -----
# Funnel plots 
# -----

# generating funnel plots for each main and subgroup
salfunnel = genFunnel(salres,saldata,"yi.c.FC","Salience")
sizefunnel = genFunnel(sizeres,sizedata,"yi.c.FC", "Surface size")
LRfunnel = genFunnel(LRres,LRdata,"yi.c.FC", "Left vs right position")
centerfunnel = genFunnel(centerres,centerdata,"yi.c.FC", "Center position")
setfunnel = genFunnel(setres,setdata,"yi.c.FC", "Set size")
setfunnel_alt = genFunnel(setmod_alt,setdata[setdata$Alt.att == "alternative"],"yi.c.FC", "Set size - altenrative")
setfunnel_att = genFunnel(setmod_att,setdata[setdata$Alt.att == "attribute"],"yi.c.FC", "Set size - attribute")
taskfunnel = genFunnel(taskres,taskdata,"yi.c.FC", "Task instructions")
taskfunnel_alt = genFunnel(taskmod_alt,taskdata[taskdata$Alt.att == "alternative"],"yi.c.FC", "Task instructions - alternative")
taskfunnel_att = genFunnel(taskmod_att,taskdata[taskdata$Alt.att == "attribute"],"yi.c.FC", "Task instructions - attribute")
preffunnel = genFunnel(prefres,prefdata,"yi.c.FC", "Preferential viewing")
preffunnel_alt = genFunnel(prefmod_alt,prefdata[prefdata$Alt.att == "alternative"],"yi.c.FC", "Preferential viewing - alternative")
preffunnel_att = genFunnel(prefmod_att,prefdata[prefdata$Alt.att == "attribute"],"yi.c.FC", "Preferential viewing - attribute")
choicefunnel = genFunnel(choiceres,choicedata,"yi.c.FC", "Choice bias")

# arrange funnel plots in panel plot - main ones
funnelpanel = plot_grid(
	salfunnel, sizefunnel, LRfunnel, centerfunnel,
	setfunnel, taskfunnel, preffunnel, choicefunnel, 
    labels = LETTERS[1:8], 
    ncol = 3
)
filename <- file.path(figsDir, "funnel_plots.pdf")
savePlots(funnelpanel, filename, fd_1c_3x2)

# arrange funnel plots in panel plot - alt vs att
funnelpanel = plot_grid(
	setfunnel_alt, setfunnel_att,
	taskfunnel_alt, taskfunnel_att,
	preffunnel_alt, preffunnel_att, 
    labels = LETTERS[1:6], 
    ncol = 2
)
filename <- file.path(figsDir, "funnel_plots_altatt.pdf")
savePlots(funnelpanel, filename, fd_SI_3x2)

