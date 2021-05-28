# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# main analysis + publication bias analyses
# all analyses are performed on fixation count


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
data = as.data.table(read_csv(
	  file.path(dataDir, "data_effect_sizes_cleaned.csv")
))
grant = as.data.table(read_csv(
    file.path(dataDir, "data_grant.csv")
))
sample = as.data.table(read_csv(
  file.path(dataDir, "data_sample.csv")
))


# ----------------------------------------------------------------------
# Processing
# ----------------------------------------------------------------------

# compute corrected (unattenuated) effect sizes 
# (performed on z transformed values)
# for psychometric meta-analysis
data$yi.c.FL = FisherZInv(FisherZ(data$fix.like.m)/data$a_acc)
data$yi.c.FC = FisherZInv(FisherZ(data$fix.count.m)/data$a_acc)


# ----------------------------------------------------------------------
# Meta analyses - visual factors
# ----------------------------------------------------------------------

# -----
# psychometric meta-analysis and trim and fill analysis of 
# visual salience
# -----

saldata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Salience"], vtype="AV"))
saldata$vi.c = saldata$vi/saldata$a_acc^2 # compute corrected variances based on artefact multiplier
salres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=saldata, method="HS")
salresrobu = robust(salres, cluster=saldata$Authors)

# top 10% most precise analysis
saldata10 = saldata[saldata$vi.c < .01]
salTrim = rma(yi.c.FC, vi.c, weights=1/vi.c, data=saldata10, method="HS")

# -----
# psychometric meta-analysis and trim and fill analysis of 
# surface size
# -----

sizedata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Size"], vtype="AV"))
sizedata$vi.c = sizedata$vi/sizedata$a_acc^2
sizeres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=sizedata, method="HS")
sizeresrobu = robust(sizeres, cluster=sizedata$Authors)

# top 10% most precise analysis
sizedata10 = sizedata[sizedata$vi.c < .008]
sizeTrim = rma(yi.c.FC, vi.c, weights=1/vi.c, data=sizedata10, method="HS")

# -----
# psychometric meta-analysis and trim and fill analysis of 
# left v right position
# -----

LRdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "LR.position"], vtype="AV"))
LRdata$vi.c = LRdata$vi/LRdata$a_acc^2
LRres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=LRdata, method="HS")
LRresrobu = robust(LRres, cluster=LRdata$Authors)

# top 10% most precise analysis
LRdata10 = LRdata[LRdata$vi.c < .025]
LRTrim = rma(yi.c.FC, vi.c, weights=1/vi.c, data=LRdata10, method="HS")

# -----
# psychometric meta-analysis and trim and fill analysis of 
# centrality position
# -----

centerdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Center.position"], vtype="AV"))
centerdata$vi.c = centerdata$vi/centerdata$a_acc^2
centerres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=centerdata, method="HS")
centerresrobu = robust(centerres, cluster=centerdata$Authors)

# top 10% most precise analysis
centerdata10 = centerdata[centerdata$vi.c < .0144]
centerTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=centerdata10, method="HS")

# -----
# psychometric meta-analysis and trim and fill analysis of 
# set size
# -----

# preparing data sets
setdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize"], vtype="AV"))
setdata$vi.c = setdata$vi/setdata$a_acc^2
setdata_att = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize" & data$Alt.att == "attribute"], vtype="AV"))
setdata_att$vi.c = setdata_att$vi/setdata_att$a_acc^2
setdata_alt = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize" & data$Alt.att == "alternative"], vtype="AV"))
setdata_alt$vi.c = setdata_alt$vi/setdata_alt$a_acc^2

# main analyses
setres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata, method="HS")
setresrobu = robust(setres, cluster=setdata$Authors)

# psychometric meta-analysis of 
# set size moderator analysis - effect of attribute vs alternatives
# moderator RVE analysis relies on clubSandwich method see:
# https://cran.r-project.org/web/packages/clubSandwich/vignettes/meta-analysis-with-CRVE.html
setmod = rma.mv(yi.c.FC ~ Alt.att, V=vi.c, random = (~ 1|Authors), data=setdata) 
setmod = coef_test(setmod, vcov = "CR2")
sandwichModTest(setmod, "moderator_setsize.tex")
setmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata_att, method="HS") 
setmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata_alt, method="HS") 
setmod_attrobu = robust(setmod_att, cluster=setdata_att$Authors)
setmod_altrobu = robust(setmod_alt, cluster=setdata_alt$Authors)

# top 10% most precise analysis using a study from each mod group (alt vs att)
setdata10 = setdata[setdata$Study == "Zuschke 2020" | setdata$Study == "Hong et al. 2016"]
setTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata10, method="HS")

# top 10% most precise analysis using alt group
setdata10alt = setdata_alt[setdata_alt$vi.c < 0.017]
setmod_altTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata10alt, method="HS")

# top 10% most precise analysis using att group
setdata10att = setdata_att[setdata_att$vi.c < 0.013]
setmod_attTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata10att, method="HS")

# ----------------------------------------------------------------------
# Meta analyses - cognitive factors
# ----------------------------------------------------------------------

# -----
# psychometric meta-analysis and trim and fill analysis of 
# task instruction 
# -----

# preparing data
taskdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task"], vtype="AV"))
taskdata$vi.c = taskdata$vi/taskdata$a_acc^2
taskdata_att = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task" & data$Alt.att == "attribute"], vtype="AV"))
taskdata_att$vi.c = taskdata_att$vi/taskdata_att$a_acc^2
taskdata_alt = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task" & data$Alt.att == "alternative"], vtype="AV"))
taskdata_alt$vi.c = taskdata_alt$vi/taskdata_alt$a_acc^2

# main analyses
taskres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata, method="HS")
taskresrobu = robust(taskres, cluster=taskdata$Authors)

# psychometric meta-analysis and trim and fill analysis of 
# task instruction moderator analysis - effect of attribute vs alternative 
taskmod = rma.mv(yi.c.FC ~ Alt.att, V=vi.c, random = (~ 1|Authors), data=taskdata) 
taskmod = coef_test(taskmod, vcov = "CR2")
sandwichModTest(taskmod, "moderator_task.tex")
taskmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata_att, method="HS")
taskmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata_alt, method="HS")
taskmod_attrobu = robust(taskmod_att, cluster=taskdata_att$Authors)
taskmod_altrobu = robust(taskmod_alt, cluster=taskdata_alt$Authors)

# top 10% most precise analysis
taskdata10 = taskdata[taskdata$vi.c < .009]
taskTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata10, method="HS")

# top 10% most precise analysis using alt group
taskdata10alt = taskdata_alt[taskdata_alt$vi.c < .036]
taskmod_altTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata10alt, method="HS")

# top 10% most precise analysis using att group
taskdata10att = taskdata_att[taskdata_att$vi.c < .014]
taskmod_attTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata10att, method="HS")

# -----
# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing
# -----

# preparing data
prefdata = data[IV == "Pref.view", list(Authors = unique(Authors), 
                                        IV = unique(IV), 
                                        Alt.att = unique(Alt.att), 
                                        fix.count.m = mean(fix.count.m), 
                                        yi.c.FC = mean(yi.c.FC), 
                                        N = unique(N), 
                                        a_acc = unique(a_acc)), "Study"]
prefdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=prefdata, vtype="AV"))
prefdata$vi.c = prefdata$vi/prefdata$a_acc^2
prefdata_att = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view" & data$Alt.att == "attribute"], vtype="AV"))
prefdata_att$vi.c = prefdata_att$vi/prefdata_att$a_acc^2
prefdata_alt = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view" & data$Alt.att == "alternative"], vtype="AV"))
prefdata_alt$vi.c = prefdata_alt$vi/prefdata_alt$a_acc^2

# main analyses
prefres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata, method="HS")
prefresrobu = robust(prefres, cluster=prefdata$Authors)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing moderator analysis - effect of alternative vs attribute 
prefmod = rma.mv(yi.c.FC ~ Alt.att, V=vi.c, random = (~ 1|Authors), data=prefdata) 
prefmod = coef_test(prefmod, vcov = "CR2")
sandwichModTest(prefmod, "moderator_pref.tex")
prefmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata_att, method="HS")
prefmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata_alt, method="HS")
prefmod_attrobu = robust(prefmod_att, cluster=prefdata_att$Authors)
prefmod_altrobu = robust(prefmod_alt, cluster=prefdata_alt$Authors)

# top 10% most precise analysis using a study from each mod group (alt vs att)
prefdata10 = prefdata[prefdata$vi.c < .006]
prefTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata10, method="HS")

# top 10% most precise analysis using alt group
prefdata10alt = prefdata_alt[prefdata_alt$vi.c < .017]
prefmod_altTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata10alt, method="HS")

# top 10% most precise analysis using att group
prefdata10att = prefdata_att[prefdata_att$vi.c < .006]
prefmod_attTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata10att, method="HS")

# -----
# psychometric meta-analysis and trim and fill analysis of 
# -----

# choice-gaze effect moderator analysis - effect of inferential vs preferential choice
choicedata_mod = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Choice.bias"], vtype="AV"))
choicedata_mod$vi.c = choicedata_mod$vi/choicedata_mod$a_acc^2
choicedata_mod$inf_prefmod = ifelse(choicedata_mod$Research.Domain == "Risky Gamble", "preferential", # create inferential vs preferential task moderator variable
                         ifelse(choicedata_mod$Research.Domain == "pref. Consumer choice", "preferential",
                         ifelse(choicedata_mod$Research.Domain == "pref. Non-consumer choice", "preferential",
                         ifelse(choicedata_mod$Research.Domain == "inf. Consumer choice", "inferential",
                         ifelse(choicedata_mod$Research.Domain == "inf. Non-consumer choice", "inferential", NA)))))
choicemod = rma.mv(yi.c.FC ~ inf_prefmod, V=vi.c, random = (~ 1|Authors), data=choicedata_mod) 
choicemod = coef_test(choicemod, vcov = "CR2")
sandwichModTest(choicemod, "moderator_choicebias.tex")

# psychometric meta-analysis and trim and fill analysis of 
# choice-gaze effect main analysis (averaging studies with more effect sizes)
choicedata = choicedata_mod[, list(yi.c.FC = mean(yi.c.FC), vi.c = mean(vi.c), IV = unique(IV), N = unique(N), Authors = unique(Authors)), by = Study]
choiceres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=choicedata, method="HS")
choiceresrobu = robust(choiceres, cluster=choicedata$Authors)

# top 10% most precise analysis using a study from each mod group (alt vs att)
choicedata10 = choicedata[choicedata$vi.c < .014]
choiceTrim  = rma(yi.c.FC, vi.c, weights=1/vi.c, data=choicedata10, method="HS")

# ----------------------------------------------------------------------
# Main and Moderator results text for manuscript
# ----------------------------------------------------------------------

# rho range
result <- paste0(
  "$\\rho=", round(salresrobu$b, 3),"$ to $\\rho=", round(choiceresrobu$b, 3), "$"
)
cat(result, file = file.path(tablesDir, "rhorange.tex"))

# I2 range
result <- paste0(
  "$I^2=", round(salresrobu$I2, 3),"$ to $I^2=", round(prefresrobu$I2, 3), "$"
)
cat(result, file = file.path(tablesDir, "I2range.tex"))


# -----
# Table with main results for manuscript
# -----

extractMain <- function(name, mainres, data) {
    res <- c(
        name,
        paste0(mainres$k),
        paste0(sum(data$N)),
        paste0(round(coef(summary(mainres))$estimate, 2)),
        paste0(round(coef(summary(mainres))$se, 2)),
        paste0(round(coef(summary(mainres))$tval, 2)),
        ifelse(coef(summary(mainres))$pval < 0.001, "<0.001",
            paste0(round(coef(summary(mainres))$pval, 3))),
        paste0(round(coef(summary(mainres))$ci.lb, 2)),
        paste0(round(coef(summary(mainres))$ci.ub, 2)),
        round(mainres$I2, 2)
    )
    return(res)
}

extractTrim <- function(trimres, data) {
  res <- c(
    NA,
    paste0("(", trimres$k, ")"),
    paste0("(", sum(data$N), ")"),
    paste0("(", round(coef(summary(trimres))$estimate, 2), ")"),
    paste0("(", round(coef(summary(trimres))$se, 2), ")"),
    paste0("(", round(coef(summary(trimres))$zval, 2), ")"),
    ifelse(coef(summary(trimres))$pval < 0.001, "(<0.001)",
           paste0("(", round(coef(summary(trimres))$pval, 3), ")")),
    paste0("(", round(coef(summary(trimres))$ci.lb, 2), ")"),
    paste0("(", round(coef(summary(trimres))$ci.ub, 2), ")"),
    paste0("(", round(trimres$I2, 2), ")")
  )
  return(res)
}

mainresults = data.frame(rbind(
  c("\\textbf{Visual factors}", rep(NA, 9)),
  extractMain("Salience",salresrobu,saldata),
  extractTrim(salTrim,saldata10),
  extractMain("Surface size",sizeresrobu,sizedata),
  extractTrim(sizeTrim,sizedata10),
  extractMain("Left vs. right position",LRresrobu,LRdata),
  extractTrim(LRTrim,LRdata10),
  extractMain("Center position",centerresrobu,centerdata),
  extractTrim(centerTrim,centerdata10),
  extractMain("Set size",setresrobu,setdata),
  extractTrim(setTrim,setdata10),
  c("\\textbf{Cognitive factors}", rep(NA, 9)),
  extractMain("Task instructions",taskresrobu,taskdata),
  extractTrim(taskTrim,taskdata10),
  extractMain("Preferential viewing",prefresrobu,prefdata),
  extractTrim(prefTrim,prefdata10),
  extractMain("Choice-gaze effect",choiceresrobu,choicedata),
  extractTrim(choiceTrim,choicedata10)
), stringsAsFactors = FALSE)

# format naming and saving csv file
setnames(mainresults, c(1:10), c("Group","$k$","$N$","$\\rho$","SE","$t$","$p$","$\\textrm{CI}^{95}_{LL}$","$\\textrm{CI}^{95}_{UL}$","$I^2$"))
write_csv(mainresults, file.path(tablesDir, "main_results.csv"))

# latex version
tab_caption <- "Main results of the meta-analysis, divided into visual and cognitive factor groups, and individual factors within them. The most important values are the corrected effect size estimate, $\\rho$, and the associated heterogeneity, $I^2$. 
                Results of the Top10 analysis are in parentheses."
tab_label <- "tab:main_results"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.95\\textwidth}}",
           "{\\scriptsize{\\textit{Note.} $k$ = number of studies; $N$ = number of participants; $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}^{95}_{LL}$ = lower limit of the 95\\% confidence interval; $\\textrm{CI}^{95}_{UL}$ = upper limit of the 95\\% confidence interval, $I^2$ = within-group heterogeneity.}} \n")
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
        setmod_altrobu,setdata[setdata$Alt.att == "alternative",]),
    extractTrim(setmod_altTrim,setdata10alt),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        setmod_attrobu,setdata[setdata$Alt.att == "attribute",]),
    extractTrim(setmod_attTrim, setdata10att),
    
    c("\\textbf{Task instruction}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        taskmod_altrobu,taskdata[taskdata$Alt.att == "alternative",]),
    extractTrim(taskmod_altTrim, taskdata10alt),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        taskmod_attrobu,taskdata[taskdata$Alt.att == "attribute",]),
    extractTrim(taskmod_attTrim, taskdata10att),

    c("\\textbf{Preferential viewing}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        prefmod_altrobu,prefdata[prefdata$Alt.att == "alternative",]),
    extractTrim(prefmod_altTrim, prefdata10alt),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        prefmod_attrobu,prefdata[prefdata$Alt.att == "attribute",]),
    extractTrim(prefmod_attTrim, prefdata10att)
), stringsAsFactors = FALSE)

# format naming
setnames(modresults, c(1:10), c("Group","$k$","$N$","$\\rho$","SE","$Z$","$p$","$\\textrm{CI}^{95}_{LL}$","$\\textrm{CI}^{95}_{UL}$","$I^2$"))

# latex version
tab_caption <- "Moderator analysis results. The most important values are the corrected effect size estimate, $\\rho$, and the associated heterogeneity, $I^2$. 
                Results of the Top10 analysis are in parentheses."
tab_label <- "tab:mod_results"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.9\\textwidth}}",
           "{\\scriptsize{\\textit{Note.} $k$ = number of studies; $N$ = number of participants; $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}^{95}_{LL}$ = lower limit of the 95\\% confidence interval; $\\textrm{CI}^{95}_{LL}$ = upper limit of the 95\\% confidence interval, $I^2$ = within-group heterogeneity.}} \n")
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
    file = file.path(tablesDir, "mod_results.tex"),
  comment = FALSE
)


# -----
# Test of task instruction vs preferential viewing vs choice-gaze effect 
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
# Forrest plots 
# -----

# generating plots for each subgroup
salplot = genForest(salresrobu,saldata,"Salience","Salience", "yi.c.FC")
sizeplot = genForest(sizeresrobu,sizedata,"Size","Surface size", "yi.c.FC")
LRplot = genForest(LRresrobu,LRdata,"LR.position","Left vs. right position", "yi.c.FC")
centerplot = genForest(centerresrobu,centerdata,"Center.position","Center position", "yi.c.FC")
setplot = genForest(setresrobu,setdata,"Setsize", "Set size", "yi.c.FC")
taskplot = genForest(taskresrobu,taskdata,"Task", "Task instructions", "yi.c.FC")
prefplot = genForest(prefresrobu,prefdata,"Pref.view", "Preferential viewing", "yi.c.FC")
choiceplot = genForest(choiceresrobu,choicedata,"Choice.bias","Choice-gaze effect", "yi.c.FC")

# att vs alt
setplot_alt = genForest(setmod_altrobu,setdata[Alt.att == "alternative"],"Setsize", "Set size - alternative", "yi.c.FC")
setplot_att = genForest(setmod_attrobu,setdata[Alt.att == "attribute"],"Setsize", "Set size - attribute", "yi.c.FC")
taskplot_alt = genForest(taskmod_altrobu,taskdata[Alt.att == "alternative"],"Task", "Task instructions - alternative", "yi.c.FC")
taskplot_att = genForest(taskmod_attrobu,taskdata[Alt.att == "attribute"],"Task", "Task instructions - attribute", "yi.c.FC")
prefplot_alt = genForest(prefmod_altrobu,prefdata[Alt.att == "alternative"],"Pref.view", "Preferential viewing - alternative", "yi.c.FC")
prefplot_att = genForest(prefmod_attrobu,prefdata[Alt.att == "attribute"],"Pref.view", "Preferential viewing - attribute", "yi.c.FC") 

# arrange all forest plots in panel plot for manuscript - visual factors
forestpanel = plot_grid(
	salplot, centerplot,
	LRplot, sizeplot,
	setplot, 
    labels = LETTERS[1:5], 
    ncol = 2
)
filename <- file.path(figsDir, "forest_plots_visual.pdf")
savePlots(forestpanel, filename, fd_1c_3x2)

# arrange all forest plots in panel plot for manuscript - visual factors
forestpanel = plot_grid(
	taskplot, prefplot, choiceplot, 
    labels = LETTERS[1:3], 
    ncol = 1
)
filename <- file.path(figsDir, "forest_plots_cognitive.pdf")
savePlots(forestpanel, filename, fd_1c_3.5x1)

# arrange all forest plots in panel plot for manuscript - alt vs att ones
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
salfunnel = genFunnel(salresrobu,saldata,"yi.c.FC","Salience")
sizefunnel = genFunnel(sizeresrobu,sizedata,"yi.c.FC", "Surface size")
LRfunnel = genFunnel(LRresrobu,LRdata,"yi.c.FC", "Left vs. right position")
centerfunnel = genFunnel(centerresrobu,centerdata,"yi.c.FC", "Center position")
setfunnel = genFunnel(setresrobu,setdata,"yi.c.FC", "Set size")
setfunnel_alt = genFunnel(setmod_altrobu,setdata[setdata$Alt.att == "alternative"],"yi.c.FC", "Set size - alternative")
setfunnel_att = genFunnel(setmod_attrobu,setdata[setdata$Alt.att == "attribute"],"yi.c.FC", "Set size - attribute")
taskfunnel = genFunnel(taskresrobu,taskdata,"yi.c.FC", "Task instructions")
taskfunnel_alt = genFunnel(taskmod_altrobu,taskdata[taskdata$Alt.att == "alternative"],"yi.c.FC", "Task instructions - alternative")
taskfunnel_att = genFunnel(taskmod_attrobu,taskdata[taskdata$Alt.att == "attribute"],"yi.c.FC", "Task instructions - attribute")
preffunnel = genFunnel(prefresrobu,prefdata,"yi.c.FC", "Preferential viewing")
preffunnel_alt = genFunnel(prefmod_altrobu,prefdata[prefdata$Alt.att == "alternative"],"yi.c.FC", "Preferential viewing - alternative")
preffunnel_att = genFunnel(prefmod_attrobu,prefdata[prefdata$Alt.att == "attribute"],"yi.c.FC", "Preferential viewing - attribute")
choicefunnel = genFunnel(choiceresrobu,choicedata,"yi.c.FC", "Choice-gaze effect")

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


# ----------------------------------------------------------------------
# Publication bias analysis
# ----------------------------------------------------------------------

# the publication analysis performs three tests for the degree of ES inflation due to publication bias
# the first inflation factor is based on the difference in ES between studies with public grants vs no pub grants
# reasoning behind the test: https://www.sciencedirect.com/science/article/abs/pii/S0895435617301348
# second inflation factor is derived from the PEESE corrected estimate relative the fixed effects uncorrected estimate
# third inflation factor is from the uncorrected ES relative to the trim fill corrected ES

# we perform analyses on fisher z transformed ES
data$fcz = FisherZ(data$fix.count.m) # effect in z
data$sdz = 1/sqrt(data$N - 3) # sd in z
data$varz = data$sdz^2 # variance in z

# analyze if public grants are associated with smaller ES
grant = grant[, list(grant = unique(grant), public = unique(public)), Study]
data = merge(data, grant, by = "Study")
pb = rma(yi=fcz, vi=varz, mods = ~ public, data = data)
pb = robust(pb, cluster = data$Author)
publicFactor = round(pb$b[1] / (pb$b[1] + pb$b[2]), digits = 3) # inflation factor due to not having public grant
cat(paste0("$", publicFactor, "$"), file = file.path(tablesDir, "publicFactor.tex"))
cat(paste0("$Q_M(1,67)=", round(pb$QM, 3),"$, $p=", round(pb$QMp, 3), "$"), file = file.path(tablesDir, "publicSig.tex"))

# PET-PEESE test
FE = lm(fcz ~ 1, weights = 1/varz, data = data) # fixed effect estimate of ES
FE = coef_test(FE, vcov = "CR2", cluster = data$Authors)
onecoefTex(FE, "FE.tex") # save coefficient to tex

PET = lm(fcz ~ sdz + a_acc, weights = 1/varz, data = data) # PET test is sig therefore perfrom PEESE
PET = coef_test(PET, vcov = "CR2", cluster = data$Authors)
oneofmanycoefTex(PET, "PETintext.tex", 1)
oneofmanycoefTex(PET, "EGGERintext.tex", 2)

PEESE = lm(fcz ~ varz + a_acc, weights = 1/varz, data = data) # PEESE estimate
PEESE = coef_test(PEESE, vcov = "CR2", cluster = data$Authors)
oneofmanycoefTex(PEESE, "PEESEintext.tex", 1)

peeseFactor = round(FE[1] / PEESE[1,1], digits = 3) # inflation factor according to PEESE
cat(paste0("$", peeseFactor, "$"), file = file.path(tablesDir, "peeseFactor.tex"))

# check inflation factor based on trim fill results
trims = c(
    salresrobu$b / salTrim$b, # extract inflation factors for each subgroup separately
    sizeresrobu$b / sizeTrim$b,
    LRresrobu$b / LRTrim$b,
    centerresrobu$b / centerTrim$b,
    setresrobu$b / setTrim$b,
    prefresrobu$b / prefTrim$b,
    taskresrobu$b / taskTrim$b,
    choiceresrobu$b / choiceTrim$b
)
trimFactor = round(mean(trims), digits = 3) # average inflation factor
cat(paste0("$", trimFactor, "$"), file = file.path(tablesDir, "trimFactor.tex"))


# -----
# Table with PET-PEESE results for manuscript
# -----

PET = round(data.frame(PET, stringsAsFactors = FALSE), digits = 3)
PEESE = round(data.frame(PEESE, stringsAsFactors = FALSE), digits = 3)
PET = cbind(Parameter=c("Intercept", "$SD$", "$A$"), PET, stringsAsFactors = FALSE)
PEESE = cbind(Parameter=c("Intercept", "$Var$", "$A$"), PEESE, stringsAsFactors = FALSE)
setnames(PET, c(4:6), c("$t$","$df$","$p$"))
setnames(PEESE, c(4:6), c("$t$","$df$","$p$"))

# latex version PET
tab_caption <- "Precision-effect test (PET) of complete data"
tab_label <- "tab:PET"
print(
  xtable(
    PET, 
    caption = tab_caption, 
    label = tab_label,
    align = "lllcccc"
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(PET)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "PET.tex"),
  comment = FALSE
)

# latex version PEESE
tab_caption <- "Precision-effect estimate test (PEESE) of complete data"
tab_label <- "tab:PEESE"
print(
  xtable(
    PEESE, 
    caption = tab_caption, 
    label = tab_label,
    align = "lllcccc"
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(PEESE)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "PEESE.tex"),
  comment = FALSE
)

# latex version PET-PEESE
PETPEESE <- data.frame(rbind(
    c("PET", rep(NA, 5)),
    PET,
    c("PEESE", rep(NA, 5)),
    PEESE,
    stringsAsFactors = FALSE
), stringsAsFactors = FALSE)
setnames(PETPEESE, c(4:6), c("$t$","$df$","$p$"))
tab_caption <- "Publication bias analysis with precision-effect test (PET) and precision-effect estimate test (PEESE) of complete data. See \\textit{Methods} for details on the tests."
tab_label <- "tab:PET-PEESE"
print(
  xtable(
    PETPEESE, 
    caption = tab_caption, 
    label = tab_label,
    align = "llccccc"
  ), 
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(PETPEESE)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "PET-PEESE.tex"),
  comment = FALSE
)

# -------------------------------------------------------------------------------------------------------
# Table with raw data for appendix
# -------------------------------------------------------------------------------------------------------

overview = data[, c("Authors", "IV", "N","a_acc","fix.count.m","Eye.tracker","Research.Domain","Alt.att")]
overview$Eye.tracker = ifelse(
    overview$Eye.tracker == "Nihon-Kohden EEG-1100", "Nihon-Kohden", 
        ifelse(overview$Eye.tracker == "SMI model unknown (acc < .5)", "SMI unknown",
            ifelse(overview$Eye.tracker == "EyeLink 1000 Plus (acc < .5)", "EyeLink 1000*",
                ifelse(overview$Eye.tracker == "EyeLink 1000 (acc = .33)", "EyeLink 1000**",
    overview$Eye.tracker)))
)
overview$Research.Domain = ifelse(
    overview$Research.Domain == "pref. Consumer choice", "Pref con", 
        ifelse(overview$Research.Domain == "pref. Non-consumer choice", "Pref non-con",
            ifelse(overview$Research.Domain == "inf. Consumer choice", "Inf con",
                ifelse(overview$Research.Domain == "inf. Non-consumer choice", "Inf non-con",
    "Lotteries")))
)
overview$Alt.att = ifelse(
    overview$Alt.att == "alternative", "Alt", 
        ifelse(overview$Alt.att == "attribute", "Att",
    NA)
)
overview$IV = ifelse(
    overview$IV == "LR.position", "LvR",
        ifelse(overview$IV == "Center.position", "Center",
            ifelse(overview$IV == "Pref.view", "Pref",
                ifelse(overview$IV == "Choice.bias", "Choice",
                    ifelse(overview$IV == "Salience", "Sal",
                        ifelse(overview$IV == "Setsize", "Set size",
    overview$IV)))))
)
overview = overview[order(Authors),]
setnames(
    overview, 
    c("N","a_acc","fix.count.m","Eye.tracker","Research.Domain","Alt.att"), 
    c("$N$","$a_a$","$r$","Eye-tracker","Domain","Alt/Att")
)

tab_caption <- "Overview of individual effect sizes: IV = independent variable (LvR = Left vs. right position, Center = Center position, Sal = Salience, Pref = Preferential viewing, Choice = Choice-gaze effect, Task = Task instructions); $N$ = number of participants; $a_a$ = artifact multiplier; $r$ = attenuated effect size correlation expressed in the fixation count metric; Domain = research domain (Pref con = Preferential consumer choice, Pref non-con = Preferential non-consumer choice, Inf con = Inferential consumer choice, Inf non-con = Inferential non-consumer choice, Lotteries = Risky gambles); Alt/Att = Alternative or attribute manipulation."
tab_label <- "tab:overviewtable"
add.to.row <- list(pos = list(0), command = NULL)
command <- paste0("\\hline\n\\endhead\n","\\hline\n","\\multicolumn{", dim(overview)[2], "}{l}","{\\footnotesize Continued on next page}\n","\\endfoot\n","\\endlastfoot\n")
add.to.row$command <- command

otab = xtable(
  overview, 
  caption = tab_caption, 
  label = tab_label, 
  align = "cp{5cm}lccclll",
  digits = c(0,0,0,0,3,3,0,0,0)
)
# align(otab) = "cp{8cm}lccclll"
print(
  otab, 
  size = "\\footnotesize",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0),
  add.to.row = add.to.row,
  tabular.environment = "longtable",
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "overview_table.tex"),
  comment = FALSE
)

# -------------------------------------------------------------------------------------------------------
# Table with participant sample characteristics
# -------------------------------------------------------------------------------------------------------

sample = merge(sample, data, by = "Study", all.x = T)
sample = sample[, list(IV = unique(IV),
                       women = unique(percent_women),
                       age = unique(mean_age),
                       ethnicity = unique(ethnicity),
                       country = unique(country)), Study]

k = sample[, list(k = NROW(Study)), IV]
k = k[order(IV), ]

age = sample[is.na(age) == F, list(age = mean(age)), IV]
age = age[order(IV), ]
ageNA = sample[is.na(age) == T, list(ageNA = NROW(age)), IV]
ageNA = ageNA[order(IV),]

women = sample[is.na(women) == F, list(women = mean(women)), IV]
women = women[order(IV),]
womenNA = sample[is.na(women) == T, list(womenNA = NROW(women)), IV]
womenNA = womenNA[order(IV),]

ethnicity = sample[is.na(ethnicity) == F, list(count = NROW(Study)), by = c("IV", "ethnicity")]
ethnicity = dcast(ethnicity, IV ~ ethnicity)
ethnicity = ethnicity[order(ethnicity$IV),]
ethnicityNA = sample[is.na(ethnicity), list(ethnicityNA = NROW(ethnicity)), IV]
ethnicityNA = merge(ethnicity, ethnicityNA, by = "IV", all.x = T)
ethnicityNA = ethnicityNA[, c(1,5)]
ethnicityNA = ethnicityNA[order(ethnicityNA$IV),]

country = sample[is.na(country) == F, list(count = NROW(Study)), by = c("IV", "country")]
country = dcast(country, IV ~ country)
country = country[order(country$IV),]
countryNA = sample[is.na(country) == T, list(countryNA = NROW(country)), IV]
countryNA = merge(country, countryNA, by = "IV", all.x = T)
countryNA = countryNA[, c(1,19)]
countryNA = countryNA[order(countryNA$IV),]

sampletable = cbind(IV=k$IV,
                    '$k$' = k$k,
                    Age=NA,
                    "\\hspace{2mm}not reported"=ageNA$ageNA,
                    "\\hspace{2mm}mean"=round(age$age, digits = 2),
                    "Gender: female"=NA,
                    "\\hspace{2mm}not reported"=womenNA$womenNA,
                    "\\hspace{2mm}percent"=round(women$women, digits = 2),
                    Ethnicity=NA,
                    "\\hspace{2mm}not reported"=ethnicityNA$ethnicityNA,
                    ethnicity[,2:4],
                    Country=NA,
                    "\\hspace{2mm}not reported"=countryNA$countryNA,
                    country[,2:18])

#sampletable[is.na(sampletable)] = 0
rownames = colnames(sampletable)
sampletable = transpose(sampletable)
sampletable = sampletable[, c(5,7,3,1,6,8,4,2)]
setnames(sampletable, c(1:8), c("Salience", "Surface size", "Left vs. right position", "Center position", "Set size", "Task instructions", "Preferential viewing", "Choice-gaze effect"))
sampletable = cbind(" "=rownames[-1],sampletable[-1,])

rows = sampletable$' '[c(10:12,15:length(sampletable$' '))]
add = rep("\\hspace{2mm}", times = length(rows))
rows = paste0(add,rows)
Dimension = paste0(sampletable$' ')
rows = c(Dimension[c(1:9)], 
         rows[c(1:3)],
         Dimension[c(13:14)],
         rows[c(4:length(rows))])
row.names(sampletable) = NULL
sampletable$' ' = rows
  

# latex version Sample Table
tab_caption <- "Participant sample characteristics grouped by visual and cognitive factors"
tab_label <- "tab:sampleTable"
print(
  xtable(
    sampletable, 
    caption = tab_caption, 
    label = tab_label,
    align = "llp{.07\\linewidth}p{.07\\linewidth}p{.09\\linewidth}p{.09\\linewidth}p{.07\\linewidth}p{.08\\linewidth}p{.1\\linewidth}p{.07\\linewidth}"
    #align = "llcccccccc"
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(sampletable)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "sampleTable.tex"),
  comment = FALSE
)
