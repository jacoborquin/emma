# ----------------------------------------------------------------------
# Information
# ----------------------------------------------------------------------

# main analysis + trim and fill analysis

# analyses of bottom up factors (salience, surface size, left vs right 
# position, central position, and set size)
# are performed using fixation likelihood. Analyses of top down factors 
# are performed on fixation count


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
# Meta analyses
# -----

# psychometric meta-analysis and trim and fill analysis of 
# visual salience
saldata = data.table(escalc(measure="COR", ri=fix.like.m, ni=N, data=data[data$IV == "Salience"], vtype="AV")[,c(1:17,19)])
saldata$vi.c = saldata$vi/saldata$a_acc^2 # compute corrected variances based on artefact multiplier
salres = rma(yi.c.FL, vi.c, weights=1/vi.c, data=saldata, method="HS")
salTrim = trimfill(salres)

# psychometric meta-analysis and trim and fill analysis of 
# surface size
sizedata = data.table(escalc(measure="COR", ri=fix.like.m, ni=N, data=data[data$IV == "Size"], vtype="AV")[,c(1:17,19)])
sizedata$vi.c = sizedata$vi/sizedata$a_acc^2
sizeres = rma(yi.c.FL, vi.c, weights=1/vi.c, data=sizedata, method="HS")
sizeTrim = trimfill(sizeres)

# psychometric meta-analysis and trim and fill analysis of 
# left v right position
LRdata = data.table(escalc(measure="COR", ri=fix.like.m, ni=N, data=data[data$IV == "LR.position"], vtype="AV")[,c(1:17,19)])
LRdata$vi.c = LRdata$vi/LRdata$a_acc^2
LRres = rma(yi.c.FL, vi.c, weights=1/vi.c, data=LRdata, method="HS")
LRTrim = trimfill(LRres)

# psychometric meta-analysis and trim and fill analysis of 
# centrality position
centerdata = data.table(escalc(measure="COR", ri=fix.like.m, ni=N, data=data[data$IV == "Center.position"], vtype="AV")[,c(1:17,19)])
centerdata$vi.c = centerdata$vi/centerdata$a_acc^2
centerres = rma(yi.c.FL, vi.c, weights=1/vi.c, data=centerdata, method="HS")
centerTrim = trimfill(centerres)

# psychometric meta-analysis and trim and fill analysis of 
# set size
setdata = data.table(escalc(measure="COR", ri=fix.like.m, ni=N, data=data[data$IV == "Setsize"], vtype="AV")[,c(1:17,19)])
setdata$vi.c = setdata$vi/setdata$a_acc^2
setres = rma(yi.c.FL, vi.c, weights=1/vi.c, data=setdata, method="HS")
setTrim = trimfill(setres)

# psychometric meta-analysis and trim and fill analysis of 
# set size moderator analysis - effect of attribute vs alternatives
setmod = rma(yi.c.FL, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=setdata, method="HS") 
setmod_att = rma(yi.c.FL, vi.c, weights=1/vi.c, data=setdata[setdata$Alt.att == "attribute",], method="HS") 
setmod_alt = rma(yi.c.FL, vi.c, weights=1/vi.c, data=setdata[setdata$Alt.att == "alternative",], method="HS") 
setmod_attTrim = trimfill(setmod_att)
setmod_altTrim = trimfill(setmod_alt)

# psychometric meta-analysis and trim and fill analysis of 
# task instruction 
taskdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task"], vtype="AV")[,c(1:17,19)])
taskdata$vi.c = taskdata$vi/taskdata$a_acc^2
taskres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata, method="HS")
taskTrim = trimfill(taskres)

# psychometric meta-analysis and trim and fill analysis of 
# task instruction moderator analysis - effect of attribute vs alternative 
taskmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=taskdata, method="HS")
taskmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata[taskdata$Alt.att == "alternative",], method="HS")
taskmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata[taskdata$Alt.att == "attribute",], method="HS")
taskmod_altTrim = trimfill(taskmod_alt)
taskmod_attTrim = trimfill(taskmod_att)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing
prefdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view"], vtype="AV")[,c(1:17,19)])
prefdata$vi.c = prefdata$vi/prefdata$a_acc^2
prefres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata, method="HS")
prefTrim = trimfill(prefres)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing moderator analysis - effect of alternative vs attribute 
prefmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=prefdata, method="HS")
prefmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata[prefdata$Alt.att == "alternative",], method="HS")
prefmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata[prefdata$Alt.att == "attribute",], method="HS")
prefmod_altTrim = trimfill(prefmod_alt)
prefmod_attTrim = trimfill(prefmod_att)

# psychometric meta-analysis and trim and fill analysis of 
# choice bias moderator analysis - effect of inferential vs preferential choice
choicedata_mod = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Choice.bias"], vtype="AV")[,c(1:17,19)])
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

mainresults = data.frame(rbind(
  c("\\textbf{Visual factors}",NA,NA,NA,NA),
  c("Salience",salres$k,sum(saldata$N),coef(summary(salres)),salres$I2),
  c("Surface size",sizeres$k,sum(sizedata$N),coef(summary(sizeres)),sizeres$I2),
  c("Left vs right position",LRres$k,sum(LRdata$N),coef(summary(LRres)),LRres$I2),
  c("Center position",centerres$k,sum(centerdata$N),coef(summary(centerres)),centerres$I2),
  c("Set size",setres$k,sum(setdata$N),coef(summary(setres)),setres$I2),
  c("\\hspace{2mm}\\textit{Alternative}",setmod_alt$k,sum(setdata$N[setdata$Alt.att == "alternative"]),coef(summary(setmod_alt)),setmod_alt$I2),
  c("\\hspace{2mm}\\textit{Attribute}",setmod_att$k,sum(setdata$N[setdata$Alt.att == "attribute"]),coef(summary(setmod_att)),setmod_att$I2),
  c("\\textbf{Cognitive factors}",NA,NA,NA,NA),
  c("Task instruction",taskres$k,sum(taskdata$N),coef(summary(taskres)),taskres$I2),
  c("\\hspace{2mm}\\textit{Alternative}",taskmod_alt$k,sum(taskdata$N[taskdata$Alt.att == "alternative"]),coef(summary(taskmod_alt)),taskmod_alt$I2),
  c("\\hspace{2mm}\\textit{Attribute}",taskmod_att$k,sum(taskdata$N[taskdata$Alt.att == "attribute"]),coef(summary(taskmod_att)),taskmod_att$I2),
  c("Preferential viewing",prefres$k,sum(prefdata$N),coef(summary(prefres)),prefres$I2),
  c("\\hspace{2mm}\\textit{Alternative}",prefmod_alt$k,sum(prefdata$N[prefdata$Alt.att == "alternative"]),coef(summary(prefmod_alt)),prefmod_alt$I2),
  c("\\hspace{2mm}\\textit{Attribute}",prefmod_att$k,sum(prefdata$N[prefdata$Alt.att == "attribute"]),coef(summary(prefmod_att)),prefmod_att$I2),
  c("Choice bias",choiceres$k,sum(choicedata$N),coef(summary(choiceres)),choiceres$I2)
))

# format main results table, e.g. rounding, variable naming etc.
mainresults[2:10] <- sapply(mainresults[2:10], as.numeric)
mainresults = mainresults %>% mutate_at(vars(4:10), round, 3)
setnames(mainresults, c(1:10), c("Group","$k$","$N$","$\\rho$","SE","$Z$","$p$","$\\textrm{CI}_{95}$ LL","$\\textrm{CI}_{95}$ UL","$I^2$"))
mainresults$Group = as.character(mainresults$Group)
write_csv(mainresults, file.path(tablesDir, "main_results.csv"))

# latex version
tab_caption <- "Main results of the meta-analysis, divided into visual and cognitive factor groups, and individual factors within them, including sub-factors used in the moderator analyses. The most important values are the corrected effect size estimate, $\\rho$, and the associated heterogeneity, $I^2$"
tab_label <- "tab:main_results"
tab_note <- paste0("\\hline \n \\multicolumn{10}{p{0.9\\textwidth}}",
           "{\\scriptsize{\\textit{Note.} $k$ = number of studies; $N$ = number of participants; $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}_{95}$ LL = lower limit of the 95\\% confidence interval; $\\textrm{CI}_{95}$ UL = upper limit of the 95\\% confidence interval, $I^2$ = within-group heterogeneity. Italicized groups are moderator subgroups.}} \n")
print(
	xtable(
		mainresults, 
		caption = tab_caption, 
		label = tab_label,
		align = "llccccccccc",
		digits = c(0,0,0,0,3,3,3,3,3,3,3)
	), 
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
# Table with trim and fill results for manuscript
# -----

trimresults = data.frame(rbind(
  c("\\textbf{Visual factors}",NA,NA),
  c("Salience",salTrim$k0,coef(summary(salTrim))),
  c("Surface size",sizeTrim$k0,coef(summary(sizeTrim))),
  c("Left vs right position",LRTrim$k0,coef(summary(LRTrim))),
  c("Center position",centerTrim$k0,coef(summary(centerTrim))),
  c("Set size",setTrim$k0,coef(summary(setTrim))),
  c("\\hspace{2mm}\\textit{Alternative}",setmod_altTrim$k0,coef(summary(setmod_altTrim))),
  c("\\hspace{2mm}\\textit{Attribute}",setmod_attTrim$k0,coef(summary(setmod_attTrim))),
  c("\\textbf{Cognitive factors}",NA,NA),
  c("Task instruction",taskTrim$k0,coef(summary(taskTrim))),
  c("\\hspace{2mm}\\textit{Alternative}",taskmod_altTrim$k0,coef(summary(taskmod_altTrim))),
  c("\\hspace{2mm}\\textit{Attribute}",taskmod_attTrim$k0,coef(summary(taskmod_attTrim))),
  c("Preferential viewing",prefTrim$k0,coef(summary(prefTrim))),
  c("\\hspace{2mm}\\textit{Alternative}",prefmod_altTrim$k0,coef(summary(prefmod_altTrim))),
  c("\\hspace{2mm}\\textit{Attribute}",prefmod_attTrim$k0,coef(summary(prefmod_attTrim))),
  c("Choice bias",choiceTrim$k0,coef(summary(choiceTrim)))
))

# format trim and fill results table, e.g. rounding, variable naming etc.
trimresults[2:8] <- sapply(trimresults[2:8],as.numeric)
trimresults = trimresults %>% mutate_at(vars(2:8), round, 3)
setnames(trimresults, c(1:8), c("Group","Studies filled","$\\rho$","SE","$Z$","$p$","$\\textrm{CI}_{95}$ LL","$\\textrm{CI}_{95}$ UL"))
trimresults$Group = as.character(trimresults$Group)
write_csv(trimresults, file.path(tablesDir, "trim_fill_results.csv"))

# latex version
tab_caption <- "Trim and fill analysis for each visual and cognitive factor, including sub-factors used in the moderator analyses."
tab_label <- "tab:trim_fill_results"
tab_note <- paste0("\\hline \n \\multicolumn{8}{p{0.9\\textwidth}}",
           "{\\scriptsize{\\textit{Note.} $\\rho$ = unattenuated effect size estimate, SE = standard error of estimate; $Z$ = Z statistic; $p$ = significance level; $\\textrm{CI}_{95}$ LL = lower limit of the 95\\% confidence interval; $\\textrm{CI}_{95}$ UL = upper limit of the 95\\% confidence interval. Italicized groups are moderator subgroups.}} \n")
print(
	xtable(
		trimresults, 
		caption = tab_caption, 
		label = tab_label,
		align = "llccccccc",
		digits = c(0,0,0,3,3,3,3,3,3)
	), 
	include.rownames = FALSE,
	caption.placement = "top", 
	hline.after = c(-1, 0),
	add.to.row = list(
		pos = list(nrow(trimresults)),
        command = tab_note
    ),
    sanitize.text.function = function(x){x},
    file = file.path(tablesDir, "trim_fill_results.tex")
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
# Forrest plots 
# -----

# generating plots for each subgroup
salplot = genForest(salres,saldata,"Salience","Salience", "yi.c.FL")
sizeplot = genForest(sizeres,sizedata,"Size","Surface size", "yi.c.FL")
LRplot = genForest(LRres,LRdata,"LR.position","Left vs right position", "yi.c.FL")
centerplot = genForest(centerres,centerdata,"Center.position","Center position", "yi.c.FL")
setplot = genForest(setres,setdata,"Setsize", "Set size", "yi.c.FL")
taskplot = genForest(taskres,taskdata,"Task", "Inferential viewing", "yi.c.FC")
prefplot = genForest(prefres,prefdata,"Pref.view", "Preferential viewing", "yi.c.FC")
choiceplot = genForest(choiceres,choicedata,"Choice.bias","Choice bias", "yi.c.FC")

# att vs alt
setplot_alt = genForest(setmod_alt,setdata[Alt.att == "alternative"],"Setsize", "Set size - alternative", "yi.c.FL")
setplot_att = genForest(setmod_att,setdata[Alt.att == "attribute"],"Setsize", "Set size - attribute", "yi.c.FL")
taskplot_alt = genForest(taskmod_alt,taskdata[Alt.att == "alternative"],"Task", "Inferential viewing - alternative", "yi.c.FC")
taskplot_att = genForest(taskmod_att,taskdata[Alt.att == "attribute"],"Task", "Inferential viewing - attribute", "yi.c.FC")
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
salfunnel = genFunnel(salres,saldata,"yi.c.FL","Salience")
sizefunnel = genFunnel(sizeres,sizedata,"yi.c.FL", "Surface size")
LRfunnel = genFunnel(LRres,LRdata,"yi.c.FL", "Left vs right position")
centerfunnel = genFunnel(centerres,centerdata,"yi.c.FL", "Center position")
setfunnel = genFunnel(setres,setdata,"yi.c.FL", "Set size")
setfunnel_alt = genFunnel(setmod_alt,setdata[setdata$Alt.att == "alternative"],"yi.c.FL", "Set size - altenrative")
setfunnel_att = genFunnel(setmod_att,setdata[setdata$Alt.att == "attribute"],"yi.c.FL", "Set size - attribute")
taskfunnel = genFunnel(taskres,taskdata,"yi.c.FC", "Inferential viewing")
taskfunnel_alt = genFunnel(taskmod_alt,taskdata[taskdata$Alt.att == "alternative"],"yi.c.FC", "Inferential viewing - alternative")
taskfunnel_att = genFunnel(taskmod_att,taskdata[taskdata$Alt.att == "attribute"],"yi.c.FC", "Inferential viewing - attribute")
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

