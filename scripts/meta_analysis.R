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
saldata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Salience"], vtype="AV"))
saldata$vi.c = saldata$vi/saldata$a_acc^2 # compute corrected variances based on artefact multiplier
salres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=saldata, method="HS")
salresrobu = robust(salres, cluster = saldata$Authors)
salTrim = trimfill(salres)

# psychometric meta-analysis and trim and fill analysis of 
# surface size
sizedata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Size"], vtype="AV"))
sizedata$vi.c = sizedata$vi/sizedata$a_acc^2
sizeres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=sizedata, method="HS")
sizeresrobu = robust(sizeres, cluster = sizedata$Authors)
sizeTrim = trimfill(sizeres)

# psychometric meta-analysis and trim and fill analysis of 
# left v right position
LRdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "LR.position"], vtype="AV"))
LRdata$vi.c = LRdata$vi/LRdata$a_acc^2
LRres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=LRdata, method="HS")
LRresrobu = robust(LRres, cluster = LRdata$Authors)
LRTrim = trimfill(LRres)

# psychometric meta-analysis and trim and fill analysis of 
# centrality position
centerdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Center.position"], vtype="AV"))
centerdata$vi.c = centerdata$vi/centerdata$a_acc^2
centerres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=centerdata, method="HS")
centerresrobu = robust(centerres, cluster = centerdata$Authors)
centerTrim = trimfill(centerres)

# psychometric meta-analysis and trim and fill analysis of 
# set size
setdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize"], vtype="AV"))
setdata$vi.c = setdata$vi/setdata$a_acc^2

setdata_att = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize" & data$Alt.att == "attribute"], vtype="AV"))
setdata_att$vi.c = setdata_att$vi/setdata_att$a_acc^2

setdata_alt = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Setsize" & data$Alt.att == "alternative"], vtype="AV"))
setdata_alt$vi.c = setdata_alt$vi/setdata_alt$a_acc^2

setres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata, method="HS")
setresrobu = robust(setres, cluster = setdata$Authors)
setTrim = trimfill(setres)

# psychometric meta-analysis and trim and fill analysis of 
# set size moderator analysis - effect of attribute vs alternatives
setmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=setdata, method="HS") 
setmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata_att, method="HS") 
setmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=setdata_alt, method="HS") 
setmod_attrobu = robust(setmod_att, cluster = setdata_att$Authors)
setmod_altrobu = robust(setmod_alt, cluster = setdata_alt$Authors)
setmod_attTrim = trimfill(setmod_att)
setmod_altTrim = trimfill(setmod_alt)

# -----
# Meta analyses - cognitive factors
# -----

# psychometric meta-analysis and trim and fill analysis of 
# task instruction 
taskdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task"], vtype="AV"))
taskdata$vi.c = taskdata$vi/taskdata$a_acc^2

taskdata_att = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task" & data$Alt.att == "attribute"], vtype="AV"))
taskdata_att$vi.c = taskdata_att$vi/taskdata_att$a_acc^2

taskdata_alt = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Task" & data$Alt.att == "alternative"], vtype="AV"))
taskdata_alt$vi.c = taskdata_alt$vi/taskdata_alt$a_acc^2

taskres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata, method="HS")
taskresrobu = robust(taskres, cluster = taskdata$Authors)
taskTrim = trimfill(taskres)

# psychometric meta-analysis and trim and fill analysis of 
# task instruction moderator analysis - effect of attribute vs alternative 
taskmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=taskdata, method="HS")
taskmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata_att, method="HS")
taskmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=taskdata_alt, method="HS")
taskmod_attrobu = robust(taskmod_att, cluster = taskdata_att$Authors)
taskmod_altrobu = robust(taskmod_alt, cluster = taskdata_alt$Authors)
taskmod_altTrim = trimfill(taskmod_alt)
taskmod_attTrim = trimfill(taskmod_att)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing
prefdata = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view"], vtype="AV"))
prefdata$vi.c = prefdata$vi/prefdata$a_acc^2

prefdata_att = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view" & data$Alt.att == "attribute"], vtype="AV"))
prefdata_att$vi.c = prefdata_att$vi/prefdata_att$a_acc^2

prefdata_alt = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Pref.view" & data$Alt.att == "alternative"], vtype="AV"))
prefdata_alt$vi.c = prefdata_alt$vi/prefdata_alt$a_acc^2

prefres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata, method="HS")
prefresrobu = robust(prefres, cluster = prefdata$Authors)
prefTrim = trimfill(prefres)

# psychometric meta-analysis and trim and fill analysis of 
# preferential viewing moderator analysis - effect of alternative vs attribute 
prefmod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ Alt.att, data=prefdata, method="HS")
prefmod_att = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata_att, method="HS")
prefmod_alt = rma(yi.c.FC, vi.c, weights=1/vi.c, data=prefdata_alt, method="HS")
prefmod_attrobu = robust(prefmod_att, cluster = prefdata_att$Authors)
prefmod_altrobu = robust(prefmod_alt, cluster = prefdata_alt$Authors)
prefmod_attTrim = trimfill(prefmod_att)
prefmod_altTrim = trimfill(prefmod_alt)

# psychometric meta-analysis and trim and fill analysis of 
# choice bias moderator analysis - effect of inferential vs preferential choice
choicedata_mod = data.table(escalc(measure="COR", ri=fix.count.m, ni=N, data=data[data$IV == "Choice.bias"], vtype="AV"))
choicedata_mod$vi.c = choicedata_mod$vi/choicedata_mod$a_acc^2
choicedata_mod$inf_prefmod = ifelse(choicedata_mod$Research.Domain == "Risky Gamble", "preferential", # create inferential vs preferential task moderator variable
                         ifelse(choicedata_mod$Research.Domain == "pref. Consumer choice", "preferential",
                         ifelse(choicedata_mod$Research.Domain == "pref. Non-consumer choice", "preferential",
                         ifelse(choicedata_mod$Research.Domain == "inf. Consumer choice", "inferential",
                         ifelse(choicedata_mod$Research.Domain == "inf. Non-consumer choice", "inferential", NA)))))
choicemod = rma(yi.c.FC, vi.c, weights=1/vi.c, mods = ~ inf_prefmod, data=choicedata_mod, method="HS")

# psychometric meta-analysis and trim and fill analysis of 
# choice bias main analysis (averaging studies with more effect sizes)
choicedata = choicedata_mod[, list(yi.c.FC = mean(yi.c.FC), vi.c = mean(vi.c), IV = unique(IV), N = unique(N), Authors = unique(Authors)), by = Study]
choiceres = rma(yi.c.FC, vi.c, weights=1/vi.c, data=choicedata, method="HS")
choiceresrobu = robust(choiceres, cluster = choicedata$Authors)
choiceTrim = trimfill(choiceres)

# -----
# Publication bias analysis
# -----

# the publicatio analysis performs three tests for the degree of ES inflation due to publication bias
# the first inflation factor is based on the difference in ES between studies with public grants vs no pub grants
# reasoning behind the test: https://www.sciencedirect.com/science/article/abs/pii/S0895435617301348
# second inflation factor is derived from the PEESE corrected estimate relative the fixed effects uncorrected estimate
# third inflation factor is from the uncorrect ES relative to the trim fill corrected ES

# we perform analyses on fisher z transformed ES
data$fcz = FisherZ(data$fix.count.m) # effect in z
data$sdz = 1/sqrt(data$N - 3) # sd in z
data$varz = data$sdz^2 # variance in z

# analyze if public grants are associated with smaller ES
grant = as.data.table(read_excel(
  file.path(dataDir, "EMMA_grant.xlsx")
))
grant = grant[, list(grant = unique(grant), public = unique(public)), Study]
data = merge(data, grant, by = "Study")
pb = rma(yi=fcz, vi=varz, mods = ~ public, data = data)
publicFactor = round(pb$b[1] / (pb$b[1] + pb$b[2]), digits = 3) # inflation factor due to not having public grant
cat(paste0("$", publicFactor, "$"), file = file.path(tablesDir, "publicFactor.tex"))
cat(paste0("$Q_M(1)=", round(pb$QM, 3),"$, $p=", round(pb$QMp, 3), "$"), file = file.path(tablesDir, "publicSig.tex"))

# PET-PEESE test
FE = lm(fcz ~ 1, weights = 1/varz, data = data) # fixed effect estimate of ES
PET = lm(fcz ~ sdz + a_acc, weights = 1/varz, data = data) # PET test is sig therefore perfrom PEESE
PEESE = lm(fcz ~ varz + a_acc, weights = 1/varz, data = data) # PEESE estimate
peeseFactor = round(summary(FE)$coef[1] / summary(PEESE)$coef[1,1], digits = 3) # inflation factor according to PEESE
cat(paste0("$", peeseFactor, "$"), file = file.path(tablesDir, "peeseFactor.tex"))

# check inflation factor based on trim fill results
trims = c(salres$b / salTrim$b, # extract inflation factors for each subgroup separately
  sizeres$b / sizeTrim$b,
  LRres$b / LRTrim$b,
  centerres$b / centerTrim$b,
  setres$b / setTrim$b,
  prefres$b / prefTrim$b,
  taskres$b / taskTrim$b,
  choiceres$b / choiceTrim$b)
trimFactor = round(mean(trims), digits = 3) # average inflation factor
cat(paste0("$", trimFactor, "$"), file = file.path(tablesDir, "trimFactor.tex"))

# -----
# Table with publication bias results for manuscript
# -----

FE = round(data.frame(summary(FE)$coef), digits = 3)
PET = round(data.frame(summary(PET)$coef), digits = 3)
PEESE = round(data.frame(summary(PEESE)$coef), digits = 3)
FE = cbind(Parameter="Intercept", FE)
PET = cbind(Parameter=c("Intercept", "$SD$", "$A$"), PET)
PEESE = cbind(Parameter=c("Intercept", "$Var$", "$A$"), PEESE)
setnames(FE, c(3:5), c("SE","$t$","$p$"))
setnames(PET, c(3:5), c("SE","$t$","$p$"))
setnames(PEESE, c(3:5), c("SE","$t$","$p$"))

# latex version FE
tab_caption <- "Fixed effects analysis of complete data"
tab_label <- "tab:FE"
print(
  xtable(
    FE, 
    caption = tab_caption, 
    label = tab_label,
    align = "lllccc"
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(FE)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "FE.tex")
)

# latex version PET
tab_caption <- "Precision-effect test (PET) of complete data"
tab_label <- "tab:PET"
print(
  xtable(
    PET, 
    caption = tab_caption, 
    label = tab_label,
    align = "lllccc"
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(PET)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "PET.tex")
)

# latex version PEESE
tab_caption <- "Precision-effect estimate test (PEESE) of complete data"
tab_label <- "tab:PEESE"
print(
  xtable(
    PEESE, 
    caption = tab_caption, 
    label = tab_label,
    align = "lllccc"
  ), 
  size = "\\small",
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0, NROW(PEESE)),
  sanitize.text.function = function(x){x},
  file = file.path(tablesDir, "PEESE.tex")
)

# -----
# Main and Moderator results text for manuscript
# -----
# rho range
result <- paste0(
  "$\\rho=", round(salres$b, 3),"$ to $\\rho=", round(choiceres$b, 3), "$"
)
cat(result, file = file.path(tablesDir, "rhorange.tex"))

# I2 range
result <- paste0(
  "$I^2=", round(salres$I2, 3),"$ to $I^2=", round(prefres$I2, 3), "$"
)
cat(result, file = file.path(tablesDir, "I2range.tex"))

# # salience summary
# result <- paste0(
# 	"($\\rho=", round(salres$b, 2), "$; 95\\% confidence interval (CI) = $[", 
# 	round(salres$ci.lb, 2), ",", round(salres$ci.ub, 2), "]$; $p=", round(salres$pval, 3), "$)"
# )
# cat(result, file = file.path(tablesDir, "saliencesummary.tex"))
# 
# # salience trim summary
# result <- paste0(
#   "$\\rho=", round(salTrim$b, 2), "$; 95\\% CI = $[", 
#   round(salTrim$ci.lb, 2), ",", round(salTrim$ci.ub, 2), "];$"
# )
# cat(result, file = file.path(tablesDir, "saliencetrimsummary.tex"))
# 
# # center position summary
# result <- paste0(
#   "($\\rho=", 
#   round(centerres$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(centerres$ci.lb, 2), 
#   ",", 
#   round(centerres$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(centerres$pval, 3) == 0, "< 0.001", paste0("p=", round(centerres$pval, 3))), 
#   "$)"
# )
# cat(result, file = file.path(tablesDir, "centersummary.tex"))
# 
# # center trim summary
# result <- paste0(
#   "$\\rho=", 
#   round(centerTrim$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(centerTrim$ci.lb, 2), 
#   ",", 
#   round(centerTrim$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(centerTrim$pval, 3) == 0, "< 0.001", paste0("p=", round(centerTrim$pval, 3))), 
#   "$;"
# )
# cat(result, file = file.path(tablesDir, "centertrimsummary.tex"))
# 
# # task instruction summary
# result <- paste0(
#   "($\\rho=", 
#   round(taskres$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(taskres$ci.lb, 2), 
#   ",", 
#   round(taskres$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(taskres$pval, 3) == 0, "< 0.001", paste0("p=", round(taskres$pval, 3))), 
#   "$)"
# )
# cat(result, file = file.path(tablesDir, "tasksummary.tex"))
# 
# # preferential viewing summary
# result <- paste0(
#   "$\\rho=", 
#   round(prefres$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(prefres$ci.lb, 2), 
#   ",", 
#   round(prefres$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(prefres$pval, 3) == 0, "< 0.001", paste0("p=", round(prefres$pval, 3))), 
#   "$;"
# )
# cat(result, file = file.path(tablesDir, "prefsummary.tex"))
# 
# # task trim summary
# result <- paste0(
#   "$\\rho=", 
#   round(taskTrim$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(taskTrim$ci.lb, 2), 
#   ",", 
#   round(taskTrim$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(taskTrim$pval, 3) == 0, "< 0.001", paste0("p=", round(taskTrim$pval, 3))), 
#   "$;"
# )
# cat(result, file = file.path(tablesDir, "tasktrimsummary.tex"))
# 
# # pref trim summary
# result <- paste0(
#   "$\\rho=", 
#   round(prefTrim$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(prefTrim$ci.lb, 2), 
#   ",", 
#   round(prefTrim$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(prefTrim$pval, 3) == 0, "< 0.001", paste0("p=", round(prefTrim$pval, 3))), 
#   "$;"
# )
# cat(result, file = file.path(tablesDir, "preftrimsummary.tex"))
# 
# # choice summary
# result <- paste0(
#   "$\\rho=", 
#   round(choiceres$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(choiceres$ci.lb, 2), 
#   ",", 
#   round(choiceres$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(choiceres$pval, 3) == 0, "< 0.001", paste0("p=", round(choiceres$pval, 3))), 
#   "$;"
# )
# cat(result, file = file.path(tablesDir, "choicesummary.tex"))
# 
# # choice trim summary
# result <- paste0(
#   "$\\rho=", 
#   round(choiceTrim$b, 2), 
#   "$; 95\\% CI = $[", 
#   round(choiceTrim$ci.lb, 2), 
#   ",", 
#   round(choiceTrim$ci.ub, 2), 
#   "]$; $p", 
#   ifelse(round(choiceTrim$pval, 3) == 0, "< 0.001", paste0("p=", round(choiceTrim$pval, 3))), 
#   "$;"
# )
# cat(result, file = file.path(tablesDir, "choicetrimsummary.tex"))

# task moderator
result <- paste0(
	"$Q_M(1)=", round(taskmod$QM, 3),"$, $p=", round(taskmod$QMp, 3), "$"
)
cat(result, file = file.path(tablesDir, "moderator_task.tex"))

# pref moderator
result <- paste0(
	"$Q_M(1)=", round(prefmod$QM, 3),"$, $p=", round(prefmod$QMp, 3), "$"
)
cat(result, file = file.path(tablesDir, "moderator_pref.tex"))

# choice moderator
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
        paste0(round(coef(summary(mainres))$tval, 2)),
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
  extractMain("Salience",salresrobu,saldata),
  extractTrim(salTrim),
  extractMain("Surface size",sizeresrobu,sizedata),
  extractTrim(sizeTrim),
  extractMain("Left vs right position",LRresrobu,LRdata),
  extractTrim(LRTrim),
  extractMain("Center position",centerresrobu,centerdata),
  extractTrim(centerTrim),
  extractMain("Set size",setresrobu,setdata),
  extractTrim(setTrim),
  c("\\textbf{Cognitive factors}", rep(NA, 9)),
  extractMain("Task instructions",taskresrobu,taskdata),
  extractTrim(taskTrim),
  extractMain("Preferential viewing",prefresrobu,prefdata),
  extractTrim(prefTrim),
  extractMain("Choice bias",choiceresrobu,choicedata),
  extractTrim(choiceTrim)
), stringsAsFactors = FALSE)

# format main results table, e.g. rounding, variable naming etc.
setnames(mainresults, c(1:10), c("Group","$k$","$N$","$\\rho$","SE","$t$","$p$","$\\textrm{CI}_{95}$ LL","$\\textrm{CI}_{95}$ UL","$I^2$"))
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
        setmod_altrobu,setdata[setdata$Alt.att == "alternative",]),
    extractTrim(setmod_altTrim),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        setmod_attrobu,setdata[setdata$Alt.att == "attribute",]),
    extractTrim(setmod_attTrim),
    
    c("\\textbf{Task instruction}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        taskmod_altrobu,taskdata[taskdata$Alt.att == "alternative",]),
    extractTrim(taskmod_altTrim),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        taskmod_attrobu,taskdata[taskdata$Alt.att == "attribute",]),
    extractTrim(taskmod_attTrim),

    c("\\textbf{Preferential viewing}", rep(NA, 9)),
    extractMain("\\hspace{2mm}\\textit{Alternative}",
        prefmod_altrobu,prefdata[prefdata$Alt.att == "alternative",]),
    extractTrim(prefmod_altTrim),
    extractMain("\\hspace{2mm}\\textit{Attribute}",
        prefmod_attrobu,prefdata[prefdata$Alt.att == "attribute",]),
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

overviewtabel = data[, c("Authors", "IV", "N","a_acc","fix.count.m","Eye.tracker")]
setnames(overviewtabel, c("N","a_acc","fix.count.m","Eye.tracker"), c("$N$","$a_a$","$r$","Eye tracker"))
overviewtabel$`Eye tracker` = ifelse(overviewtabel$`Eye tracker` == "Nihon-Kohden EEG-1100", "Nihon-Kohden", overviewtabel$`Eye tracker`)
overviewtabel$IV = ifelse(overviewtabel$IV == "LR.position", "LvR",
                          ifelse(overviewtabel$IV == "Center.position", "Center",
                                 ifelse(overviewtabel$IV == "Pref.view", "Pref",
                                        ifelse(overviewtabel$IV == "Choice.bias", "Choice",
                                               ifelse(overviewtabel$IV == "Salience", "Sal",overviewtabel$IV)))))
overviewtabel = overviewtabel[order(Authors),]

tab_caption <- "Overview of individual effect sizes: IV = independent variable; $N$ = number of participants; $a_a$ = artifact multiplier; $r$ = attenuated effect size correlation expressed in the fixation count metric."
tab_label <- "tab:overviewtable"
add.to.row <- list(pos = list(0), command = NULL)
command <- paste0("\\hline\n\\endhead\n","\\hline\n","\\multicolumn{", dim(overviewtabel)[2], "}{l}","{\\footnotesize Continued on next page}\n","\\endfoot\n","\\endlastfoot\n")
add.to.row$command <- command

otab = xtable(
  overviewtabel, 
  caption = tab_caption, 
  label = tab_label, 
  digits = c(0,0,0,0,3,3,0)
)
align(otab) = "cp{8cm}lcccl"
print(otab
, 
  include.rownames = FALSE,
  caption.placement = "top", 
  hline.after = c(-1, 0),
  add.to.row = add.to.row,
  tabular.environment = "longtable",
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

