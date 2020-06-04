library(lme4)
library(ggplot2)
library(data.table)
library(lmerTest)
library(metafor)
library(tidyr)
library(dplyr)
library(readxl)
library(DescTools)
library(readxl)
library(cowplot)
library(jpeg)
library(gridExtra)
library(irr)

# set paths for retrieving data and storing plots
data_path = "/Users/au161118/Dropbox/ASB/Admin stuff/Posters & Papers/PAPERS/EMMA/data"
plot_path = "/Users/au161118/Dropbox/ASB/Admin stuff/Posters & Papers/PAPERS/EMMA/plots"
setwd(data_path)

#------------------------------------------------------------------
# compute artefact multiplier based on accuracy ####
#------------------------------------------------------------------

# read in data from excel coding sheet
data = as.data.table(read_excel("EMMA_ES_data.xlsx"))

# impute eye tracker "Unknown" when eye tracker name not identified
data$Eye.tracker[is.na(data$Eye.tracker)] = "Unknown"

# read in eye tracking accuracy and precision to merge with data file for precision values of ET
ET_specs = as.data.table(read.csv("EMMA_ET_specs_data.csv", header=T ,sep=";",dec=","))
ET_specs$Eye.tracker = as.character(ET_specs$Eye.tracker)
data = merge(data, ET_specs, by = c("Eye.tracker"), all.x = T)

# make long format for analysis of accuracy and precision on effect size attenuation
data_long = as.data.table(gather(data, DV, R, fix.count:dwell.count, factor_key = T))
data_long = data_long[data_long$R != "NA"]

# test whether eye tracker accuracy or precision attenuate effect sizes 
data_long$R_abs = abs(data_long$R) # the models are computed on the absolute effect sizes
data_long = data_long[, list(R_abs = mean(R_abs), Precision = unique(Precision), # averaging studies with multiple effect sizes
                             Accuracy = unique(Accuracy), Eye.tracker = unique(Eye.tracker), IV = unique(IV)), by = c("Study", "IV")]
table = table(data_long$IV, data_long$Eye.tracker)
chisq.test(table)
model1 <- lmer(R_abs~ + (1|Eye.tracker), data_long)
model2 <- lmer(R_abs~ + (1|Study), data_long)
model3 <- lmer(R_abs~ + (1|Eye.tracker) + (1|Study), data_long)
model21 <- lmer(R_abs~ Accuracy + (1|Study), data_long)
model212 <- lmer(R_abs~ Accuracy + IV + (1|Study), data_long)
model22 <- lmer(R_abs~ Precision + (1|Study), data_long)
model222 <- lmer(R_abs~ Precision + IV + (1|Study), data_long)
model23 <- lmer(R_abs~ Precision + Accuracy + (1|Study), data_long)
anova(model1,model2,model3,model21,model22,model23,model212,model222) # winning model includes accuracy

# make accuracy results table for manuscript
acc_table = data.table(summary(model21)$coef) 
obs = data.table("Number of observations =",NROW(data_long),NA,NA,NA)
Ll = data.table("Log likelihood =",as.numeric(summary(model21)$logLik),NA,NA,NA)
AIC = data.table("AIC =",AIC(model21),NA,NA,NA)
BIC = data.table("BIC =",BIC(model21),NA,NA,NA)
Raneff = data.table("SD(study) =",data.frame(summary(model21)$varcor)[1,5],NA,NA,NA)
l = list(acc_table,obs,Ll,AIC,BIC,Raneff)
acc_table = rbindlist(l)
acc_table = acc_table %>% mutate_at(vars(2:5), round, 3)
setnames(acc_table, c(2:5), c("SE", "df","t","p"))
write.csv(acc_table, "accuracytable.csv", row.names = F)

# compute artefact multiplier based on accuracy for psychometric meta-analysis
b0 = coef(summary(model21))[1,1] 
b1 = coef(summary(model21))[2,1]
data$a_acc = (b0 + b1*data$Accuracy)/b0 # see eq. 1 in manuscript

# make eye tracker specifications table for manuscript appendix
ET_specs_final = data[, list(A = unique(a_acc), accuracy = unique(Accuracy), precision = unique(Precision)), by= c("Eye.tracker")]
ET_specs_final$A <- round(ET_specs_final$A, 4)
write.csv(ET_specs_final, "ET_specstable.csv", row.names = F)

#------------------------------------------------------------------
# DV correction factors ####
#------------------------------------------------------------------

# compute correction factor am for transforming effect sizes in different metrics
# correction factor ratio = M_to / M_from
# M_. = FisherZ(sum(Fix*N)/N) see eq. 2 in manuscript

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
DC_to_FC = DC_to_TFD*TFD_to_FC # no cases for computing DC to FC or FL hence the current approach
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

#------------------------------------------------------------------
# HS MA main analysis + trim and fill analysis ####
#------------------------------------------------------------------

# analyses of bottom up factors (salience, surface size, left vs right position, central position, and set size)
# are performed using fixation likelihood. Analyses of top down factors are performed on fixation count

# compute corrected (unattenuated) effect sizes (performed on z transformed values)
# for psychometric meta-analysis
data$yi.c.FL = FisherZInv(FisherZ(data$fix.like.m)/data$a_acc)
data$yi.c.FC = FisherZInv(FisherZ(data$fix.count.m)/data$a_acc)

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

# paste a summary of moderator test results for manuscript
paste(
  "set size, QM(1) = ", round(setmod$QM, 3),", p = ", round(setmod$QMp, 3),", ",
  "task instruction, QM(1) = ", round(taskmod$QM, 3),", p = ", round(taskmod$QMp, 3),", ",
  "preferential viewing, QM(1) = ", round(prefmod$QM, 3),", p = ", round(prefmod$QMp, 3),", ",
  "choice bias, QM(1) = ", round(choicemod$QM, 3),", p = ", round(choicemod$QMp, 3), 
sep = "")

# combine meta-analysis results for the main results table in manuscript
mainresults = data.frame(rbind(
  c("Salience",salres$k,sum(saldata$N),coef(summary(salres)),salres$I2),
  c("Surface size",sizeres$k,sum(sizedata$N),coef(summary(sizeres)),sizeres$I2),
  c("Left vs right position",LRres$k,sum(LRdata$N),coef(summary(LRres)),LRres$I2),
  c("Center position",centerres$k,sum(centerdata$N),coef(summary(centerres)),centerres$I2),
  c("Set size",setres$k,sum(setdata$N),coef(summary(setres)),setres$I2),
  c("Alternative",setmod_alt$k,sum(setdata$N[setdata$Alt.att == "alternative"]),coef(summary(setmod_alt)),setmod_alt$I2),
  c("Attribute",setmod_att$k,sum(setdata$N[setdata$Alt.att == "attribute"]),coef(summary(setmod_att)),setmod_att$I2),
  c("Task instruction",taskres$k,sum(taskdata$N),coef(summary(taskres)),taskres$I2),
  c("Alternative",taskmod_alt$k,sum(taskdata$N[taskdata$Alt.att == "alternative"]),coef(summary(taskmod_alt)),taskmod_alt$I2),
  c("Attribute",taskmod_att$k,sum(taskdata$N[taskdata$Alt.att == "attribute"]),coef(summary(taskmod_att)),taskmod_att$I2),
  c("Preferential viewing",prefres$k,sum(prefdata$N),coef(summary(prefres)),prefres$I2),
  c("Alternative",prefmod_alt$k,sum(prefdata$N[prefdata$Alt.att == "alternative"]),coef(summary(prefmod_alt)),prefmod_alt$I2),
  c("Attribute",prefmod_att$k,sum(prefdata$N[prefdata$Alt.att == "attribute"]),coef(summary(prefmod_att)),prefmod_att$I2),
  c("Choice bias",choiceres$k,sum(choicedata$N),coef(summary(choiceres)),choiceres$I2)
  ))

# format main results table, e.g. rounding, variable naming etc.
mainresults[2:10] <- sapply(mainresults[2:10],as.numeric)
mainresults = mainresults %>% mutate_at(vars(4:10), round, 3)
setnames(mainresults, c(1:10), c("Group","k","N","Estimate","SE","Z","p","CI UL","CI LL","I2"))
mainresults$Group = as.character(mainresults$Group)
write.csv(mainresults, "mainresultstable.csv", row.names = F)

# combine trim and fill results for table in manuscript
trimresults = data.frame(rbind(
  c("Salience",salTrim$k0,coef(summary(salTrim))),
  c("Surface size",sizeTrim$k0,coef(summary(sizeTrim))),
  c("Left vs right position",LRTrim$k0,coef(summary(LRTrim))),
  c("Center position",centerTrim$k0,coef(summary(centerTrim))),
  c("Set size",setTrim$k0,coef(summary(setTrim))),
  c("Alternative",setmod_altTrim$k0,coef(summary(setmod_altTrim))),
  c("Attribute",setmod_attTrim$k0,coef(summary(setmod_attTrim))),
  c("Task instruction",taskTrim$k0,coef(summary(taskTrim))),
  c("Alternative",taskmod_altTrim$k0,coef(summary(taskmod_altTrim))),
  c("Attribute",taskmod_attTrim$k0,coef(summary(taskmod_attTrim))),
  c("Preferential viewing",prefTrim$k0,coef(summary(prefTrim))),
  c("Alternative",prefmod_altTrim$k0,coef(summary(prefmod_altTrim))),
  c("Attribute",prefmod_attTrim$k0,coef(summary(prefmod_attTrim))),
  c("Choice bias",choiceTrim$k0,coef(summary(choiceTrim)))
  ))

# format trim and fill results table, e.g. rounding, variable naming etc.
trimresults[2:8] <- sapply(trimresults[2:8],as.numeric)
trimresults = trimresults %>% mutate_at(vars(2:8), round, 3)
setnames(trimresults, c(1:8), c("Group","Studies filled","Estimate","SE","Z","p","CI UL","CI LL"))
trimresults$Group = as.character(trimresults$Group)
write.csv(trimresults, "trimresultstable.csv", row.names = F)

#------------------------------------------------------------------
# test of task instruction vs preferential viewing vs choice bias ####
#------------------------------------------------------------------

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
paste("task instruction vs preferential viewing, ", "z = ", round(task_vs_pref_z, 3), ", p = ", round(task_vs_pref_p, 3), sep = "")
paste("task instruction vs choice bias, ", "z = ", round(task_vs_choice_z, 3), ", p = ", round(task_vs_choice_p, 3), sep = "")
paste("preferential viewing vs choice bias, ", "z = ", round(pref_vs_choice_z, 3), ", p = ", round(pref_vs_choice_p, 3), sep = "")

#------------------------------------------------------------------
# forrest plots ####
#------------------------------------------------------------------

# forest plot function 
forestfunc = function(rmaobj, data, IV_level, label, varname) {
  fig = 
    data %>%
    select(Study, DV = varname, Domain = IV, vi.c) %>%
    mutate(SE = sqrt(vi.c)) %>%
    select(-(vi.c))
  Study = c("                               Summary effect")
  Domain = c("")
  DV = rmaobj[1]
  SE = rmaobj$se
  temp = data.table(Study, DV, SE, Domain)
  fig = rbind(fig,temp)
  fig$DV = as.numeric(fig$DV)
  fig$UL = FisherZInv(FisherZ(fig$DV) + 1.96*fig$SE)
  fig$LL = FisherZInv(FisherZ(fig$DV) - 1.96*fig$SE)
  fig$Domain = factor(fig$Domain, levels = c(IV_level,""), labels = c(label,""))

  plot = ggplot(fig, aes(x=reorder(Study,-DV),y=DV))+
    geom_point()+
    geom_errorbar(aes(ymin=LL, ymax=UL), width=.2)+
    coord_flip()+
    facet_grid(Domain ~., scales = "free_y",space = "free")+
    theme_bw()+
    ylim(-1, 1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("")+
    ylab("")+ 
    scale_color_manual(values="black")+
    geom_hline(yintercept = 0) +
    theme(legend.position="none")
  plot
}

# generating plots for each subgroup
salplot = forestfunc(salres,saldata,"Salience","Salience", "yi.c.FL")
sizeplot = forestfunc(sizeres,sizedata,"Size","Surface size", "yi.c.FL")
LRplot = forestfunc(LRres,LRdata,"LR.position","Left vs right position", "yi.c.FL")
centerplot = forestfunc(centerres,centerdata,"Center.position","Center position", "yi.c.FL")
setplot_alt = forestfunc(setmod_alt,setdata[Alt.att == "alternative"],"Setsize", "Set size - alternative", "yi.c.FL")
setplot_att = forestfunc(setmod_att,setdata[Alt.att == "attribute"],"Setsize", "Set size - attribute", "yi.c.FL")
taskplot_alt = forestfunc(taskmod_alt,taskdata[Alt.att == "alternative"],"Task", "Inferential viewing - alternative", "yi.c.FC")
taskplot_att = forestfunc(taskmod_att,taskdata[Alt.att == "attribute"],"Task", "Inferential viewing - attribute", "yi.c.FC")
prefplot_alt = forestfunc(prefmod_alt,prefdata[Alt.att == "alternative"],"Pref.view", "Preferential viewing - alternative", "yi.c.FC")
prefplot_att = forestfunc(prefmod_att,prefdata[Alt.att == "attribute"],"Pref.view", "Preferential viewing - attribute", "yi.c.FC") 
choiceplot = forestfunc(choiceres,choicedata,"Choice.bias","Choice bias", "yi.c.FC")

# arange all forest plots in panel plot for manuscript 
BUleftside = plot_grid(salplot,LRplot,setplot_alt,labels = c("A","C","E"), ncol = 1, rel_heights = c(9,5,7))
BUrightside = plot_grid(sizeplot,centerplot,setplot_att,labels = c("B","D","F"), ncol = 1, rel_heights = c(7,11,6))
BUplot = plot_grid(BUleftside, BUrightside, ncol = 2)
TDleftside = plot_grid(taskplot_alt,prefplot_alt,choiceplot,labels = c("G","I","K"), ncol = 1, rel_heights = c(8,6,11))
TDrightside = plot_grid(taskplot_att,prefplot_att,labels = c("H","J","F"), nrow = 3, rel_heights = c(11,12,9))
TDplot = plot_grid(TDleftside,TDrightside, ncol = 2)
completeplot = plot_grid(BUplot,TDplot, nrow = 2, rel_heights = c(2,3))
ggsave(paste(plot_path,"/completeplot.jpg", sep = ""), completeplot, width = 18, height = 23)

#------------------------------------------------------------------
# funnel plots ####
#------------------------------------------------------------------

# funnel plot function
funnelfunc = function(rmaobj, data, varname){ 
  fig = 
    data %>%
    select(DV = varname, vi.c) %>%
    mutate(SE = sqrt(vi.c)) %>%
    select(-(vi.c))  
  estimate = FisherZ(as.numeric(rmaobj[2]))  # rma estimate 
  se.seq=seq(0, sqrt(max(data$vi.c)), 0.001) # make SE vector from 0 to max SE
  ll95 = estimate-(1.96*se.seq) 
  ul95 = estimate+(1.96*se.seq)
  dfCI = data.frame(ll95, ul95, se.seq, estimate)
  ggplot(aes(x = SE, y = DV), data = fig) +
    geom_point(shape = 16) +
    xlab('Standard Error') + 
    ylab('z')+
    geom_line(aes(x = se.seq, y = ll95), data = dfCI) +
    geom_line(aes(x = se.seq, y = ul95), data = dfCI) +
    geom_segment(aes(x=sqrt(max(data$vi.c)),xend=0,y=estimate,yend=estimate))+
    scale_x_reverse()+
    coord_flip()+
    theme_bw()+                                      
    theme(panel.grid.major = element_blank(),        
          panel.grid.minor = element_blank(),        
          axis.title.x = element_text(face = "italic"))         
}

# generating funnel plots for each main and subgroup
salfunnel = funnelfunc(salres,saldata,"yi.c.FL")
sizefunnel = funnelfunc(sizeres,sizedata,"yi.c.FL")
LRfunnel = funnelfunc(LRres,LRdata,"yi.c.FL")
centerfunnel = funnelfunc(centerres,centerdata,"yi.c.FL")
setfunnel = funnelfunc(setres,setdata,"yi.c.FL")
setfunnel_alt = funnelfunc(setmod_alt,setdata[setdata$Alt.att == "alternative"],"yi.c.FL")
setfunnel_att = funnelfunc(setmod_att,setdata[setdata$Alt.att == "attribute"],"yi.c.FL")
taskfunnel = funnelfunc(taskres,taskdata,"yi.c.FC")
taskfunnel_alt = funnelfunc(taskmod_alt,taskdata[taskdata$Alt.att == "alternative"],"yi.c.FC")
taskfunnel_att = funnelfunc(taskmod_att,taskdata[taskdata$Alt.att == "attribute"],"yi.c.FC")
preffunnel = funnelfunc(prefres,prefdata,"yi.c.FC")
preffunnel_alt = funnelfunc(prefmod_alt,prefdata[prefdata$Alt.att == "alternative"],"yi.c.FC")
preffunnel_att = funnelfunc(prefmod_att,prefdata[prefdata$Alt.att == "attribute"],"yi.c.FC")
choicefunnel = funnelfunc(choiceres,choicedata,"yi.c.FC")

# arrange funnel plots in panel plot for manuscript 
funnelpanel = plot_grid(salfunnel, sizefunnel,LRfunnel,centerfunnel,
          setfunnel,setfunnel_alt,setfunnel_att,
          taskfunnel,taskfunnel_alt,taskfunnel_att,
          preffunnel,preffunnel_alt,preffunnel_att,
          choicefunnel, labels = LETTERS[1:14], ncol = 4)
ggsave(paste(plot_path,"/funnelplot.jpg", sep = ""), funnelpanel, width = 10, height = 9)

#------------------------------------------------------------------
# accuracy + metric conversion plots ####
#------------------------------------------------------------------

# plot of accuracy on effect size
accuracy_plot = ggplot(data=data_long, aes(Accuracy,R_abs))+
  geom_point(alpha=.5)+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(0,1,by=.2))+
  xlab("Accuracy")+
  ylab("Observed effect size")

# fix like to fix count plot
FLFCfig = data[, fix.count:N]
FLFCfig = FLFCfig[FLFCfig$fix.count != "NA" & FLFCfig$fix.like != "NA"]
FLFC_plot = ggplot(data=FLFCfig, aes(fix.count,fix.like))+
  geom_point(alpha=.5)+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(0,1,by=.2))+
  xlab("Fixation count")+
  ylab("Fixation likelihood")

# fix like to TFD plot
FLTFDfig = data[, fix.count:N]
FLTFDfig = FLTFDfig[FLTFDfig$TFD != "NA" & FLTFDfig$fix.like != "NA"]
FLTFD_plot = ggplot(data=FLTFDfig, aes(TFD,fix.like))+
  geom_point(alpha=.5)+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(0,1,by=.2))+
  xlab("Total fixation duration")+
  ylab("Fixation likelihood")

# fix count to TFD plot
FCTFDfig = data[, fix.count:N]
FCTFDfig = FCTFDfig[FCTFDfig$TFD != "NA" & FCTFDfig$fix.count != "NA"]
FCTFD_plot = ggplot(data=FCTFDfig, aes(TFD,fix.count))+
  geom_point(alpha=.5)+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(0,1,by=.2))+
  xlab("Total fixation duration")+
  ylab("Fixation count")

# dwell count to TFD plot
DCTFDfig = data[, fix.count:N]
DCTFDfig = DCTFDfig[DCTFDfig$TFD != "NA" & DCTFDfig$dwell.count != "NA"]
DCTFD_plot = ggplot(data=DCTFDfig, aes(TFD,dwell.count))+
  geom_point(alpha=.5)+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = seq(0,1,by=.2))+
  xlab("Total fixation duration")+
  ylab("Dwell count")

# arrange plots in panel for manuscript
acc_metric_plot = plot_grid(FLFC_plot, FLTFD_plot, FCTFD_plot, DCTFD_plot, accuracy_plot, labels = LETTERS[1:5], ncol = 2)
ggsave(paste(plot_path,"/acc_metric_plot.jpg", sep = ""), acc_metric_plot, width = 7, height = 8)

#------------------------------------------------------------------
# intercoder reliability ####
#------------------------------------------------------------------

# read in coder realiability data
coder1 = read_excel("EMMA_intercoder_reliability_data.xlsx", sheet= 1, col_names = T)
coder2 = read_excel("EMMA_intercoder_reliability_data.xlsx", sheet= 2, col_names = T)
ICCdata = cbind(coder1,coder2)

# kappa for categorical variables
DV_kappa = kappa2(ICCdata[,c(2,8)], "unweighted")
IV_kappa = kappa2(ICCdata[,c(3,9)], "unweighted")
ET_kappa = kappa2(ICCdata[,c(4,10)], "unweighted")
domain_kappa = kappa2(ICCdata[,c(7,13)], "unweighted")

# ICC for continuous variables
ES_ICC = icc(ICCdata[,c(5,11)], model="oneway", type="agreement")
N_ICC = icc(ICCdata[,c(6,12)], model="oneway", type="agreement")

# paste summary of inter coder reliability for manuscript
paste("effect size = ", round(ES_ICC$value, 3),",",
      " sample size = ", round(N_ICC$value, 3),",",
      " research domain = ", round(domain_kappa$value, 3),",",
      " eye tracker model = ", round(ET_kappa$value, 3),",",
      " dependent variable = ", round(DV_kappa$value, 3),",", 
      " independent variable = ", round(IV_kappa$value, 3),",",
      sep = "")
