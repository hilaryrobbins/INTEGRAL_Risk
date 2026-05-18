# ###############################################
# Script necessary to reproduce all tables and figures from the Zahed et al. JAMA 2026
# Paper title: Biomarker-based eligibility for lung cancer screening
#              Validation of the protein-based INTEGRAL-Risk model 
# Authors: Mainly: Hana Zahed; Secondary: Xiaoshuang Feng
# Senior authors: Hilary Robbins, Mattias Johansson
# Contact email: fengx@iarc.who.int; robbinsh@iarc.who.int; johanssonm@iarc.who.int
# #################################################


# #####################
### Set up ####
rm(list=ls())
path_to_save <- "your path of saving results"

# 1. Load packages
packages <- c("tidyverse","dplyr","arsenal","gmodels","survival","splines","glmnet","WeightedROC","survey", 
              "rstpm2", "stringr","grplasso", "ggpubr", "ggtext", "gridExtra", "grid", "ggpmisc", "ggforce")

lapply(packages, require, c = T)
source("https://github.com/hilaryrobbins/INTEGRAL_Risk/INTEGRAL_function_validation.R")

##2. Load models
#Please download models INTEGRAL RISK model for the Github
#INTEGRAL_RISK_AIC: The risk model identified as INTEGRAL-Risk
INTEGRAL_RISK_AIC<-read_rds("INTEGRAL_RISK_AIC.rds")

##3. Load datasets:
# data_train: Training dataset
# testing_all_imp: Testing dataset for validation, For paper, we used 15 imputed datasets pooled together 
# testing_1_imp: Random 1 of 15 imputed testing datasets

##4. Key variables:
#LC_event: lung cancer event
#t0: enrollment time
#tevent: follow-up time until  lung cancer diagnosis, death, or end of follow-up (3 years) 


# ###########################################################################
#########Table 1
##Characteristics of participants included in the training and testing sets 
##for development of the INTEGRAL-Risk model.
# ###########################################################################
render.continuous.onlyIQR <- function(x, name, ...) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  m   <- median(x, na.rm = TRUE)
  c(`Median [IQR]` = sprintf("%0.2f [%0.2f, %0.2f]", m, q1, q3))}

Table1_test<-table1::table1(~age+factor(sex)+factor(race)+factor(smoke_status)+intensity+years_smoked+packyears+factor(history_cancer)+factor(copd)+cohort|factor(LC_event),
                            data=testing_1_imp,render.continuous = render.continuous.onlyIQR)
summary(testing_1_imp[testing_1_imp$smoke_status==2&testing_1_imp$LC_event==1,]$quit_years)
summary(testing_1_imp[testing_1_imp$smoke_status==2&testing_1_imp$LC_event==0,]$quit_years)

Table1_train<-table1::table1(~age+factor(sex)+factor(race)+factor(smoke_status)+intensity+years_smoked+packyears+factor(history_cancer)+factor(copd)+cohort|factor(LC_event),
                             data=data_train,render.continuous = render.continuous.onlyIQR)
summary(data_train[data_train$smoke_status==2&data_train$LC_event==1,]$quit_years)
summary(data_train[data_train$smoke_status==2&data_train$LC_event==0,]$quit_years)

write.csv(Table1_test,paste0(path_to_save,"Table1_test.csv"),row.names = FALSE)
write.csv(Table1_train,paste0(path_to_save,"Table1_train.csv"),row.names = FALSE)

###Weighted table1
## ---- function to build weighted Table 1 for one dataset ----
make_table1 <- function(data_name, label) {
  des <- svydesign(ids = ~1, weights = ~U19_weight, data = data_name)
  cont_vars <- c("age", "intensity", "years_smoked","quit_years","packyears")
  cat_vars  <- c("sex", "race","smoke_status","history_cancer","copd","cohort")
  cont <- lapply(cont_vars, function(v) {
    q <- svyquantile(as.formula(paste0("~", v)),des,c(0.25, 0.5, 0.75),na.rm = TRUE)
    qv <- coef(q)
    tibble(variable = v,level = NA,
           !!label := sprintf("%.1f (%.1f–%.1f)",qv[2],qv[1],qv[3]))}) %>% bind_rows()
  
  cat <-lapply(cat_vars, function(v) {
    data_name %>% group_by(level = as.character(.data[[v]])) %>%
      summarise(wt_n = sum(U19_weight),.groups = "drop") %>%
      mutate(pct = 100 * wt_n / sum(wt_n),variable = v,
             !!label := sprintf("%.0f (%.1f%%)", wt_n, pct)) %>%
      select(variable, level, !!label)}) %>% bind_rows()
  
  bind_rows(cont, cat)}

## ---- build tables ----
t1 <- make_table1(data_train %>% filter(LC_event==0), "Training")
t2 <- make_table1(testing_1_imp%>% filter(LC_event==0), "Testing")

t1_quit <- make_table1(data_train %>% filter(LC_event==0&smoke_status==2), "Training")
t2_quit <- make_table1(testing_1_imp%>% filter(LC_event==0&smoke_status==2), "Testing")

write.csv(t1,paste0(path_to_save,"Table1_train_weight.csv"),row.names = FALSE)
write.csv(t2,paste0(path_to_save,"Table1_test_weight.csv"),row.names = FALSE)

write.csv(t1_quit,paste0(path_to_save,"Table1_train_weight_former.csv"),row.names = FALSE)
write.csv(t2_quit,paste0(path_to_save,"Table1_test_weight_former.csv"),row.names = FALSE)


# ###########################################################################
### Figure 2 ####
# Performance of the INTEGRAL-Risk model for 1-, 2- and 3-year lung cancer risk in the testing set. 
# ###########################################################################

## add PLCOm2012 score #### 
testing_all_imp <- testing_all_imp %>% mutate( plcom2012_logit = -4.532506 + 0.0778868*(age-62) + 0*(race_plco==1 | race_plco==6) +
                                                0.3944778*(race_plco==2) + (-0.7434744)*(race_plco==3) + (-0.466585)*(race_plco==4) + 
                                                1.027152*(race_plco==5) +
                                                (-0.0812744)*(education-4) + (-0.0274194)*(bmi-27) + 0.3553063*(copd==1) + 0.4589971*(history_cancer==1) +
                                                0.587185*(family_history_binary==1) + 0.2597431*(smoke_status==3) + 
                                                (-1.822606)*((1/(intensity/10))-0.4021541613) +
                                                0.0317321*(years_smoked-27) + (-0.0308572)*(quit_years-10),
                                              plcom2012 = as.numeric(exp(plcom2012_logit) / (1+exp(plcom2012_logit))))
## add USPSTF eligibility ####
testing_all_imp <- testing_all_imp%>% 
  mutate(USPSTF2021=ifelse(age>=50& age<=80&packyears>=20 &quit_years<=15,1,0),
         USPSTF2013=ifelse(age>=55& age<=80&packyears>=30 &quit_years<=15,1,0))

# add predicted model risks in the Testing set ####
testing_all_imp <- calculate_modelrisk_time(testing_all_imp, c("INTEGRAL_RISK_AIC"), c("proteinmodelrisk"))

bootstraping_res <- data.frame(NULL)
modelstest_vec <- c("proteinmodelrisk", "plcom2012") 
for(imp in 1:15){
  print(imp)
  dataimp <- testing_all_imp %>% filter(name==imp)###imp represents the number of imputed datasets, 1,2,3,,,,15
  for(t in 1:3){
    print(t)
    LC_event_t <- paste0("LC_event", t,"y")
    dataimp_t <- dataimp %>% filter(!is.na(!!sym(LC_event_t)))
    # Models ROC
    for(mod in modelstest_vec){ #ROC by model
      mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
      boot_res_iter <- bootstrapped_AUC_cal(dataimp_t, mod_t, LC_event_t) %>% 
        mutate(model=mod, time=t, iter=imp)
      ### bind all bootstrap
      bootstraping_res <- bind_rows(bootstraping_res, boot_res_iter) # within each iteration of imputation: Bootstrap inference when using multiple imputation
    }
  }
}
bootstraping_res %>% head

with(bootstraping_res, table(model,time)) #15 imputed datasets, 1000 bootstrapped iter, 3 time, 3 models

pooled_discriminatory_perf_main <- summary_res_pooled_ROC(bootstraping_res, c("plcom2012", "proteinmodelrisk"))

# Discrimination assessment in the testing set ####
USPSTF_perf2013 <- data.frame(time=NULL)
for (t in 1:3) {
  perf <- perf_uspstf(testing_1_imp, "USPSTF2013",paste0("LC_event",t,"y"))
  names(perf) <- paste0(names(perf),"_2013")
  perf <- perf %>% mutate(time=t)
  USPSTF_perf2013 <- bind_rows(USPSTF_perf2013,perf)
}

USPSTF_perf2021 <- data.frame(time=NULL)
for (t in 1:3) {
  perf <- perf_uspstf(testing_1_imp, "USPSTF2021",paste0("LC_event",t,"y"))
  names(perf) <- paste0(names(perf),"_2021")
  perf <- perf %>% mutate(time=t)
  USPSTF_perf2021 <- bind_rows(USPSTF_perf2021,perf)
}

# ROCs FPR and TPR for 1 iteration(data) of the imputation 
ROC_mainplot <- data.frame(NULL)
for(t in 1:3){
  LC_event_t <- paste0("LC_event", t,"y")
  dataimp_t <- testing_1_imp %>% filter(!is.na(!!sym(LC_event_t)))
  for (mod in c("plcom2012", "proteinmodelrisk")){
    mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
    ROC_data_interim <- with(dataimp_t, WeightedROC(guess=get(mod_t), label=get(LC_event_t), 
                                                    weight=U19_weight)) %>% mutate(time=t, model=mod)
    ROC_mainplot <- bind_rows(ROC_mainplot, ROC_data_interim)
  }
}
dim(ROC_mainplot)
ROC_mainplot <- ROC_mainplot %>% merge(dplyr::select(pooled_discriminatory_perf_main, time, model, AUC_CI_pooled,
                                                     pvalue_ref), 
                                       by=c("time","model"), all.x = T, all.y = F)
dim(ROC_mainplot) 

ROC_mainplot <- ROC_mainplot%>% 
  mutate(modelraw=case_when(model=="plcom2012"~"PLCOm2012", 
                            model=="proteinmodelrisk"~"INTEGRAL Risk model", 
                            T~model),
         model=paste0(modelraw,"\nAUC=", AUC_CI_pooled))


# ###
#### Sensitivity, True positive and NNS for INTEGRAL-RISK and PLCOm2012 at the same Specificity as USPSTF-2013/21
#    15 imputed datasets with confidence interval 
# ##
bootstrapped_sensitivity <- function(datax, uspstfy, LC_event_time,mod,t){
  results_all <- list()
  datax <- datax %>% filter(!is.na(!!sym(LC_event_time)))
  for(i in 1:1000){
    set.seed(123*i)
    ## bootstrapping data: resampling data similar to original pop/ full data will not always have same IR, that is 
    bootstrapped_data_m3 <- datax %>% group_by(cohort) %>%
      sample_n(size=n(), replace = T) %>% ungroup()
    
    ###USPSTF
    TP_uspstf <- bootstrapped_data_m3 %>% filter(!!sym(LC_event_time)==1 & !!sym(uspstfy)==1) %>% nrow()
    FN_uspstf <- bootstrapped_data_m3 %>% filter(!!sym(LC_event_time)==1 & !!sym(uspstfy)==0) %>% nrow() 
    FP_uspstf <- bootstrapped_data_m3 %>% filter(!!sym(LC_event_time)==0 & !!sym(uspstfy)==1) %>% pull(U19_weight) %>% sum() 
    TN_uspstf <- bootstrapped_data_m3 %>% filter(!!sym(LC_event_time)==0 & !!sym(uspstfy)==0) %>% pull(U19_weight) %>% sum()  
    
    TPR_uspstf  <- TP_uspstf/(TP_uspstf+FN_uspstf )
    FPR_uspstf  <- FP_uspstf/(FP_uspstf+TN_uspstf )
    TNR_uspstf <- TN_uspstf / (TN_uspstf + FP_uspstf)  # specificity (weighted)

    ####Model
    for (mod in c("plcom2012", "proteinmodelrisk")){
      mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
      ROC_boot <- with(bootstrapped_data_m3, WeightedROC(guess=get(mod_t), label=get(LC_event_time), 
                                                         weight=U19_weight)) %>%as.data.frame() %>%mutate(TNR = 1 - FPR)
      
      # Match model specificity to USPSTF specificity (closest FPR)
      matched_row <- ROC_boot %>%
        mutate(fpr_diff = abs(FPR - FPR_uspstf)) %>%
        slice_min(fpr_diff, n = 1, with_ties = FALSE)
      
      results_all[[length(results_all) + 1]] <- data.frame(
        boot_i           = i,
        model            = mod,
        time             = t,
        uspstf_criterion = uspstfy,
        sens_model       = matched_row$TPR,
        spec_model       = matched_row$TNR,
        TP_model         = matched_row$TPR*(TP_uspstf+FN_uspstf),
        TN_model         = matched_row$TNR*(FP_uspstf+TN_uspstf),
        FP_model         = matched_row$FP,
        FN_model         = matched_row$FN,
        thredshold_model       = matched_row$threshold,
        sens_uspstf      = TPR_uspstf,
        spec_uspstf      = TNR_uspstf,
        TP_uspstf        =TP_uspstf,
        FP_uspstf        =FP_uspstf,
        TN_uspstf        =TN_uspstf)
    }}
  bind_rows(results_all)
}


# ── Run across 15 imputed datasets & both USPSTF criteria ────────────────────
all_results_1y <- list()
for (imp_id in 1:15) {
  datax <- testing_all_imp %>% filter(name==imp_id)  # replace with your indexing
  for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
    res <- bootstrapped_sensitivity(
      datax         = datax,uspstfy       = uspstfy,
      LC_event_time = "LC_event1y",t= 1)
    res$imp_id <- imp_id
    all_results_1y[[length(all_results_1y) + 1]] <- res
  }}
final_results_1y <- bind_rows(all_results_1y)%>% mutate(time=1)

all_results_2y <- list()
for (imp_id in 1:15) {
  datax <- testing_all_imp%>% filter(name==imp_id)  # replace with your indexing
  for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
    res <- bootstrapped_sensitivity(
      datax         = datax,uspstfy       = uspstfy,
      LC_event_time = "LC_event2y",t= 2)
    res$imp_id <- imp_id
    all_results_2y[[length(all_results_2y) + 1]] <- res
  }}
final_results_2y <- bind_rows(all_results_2y)%>% mutate(time=2)

all_results_3y <- list()
for (imp_id in 1:15) {
  datax <- testing_all_imp%>% filter(name==imp_id)   # replace with your indexing
  for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
    res <- bootstrapped_sensitivity(
      datax  = datax,uspstfy       = uspstfy,
      LC_event_time = "LC_event3y",t= 3)
    res$imp_id <- imp_id
    all_results_3y[[length(all_results_3y) + 1]] <- res
  }}
final_results_3y <- bind_rows(all_results_3y) %>% mutate(time=3)

final_results<-bind_rows(final_results_1y,final_results_2y,final_results_3y)

final_results<-final_results %>% mutate( NNS_model=(TP_model+FP_model)/TP_model,
                                        NNS_uspstf=(TP_uspstf+FP_uspstf)/TP_uspstf)

write.csv(final_results,paste0(path_to_save,"Figure2_SE_SP_bootstrap_results.csv"),row.names = FALSE)

# ── Summarise: mean + 95% CI across 15,000 bootstraps ────────────────────────
summary_results <- final_results %>%
  group_by(time,uspstf_criterion,model) %>%
  summarise(
    mean_threshold=mean(thredshold_model,       na.rm = TRUE),
    median_threshold=median(thredshold_model,       na.rm = TRUE),
    
    mean_sens  = mean(sens_model,       na.rm = TRUE),
    lower_sens = quantile(sens_model,   0.025, na.rm = TRUE),
    upper_sens = quantile(sens_model,   0.975, na.rm = TRUE),
    mean_TP=mean(TP_model,na.rm = TRUE ),
    mean_spec  = mean(spec_model,       na.rm = TRUE),
    
    mean_NNS  = mean(NNS_model,       na.rm = TRUE),
    lower_NNS = quantile(NNS_model,   0.025, na.rm = TRUE),
    upper_NNS = quantile(NNS_model,   0.975, na.rm = TRUE),
    
    # USPSTF reference
    mean_sens_uspstf = mean(sens_uspstf, na.rm = TRUE),
    lower_sens_uspstf = quantile(sens_uspstf,   0.025, na.rm = TRUE),
    upper_sens_uspstf = quantile(sens_uspstf,   0.975, na.rm = TRUE),
    mean_spec_uspstf = mean(spec_uspstf, na.rm = TRUE),
    mean_TP_uspstf=mean(TP_uspstf ,na.rm = TRUE ),
    mean_TN_uspstf=mean(TN_uspstf ,na.rm = TRUE ),
    
    mean_NNS_uspstf  = mean(NNS_uspstf,       na.rm = TRUE),
    lower_NNS_uspstf = quantile(NNS_uspstf,   0.025, na.rm = TRUE),
    upper_NNS_uspstf = quantile(NNS_uspstf,   0.975, na.rm = TRUE),
    
    n_boots    = n(),
    .groups    = "drop"
  )

print(summary_results)
write.csv(summary_results,paste0(path_to_save,"Figure2_SE_SP_bootstrap_summary_results.csv"),row.names = FALSE)
summary_results %>% head()

####Plot ROC curves
ROC_time_plot_update <- function(data_plot, title_p,  linesize=0.75, legendsize, legendpos,
                                 axissize=12,axistitlesize=15){
  pts <- data_plot %>%
    dplyr::slice(1) %>%
    dplyr::select(FPR_2013, TPR_2013, FPR_2021, TPR_2021)
  plot <-data_plot %>%ggplot()+
    geom_path(aes(x=FPR, y=TPR, color=model), linewidth = linesize)+
    ggtitle(title_p)+
    labs(color="")+coord_equal()+theme_bw()+ggsci::scale_color_nejm()+ 
    geom_abline(slope=1, intercept=0, color="#e4e4e4")+
    labs(x = "False positive rate", y = "True positive rate")+
    theme(legend.text = element_text(margin = margin(t = 10), size=legendsize), 
          legend.key = element_rect(colour = NA, fill = NA),
          legend.position = legendpos,
          axis.text = element_text(size= axissize), 
          axis.title=element_text(size=axistitlesize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_continuous(labels = scales::percent, expand = c(0, 0),limits = c(0,NA),
                       breaks=c(0,0.25,0.50,0.75,1))+
    scale_y_continuous(labels = scales::percent, expand = c(0, 0),limits = c(0,NA),
                       breaks=c(0.25,0.50,0.75,1))
  plot <- plot+
    geom_point(data = pts, aes(x = FPR_2013, y = TPR_2013),shape = 23, size=5)+
    geom_text(data = pts, aes(x = FPR_2013+0.03, y = TPR_2013, label = "USPSTF-2013"),  fontface = "plain", hjust=0,size=5)+
    geom_point(data = pts, aes(x = FPR_2021, y = TPR_2021),shape = 1, size=5)+
    geom_text(data = pts, aes(x = FPR_2021+0.03, y = TPR_2021, label = "USPSTF-2021"), fontface = "plain", hjust=0,size=5)
  return(plot)
}

for(t in 1:3){
  title <- paste0("Discrimination at ", t,"y")
  data_f <- filter(ROC_mainplot, time==t)
  data_f <- data_f %>% merge(USPSTF_perf2021, by="time", all.x = TRUE, all.y = FALSE)
  data_f <- data_f %>% merge(USPSTF_perf2013, by="time", all.x = TRUE, all.y = FALSE)
  pvalue <- formatC(unique(data_f$pvalue), digits=3, format="f")
  pvalue_annot <- ifelse(pvalue<0.001, "pdiff<0.001", paste0("pdiff=",pvalue))
  plot <- ROC_time_plot_update(data_f, title, 1.5, 15, c(.7, 0.35),15,15)+
    annotate("text", y=0.15, x=0.60, label=pvalue_annot, 
             fontface="italic", size=6, hjust=1, vjust=1)
  assign(paste0("ROC_main_", t,"y"), plot)
  print(plot)
}

#########Sensitivity, True positive, NNS plot
summary_results_plot<-summary_results %>% 
  mutate(`Sensitity, % (95CI)`=paste0(sprintf("%.1f", mean_sens*100)," (",sprintf("%.1f", lower_sens*100),"-",sprintf("%.1f", upper_sens*100),")"),
         `NNS, % (95CI)`=paste0(sprintf("%.0f", mean_NNS)," (",sprintf("%.0f", lower_NNS),"-",sprintf("%.0f", upper_NNS),")"),
         `Sensitity USPSTF, % (95CI)`=paste0(sprintf("%.1f", mean_sens_uspstf*100)," (",sprintf("%.1f", lower_sens_uspstf*100),"-",sprintf("%.1f", upper_sens_uspstf*100),")"),
         `NNS USPSTF, % (95CI)`=paste0(sprintf("%.0f", mean_NNS_uspstf)," (",sprintf("%.0f", lower_NNS_uspstf),"-",sprintf("%.0f", upper_NNS_uspstf),")"),
         mean_TP_uspstf=round(mean_TP_uspstf,0),
         mean_TP=round(mean_TP,0)) %>% 
  dplyr::select(- c("mean_sens","lower_sens", "upper_sens","mean_spec","n_boots","mean_sens_uspstf","lower_sens_uspstf","upper_sens_uspstf","mean_TN_uspstf",
                    "mean_NNS_uspstf","lower_NNS_uspstf","upper_NNS_uspstf", "mean_NNS","lower_NNS","upper_NNS")) 

data_figure2def <- data.frame(NULL)
for(t in 1:3){
  data_figure2def_sub<-data.frame(
    Metric = c("USPSTF-2013/2021", "PLCOm2012 model","INTEGRAL Risk model"),
    time=t,
    us13_TP   = c(summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="plcom2012",] $mean_TP_uspstf,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="plcom2012",] $mean_TP,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="proteinmodelrisk",] $mean_TP),
    us13_SE   = c(summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="plcom2012",] $`Sensitity USPSTF, % (95CI)`,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="plcom2012",] $`Sensitity, % (95CI)`,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="proteinmodelrisk",] $`Sensitity, % (95CI)`),
    us13_NNS   = c(summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="plcom2012",] $`NNS USPSTF, % (95CI)`,
                   summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="plcom2012",] $`NNS, % (95CI)`,
                   summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2013"&summary_results_plot$model=="proteinmodelrisk",] $`NNS, % (95CI)`),
    
    us21_TP   = c(summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="plcom2012",] $mean_TP_uspstf,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="plcom2012",] $mean_TP,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="proteinmodelrisk",] $mean_TP),
    us21_SE   = c(summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="plcom2012",] $`Sensitity USPSTF, % (95CI)`,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="plcom2012",] $`Sensitity, % (95CI)`,
                  summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="proteinmodelrisk",] $`Sensitity, % (95CI)`),
    
    us21_NNS   = c(summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="plcom2012",] $`NNS USPSTF, % (95CI)`,
                   summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="plcom2012",] $`NNS, % (95CI)`,
                   summary_results_plot[summary_results_plot$time==t&summary_results_plot$uspstf_criterion=="USPSTF2021"&summary_results_plot$model=="proteinmodelrisk",] $`NNS, % (95CI)`))
  
  data_figure2def<-bind_rows(data_figure2def,data_figure2def_sub)
}


data_figure2def_long <- data_figure2def %>%
  mutate(across(c(us13_TP, us21_TP), as.character)) %>%  # make all columns character
  pivot_longer(
    cols = c(us13_TP, us13_SE,us13_NNS, us21_TP, us21_SE,us21_NNS),
    names_to = c("group", "metric"),
    names_sep = "_"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = value
  ) %>%
  mutate(metric = case_when(metric == "TP"~ "TP, N", 
                            metric =="SE"~"Sensitivity, %",
                            metric =="NNS"~"quasi-NNS"))
data_figure2def_long$time<-as.numeric(data_figure2def_long$time)

library(ggplotify)
make_table <- function(data, group_name) {
  data<-data %>% dplyr::select(-time)%>%
    group_by(Metric) %>%
    mutate(Label = if_else(
      row_number() == 1,
      paste0(Metric, "\n  ", metric),
      paste0("  ", metric))) %>%
    ungroup()%>%
    dplyr::select(
      Label,us13,us21)
  header_bottom <- c("Model", "USPSTF-2013\n Specificity at 87.2%", "USPSTF-2021\n Specificity at 75.1%")
  mat <- rbind(header_bottom, as.matrix(data))
  tg <- tableGrob(mat, rows = NULL, cols = NULL,
                  theme = ttheme_minimal(
                    core    = list(fg_params = list(hjust = 0.5, x = 0.5, fontsize = 14),
                                   bg_params = list(fill = NA, col = NA))))
  hline <- function(tg, r, lwd = 0.8) {
    gtable_add_grob(tg,segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"),
      y0 = unit(0, "npc"), y1 = unit(0, "npc"),
      gp = gpar(lwd = lwd)),
      t = r, b = r, l = 1, r = ncol(tg), z = Inf)}
  n_rows <- nrow(mat)
  # top line
  tg <- hline(tg, r = 1, lwd = 1.5)
  # between-group lines (every 2 data rows, skip within-group)
  group_rows <- seq(4, n_rows - 2, by = 3)  # bottom of first row of each group
  for (r in group_rows) tg <- hline(tg, r = r, lwd = 0.8)
  # bottom line
  tg <- hline(tg, r = n_rows, lwd = 1.5)
  as.ggplot(tg)
}

table_plot_list <- list()
for(t in 1:3){
  title <- paste0("Sensitivity at ", t,"y")
  data_f <- filter(data_figure2def_long, time==t)
  plot <- make_table(data_f, title)
  table_plot_list[[t]] <- plot
  assign(paste0("SE_main_", t,"y"), table_plot_list[[t]] )
  print(table_plot_list)
}

ggsave(paste0(path_to_save, "Figure2_modelperformance.eps"), 
       cowplot::plot_grid(ROC_main_1y+ggtitle("A) Discrimination over 1 year")+theme(plot.title = element_text(size = 17)),
                          ROC_main_2y+ggtitle("B) Discrimination over 2 years")+theme(plot.title = element_text(size = 17)),
                          ROC_main_3y+ggtitle("C) Discrimination over 3 years")+theme(plot.title = element_text(size = 17)),
                          
                          SE_main_1y+ggtitle("")+theme(plot.title = element_text(size = 17)),
                          SE_main_2y+ggtitle("")+theme(plot.title = element_text(size = 17)),
                          SE_main_3y+ggtitle("")+theme(plot.title = element_text(size = 17)),
                          ncol=3, nrow=2, align = "hv"),
       width = 25*0.8, height = 15*0.8,  bg = "white",
       device = "eps")


#######################################
## model threshold calculation:Full testing dataset
results_all <- list()
for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
  LC_event_time<-"LC_event1y"
  t<-1
  datax <- testing_1_imp %>% filter(!is.na(!!sym(LC_event_time)))
  ###USPSTF
  TP_uspstf <- datax %>% filter(!!sym(LC_event_time)==1 & !!sym(uspstfy)==1) %>% nrow()
  FN_uspstf <- datax %>% filter(!!sym(LC_event_time)==1 & !!sym(uspstfy)==0) %>% nrow() 
  FP_uspstf <- datax %>% filter(!!sym(LC_event_time)==0 & !!sym(uspstfy)==1) %>% pull(U19_weight) %>% sum() 
  TN_uspstf <- datax %>% filter(!!sym(LC_event_time)==0 & !!sym(uspstfy)==0) %>% pull(U19_weight) %>% sum()  
  
  TPR_uspstf  <- TP_uspstf/(TP_uspstf+FN_uspstf )
  FPR_uspstf  <- FP_uspstf/(FP_uspstf+TN_uspstf )
  TNR_uspstf <- TN_uspstf / (TN_uspstf + FP_uspstf)  # specificity (weighted)

  ####Model
  for (mod in c("plcom2012", "Integral_risk")){
    mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
    ROC_boot <- with(datax, WeightedROC(guess=get(mod_t), label=get(LC_event_time), 
                                        weight=U19_weight)) %>%as.data.frame() %>%mutate(TNR = 1 - FPR)
    
    # Match model specificity to USPSTF specificity (closest FPR)
    matched_row <- ROC_boot %>%
      mutate(fpr_diff = abs(FPR - FPR_uspstf)) %>%
      slice_min(fpr_diff, n = 1, with_ties = FALSE)
    
    results_all[[length(results_all) + 1]] <- data.frame(
      model            = mod,
      time             = t,
      uspstf_criterion = uspstfy,
      sens_model       = matched_row$TPR,
      spec_model       = matched_row$TNR,
      TP_model         = matched_row$TPR*(TP_uspstf+FN_uspstf),
      TN_model         = matched_row$TNR*(FP_uspstf+TN_uspstf),
      thredshold_model       = matched_row$threshold,
      sens_uspstf      = TPR_uspstf,
      spec_uspstf      = TNR_uspstf)
  } 
  bind_rows(results_all)
}
results_all_onetime <- bind_rows(results_all) %>% mutate(time=1)
write.csv(results_all_onetime,paste0(path_to_save,"Figure2_thresholdvalue.csv"),row.names = FALSE)

# ###########################################################################
#########Table 2. Discrimination and calibration of the INTEGRAL-Risk model for 1-year 
#lung cancer risk overall and across subgroups in the testing set.
#########eTable 8. Performance of the INTEGRAL-Risk model (over 1 year) stratified by 
#cohort in the testing set
#########eTable 10. Performance of the INTEGRAL-Risk model in the testing set for 
#predicting lung cancer occurring within 2 and 3 years, respectively, stratified 
#by participant characteristics and disease stage
# ###########################################################################
subgroup_AUC_function<- function(data){
  bootstraping_res_subgroups <- data.frame(NULL)
  modelstest_vec <- c("proteinmodelrisk", "plcom2012")
  for(imp in 1:15){
    print(imp)
    dataimp <- data %>% filter(name==imp) 
    for(t in 1:3){
      print(t)
      LC_event_t <- paste0("LC_event", t,"y")
      dataimp_t <- dataimp %>% filter(!is.na(!!sym(LC_event_t)))
      # Models ROC
      for(mod in modelstest_vec){ #ROC by model
        print(mod)
        mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
        
        boot_res_iter <- bootstrapped_AUC_cal(dataimp_t, mod_t, LC_event_t) %>% 
          mutate(model=mod, time=t, iter=imp)
        ### bind all bootstrap
        bootstraping_res_subgroups <- bind_rows(bootstraping_res_subgroups, boot_res_iter) # within each iteration of imputation: Bootstrap inference when using multiple imputation
      }}}
  return(bootstraping_res_subgroups)}

result_age_1<-subgroup_AUC_function( testing_all_imp%>% filter(age<=55)) %>% mutate(subgroup="age<=55")
result_age_2<-subgroup_AUC_function( testing_all_imp%>% filter(age>55&age<=65))%>% mutate(subgroup="age 56-65")
result_age_3<-subgroup_AUC_function( testing_all_imp%>% filter(age>65))%>% mutate(subgroup="age 65+")

result_race_1<-subgroup_AUC_function( testing_all_imp%>% filter(race==1))%>% mutate(subgroup="White")
result_race_2<-subgroup_AUC_function( testing_all_imp%>% filter(race==2))%>% mutate(subgroup="Black")
result_race_3<-subgroup_AUC_function( testing_all_imp%>% filter(race==4))%>% mutate(subgroup="Asia")

result_sex_1<-subgroup_AUC_function( testing_all_imp%>% filter(sex==1))%>% mutate(subgroup="Men")
result_sex_2<-subgroup_AUC_function( testing_all_imp%>% filter(sex==2))%>% mutate(subgroup="Women")

result_smoke_2<-subgroup_AUC_function( testing_all_imp%>% filter(smoke_status==2))%>% mutate(subgroup="Former")
result_smoke_3<-subgroup_AUC_function( testing_all_imp%>% filter(smoke_status==3))%>% mutate(subgroup="Current")

result_copd_1<-subgroup_AUC_function( testing_all_imp%>% filter(copd==1))%>% mutate(subgroup="COPD YES")
result_copd_0<-subgroup_AUC_function( testing_all_imp%>% filter(copd==0))%>% mutate(subgroup="COPD NO")

result_uspstf_1<-subgroup_AUC_function( testing_all_imp%>% filter(USPSTF2021==1))%>% mutate(subgroup="uspstf eligible")
result_uspstf_0<-subgroup_AUC_function( testing_all_imp%>% filter(USPSTF2021==0))%>% mutate(subgroup="uspstf ineligible")

result_plco_1<-subgroup_AUC_function( testing_all_imp%>% filter(plcom2012<0.0151))%>% mutate(subgroup="plcom2012<1.51%")
result_plco_2<-subgroup_AUC_function( testing_all_imp%>% filter(plcom2012<0.03&plcom2012>=0.0151))%>% mutate(subgroup="plcom2012 1.51%-2.9%")
result_plco_3<-subgroup_AUC_function( testing_all_imp%>% filter(plcom2012>=0.03))%>% mutate(subgroup="plcom2012>=3%")

result_GCS<-subgroup_AUC_function( testing_all_imp%>% filter(cohort=="Golestan"))%>% mutate(subgroup="Cohort-Golestan")
result_NYUWHS<-subgroup_AUC_function( testing_all_imp%>% filter(cohort=="NYUWHS"))%>% mutate(subgroup="Cohort-NYUWHS")
result_SCCS<-subgroup_AUC_function( testing_all_imp%>% filter(cohort=="SCCS"))%>% mutate(subgroup="Cohort-SCCS")
result_SCS<-subgroup_AUC_function( testing_all_imp%>% filter(cohort=="SCS"))%>% mutate(subgroup="Cohort-SCS")
result_SMHS<-subgroup_AUC_function( testing_all_imp%>% filter(cohort=="SMHS"))%>% mutate(subgroup="Cohort-SMHS")
result_WHI_FUP <-subgroup_AUC_function( testing_all_imp%>% filter(cohort="WHI_FUP"))%>% mutate(subgroup="Cohort-WHI_FUP")
result_WHS <-subgroup_AUC_function( testing_all_imp%>% filter(cohort=="WHS"))%>% mutate(subgroup="Cohort-WHS")

result_stage_1<-subgroup_AUC_function( testing_all_imp%>% filter(lc_early_late=="early"|LC_event==0))%>% mutate(subgroup="Early stage")
result_stage_0<-subgroup_AUC_function( testing_all_imp%>% filter(lc_early_late=="late"|LC_event==0))%>% mutate(subgroup="Late stage")

result_survival_1<-subgroup_AUC_function( testing_all_imp%>% filter(survival_time=="<0.5y"|LC_event==0))%>% mutate(subgroup="Survival<0.5y")
result_survival_2<-subgroup_AUC_function( testing_all_imp%>% filter(survival_time=="0.5-2y"|LC_event==0))%>% mutate(subgroup="Survival 0.5-2y")
result_survival_3<-subgroup_AUC_function( testing_all_imp%>% filter(survival_time=="2y+"|LC_event==0))%>% mutate(subgroup="Survival>=2y")


results_subgroup_all<-bind_rows(result_age_1,result_age_2,result_age_3,result_race_1,result_race_2,result_race_3,
                                result_sex_1,result_sex_2,result_smoke_2,result_smoke_3,result_copd_1,result_copd_0,result_uspstf_1,result_uspstf_0,
                                result_plco_1,result_plco_2,result_plco_3,result_GCS,result_NYUWHS,result_SCCS,result_SCS,result_SMHS,result_WHI_FUP,result_WHS,
                                result_stage_1,result_stage_0,
                                result_survival_1,result_survival_2,result_survival_3)

rm(result_age_1,result_age_2,result_age_3,result_race_1,result_race_2,result_race_3,
   result_sex_1,result_sex_2,result_smoke_2,result_smoke_3,result_copd_1,result_copd_0,result_uspstf_1,result_uspstf_0,
   result_plco_1,result_plco_2,result_plco_3,result_GCS,result_NYUWHS,result_SCCS,result_SCS,result_SMHS,result_WHI_FUP,result_WHS,
   result_stage_1,result_stage_0,result_survival_1,result_survival_2,result_survival_3)

###Confidence interval calculation
mean_ci <- function(col){
  mean <- mean(col, na.rm=TRUE)
  lowerCIpooled <- quantile(col, probs=0.025, na.rm=TRUE)
  upperCIpooled <- quantile(col, probs=0.975, na.rm=TRUE)
  meanci <- paste0(round(mean,3), " (", round(lowerCIpooled,3), "-", round(upperCIpooled,3), ")")
  return(meanci)
}
results_subgroup_all %>% group_by(subgroup, time,model) %>% 
  dplyr::summarise(across(c("ncase_overall", "ncases_boot",
                            "pop_overall","totpop_boot",
                            "bootstrappedAUC_vec" , "bootstrappedcal_vec"), 
                          ~mean_ci(.x))) 


# Make AUC summary table for stratified performance  
strat_group <- unique(results_subgroup_all$subgroup)
summary_perf_allstrat <- data.frame(NULL)
for(g in strat_group){
  interim <- summary_res_pooled_ROC(filter(results_subgroup_all, subgroup==g),c("proteinmodelrisk", "plcom2012")) %>% 
    mutate(subgroup=g)
  summary_perf_allstrat <- dplyr::bind_rows(summary_perf_allstrat, interim)
}

# Make Calibration summary table
cal_stratified_clean <- results_subgroup_all %>%  group_by(model, time,subgroup) %>% 
  filter(model=="proteinmodelrisk") %>% 
  dplyr::summarise(mean_EO=mean(as.numeric(cal_overall)),
                   lowerCI_pooled=quantile(bootstrappedcal_vec, prob=c(0.025), na.rm = TRUE)[[1]],
                   upperCI_pooled=quantile(bootstrappedcal_vec, prob=c(0.975), na.rm = TRUE)[[1]]) %>% 
  mutate(cal_CI=paste0(formatC(mean_EO, digits=2, format = "f")," (", 
                       formatC(lowerCI_pooled, digits = 2, format = "f"), "-",
                       formatC(upperCI_pooled, digits=2, format="f"), ")")) %>% ungroup() %>% 
  dplyr::select(model, time, subgroup, cal_CI)
summary_perf_allstrat <- summary_perf_allstrat %>% merge(cal_stratified_clean, by=c("model", "time", "subgroup"), all = T)

# add number of cases observed and number of cases bootstrapped
data_var_IR_boot <- results_subgroup_all %>% group_by(subgroup, time,model) %>% 
  dplyr::summarise(across(c("ncase_overall", "ncases_boot",
                            "pop_overall","totpop_boot"), 
                          ~mean_ci(.x))) 

summary_perf_allstrat <- summary_perf_allstrat %>% merge(data_var_IR_boot, by=c("subgroup", "time", "model"), all = T)

summary_perf_allstrat %>% head()
# wider for save 
summary_perf_allstrat_clean <- summary_perf_allstrat %>% 
  dplyr::select(model,ncase_overall,pop_overall, ncases_boot, totpop_boot,
                subgroup,time,AUC_CI_pooled, pvalue, pvalue_ref, cal_CI) %>% 
  pivot_wider(id_cols=c(subgroup), names_from = c("time", "model"), 
              values_from = c("ncase_overall","pop_overall", "ncases_boot", "totpop_boot",
                              "AUC_CI_pooled", "pvalue", "pvalue_ref","cal_CI")) 

summary_perf_allstrat_clean%>% write.csv(file=paste0(path_to_save,"Table2_stratifiedperformance.csv" ), row.names = F)

# ###########################################################################
### Figure 3 ####
# Reclassification of screening eligibility using the INTEGRAL-Risk model 
# compared to USPSTF-2021 and PLCOm2012.
# ###########################################################################
library(forcats)
library(ggsankey)
library(ggalluvial)

thresholds_1y <- getdata_lorezplot(testing_1_imp, "proteinmodelrisk_1y", "LC_event1y",1)
threshold_est1y <- thresholds_1y[[3]]
#### total population
#########USPSTF---INTEGRAL
testing_1_imp <- testing_1_imp %>% mutate(elig_USPSTF2021=ifelse(USPSTF2021==1,"EligUSPSTF", "Non elig"))

summary_alluvial_baselinepop_a <- testing_1_imp %>% filter(LC_event==0)  %>% 
  group_by(elig_USPSTF2021) %>% 
  dplyr::summarise(`Eligible by INTEGRAL risk model`=sum(U19_weight[proteinmodelrisk_1y >= thresholds_1y[[3]]], na.rm = TRUE),
                   `Non eligible by INTEGRAL risk model`=sum(U19_weight[proteinmodelrisk_1y< thresholds_1y[[3]]], na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_longer(c("Eligible by INTEGRAL risk model", "Non eligible by INTEGRAL risk model"), names_to = "integral_elig", values_to = "freq")

senkey_baselinepop_a <- summary_alluvial_baselinepop_a %>%
  mutate(
    color_s1=case_when(grepl("Non elig",elig_USPSTF2021) & grepl("Non elig", integral_elig)~"Inelig_s1",
                       !grepl("Non elig",elig_USPSTF2021) & grepl("Non elig", integral_elig)~"downgrade_s1",
                       grepl("Non elig",elig_USPSTF2021) & !grepl("Non elig", integral_elig)~"upgrade_s1",
                       !grepl("Non elig",elig_USPSTF2021) & !grepl("Non elig", integral_elig)~"elig_s1"),
    elig_USPSTF2021=ifelse(elig_USPSTF2021=="Non elig","Not\neligible", "25%\neligible"),
    elig_USPSTF2021=factor(elig_USPSTF2021, levels=c("25%\neligible","Not\neligible")),
    integral_elig=ifelse(grepl("Non elig", integral_elig), "  Not  \n  eligible  ", "  25%  \n  eligible  "), 
    integral_elig=factor(integral_elig, levels=c("  25%  \n  eligible  ","  Not  \n  eligible  "))
  ) %>%
  ggplot(aes(y = freq,axis=elig_USPSTF2021, axis2 =  integral_elig)) +
  geom_alluvium(aes(fill=color_s1), width = 1/5, curve_type = "cubic") + 
  geom_stratum(aes(fill=after_stat(stratum)), width = 1/3, color = "black") + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            color="black", size=6,
            hjust=0.5)+
  scale_fill_manual(values=c("elig_s1" = "#5DA5DA", "Inelig_s1" = "gray", 
                             "upgrade_s1" = "#60BD68", "downgrade_s1" = "#FAA43A",
                             "25%\neligible"="#e18727ff", "Not\neligible"="gray",
                             "  25%  \n  eligible  "="#bc3c29ff", "  Not  \n  eligible  "="gray"))+
  scale_x_discrete(limits = c("USPSTF2021","INTEGRAL\nRisk Model"), 
                   expand = c(0.25,0.25))+
  theme_alluvial()+theme(axis.text.y  = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "none", 
                         axis.text.x = element_text(size=15, face="bold", color="black"),
                         plot.title = element_text(hjust = 0.5, vjust=0,face="bold", size=20)) +
  labs(title="A) Screening eligibility for the baseline-population based on\nUSPSTF-2021 and the INTEGRAL-Risk model")
senkey_baselinepop_a
###Generate as Table
summary_alluvial_baselinepop_a_wide<-summary_alluvial_baselinepop_a%>% pivot_wider(
  names_from = integral_elig  ,values_from = c(freq))
write.csv(summary_alluvial_baselinepop_a_wide,paste0(path_to_save,"summary_alluvial_baselinepop_a_wide.csv"),row.names = FALSE)


#########PLCO---INTEGRAL
summary_alluvial_baselinepop_b <- testing_1_imp %>% filter(LC_event==0) %>% 
  mutate(elig_plco=ifelse(plcom2012<thresholds_1y[[4]], "Non eligible by PLCOm2012", "Eligible by PLCOm2012"), 
         plco_riskgroups=case_when(plcom2012<thresholds_1y[[4]]~"<0.92%",
                                   plcom2012>=thresholds_1y[[4]] ~"0.92%+"), 
         plco_riskgroups=factor(plco_riskgroups, levels=c("<0.92%","0.92%+") )) %>% 
  group_by(plco_riskgroups , elig_plco) %>% 
  dplyr::summarise(`Eligible by INTEGRAL risk model`=sum(U19_weight[proteinmodelrisk_1y >= thresholds_1y[[3]]], na.rm = TRUE),
                   `Non eligible by INTEGRAL risk model`=sum(U19_weight[proteinmodelrisk_1y< thresholds_1y[[3]]], na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_longer(c("Eligible by INTEGRAL risk model", "Non eligible by INTEGRAL risk model"), names_to = "integral_elig", values_to = "freq")


senkey_baselinepop_b <- summary_alluvial_baselinepop_b %>%
  mutate(plco_riskgroups=fct_rev(plco_riskgroups), 
         color_s2=case_when(grepl("<0.92%",plco_riskgroups) & grepl("Non elig", integral_elig)~"Inelig_s2", 
                            !grepl("<0.92%",plco_riskgroups) & grepl("Non elig", integral_elig)~"downgrade_s2", 
                            grepl("<0.92%",plco_riskgroups) & !grepl("Non elig", integral_elig)~"upgrade_s2", 
                            !grepl("<0.92%",plco_riskgroups) & !grepl("Non elig", integral_elig)~"elig_s2",
                            T~"error"), 
         plco_riskgroups=ifelse(plco_riskgroups=="<0.92%", " Not \n eligible ", " 25% \n eligible "),
         plco_riskgroups=factor(plco_riskgroups, levels=c(" 25% \n eligible "," Not \n eligible ")),
         integral_elig=ifelse(grepl("Non elig", integral_elig), "  Not  \n  eligible  ", "  25%  \n  eligible  "), 
         integral_elig=factor(integral_elig, levels=c("  25%  \n  eligible  ","  Not  \n  eligible  "))
  ) %>%
  ggplot(aes(y = freq,axis=plco_riskgroups ,axis2 = integral_elig)) +
  geom_alluvium(aes(fill=color_s2), width = 1/5, curve_type = "cubic") + 
  geom_stratum(aes(fill=after_stat(stratum)), width = 1/3, color = "black") + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            color="black", size=6,
            hjust=0.5)+
  scale_fill_manual(values=c("elig_s2" = "#5DA5DA", "Inelig_s2" = "gray", 
                             "upgrade_s2" = "#60BD68", "downgrade_s2" = "#FAA43A",
                             " 25% \n eligible "="#0072b5ff"," Not \n eligible "="gray",
                             "  25%  \n  eligible  "="#bc3c29ff", "  Not  \n  eligible  "="gray"))+
  scale_x_discrete(limits = c("PLCOm2012", "INTEGRAL\nRisk Model"), 
                   expand = c(0.25,0.25))+
  theme_alluvial()+theme(axis.text.y  = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "none", 
                         axis.text.x = element_text(size=15, face="bold", color="black"),
                         plot.title = element_text(hjust = 0.5, vjust=0,face="bold", size=20))+
  labs(title="C) Screening eligibility for the baseline-population based on\nthe PLCOm2012  and the INTEGRAL-Risk models")

senkey_baselinepop_b

###Generate as Table
summary_alluvial_baselinepop_b_wide<-summary_alluvial_baselinepop_b %>% dplyr::select(-c("plco_riskgroups"))%>%  pivot_wider(
  names_from = integral_elig  ,values_from = c(freq))
write.csv(summary_alluvial_baselinepop_b_wide,paste0(path_to_save,"summary_alluvial_baselinepop_b_wide.csv"),row.names = FALSE)


#######1-year lung cancer case
#########USPSTF----INTEGRAL

summary_alluvial_a <- testing_1_imp %>% filter(LC_event1y==1) %>% 
  group_by(elig_USPSTF2021) %>% 
  dplyr::summarise(`Eligible by INTEGRAL risk model`=sum(proteinmodelrisk_1y>=thresholds_1y[[3]]),
                   `Non eligible by INTEGRAL risk model`=sum(proteinmodelrisk_1y<thresholds_1y[[3]])) %>% 
  ungroup() %>% 
  pivot_longer(c("Eligible by INTEGRAL risk model", "Non eligible by INTEGRAL risk model"), names_to = "integral_elig", values_to = "freq")

senkey_cases_a <- summary_alluvial_a %>%
  mutate(
    color_s1=case_when(grepl("Non elig",elig_USPSTF2021) & grepl("Non elig", integral_elig)~"Inelig_s1",
                       !grepl("Non elig",elig_USPSTF2021) & grepl("Non elig", integral_elig)~"downgrade_s1",
                       grepl("Non elig",elig_USPSTF2021) & !grepl("Non elig", integral_elig)~"upgrade_s1",
                       !grepl("Non elig",elig_USPSTF2021) & !grepl("Non elig", integral_elig)~"elig_s1"), 
    elig_USPSTF2021=ifelse(elig_USPSTF2021=="Non elig","Not\neligible", "63%\neligible"),
    elig_USPSTF2021=factor(elig_USPSTF2021, levels=c("63%\neligible","Not\neligible")),
    integral_elig=ifelse(grepl("Non elig", integral_elig), "  Not  \n  eligible  ", "85%\neligible"), 
    integral_elig=factor(integral_elig, levels=c("85%\neligible","  Not  \n  eligible  "))
  ) %>% 
  ggplot(aes(y = freq,axis=elig_USPSTF2021, axis2 = integral_elig)) +
  geom_alluvium(aes(fill=color_s1), width = 1/5, curve_type = "cubic") + 
  geom_stratum(aes(fill=after_stat(stratum)), width =  1/3, color = "black") + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            color="black", size=6 , hjust=0.5)+
  scale_fill_manual(values=c("elig_s1" = "#5DA5DA", "Inelig_s1" = "gray", 
                             "upgrade_s1" = "#60BD68", "downgrade_s1" = "#FAA43A", 
                             "Not\neligible"="gray",
                             "  Not  \n  eligible  "="gray",
                             "63%\neligible"="#e18727ff", 
                             "85%\neligible"="#bc3c29ff"))+
  scale_x_discrete(limits = c("USPSTF2021", "INTEGRAL\nRisk Model"), 
                   expand = c(0.25,0.25))+theme_alluvial()+theme(axis.text.y  = element_blank(),
                                                                 axis.ticks = element_blank(),
                                                                 axis.title.y = element_blank(),
                                                                 legend.position = "none", 
                                                                 axis.text.x = element_text(size=15, face="bold", color="black"),
                                                                 plot.title = element_text(hjust = 0.5, vjust=0,
                                                                                           face="bold", size=20))+
  labs(title="B) Screening eligibility for incident lung cancer cases based on\nUSPSTF-2021 and the INTEGRAL-Risk model")

senkey_cases_a
###Generate as Table
summary_alluvial_a_wide<-summary_alluvial_a %>%  pivot_wider(
  names_from = integral_elig ,values_from = c(freq))
write.csv(summary_alluvial_a_wide,paste0(path_to_save,"summary_alluvial_a_wide.csv"),row.names = FALSE)

#########PLCO----INTEGRAL
summary_alluvial_b <- testing_1_imp %>% filter(LC_event1y==1) %>% 
  mutate(elig_plco=ifelse(plcom2012<thresholds_1y[[4]], "Non eligible by PLCOm2012", "Eligible by PLCO 2012"), 
         plco_riskgroups=case_when(plcom2012<thresholds_1y[[4]]~"<0.92%",
                                   plcom2012>=thresholds_1y[[4]] ~"0.92%+"), 
         plco_riskgroups=factor(plco_riskgroups, levels=c("<0.92%","0.92%+") )) %>% 
  group_by(plco_riskgroups,elig_plco) %>% 
  dplyr::summarise(`Eligible by INTEGRAL risk model`=sum(proteinmodelrisk_1y>=thresholds_1y[[3]]),
                   `Non eligible by INTEGRAL risk model`=sum(proteinmodelrisk_1y<thresholds_1y[[3]])) %>% 
  ungroup() %>% 
  pivot_longer(c("Eligible by INTEGRAL risk model", "Non eligible by INTEGRAL risk model"), names_to = "integral_elig", values_to = "freq")


senkey_cases_b <- summary_alluvial_b %>%
  mutate(color_s2=case_when(grepl("<0.92%",plco_riskgroups) & grepl("Non elig", integral_elig)~"Inelig_s2", 
                            !grepl("<0.92%",plco_riskgroups) & grepl("Non elig", integral_elig)~"downgrade_s2", 
                            grepl("<0.92%",plco_riskgroups) & !grepl("Non elig", integral_elig)~"upgrade_s2", 
                            !grepl("<0.92%",plco_riskgroups) & !grepl("Non elig", integral_elig)~"elig_s2",
                            T~"error"), 
         plco_riskgroups=ifelse(plco_riskgroups=="<0.92%", " Not \n eligible ", "70%\neligible"),
         plco_riskgroups=factor(plco_riskgroups, levels=c("70%\neligible"," Not \n eligible ")),
         integral_elig=ifelse(grepl("Non elig", integral_elig), "  Not  \n  eligible  ", "85%\neligible"), 
         integral_elig=factor(integral_elig, levels=c("85%\neligible","  Not  \n  eligible  "))
  ) %>% 
  ggplot(aes(y = freq,axis=plco_riskgroups, axis2 = integral_elig)) +
  geom_alluvium(aes(fill=color_s2), width = 1/5, curve_type = "cubic") + 
  geom_stratum(aes(fill=after_stat(stratum)), width = 1/3, color = "black") + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            color="black", size=6 , hjust=0.5)+
  scale_fill_manual(values=c("elig_s2" = "#5DA5DA", "Inelig_s2" = "gray", 
                             "upgrade_s2" = "#60BD68", "downgrade_s2" = "#FAA43A", 
                             " Not \n eligible "="gray",
                             "  Not  \n  eligible  "="gray",
                             "70%\neligible"="#0072b5ff",
                             "85%\neligible"="#bc3c29ff"))+
  scale_x_discrete(limits = c("PLCOm2012", "INTEGRAL\nRisk Model"), 
                   expand = c(0.25,0.25))+theme_alluvial()+theme(axis.text.y  = element_blank(),
                                                                 axis.ticks = element_blank(),
                                                                 axis.title.y = element_blank(),
                                                                 legend.position = "none", 
                                                                 axis.text.x = element_text(size=15, face="bold", color="black"),
                                                                 plot.title = element_text(hjust = 0.5, vjust=0,
                                                                                           face="bold", size=20))+
  labs(title="D) Screening eligibility for incident lung cancer cases based on\nthe PLCOm2012 and the INTEGRAL-Risk models")

senkey_cases_b

###Generate as Table
summary_alluvial_b_wide<-summary_alluvial_b %>%dplyr::select(-plco_riskgroups) %>%  pivot_wider(
  names_from = integral_elig ,values_from = c(freq))
write.csv(summary_alluvial_b_wide,paste0(path_to_save,"summary_alluvial_b_wide.csv"),row.names = FALSE)

cowplot::plot_grid(senkey_baselinepop_a, senkey_cases_a, 
                   NULL, NULL,   # ← spacer row
                   senkey_baselinepop_b, senkey_cases_b, 
                   NULL, NULL,   
                   ncol = 2, 
                   rel_heights = c(1, 0.4, 1,0.5))

####################################################
##Add table in the plot
##This step (line 884-1000) is pure visualization using the results generated above.
summary_alluvial_baselinepop_a_wide
summary_alluvial_a_wide
data_1<-data.frame("USPSTF-2021"=c("Eligible","Eligible","Eligible",
                                   "Not eligible","Not eligible","Not eligible",
                                   "Total","Total","Total"),
                   "INTEGRAL-Risk model"=c("Eligible","Not eligible","Total",
                                           "Eligible","Not eligible","Total",
                                           "Eligible","Not eligible","Total"),
                   value_population=c("16998 (13%)","15526 (12%)","35224 (25%)",
                                      "15509 (12%)","82743 (63%)","98252 (75%)",
                                      "35207 (25%)","98268 (75%)","130776 (100%)"),
                   value_case=c("104 (57%)","10 (6%)","114 (63%)",
                                "50 (28%)","17 (9%)","67 (37%)",
                                "154 (85%)","27 (15%)","181 (100%)"),check.names = FALSE)

data_1$fill_group <- c(
  "blue","orange","dark_orange",
  "green","grey","grey",
  "red","grey","grey")

data_1$`USPSTF-2021`<-factor(data_1$`USPSTF-2021`,levels=c("Total","Not eligible","Eligible"))

p_population_uspstf<-ggplot(data_1, aes(x = `INTEGRAL-Risk model` , y =`USPSTF-2021`   ))+
  geom_tile(aes(fill = fill_group), color = "black")+
  geom_text(aes(label = value_population), size = 4) +
  scale_fill_manual(values = c("blue" = "#5DA5DA","orange" = "#FAA43A","dark_orange" = "#e18727ff",
                               "green" = "#60BD68","red" = "#bc3c29ff","grey" = "gray")) +
  scale_x_discrete(position = "top", expand = c(0,0))+
  theme_minimal() +
  theme(legend.position = "none",  aspect.ratio = 0.2 ,
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12))

p_case_uspstf<-ggplot(data_1, aes(x = `INTEGRAL-Risk model` , y =`USPSTF-2021`   ))+
  geom_tile(aes(fill = fill_group), color = "black")+
  geom_text(aes(label = value_case), size = 4) +
  scale_fill_manual(values = c("blue" = "#5DA5DA","orange" = "#FAA43A","dark_orange" = "#e18727ff",
                               "green" = "#60BD68","red" = "#bc3c29ff","grey" = "gray")) +
  scale_x_discrete(position = "top", expand = c(0,0))+
  theme_minimal() +
  theme(legend.position = "none",  aspect.ratio = 0.2 ,
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12))

###PLCOm2012
summary_alluvial_baselinepop_b_wide
summary_alluvial_b_wide
data_2<-data.frame("PLCOm2012"=c("Eligible","Eligible","Eligible",
                                 "Not eligible","Not eligible","Not eligible",
                                 "Total","Total","Total"),
                   "INTEGRAL-Risk model"=c("Eligible","Not eligible","Total",
                                           "Eligible","Not eligible","Total",
                                           "Eligible","Not eligible","Total"),
                   value_population=c("19565 (15%)","12934 (10%)","32499 (25%)",
                                      "12943 (10%)","85334 (65%)","98277 (75%)",
                                      "35207 (25%)","98268 (75%)","130776 (100%)"),
                   value_case=c("116 (64%)","11 (6%)","127 (70%)",
                                "38 (21%)","16 (9%)","54 (30%)",
                                "154 (85%)","27 (15%)","181 (100%)"),check.names = FALSE)

data_2$fill_group <- c(
  "blue","orange","dark_blue",
  "green","grey","grey",
  "red","grey","grey")

data_2$PLCOm2012<-factor(data_2$PLCOm2012,levels=c("Total","Not eligible","Eligible"))

p_population_plco<-ggplot(data_2, aes(x = `INTEGRAL-Risk model` , y =PLCOm2012))+
  geom_tile(aes(fill = fill_group), color = "black")+
  geom_text(aes(label = value_population), size = 4) +
  scale_fill_manual(values = c("blue" = "#5DA5DA","orange" = "#FAA43A","dark_blue" = "#0072b5ff",
                               "green" = "#60BD68","red" = "#bc3c29ff","grey" = "gray")) +
  scale_x_discrete(position = "top", expand = c(0,0))+
  theme_minimal() +
  theme(legend.position = "none",  aspect.ratio = 0.2 ,
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12))

p_case_plco<-ggplot(data_2, aes(x = `INTEGRAL-Risk model` , y =PLCOm2012   ))+
  geom_tile(aes(fill = fill_group), color = "black")+
  geom_text(aes(label = value_case), size = 4) +
  scale_fill_manual(values = c("blue" = "#5DA5DA","orange" = "#FAA43A","dark_blue" = "#0072b5ff",
                               "green" = "#60BD68","red" = "#bc3c29ff","grey" = "gray")) +
  scale_x_discrete(position = "top", expand = c(0,0))+
  theme_minimal() +
  theme(legend.position = "none",  aspect.ratio = 0.2 ,
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12))

cowplot::plot_grid(senkey_baselinepop_a, senkey_cases_a, 
                   p_population_uspstf, p_case_uspstf,   # ← spacer row
                   senkey_baselinepop_b, senkey_cases_b, 
                   p_population_plco, p_case_plco,   
                   ncol = 2,
                   rel_heights = c(1, 0.5, 1,0.5))

ggsave(paste0(path_to_save, "Figure3_panel_reclassification.eps"), 
       cowplot::plot_grid(senkey_baselinepop_a, senkey_cases_a, 
                          p_population_uspstf, p_case_uspstf,   # ← spacer row
                          senkey_baselinepop_b, senkey_cases_b, 
                          p_population_plco, p_case_plco,   
                          ncol = 2,
                          rel_heights = c(1, 0.5, 1,0.5)),
       width = 20*0.88, height = 18*0.9,  bg = "white",
       device = "eps")

# ##############################################################################
# ##############################################################################
# #########################Supplemental Online Content
# ##############################################################################
# ##############################################################################
#Check model parameters 
# ###########################################################################
#eTable 6. Model parameters (Beta estimates) for the INTEGRAL-Risk model (M1).
# ###########################################################################
summary(INTEGRAL_RISK_AIC)

# ###########################################################################
#eTable 7. Model parameters for additional models evaluated in the testing set
# ###########################################################################
#summary(INTEGRAL_RISK_AUC)
#summary(INTEGRAL_RISK_AICAsia)
#summary(INTEGRAL_RISK_AUCAsia)

# ###########################################################################
# eTable 8 Performance of the INTEGRAL-Risk model (over 1 year) stratified by cohort in the testing set
# eTable 10. Performance of the INTEGRAL-Risk model in the testing set for predicting lung cancer
#            occurring within 2 and 3 years, respectively, stratified by participant characteristics and disease stage
# Please see script for Table 2.
# ###########################################################################


# ###########################################################################
#eTable 11. Characteristics of study participants selected as screening eligible in the baseline
#           population and in lung cancer cases diagnosed over one, two or three years of follow-up. Risk
#           thresholds for the PLCOm2012 and INTEGRAL-Risk model were set to select the same baseline
#           population size as the USPSTF2021 screening criteria
# ###########################################################################
testing_1_imp <- testing_1_imp %>% mutate(elig_USPSTF2021=ifelse(USPSTF2021==1,"EligUSPSTF", "Non elig"),
                                          PLCO_08=ifelse(plcom2012>=thresholds_1y[[4]], "eligPLCO_08", "Non elig"),
                                          elig_prot1y=ifelse(proteinmodelrisk_1y>=thresholds_1y[[3]], "Elig_prot1y", "Non elig"),
                                          elig_prot2y=ifelse(proteinmodelrisk_1y>=thresholds_1y[[3]], "Elig_prot2y", "Non elig"),
                                          elig_prot3y=ifelse(proteinmodelrisk_1y>=thresholds_1y[[3]], "Elig_prot3y", "Non elig"))

eTable11_fulleligible_description <- 
  bind_rows(baselineelig_1=summary_eligibility_charac(testing_1_imp,1,"baseline",0, T),
            Cases_1y=summary_eligibility_charac(testing_1_imp,1,"cases",0,T),
            Cases_2y=summary_eligibility_charac(testing_1_imp,2,"cases",0,T),
            Cases_3y=summary_eligibility_charac(testing_1_imp,3,"cases",0,T),
            .id="Set") %>% 
  mutate(elig_by=factor(elig_by, levels=c("elig_USPSTF2021", "PLCO_08", "elig_prot1y","elig_prot2y","elig_prot3y")),
         Set=factor(Set, levels=c("baselineelig_1","Cases_1y","Cases_2y","Cases_3y"))) %>% 
  arrange(Set, elig_by)%>% 
  dplyr::select(Set,elig_by,N_cases_screened_per , Age, CPD, `Years smoked`,`per_former`,quit_years_former) %>% 
  dplyr::rename(`Screened by`=elig_by,`Former smokers (%)`=`per_former`, `Quit years\n(Former)`=quit_years_former)
write.csv(eTable11_fulleligible_description,paste0(path_to_save,"eTable11_eligiblepopulation_description.csv"),row.names = FALSE)

# ###########################################################################
##eFigure 2. Depiction of the weights that each predictor exerts in 
##the multivariate INTEGRAL-Risk model and include hazard ratios after 
##refitting the model with each parameter scaled to have a standard deviation 
##of one in the training set
# ###########################################################################
###1_ centering scaled,  no weighting
scale_this_center <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

protein_names_log <- c("MMP12_log", "TRAILR2_log", "SCF_log", "LAMP3_log", "MUC16_log", "KRT19_log", "SYND1_log", "CEACAM5_log", 
                       "FASLG_log", "IL6_log", "CLEC4D_log", "WFDC2_log", "LPL_log", "CDCP1_log", "SCGB3A2_log",
                       "CXCL13_log", "ALPP_log", "TFPI2_log", "CXL17_log", "CHI3L1_log", "CXCL9_log")

###4 no centering scaled,  subcohort SD, with weighting
data_scaled_subcohort_weight <- data_train %>%
  mutate(across(all_of(protein_names_log),
                ~ {s <- sqrt(Hmisc::wtd.var(.x[subcohort_status == 5], 
                                            weights = U19_weight[subcohort_status == 5], na.rm = TRUE))
                .x / s})) %>% ungroup()

### age, intensity, years_smoked scale as well
other_variables<-c("age","intensity","years_smoked")
data_scaled_subcohort_weight <- data_scaled_subcohort_weight %>%
  mutate(across(all_of(other_variables),
                ~ {s <- sqrt(Hmisc::wtd.var(.x[subcohort_status == 5], 
                                            weights = U19_weight[subcohort_status == 5], na.rm = TRUE))
                .x / s})) %>% ungroup()
weighted_var <- function(x, w) {
  wm <- weighted.mean(x, w)
  sum(w * (x - wm)^2) / sum(w)
}

model_3y<- rstpm2::stpm2(Surv(t0,tevent, LC_event)~
                           age+intensity+years_smoked+TRAILR2_log+WFDC2_log+CEACAM5_log+LAMP3_log+MMP12_log+LPL_log+SCF_log+
                           CXCL9_log+CXL17_log+FASLG_log+SYND1_log+CXCL13_log+IL6_log, 
                         weights=U19_weight, tvc=c(IL6_log=2),df=3, 
                         control = list(optimiser = "NelderMead"),robust = T, 
                         data=data_scaled_subcohort_weight)
##IL6 harzard spline
plot(model_3y,newdata =data.frame( IL6_log =0,TRAILR2_log=0,WFDC2_log=0,CEACAM5_log=0,LAMP3_log=0,MMP12_log=0,LPL_log=0,SCF_log=0,
                                   CXCL9_log=0,CXL17_log=0,FASLG_log=0,SYND1_log=0,CXCL13_log=0,age=0,intensity=0,years_smoked=0),
     type = 'hr',var ='IL6_log',ci = TRUE, rug = FALSE,  main="IL6_log")
abline(h=1, lty=2)    

hr_IL6_123 <- predict(model_3y,type = "hr", var = "IL6_log",
                      newdata = data.frame(t0=0, tevent    = c(1, 2, 3),
                                           IL6_log = 0,TRAILR2_log=0,WFDC2_log=0,CEACAM5_log=0,LAMP3_log=0,MMP12_log=0,LPL_log=0,SCF_log=0,
                                           CXCL9_log=0,CXL17_log=0,FASLG_log=0,SYND1_log=0,CXCL13_log=0,age=0,intensity=0,years_smoked=0),se.fit  = TRUE)
hr_IL6_123
hr_IL6_123$Terms<-c("IL6-1y","IL6-2y","IL6-3y")
hr_IL6_123<-hr_IL6_123 %>% rename(HR=Estimate)

####Get dataset for forest plot
output_model_scale<-summary(model_3y)@coef
forest_data_scale_update<- as.data.frame(output_model_scale)
forest_data_scale_update$Terms<-row.names(forest_data_scale_update)
row.names(forest_data_scale_update)<-NULL
forest_data_scale_update<-forest_data_scale_update %>% filter(Terms %in%c(other_variables,protein_names_log)) %>% filter(Terms!="IL6_log")
forest_data_scale_update<-forest_data_scale_update %>% mutate(HR=exp(Estimate),lower = exp(Estimate - 1.96 * `Std. Error`),
                                                              upper = exp(Estimate + 1.96 *`Std. Error`)) %>% dplyr::select(Terms,HR, lower, upper)
###Add IL-6 HR
forest_data_scale_update<-bind_rows(forest_data_scale_update,hr_IL6_123)

###Weighted mean and SD
variable_figure2<-c("age","intensity","years_smoked",protein_names_log)
weighted_stats <- lapply(variable_figure2, function(var) {
  data_train %>%summarise(Terms= var, 
                          wmean= weighted.mean(.data[[var]], U19_weight),
                          wsd  = sqrt(weighted_var(.data[[var]], U19_weight)))}) 

weighted_stats<-weighted_stats%>%bind_rows()  %>%mutate(mean_sd_label = sprintf("%.1f (%.1f)", wmean, wsd))
###Change the name of IL6 for the merging purpose
weighted_stats<-weighted_stats %>% mutate(Terms=ifelse(Terms=="IL6_log","IL6-1y",Terms))

# --- Merge the weighted mean (SD) into forest plot data ---
forest_data_scale_update <- forest_data_scale_update %>%
  left_join(weighted_stats, by = "Terms")

forest_data_scale_update <- forest_data_scale_update %>%
  mutate(hr_ci_label = sprintf("%.2f (%.2f-%.2f)", HR, lower, upper))

##Change the name for better visulizaion
forest_data_scale_update <- forest_data_scale_update %>%mutate(Terms = sub("_log$", "", Terms))
forest_data_scale_update <- forest_data_scale_update %>%mutate(Terms=ifelse(Terms=="age","Age",Terms))
forest_data_scale_update <- forest_data_scale_update %>%mutate(Terms=ifelse(Terms=="intensity","CPD",Terms))
forest_data_scale_update <- forest_data_scale_update %>%mutate(Terms=ifelse(Terms=="years_smoked","Years smoked",Terms))

###Arrange the order of proteins:
forest_data_scale_update<-forest_data_scale_update %>% arrange(desc(HR))
rev(unique(forest_data_scale_update$Terms))
forest_data_scale_update<-forest_data_scale_update%>%
  mutate(Terms = factor(Terms, levels = c("TRAILR2", "LAMP3", "LPL", "SCF" ,  "FASLG","CXCL13","SYND1","IL6-3y","IL6-2y","IL6-1y","CXL17", "CXCL9",
                                          "MMP12","WFDC2", "CEACAM5","CPD","Years smoked","Age" )))
# --- Build plot ---
forest_data_scale_update<-forest_data_scale_update %>% arrange(desc(HR))
rev(unique(forest_data_scale_update$Terms))

efigure2_HR_scale<-ggplot(forest_data_scale_update, aes(y = Terms, x = HR)) +
  annotate("text", x = 0.3, y = Inf,label = "Mean (SD)", fontface = "bold", size = 4, vjust = -0.5, hjust = 0.5) +
  annotate("text", x = 4, y = Inf,label = "HR (95% CI)", fontface = "bold", size = 4, vjust = -0.5, hjust = 0.5) +
  # Mean (SD) labels on the left
  geom_text(aes(x = 0.3, label = mean_sd_label),
            hjust = 0.5, size = 3.8, color = "black") +
  # HR (95% CI) labels on the right
  geom_text(aes(x = 4, label = hr_ci_label),
            hjust = 0.5, size = 3.8, color = "black") +
  # --- Variable name labels (replaces y-axis text) ---
  geom_text(aes(x = 0.2, label = Terms),
            hjust = 1, size = 3.8, color = "black") +
  # Forest plot elements
  geom_point(size = 3, color = "darkred") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, color = "darkred") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed",colour = "gray", linewidth = 0.3) +
  geom_vline(xintercept = 1.5, linetype = "dashed",colour = "gray", linewidth = 0.3) +
  geom_vline(xintercept = 2, linetype = "dashed",colour = "gray", linewidth = 0.3) +
  scale_x_log10(name   = "Hazard Ratio (log scale) for developing lung cancer over three years",
                breaks = c(0.5, 1, 1.5, 2)) +
  coord_cartesian(clip = "off") +  
  ylab(NULL) +
  theme_minimal(base_size = 16) +
  theme( panel.grid.minor  = element_blank(),
         panel.grid.major  = element_blank(),
         axis.ticks.y      = element_blank(),
         axis.text.y       = element_blank(),  
         axis.line.x.bottom = element_line(colour = "gray", linewidth = 0.6), 
         plot.margin       = margin(t = 20, r = 160, b = 10, l = 120))
ggsave(efigure2_HR_scale,filename = paste0(path_to_save,"Supple_Figure2_HR_scale.pdf"),height = 6, width=10)


####################################################################################
#######eFigure 3. Risk-discriminative performance and calibration of the 
##INTEGRAL-Risk model in the LC3 testing set for predicting lung cancer occurring 
##within 0 and 1 year of follow-up, within 1 and 2 years of follow-up, and within 
##2 and 3 years of follow-up
####################################################################################
testing_all_imp_exclusive<-testing_all_imp %>% mutate(LC_event3y=ifelse((LC_event==1 & tevent>3)|(LC_event==1 & tevent<=2), NA, LC_event),   
                                                      LC_event2y=ifelse((LC_event==1 & tevent>2)|(LC_event==1 & tevent<=1), NA, LC_event))
testing_1_imp_exclusive <- testing_all_imp_exclusive %>% filter(name==1)

table (testing_1_imp_exclusive$LC_event3y,exclude = NULL)
table (testing_1_imp_exclusive$LC_event2y,exclude = NULL)
###Good
# Models performance test ####
bootstraping_res_exclusive <- data.frame(NULL)
modelstest_vec <- c("proteinmodelrisk", "plcom2012")
for(imp in 1:15){
  print(imp)
  dataimp <- testing_all_imp_exclusive %>% filter(name==imp)
  for(t in 1:3){
    print(t)
    LC_event_t <- paste0("LC_event", t,"y")
    dataimp_t <- dataimp %>% filter(!is.na(!!sym(LC_event_t)))
    # Models ROC
    for(mod in modelstest_vec){ #ROC by model
      print(mod)
      mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
      
      boot_res_iter <- bootstrapped_AUC_cal(dataimp_t, mod_t, LC_event_t) %>% 
        mutate(model=mod, time=t, iter=imp)
      ### bind all bootstrap
      bootstraping_res_exclusive <- bind_rows(bootstraping_res_exclusive, boot_res_iter) # within each iteration of imputation: Bootstrap inference when using multiple imputation
    }
  }
}
bootstraping_res_exclusive %>% head
with(bootstraping_res_exclusive, table(model,time)) #15 imputed data set, 1000 bootstrapped iter, 3 time, 2model s

#random manual checks overall AUCs by leadtime and model 
bootstraping_res_exclusive %>% group_by(model,time, iter) %>% distinct(AUC_overall) %>% arrange(model, time)
bootstraping_res_exclusive %>% group_by(model,time) %>% summarise(AUC_mean=mean(AUC_overall), 
                                                                  lowerci=quantile(bootstrappedAUC_vec, probs=c(0.025)))

# MANUAL CALCULATION TO VERIFY
with(filter(testing_1_imp_exclusive, name==1, !is.na(LC_event2y)), 
     WeightedAUC(WeightedROC(Integral_risk_2y, LC_event2y, U19_weight))) 

# ROC and calibration plots
# Pool AUCs, CI and get pvalues
pooled_discriminatory_perf_main_exclusive <- summary_res_pooled_ROC(bootstraping_res_exclusive, c("plcom2012", "proteinmodelrisk"))

# ROCs FPR and TPR for 1 iteration(data) of the imputation 
ROC_mainplot_exclusive <- data.frame(NULL)
for(t in 1:3){
  LC_event_t <- paste0("LC_event", t,"y")
  dataimp_t <- testing_1_imp_exclusive %>% filter(!is.na(!!sym(LC_event_t)))
  for (mod in c("plcom2012", "proteinmodelrisk")){
    mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
    ROC_data_interim <- with(dataimp_t, WeightedROC(guess=get(mod_t), label=get(LC_event_t), 
                                                    weight=U19_weight)) %>% mutate(time=t, model=mod)
    ROC_mainplot_exclusive <- bind_rows(ROC_mainplot_exclusive, ROC_data_interim)
  }
}
dim(ROC_mainplot_exclusive) 
ROC_mainplot_exclusive <- ROC_mainplot_exclusive %>% merge(dplyr::select(pooled_discriminatory_perf_main_exclusive, time, model, AUC_CI_pooled,
                                                                         pvalue_ref), 
                                                           by=c("time","model"), all.x = T, all.y = F)
dim(ROC_mainplot_exclusive)

##### ROC ####
ROC_mainplot_exclusive <- ROC_mainplot_exclusive%>% 
  mutate(modelraw=case_when(model=="plcom2012"~"PLCOm2012", 
                            model=="proteinmodelrisk_"~"INTEGRAL Risk model", 
                            T~model),
         model=paste0(modelraw,"\nAUC=", AUC_CI_pooled))

head(ROC_mainplot_exclusive)
for(t in 1:3){
  title <- paste0("Discrimination at ", t,"y")
  data_f <- filter(ROC_mainplot_exclusive, time==t)
  pvalue <- formatC(unique(data_f$pvalue), digits=2, format="f")
  pvalue_annot <- ifelse(pvalue<0.01, "pdiff<0.01", paste0("pdiff=",pvalue))
  plot <- ROC_time_plot(data_f, title, "no",1.5, 17, c(.6, 0.35),17,17)+
    annotate("text", y=0.15, x=0.60, label=pvalue_annot, 
             fontface="italic", size=6, hjust=1, vjust=1)
  assign(paste0("ROC_main_", t,"y"), plot)
  print(plot)
}


# ###
#### Sensitivity, specificity, for 15 imputed datasets with confidence interval 
# ###
# ── Run across 15 imputed datasets & both USPSTF criteria ────────────────────
all_results_1y <- list()
for (imp_id in 1:15) {
  datax <- testing_all_imp_exclusive %>% filter(name==imp_id)  # replace with your indexing
  for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
    res <- bootstrapped_sensitivity(
      datax         = datax,uspstfy       = uspstfy,
      LC_event_time = "LC_event1y",t= 1)
    res$imp_id <- imp_id
    all_results_1y[[length(all_results_1y) + 1]] <- res
  }}
final_results_1y <- bind_rows(all_results_1y)%>% mutate(time=1)

all_results_2y <- list()
for (imp_id in 1:15) {
  datax <- testing_all_imp_exclusive%>% filter(name==imp_id)  # replace with your indexing
  for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
    res <- bootstrapped_sensitivity(
      datax         = datax,uspstfy       = uspstfy,
      LC_event_time = "LC_event2y",t= 2)
    res$imp_id <- imp_id
    all_results_2y[[length(all_results_2y) + 1]] <- res
  }}
final_results_2y <- bind_rows(all_results_2y)%>% mutate(time=2)

all_results_3y <- list()
for (imp_id in 1:15) {
  datax <- testing_all_imp_exclusive%>% filter(name==imp_id)   # replace with your indexing
  for (uspstfy in c("USPSTF2013", "USPSTF2021")) {
    res <- bootstrapped_sensitivity(
      datax  = datax,uspstfy       = uspstfy,
      LC_event_time = "LC_event3y",t= 3)
    res$imp_id <- imp_id
    all_results_3y[[length(all_results_3y) + 1]] <- res
  }}
final_results_3y <- bind_rows(all_results_3y) %>% mutate(time=3)

final_results_exclusive<-bind_rows(final_results_1y,final_results_2y,final_results_3y)

final_results_exclusive<-final_results_exclusive %>% mutate( NNS_model=(TP_model+FP_model)/TP_model,
                                                             NNS_uspstf=(TP_uspstf+FP_uspstf )/TP_uspstf)

# ── Summarise: mean + 95% CI across 15,000 bootstraps ────────────────────────
summary_results_exclusive <- final_results_exclusive %>%
  group_by(time,uspstf_criterion,model) %>%
  summarise(
    mean_sens  = mean(sens_model,       na.rm = TRUE),
    lower_sens = quantile(sens_model,   0.025, na.rm = TRUE),
    upper_sens = quantile(sens_model,   0.975, na.rm = TRUE),
    mean_TP=mean(TP_model,na.rm = TRUE ),
    mean_spec  = mean(spec_model,       na.rm = TRUE),

    mean_NNS  = mean(NNS_model,       na.rm = TRUE),
    lower_NNS = quantile(NNS_model,   0.025, na.rm = TRUE),
    upper_NNS = quantile(NNS_model,   0.975, na.rm = TRUE),
    
    # USPSTF reference
    mean_sens_uspstf = mean(sens_uspstf, na.rm = TRUE),
    lower_sens_uspstf = quantile(sens_uspstf,   0.025, na.rm = TRUE),
    upper_sens_uspstf = quantile(sens_uspstf,   0.975, na.rm = TRUE),
    mean_spec_uspstf = mean(spec_uspstf, na.rm = TRUE),
    mean_TP_uspstf=mean(TP_uspstf ,na.rm = TRUE ),
    mean_TN_uspstf=mean(TN_uspstf ,na.rm = TRUE ),
    
    mean_NNS_uspstf  = mean(NNS_uspstf,       na.rm = TRUE),
    lower_NNS_uspstf = quantile(NNS_uspstf,   0.025, na.rm = TRUE),
    upper_NNS_uspstf = quantile(NNS_uspstf,   0.975, na.rm = TRUE),
    
    n_boots    = n(),
    .groups    = "drop"
  )
summary_results_exclusive %>% head()


####Plot ROC curves
USPSTF_perf2013_exclusive <- data.frame(time=NULL)
for (t in 1:3) {
  perf <- perf_uspstf(testing_1_imp_exclusive, "USPSTF2013",paste0("LC_event",t,"y"))
  names(perf) <- paste0(names(perf),"_2013")
  perf <- perf %>% mutate(time=t)
  USPSTF_perf2013_exclusive <- bind_rows(USPSTF_perf2013_exclusive,perf)
}
USPSTF_perf2013_exclusive
USPSTF_perf2021_exclusive <- data.frame(time=NULL)
for (t in 1:3) {
  perf <- perf_uspstf(testing_1_imp_exclusive, "USPSTF2021",paste0("LC_event",t,"y"))
  names(perf) <- paste0(names(perf),"_2021")
  perf <- perf %>% mutate(time=t)
  USPSTF_perf2021_exclusive <- bind_rows(USPSTF_perf2021_exclusive,perf)
}
USPSTF_perf2021_exclusive
for(t in 1:3){
  title <- paste0("Discrimination at ", t,"y")
  data_f <- filter(ROC_mainplot_exclusive, time==t)
  data_f <- data_f %>% merge(USPSTF_perf2021_exclusive, by="time", all.x = TRUE, all.y = FALSE)
  data_f <- data_f %>% merge(USPSTF_perf2013_exclusive, by="time", all.x = TRUE, all.y = FALSE)
  pvalue <- formatC(unique(data_f$pvalue), digits=2, format="f")
  pvalue_annot <- ifelse(pvalue<0.01, "pdiff<0.01", paste0("pdiff=",pvalue))
  plot <- ROC_time_plot_update(data_f, title, 1.5, 15, c(.7, 0.35),15,15)+
    annotate("text", y=0.15, x=0.60, label=pvalue_annot, 
             fontface="italic", size=6, hjust=1, vjust=1)
  assign(paste0("ROC_main_exclusive_", t,"y"), plot)
  print(plot)
}

#########Sensitivity plot
summary_results_plot_exclusive<-summary_results_exclusive %>% 
  mutate(`Sensitity, % (95CI)`=paste0(sprintf("%.1f", mean_sens*100)," (",sprintf("%.1f", lower_sens*100),"-",sprintf("%.1f", upper_sens*100),")"),
         `NNS, % (95CI)`=paste0(sprintf("%.0f", mean_NNS)," (",sprintf("%.0f", lower_NNS),"-",sprintf("%.0f", upper_NNS),")"),
         `Sensitity USPSTF, % (95CI)`=paste0(sprintf("%.1f", mean_sens_uspstf*100)," (",sprintf("%.1f", lower_sens_uspstf*100),"-",sprintf("%.1f", upper_sens_uspstf*100),")"),
         `NNS USPSTF, % (95CI)`=paste0(sprintf("%.0f", mean_NNS_uspstf)," (",sprintf("%.0f", lower_NNS_uspstf),"-",sprintf("%.0f", upper_NNS_uspstf),")"),
         mean_TP_uspstf=round(mean_TP_uspstf,0),
         mean_TP=round(mean_TP,0)) %>% 
  dplyr::select(- c("mean_sens","lower_sens", "upper_sens","mean_spec","n_boots","mean_sens_uspstf","lower_sens_uspstf","upper_sens_uspstf","mean_TN_uspstf",
                    "mean_NNS_uspstf","lower_NNS_uspstf","upper_NNS_uspstf", "mean_NNS","lower_NNS","upper_NNS")) 

data_figure3def_exclusive <- data.frame(NULL)
for(t in 1:3){
  data_figure3def_sub<-data.frame(
    Metric = c("USPSTF-2013/2021", "PLCOm2012 model","INTEGRAL Risk model"),
    time=t,
    us13_TP   = c(summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="plcom2012",] $mean_TP_uspstf,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="plcom2012",] $mean_TP,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="Integral_risk",] $mean_TP),
    us13_SE   = c(summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="plcom2012",] $`Sensitity USPSTF, % (95CI)`,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="plcom2012",] $`Sensitity, % (95CI)`,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="Integral_risk",] $`Sensitity, % (95CI)`),
    us13_NNS   = c(summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="plcom2012",] $`NNS USPSTF, % (95CI)`,
                   summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="plcom2012",] $`NNS, % (95CI)`,
                   summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2013"&summary_results_plot_exclusive$model=="Integral_risk",] $`NNS, % (95CI)`),
    
    us21_TP   = c(summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="plcom2012",] $mean_TP_uspstf,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="plcom2012",] $mean_TP,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="Integral_risk",] $mean_TP),
    us21_SE   = c(summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="plcom2012",] $`Sensitity USPSTF, % (95CI)`,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="plcom2012",] $`Sensitity, % (95CI)`,
                  summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="Integral_risk",] $`Sensitity, % (95CI)`),
    
    us21_NNS   = c(summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="plcom2012",] $`NNS USPSTF, % (95CI)`,
                   summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="plcom2012",] $`NNS, % (95CI)`,
                   summary_results_plot_exclusive[summary_results_plot_exclusive$time==t&summary_results_plot_exclusive$uspstf_criterion=="USPSTF2021"&summary_results_plot_exclusive$model=="Integral_risk",] $`NNS, % (95CI)`))
  
  data_figure3def_exclusive<-bind_rows(data_figure3def_exclusive,data_figure3def_sub)
}

data_figure3def_long_exclusive <- data_figure3def_exclusive %>%
  mutate(across(c(us13_TP, us21_TP), as.character)) %>%  # make all columns character
  pivot_longer(
    cols = c(us13_TP, us13_SE,us13_NNS, us21_TP, us21_SE,us21_NNS),
    names_to = c("group", "metric"),
    names_sep = "_"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = value
  ) %>%
  mutate(metric = case_when(metric == "TP"~ "TP, N", 
                            metric =="SE"~"Sensitivity, %",
                            metric =="NNS"~"quasi-NNS"))
data_figure3def_long_exclusive$time<-as.numeric(data_figure3def_long_exclusive$time)

USPSTF_perf2021_exclusive
USPSTF_perf2013_exclusive
make_table2 <- function(data, group_name) {
  data<-data %>% dplyr::select(-time)%>%
    group_by(Metric) %>%
    mutate(Label = if_else(
      row_number() == 1,
      paste0(Metric, "\n  ", metric),
      paste0("  ", metric))) %>%
    ungroup()%>%
    dplyr::select(
      Label,us13,us21)
  header_bottom <- c("Model", "USPSTF-2013\n Specificity at 87.2%", "USPSTF-2021\n Specificity at 75.1%")
  mat <- rbind(header_bottom, as.matrix(data))
  tg <- tableGrob(mat, rows = NULL, cols = NULL,
                  theme = ttheme_minimal(
                    core    = list(fg_params = list(hjust = 0.5, x = 0.5, fontsize = 14),
                                   bg_params = list(fill = NA, col = NA))))
  # helper: add a horizontal line at the BOTTOM of row r
  hline <- function(tg, r, lwd = 0.8) {
    gtable_add_grob(tg,segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"),
      y0 = unit(0, "npc"), y1 = unit(0, "npc"),
      gp = gpar(lwd = lwd)),
      t = r, b = r, l = 1, r = ncol(tg), z = Inf)}
  n_rows <- nrow(mat)
  # top line
  tg <- hline(tg, r = 1, lwd = 1.5)
  # between-group lines (every 2 data rows, skip within-group)
  group_rows <- seq(4, n_rows - 2, by = 3)  # bottom of first row of each group
  for (r in group_rows) tg <- hline(tg, r = r, lwd = 0.8)
  # bottom line
  tg <- hline(tg, r = n_rows, lwd = 1.5)
  as.ggplot(tg)
}


table_plot_list <- list()
for(t in 1:3){
  title <- paste0("Sensitivity at ", t,"y")
  data_f <- filter(data_figure3def_long_exclusive , time==t)
  plot <- make_table2(data_f, title)
  table_plot_list[[t]] <- plot
  assign(paste0("SE_main_", t,"y"), table_plot_list[[t]] )
  print(table_plot_list)
}

cowplot::plot_grid(ROC_main_1y+ggtitle("A) Discrimination over 0-1 year")+theme(plot.title = element_text(size = 17)),
                   ROC_main_2y+ggtitle("B) Discrimination over 1-2 years")+theme(plot.title = element_text(size = 17)),
                   ROC_main_3y+ggtitle("C) Discrimination over 2-3 years")+theme(plot.title = element_text(size = 17)),
                   
                   SE_main_1y+ggtitle("")+theme(plot.title = element_text(size = 17)),
                   SE_main_2y+ggtitle("")+theme(plot.title = element_text(size = 17)),
                   SE_main_3y+ggtitle("")+theme(plot.title = element_text(size = 17)),
                   ncol=3, nrow=2, align = "hv")

ggsave(filename="Supple_Figure3_modelperformance_exclusive.pdf", path = paste0(path_to_save), 
       height = 15, width = 25, units="in", device = cairo_pdf, scale=0.8)



# #############################################################
########eFigure 4. Calibration of the INTEGRAL-Risk model in 
##the LC3 testing set for predicting lung cancer occurring within 
##1-, 2-, and 3-years of follow-up
###############################################################

#########Hosmer-Lemehow calibration test######
library(ResourceSelection)  # for hoslem.test

###Use ~10 categories that divides the cases more equally
hl_results <- data.frame()
hl_plot_data <- data.frame()

for (t in 1:3) {
  LC_event_time <- paste0("LC_event", t, "y")
  risk_predtime <- paste0("proteinmodelrisk_", t, "y")
  
  dat_t <- testing_1_imp %>%
    mutate(U19_weight = ifelse(is.na(.data[[LC_event_time]]), 0, U19_weight)) %>%
    ungroup() %>%
    filter(!is.na(.data[[LC_event_time]]),
           !is.na(.data[[risk_predtime]]))
  
  obs  <- dat_t[[LC_event_time]]
  pred <- dat_t[[risk_predtime]]
  wt   <- dat_t$U19_weight
  pred_case <- pred[obs == 1]
  g <- 10
  
  breaks <- quantile(pred_case, probs = seq(0, 1, length.out = g + 1),na.rm = TRUE)
  breaks[1]   <- -Inf
  breaks[g+1] <- Inf
  
  ## avoid duplicated cut points
  breaks <- unique(breaks)
  group <- cut(pred, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  dat_t <- dat_t %>% mutate(group = cut(.data[[risk_predtime]],
                                        breaks = breaks, labels = FALSE, include.lowest = TRUE))
  
  hl_grouped <- dat_t %>% group_by(group) %>%
    summarise(
      O = sum(.data[[LC_event_time]]),
      E = sum(.data[[risk_predtime]]*U19_weight),
      N = sum(U19_weight),               # weighted group size
      OE_factor = ifelse(E > 0, O / E, NA_real_),  
      obs_risk = O / N,
      lowerCI=obs_risk-1.96*sqrt((obs_risk*(1-obs_risk))/N ),
      upperCI=obs_risk+1.96*sqrt((obs_risk*(1-obs_risk))/N ),
      exp_risk = E / N,
      lower = breaks[first(group)],
      upper = breaks[first(group) + 1],
      .groups = "drop") %>%
    mutate( threshold = paste0( ifelse(is.infinite(lower), "Min", sprintf("%.4f", lower)),
                                " to ", ifelse(is.infinite(upper), "Max", sprintf("%.4f", upper))),
            time = t )
  
  #### Current calibration test
  hl_stat <- with(hl_grouped, sum((O - E)^2 / (E * (1 - E / N))))
  df_hl   <- nrow(hl_grouped) - 2
  p_val   <- pchisq(hl_stat, df = df_hl, lower.tail = FALSE)
  
  hl_results <- bind_rows(
    hl_results,
    data.frame(
      time = t,
      HL_stat = round(hl_stat, 3),
      df = df_hl,
      p_value = ifelse(p_val<0.01, "<0.01", paste0(round(p_val, 4))),
      n_groups = nrow(hl_grouped)))
  
  hl_plot_data  <- bind_rows(hl_plot_data , hl_grouped)
}

hl_plot_data %>%
  select(time, group, threshold) 


plot_list <- vector("list", length = 3)
for (i in 1:3) {
  time_i <-i
  dat_i <- hl_plot_data %>% filter(time == time_i)
  res_i <- hl_results  %>% filter(time == time_i)
  max_y <- max(dat_i$upperCI)+0.001
  max_x <- max(dat_i$exp_risk)+0.001
  max_eq <- max(max_y, max_x)
  
  plot_list[[i]] <- ggplot(dat_i, aes(x = exp_risk, y = obs_risk)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8) +
    geom_point(size = 3) +
    geom_pointrange(aes(ymin=lowerCI, ymax=upperCI))+
    geom_line(aes(group = 1), linewidth = 0.8) +
    annotate("text",
             x =  0.0002,
             y =  0.1 ,
             label = paste0(
               "HL = ", res_i$HL_stat,
               "\ndf = ", res_i$df,
               "\np  ", res_i$p_value),
             hjust = 0,size = 4.5) +
    labs( title = paste0(LETTERS[i], ") ", time_i, "-year calibration"),
          x = "Expected risk",
          y = "Observed risk") +
    theme_bw() +
    scale_x_log10(breaks = c(0.0001,0.001,0.005,0.01,0.1),
                  labels = c("0.01%", "0.1%", "0.5%","1.0%",  "10.0%"))+
    scale_y_log10(breaks = c(0.0001,0.001,0.005,0.01,0.1),
                  labels = c("0.01%", "0.1%", "0.5%","1.0%","10.0%"))+
    coord_equal(                                       
      xlim = c(0.0001,  0.35),
      ylim = c(0.0001,  0.35)) +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

cowplot::plot_grid(
  plotlist = plot_list,
  ncol = 3)

ggsave(filename="Supple_Figure4_calibration_10groups.pdf", path = paste0(path_to_save), 
       height = 13, width = 20, units="in", device = cairo_pdf, scale=0.8)

###Get the values for group
data_eo_10group <- data.frame(NULL)
for(t in 1:3){
  LC_event_time <- paste0("LC_event",t,"y")
  risk_predtime <- paste0("Integral_risk","_",t,"y")
  for(i in 1:1000){
    set.seed(123*i)
    ## bootstrapping data: resampling data similar to original pop/ full data wont always have same IR, that's 
    #why we're allowing change in the IR , stratified resampling only by cohort
    bootstrapped_data_m3 <- testing_all_imp %>% group_by(cohort) %>%
      sample_n(size=n(), replace = T) %>% ungroup()
    
    obs  <- bootstrapped_data_m3[[LC_event_time]]
    pred <- bootstrapped_data_m3[[risk_predtime]]
    pred_case <- pred[obs == 1]
    g <- 10
    
    breaks <- quantile(pred_case, probs = seq(0, 1, length.out = g + 1),na.rm = TRUE)
    breaks[1]   <- -Inf
    breaks[g+1] <- Inf
    
    ## avoid duplicated cut points
    breaks <- unique(breaks)
    group <- cut(pred, breaks = breaks, labels = FALSE, include.lowest = TRUE)
    
    bootstrapped_data_m3 <- bootstrapped_data_m3 %>% mutate(group = cut(.data[[risk_predtime]],
                                                                        breaks = breaks, labels = FALSE, include.lowest = TRUE))
    
    data_plot_t <- bootstrapped_data_m3 %>% 
      group_by(name) %>% 
      mutate(U19_weight=ifelse(is.na(!!sym(LC_event_time)), 0, U19_weight)) %>% ungroup()
    
    data_plot_t <- data_plot_t %>% group_by(name,group) %>% 
      summarise(n_pop=sum(U19_weight, na.rm = TRUE),
                predicted_event=sum(get(risk_predtime)*U19_weight),
                observed_event=sum(get(LC_event_time), na.rm = TRUE),
                observed_risk=observed_event/n_pop,
                predicted_risk=predicted_event/n_pop,
                .groups         = "drop")%>% 
      mutate(E_O= predicted_risk/observed_risk, 
             time=t,boot=i) 
    data_eo_10group <- bind_rows(data_eo_10group, data_plot_t)
  }}


data_eo_10group<-data_eo_10group%>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.), NA, .)))


summary_EO_GMs_10sub <- data_eo_10group %>%  group_by(time,group) %>%
  summarise(n_population=round(mean(as.numeric(n_pop)),0),
            mean_Exp=mean(as.numeric(predicted_event)), lowerCI_pooled_Exp=quantile(predicted_event, prob=c(0.025))[[1]],upperCI_pooled_Exp=quantile(predicted_event, prob=c(0.975))[[1]],
            mean_Obs=mean(as.numeric(observed_event)),lowerCI_pooled_Obs=quantile(observed_event, prob=c(0.025))[[1]],upperCI_pooled_Obs=quantile(observed_event, prob=c(0.975))[[1]],
            mean_E_O=mean(as.numeric(E_O), na.rm = TRUE),
            lowerCI_pooled_E_O=quantile(E_O, prob=c(0.025), na.rm = TRUE)[[1]],upperCI_pooled_E_O=quantile(E_O, prob=c(0.975), na.rm = TRUE)[[1]],
            mean_Exprisk=100*mean(as.numeric(predicted_risk)),lowerCI_pooled_Exprisk=100*quantile(predicted_risk, prob=c(0.025))[[1]],upperCI_pooled_Exprisk=100*quantile(predicted_risk, prob=c(0.975))[[1]],
            mean_Obsrisk=100*mean(as.numeric(observed_risk)), lowerCI_pooled_Obsrisk=100*quantile(observed_risk, prob=c(0.025))[[1]],upperCI_pooled_Obsrisk=100*quantile(observed_risk, prob=c(0.975))[[1]]) %>% 
  mutate(EO_CI=paste0(formatC(mean_E_O, digits=2, format = "f")," (", 
                      formatC(lowerCI_pooled_E_O, digits = 2, format = "f"), "-",
                      formatC(upperCI_pooled_E_O, digits=2, format="f"), ")"),
         Exp_CI=paste0(formatC(mean_Exp, digits=0, format = "f")," (", 
                       formatC(lowerCI_pooled_Exp, digits = 0, format = "f"), "-",
                       formatC(upperCI_pooled_Exp, digits=0, format="f"), ")"),
         Obs_CI=paste0(formatC(mean_Obs, digits=0, format = "f")," (", 
                       formatC(lowerCI_pooled_Obs, digits = 0, format = "f"), "-",
                       formatC(upperCI_pooled_Obs, digits=0, format="f"), ")"),
         Exprisk_CI=paste0(formatC(mean_Exprisk, digits=2, format = "f")," (", 
                           formatC(lowerCI_pooled_Exprisk, digits = 2, format = "f"), "-",
                           formatC(upperCI_pooled_Exprisk, digits=2, format="f"), ")"),
         Obsrisk_CI=paste0(formatC(mean_Obsrisk, digits=2, format = "f")," (", 
                           formatC(lowerCI_pooled_Obsrisk, digits = 2, format = "f"), "-",
                           formatC(upperCI_pooled_Obsrisk, digits=2, format="f"), ")"),) %>% ungroup() %>% 
  dplyr::select(time,group,n_population,Exprisk_CI,Obsrisk_CI,EO_CI,Exp_CI,Obs_CI)
summary_EO_GMs_10sub

write.csv(summary_EO_GMs_10sub,paste0(path_to_save,"summary_EO_GMs_10sub.csv"),row.names = FALSE)

###Overall 
summary_EO_GMs_10sub_overall<-bootstraping_res %>%  group_by(model, time) %>%summarise(mean_EO=mean(as.numeric(cal_overall)),
                                                                                       lowerCI_pooled=quantile(bootstrappedcal_vec, prob=c(0.025))[[1]],
                                                                                       upperCI_pooled=quantile(bootstrappedcal_vec, prob=c(0.975))[[1]],
                                                                                       mean_N=mean(as.numeric(totpop_boot )),
                                                                                       mean_O=mean(ncases_boot),
                                                                                       lowerCI_pooled_O=quantile(ncases_boot, prob=c(0.025))[[1]],
                                                                                       upperCI_pooled_O=quantile(ncases_boot, prob=c(0.975))[[1]],
                                                                                       mean_E=mean(ncases_boot*bootstrappedcal_vec),
                                                                                       lowerCI_pooled_E=quantile(ncases_boot*bootstrappedcal_vec, prob=c(0.025))[[1]],
                                                                                       upperCI_pooled_E=quantile(ncases_boot*bootstrappedcal_vec, prob=c(0.975))[[1]]) %>% 
  mutate(cal_CI=paste0(formatC(mean_EO, digits=2, format = "f")," (", 
                       formatC(lowerCI_pooled, digits = 2, format = "f"), "-",
                       formatC(upperCI_pooled, digits=2, format="f"), ")"),
         O_CI=paste0(formatC(mean_O, digits=2, format = "f")," (", 
                     formatC(lowerCI_pooled_O, digits = 2, format = "f"), "-",
                     formatC(upperCI_pooled_O, digits=2, format="f"), ")"),
         E_CI=paste0(formatC(mean_E, digits=2, format = "f")," (", 
                     formatC(lowerCI_pooled_E, digits = 2, format = "f"), "-",
                     formatC(upperCI_pooled_E, digits=2, format="f"), ")"),
         N=mean_N) %>% ungroup()
write.csv(summary_EO_GMs_10sub_overall,paste0(path_to_save,"summary_EO_GMs_10sub_overall.csv"),row.names = FALSE)


# ####################################################################################################
########eFigure 5. Clinical impact of using the INTEGRAL-Risk model to assess screening eligibility, 
##compared to the USPSTF 2021 criteria and the PLCOm2012 model. The Lorenz curve depicts the fraction 
##of 1-year incident lung cancer cases identified as screening eligible as a function of the size of 
##the overall eligible population (expressed as the fraction of eligible participants with a 
##smoking history)
######################################################################################################
lorenz_1y <- lorenz_plot_function(testing_1_imp, "Integral_risk_1y", "LC_event1y",1,4.5)+theme(legend.position = "none", axis.text  = element_text(size=13),
                                                                                               plot.margin = margin(0.5,1,0.5,0.5, "cm"))

ggsave(filename="Supple_Figure5_lorenz_1y.pdf",
       path = paste0(path_to_save), 
       height = 9, width = 9, units="in",scale=0.88, device = cairo_pdf)

# ####################################################################################################
########eFigure 6. The Venn-diagram depicts the characteristics for research participants identified 
##by the USPSTF2021 criteria, the PLCOm2012 and INTEGRAL-Risk models. Risk thresholds for the PLCOm2012 
##and INTEGRAL-Risk model were set to select the same baseline population as the USPSTF2021 screening 
##criteria
######################################################################################################
threshold_plco_USPSTF <- thresholds_1y[[4]]
threshold_protmodel1y_USPSTF <- thresholds_1y[[3]]
# include training and testing set since we're doing mutually exclusive groups and the numbers are not that big 
data_train <- data_train %>% mutate(USPSTF2021=ifelse(age>=50& age<=80&packyears>=20 &quit_years<=15,1,0))
data_train <- calculate_modelrisk_time(data_train, c("model_fit_AIC"), c("proteinmodelrisk")) 

data_venn_description <- bind_rows(dplyr::select(testing_1_imp, tevent, cohort, LC_event,LC_event1y,LC_event2y,LC_event3y, U19_weight, age, years_smoked,quit_years,
                                                 smoke_status,cohort, plcom2012, USPSTF2021, proteinmodelrisk_1y,proteinmodelrisk_2y, proteinmodelrisk_3y),
                                   dplyr::select(data_train,tevent, cohort, LC_event, LC_event1y,LC_event2y,LC_event3y,U19_weight, age, years_smoked,quit_years,
                                                 smoke_status, plcom2012, USPSTF2021, proteinmodelrisk_1y,proteinmodelrisk_2y, proteinmodelrisk_3y)) %>%
  mutate(quit_years=ifelse(smoke_status==3, NA, quit_years),
         elig_USPSTF2021=ifelse(USPSTF2021==1,"EligUSPSTF", "Non elig"),
         PLCO_08=ifelse(plcom2012>=threshold_plco_USPSTF, "eligPLCO_08", "Non elig"),
         PLCO_1=ifelse(plcom2012>=0.0151, "eligPLCO_1", "Non elig"),
         elig_prot=ifelse(proteinmodelrisk_1y>=threshold_protmodel1y_USPSTF, "Elig_prot1y", "Non elig"),
         elig_comp_3g=case_when(PLCO_08=="eligPLCO_08" & elig_prot=="Non elig" &elig_USPSTF2021=="Non elig" ~"PLCO_08only",
                                PLCO_08=="Non elig" & elig_prot=="Elig_prot1y" &elig_USPSTF2021=="Non elig" ~"INTEGRAL_only",
                                PLCO_08=="Non elig" & elig_prot=="Non elig" &elig_USPSTF2021=="EligUSPSTF" ~"USPSTF_only",
                                PLCO_08=="eligPLCO_08" & elig_prot=="Elig_prot1y" &elig_USPSTF2021=="EligUSPSTF" ~"elig_all",
                                PLCO_08=="eligPLCO_08" & elig_prot=="Non elig" &elig_USPSTF2021=="EligUSPSTF" ~"PLCO_USPSTF",
                                PLCO_08=="Non elig" & elig_prot=="Elig_prot1y" &elig_USPSTF2021=="EligUSPSTF" ~"INTEGRAL_USPSTF",
                                PLCO_08=="eligPLCO_08" & elig_prot=="Elig_prot1y" &elig_USPSTF2021=="Non elig" ~"INTEGRAL_plco",
                                PLCO_08=="Non elig" & elig_prot=="Non elig" &elig_USPSTF2021=="Non elig" ~"Not eligible"))

data_venn_y <- function(datax, LC_event_time, y){
  data_venn_3g_y <- datax %>% 
    filter(!is.na(!!sym(LC_event_time))) %>% 
    mutate(weight_baseline=ifelse(LC_event==0, U19_weight, 0),
           total_pop=sum(weight_baseline),
           total_cases=sum(!!sym(LC_event_time))) %>% 
    group_by(elig_comp_3g) %>% 
    summarise(total_pop=unique(total_pop) ,
              total_cases=unique(total_cases),
              N_screened=sum(weight_baseline), 
              N_obs=n(),
              per_screened= 100*N_screened/total_pop,
              mean_age=formatC(Hmisc::wtd.mean(age, weight_baseline ,na.rm = TRUE), digits=0, format="f"),
              age=wtd_distributions(age, weight_baseline), 
              mean_quity=formatC(Hmisc::wtd.mean(quit_years, weight_baseline ,na.rm = TRUE), digits=0, format="f"),
              quit_years=wtd_distributions(quit_years, weight_baseline), 
              per_former=100*sum((smoke_status==2)*weight_baseline)/N_screened,
              N_cases=sum(!!sym(LC_event_time)), 
              per_cases_screened=100*N_cases/total_cases,
              cumulative_risk=100*N_cases/N_screened) %>% 
    mutate(annot_old=paste0("N=",ceiling(N_screened), " people (", formatC(per_screened, format = "f", digits=0),"%)\n",
                            mean_age,"y | ",mean_quity,"y cessation\n",
                            y,"y risk=",formatC(cumulative_risk, format="f", digits=2), "t %"),
           annot=paste0("N=",ceiling(N_screened), " people (", formatC(per_screened, format = "f", digits=0),"%)\n",
                        mean_age,"y | ",mean_quity,"y cessation\n",
                        "Cases ",y, "y=",N_cases," (", formatC(per_cases_screened, format="f", digits=0), "%)"))
  
  annot_venn_y_3g <- c("INTEGRAL_USPSTF"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="INTEGRAL_USPSTF"][1],
                       "INTEGRAL_only"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="INTEGRAL_only"][1],
                       "INTEGRAL_plco"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="INTEGRAL_plco"][1],
                       "PLCO_08only"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="PLCO_08only"][1],
                       "PLCO_USPSTF"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="PLCO_USPSTF"][1],
                       "USPSTF_only"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="USPSTF_only"][1],
                       "elig_all"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="elig_all"][1],
                       "Not eligible"=data_venn_3g_y$annot[data_venn_3g_y$elig_comp_3g=="Not eligible"][1])
  
  
  return(list(data_venn_3g_y,annot_venn_y_3g))
}

plot_venn <- function(annot_y){
  plot <- data.frame(tool=c("USPSTF2021", "PLCO", "INTEGRAL"),
                     x0=c(2.5,6.5,4.5),#c(2.5,6.5,4.5),
                     y0=c(2,2,5.5)) %>% 
    ggplot(aes(fill=tool))+
    geom_circle(aes(x0 = x0, y0 = y0, r=4), alpha=0.5)+
    geom_circle(aes(x0=9.5,y0=10, r=2, fill="Not eligible"))+# a = 5, b = 5, angle = 0 / for ellipse
    annotate("text",
             y=c(4.5, #Integral-uspstf
                 7.5, #Integral oinly
                 4.5, #INTEGRAL_plco
                 1.5, #PLCO_08only
                 0.2, #PLCO_USPSTF
                 1.5, #USPSTF_only
                 3,#elig all
                 10), #not eligble
             x=c(1.2,#Integral-uspstg
                 3.5, #integral only
                 6.1, #INTEGRAL_plco
                 7.5, #PLCO_08only
                 3.5, #PLCO_USPSTF
                 -0.5, #USPSTF_only
                 3.5,# elig all
                 8.5), #not eligible
             label=annot_y,
             #fontface="italic",
             size=5, hjust=0, vjust=0.5)+
    annotate("text", y=c(-2.5,-2.5,10, 12.5), x=c(1,7.2, 3,9), 
             label=c("USPSTF2021", "PLCOm2012", "INTEGRAL Risk Model", "Not eligible"), 
             fontface="bold", size=6, hjust=0, colour=c("#e18727ff", "#0072b5ff", "#bc3c29ff","gray30"))+
    scale_fill_manual(values=c("USPSTF2021" = "#f3cfa8", "PLCO" = "#99c6e1", "INTEGRAL" = "#e4b1a9","Not eligible"= "gray90"))+
    theme_void()+ 
    theme(legend.position = "none")
  return(plot)
}

## Venn 1y
annot_venn_1y_3g <- data_venn_y(data_venn_description,"LC_event1y", "1")[[2]]
Venn_1y_3g <- plot_venn(annot_venn_1y_3g)

## Venn 2y
annot_venn_2y_3g <- data_venn_y(data_venn_description,"LC_event2y", "2")[[2]]
Venn_2y_3g <- plot_venn(annot_venn_2y_3g)

## Venn 3y
annot_venn_3y_3g <- data_venn_y(data_venn_description,"LC_event3y", "3")[[2]]
Venn_3y_3g <- plot_venn(annot_venn_3y_3g)


ggarrange(Venn_1y_3g+ggtitle("A) Characteristics of individuals and cases within 1y selected")+
            theme(plot.title = element_text(size = 20, face = "bold")),
          Venn_2y_3g+ggtitle("B) Characteristics of individuals and cases within 2y selected")+
            theme(plot.title = element_text(size = 20, face = "bold")),
          Venn_3y_3g+ggtitle("C) Characteristics of individuals and cases within 3y selected")+
            theme(plot.title = element_text(size = 20, face = "bold")), ncol=1)

ggsave(filename="Supple_Figure6_Venn_diagram.pdf", path = paste0(path_save, "supp"), 
       height = 19, width = 16, scale=0.9 , units="in", device = cairo_pdf)




