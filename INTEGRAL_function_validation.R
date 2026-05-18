# ###############################################
# Functions necessary to assist the INTEGRAL_validationscript.R script from the Zahed et al. JAMA 2026
# Paper title: Biomarker-based eligibility for lung cancer screening
#              Validation of the protein-based INTEGRAL-Risk model 
# Authors: Mainly: Hana Zahed; Secondary: Xiaoshuang Feng
# Senior authors: Hilary Robbins, Mattias Johansson
# Contact email: fengx@iarc.who.int; robbinsh@iarc.who.int; johanssonm@iarc.who.int
# #################################################

## Function to calculate absolute risk values in the validation dataset
calculate_modelrisk_time <- function(datax, modelobjectname, column_name){
  for(mod in 1:length(modelobjectname)){ #or maybe create a vector with names here might 
    mod_name <- modelobjectname[mod]
    col_name <- column_name[mod]
    for(t in c(0.5, 1:3)){
      LC_event_y <- paste0("LC_event",t,"y")
      col_name_t <- paste0(col_name,"_", t,"y")
      #print(col_name_t)
      datax <- datax %>% 
        mutate(!!sym(col_name_t):=1-predict(get(mod_name), newdata=mutate(datax, t0=0, tevent=t),type="surv"),
               # create LC event by lead time
               !!sym(LC_event_y):=ifelse(LC_event==1 & tevent>t, NA, LC_event))  
    }
  }
  return(datax)
}

# ###
bootstrapped_AUC_cal <- function(datax, risk_predtime, LC_event_time){
  bootstrapAUC <- c(NULL)
  bootstrappedcal <- c(NULL)
  cases_observed_boot <- c(NULL)
  tot_pop_boot<- c(NULL)
  
  datax <- datax %>% filter(!is.na(!!sym(LC_event_time)))
  totpop_overall <- sum(datax$U19_weight, na.rm = TRUE)
  # AUC
  AUC_overall <- with(datax, 
                      WeightedROC::WeightedAUC(WeightedROC(guess=get(risk_predtime), label=get(LC_event_time), 
                                                           weight=U19_weight)))
  #AUC_overall <- formatC(AUC_overall, digits = 2, format = "f")
  # Cal 
  expected_test <- sum(datax[,risk_predtime, drop=T]*datax$U19_weight, na.rm = TRUE)
  observed_test <- sum(datax[,LC_event_time, drop=T], na.rm = TRUE)
  Calibration_overall <-ifelse(observed_test!=0, expected_test/observed_test, expected_test)
  
  #Calibration_overall <- formatC(Calibration_overall, digits=2, format="f")
  iter_seq <- c(1:1000)
  for(i in iter_seq){
    set.seed(123*i)
    ## bootstrapping data: resampling data similar to original pop/ full data wont always have same IR, that's 
    #why we're allowing change in the IR , stratified resampling only by cohort
    bootstrapped_data_m3 <- datax %>% group_by(cohort) %>%
      sample_n(size=n(), replace = T) %>% ungroup()
    tot_pop_i <- sum(bootstrapped_data_m3$U19_weight, na.rm=TRUE)
    tot_pop_boot <- c(tot_pop_boot, tot_pop_i)
  
    AUC_i <- tryCatch(with(bootstrapped_data_m3, 
                           WeightedROC::WeightedAUC(WeightedROC(guess=get(risk_predtime), label=get(LC_event_time), 
                                                                weight=U19_weight))),
                      error=function(e){return(NA)})
    
    bootstrapAUC <- c(bootstrapAUC,AUC_i )
    #cal
    expected_testi <- sum(bootstrapped_data_m3[,risk_predtime, drop=T]*bootstrapped_data_m3$U19_weight, na.rm = TRUE)
    observed_testi <- sum(bootstrapped_data_m3[,LC_event_time, drop=T], na.rm=TRUE)
    Calibration_test <-ifelse(observed_testi!=0,expected_testi/observed_testi, NA)
    bootstrappedcal <- c(bootstrappedcal,Calibration_test)
    cases_observed_boot <- c(cases_observed_boot,observed_testi)
    
  }
  summary_res <- data.frame(iter=iter_seq,
                            bootstrappedAUC_vec=bootstrapAUC, bootstrappedcal_vec=bootstrappedcal,
                            ncases_boot=cases_observed_boot,totpop_boot=tot_pop_boot) %>% 
    mutate(cal_overall=Calibration_overall, AUC_overall=AUC_overall, ncase_overall=observed_test, pop_overall=totpop_overall)
  return(summary_res)
} #IN ONE ITERATION 


# Function to pull all vectors of bootstrap across 15 datasets
run_bootstrap_t <- function(datax_1imp,modelstest_vec){
  bootstraping_res <- data.frame(NULL)
  for(t in 1:3){
    print(t)
    LC_event_t <- paste0("LC_event", t,"y")
    dataimp_t <- datax_1imp %>% filter(!is.na(!!sym(LC_event_t)))
    # Models ROC
    for(mod in modelstest_vec){ #ROC BY MOD
      #print(mod)
      mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
      boot_res_iter <- bootstrapped_AUC_cal(dataimp_t, mod_t, LC_event_t) %>% 
        mutate(model=mod, time=t)
      ### bind all bootstrap
      bootstraping_res <- bind_rows(bootstraping_res, boot_res_iter) # within each iteration of imputation: Bootstrap inference when using multiple imputation
    }
  }
  return(bootstraping_res)
} 

##Function to calculate P value for the AUC difference
pvalue_bootstrapAUC <- function(bootstrapAUCvec1, bootstrapAUCvec2, maindiff_integer){
  # Get pvalues for differences ###
  # 1. Get distributions of differences 
  diff_AUC_vec <- bootstrapAUCvec1-bootstrapAUCvec2
  
  #2. Center the differences of bootstrapped distribution
  diff_AUC_veccenter <- diff_AUC_vec-mean(diff_AUC_vec)
  # 3. Compare if the main difference comes from the distribution of the centered disribution and calculate pvalues 
  # get pvalues:
  #print(pnorm(maindiff_integer, mean=mean(diff_AUC_veccenter), sd=sd(diff_AUC_veccenter), lower.tail = F))
  pvalue_diff <- 2*pnorm(maindiff_integer, mean=mean(diff_AUC_veccenter), sd=sd(diff_AUC_veccenter), lower.tail = F)
  
  ###
  N_vector <- length(diff_AUC_vec)
  
  pvalue_diff2 <- length(which(-abs(maindiff_integer)>diff_AUC_veccenter))/N_vector+
    length(which(abs(maindiff_integer)<diff_AUC_veccenter))/N_vector
  return(c(pvalue_diff,pvalue_diff2))
}

##Function to summarize AUC, 95%CI and p value for the difference between the studied model compared to PLCOm2012
summary_res_pooled_ROC <- function(bootstrapdatax,modelstest_vec){
  dataframe_summmaryROC <- data.frame(NULL)
  for(t in 1:3){
    bootstraping_res_pooled_t <- bootstrapdatax %>% filter(time==t, model %in% modelstest_vec)
    AUC_models_pooled_t <- bootstrapdatax %>% filter(time==t, model %in% modelstest_vec) %>% 
      group_by(model) %>% summarise(AUCpooled=mean(AUC_overall), 
                                    lowerCIpooled=quantile(bootstrappedAUC_vec, probs=0.025, na.rm=TRUE), 
                                    upperCIpooled=quantile(bootstrappedAUC_vec, probs=0.975, na.rm=TRUE)) %>% 
      ungroup()
    models_notplco <- modelstest_vec[modelstest_vec!="plcom2012"]
    for(n in models_notplco){
      bootstrapAUCvec1 <- bootstraping_res_pooled_t %>% filter(model=="plcom2012") %>% 
        pull(bootstrappedAUC_vec)
      bootstrapAUCvec2 <- bootstraping_res_pooled_t %>% filter(model==n) %>% 
        pull(bootstrappedAUC_vec)
      name_mod_compare <- paste0("plcom2012-",n)
      # Main AUC diff 
      AUC_mod1 <- AUC_models_pooled_t$AUCpooled[AUC_models_pooled_t$model=="plcom2012"]
      AUC_mod2 <- AUC_models_pooled_t$AUCpooled[AUC_models_pooled_t$model==n]
      maindiff_integer <- abs(AUC_mod1-AUC_mod2) ## SHOULD IT BE ABS or doesnt matter?
      # pvalue computing
      p_valuediff <- pvalue_bootstrapAUC(bootstrapAUCvec1, bootstrapAUCvec2, maindiff_integer)[2] #chose second method like xiaoshuang's (both gave almost similar results)
      p_valuediff_pm <- pvalue_bootstrapAUC(bootstrapAUCvec1, bootstrapAUCvec2, maindiff_integer)[1]
      data_inter_p <- AUC_models_pooled_t %>% mutate(pvalue=p_valuediff, pvalue_ref=p_valuediff_pm,
                                                     time=t, 
                                                     AUC_CI_pooled=paste0(formatC(AUCpooled,digits=2, format = "f"), 
                                                                          " (", 
                                                                          formatC(lowerCIpooled, digits=2, format="f"), "-", 
                                                                          formatC(upperCIpooled, digits=2, format = "f"), ")"),
                                                     pvaluecomp=name_mod_compare) 
      dataframe_summmaryROC <- dplyr::bind_rows(dataframe_summmaryROC, data_inter_p)
    }
  }
  if(length(modelstest_vec)>2){
    dataframe_summmaryROC <- dataframe_summmaryROC %>% pivot_wider(names_from = pvaluecomp, values_from = pvalue_ref) #use ref methods published in the protol to get pvalues
  }
  return(dataframe_summmaryROC)
}

## ROC curves for one dataset ###
getDATA_ROCplot <- function(datax_1imp,modelstest_vec, pooledsummary_data){
  pooledsummary_data <- pooledsummary_data %>% filter(model %in%modelstest_vec )
  ROC_mainplot <- data.frame(NULL)
  for(t in 1:3){
    LC_event_t <- paste0("LC_event", t,"y")
    dataimp_t <- datax_1imp %>% filter(!is.na(!!sym(LC_event_t)))
    for (mod in modelstest_vec){
      mod_t <- ifelse(mod!="plcom2012",paste0(mod, "_", t,"y"), mod)
      ROC_data_interim <- with(dataimp_t, WeightedROC(guess=get(mod_t), label=get(LC_event_t), 
                                                      weight=U19_weight)) %>% mutate(time=t, model=mod)
      ROC_mainplot <- bind_rows(ROC_mainplot, ROC_data_interim)
    }
  }
  
  ROC_mainplot <- ROC_mainplot %>% merge(pooledsummary_data, 
                                         by=c("time","model"), all.x = T, all.y = F) 
  return(ROC_mainplot)
  
}

ROC_time_plot <- function(data_plot, title_p, uspstfy="USPSTF2021", linesize=0.75, legendsize, legendpos,
                          axissize=12,axistitlesize=15){
  if(uspstfy!="no"){
    TPR_USPSTF_f <- ifelse(uspstfy=="USPSTF2013", "TPR_2013", "TPR_2021")
    FPR_USPSTF_f <- ifelse(uspstfy=="USPSTF2013", "FPR_2013", "FPR_2021")
  }
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

   
  
  if(uspstfy!="no"){
    plot <- plot+#USPSTF20
      # Add locations
      geom_segment(aes(x=0.08,xend=get(FPR_USPSTF_f),y=get(TPR_USPSTF_f),yend=get(TPR_USPSTF_f)),
                   linetype="dotted", color="#E18727FF")+
      geom_segment(aes(x=get(FPR_USPSTF_f),xend=get(FPR_USPSTF_f),y=0.05,yend=1),linetype="dotted", 
                   color="#E18727FF")+
      geom_point(aes(x = get(FPR_USPSTF_f), y = get(TPR_USPSTF_f)),shape = 23,color="#E18727FF",fill="#E18727FF", size=5)+
      # add text annotation
      geom_text(data=. %>% distinct(get(TPR_USPSTF_f), .keep_all = T) %>% mutate(USPSTF=uspstfy), 
                aes(x = get(FPR_USPSTF_f)+0.03, y = get(TPR_USPSTF_f), label = USPSTF), color="#E18727FF", hjust=0,
                fontface = "bold")+
      geom_text(data=. %>% distinct(get(TPR_USPSTF_f), .keep_all = T),
                aes(x=0.04, y=get(TPR_USPSTF_f),
                    label=formatC(get(TPR_USPSTF_f),digits=2, format="f")), color="#E18727FF")+
      geom_text(data=. %>% distinct(get(TPR_USPSTF_f), .keep_all = T),
                aes(x=get(FPR_USPSTF_f), y=0.02,
                    label=formatC(get(FPR_USPSTF_f),2, digits=2, format="f")), color="#E18727FF")+
      # fixed spe and add sens for plco and prot model #aes(x=FPR, y=TPR, color=model)
      geom_segment(data=. %>% mutate(closesetFPR=abs(FPR-get(FPR_USPSTF_f)))%>%group_by(time,model) %>%
                     slice_min(closesetFPR) %>% slice_max(TPR),
                   aes(x=0.08,xend=get(FPR_USPSTF_f),y=TPR,yend=TPR, color=model),linetype="dotted")+
      geom_text(data=. %>% mutate(closesetFPR=abs(FPR-get(FPR_USPSTF_f)))%>%group_by(time,model) %>%
                  slice_min(closesetFPR) %>% slice_max(TPR),
                aes(x=0.04, y=TPR,label=formatC(TPR,digits=2, format = "f"), color=model),show.legend = FALSE)
  }
  print(plot)
  return(plot)
}


## calibration plots ###
getdata_calplot <- function(datax_1imp, modelriskpred, summary_EO_res){ #proteinmodelrisk_
  data_plot <- data.frame(NULL)
  for(t in 1:3){
    summary_EO_data <- summary_EO_res %>% filter(model==modelriskpred, time==t)
    LC_event_time <- paste0("LC_event",t,"y")
    risk_predtime <- paste0(modelriskpred,"_",t,"y")
    data_plot_t <- datax_1imp %>% 
      mutate(U19_weight=ifelse(is.na(!!sym(LC_event_time)), 0, U19_weight), #so we don't sum cases after the event that are NA
             ntile_predrisk=dineq::ntiles.wtd(get(risk_predtime),5, U19_weight)) 
    data_plot_t <- data_plot_t %>% group_by(ntile_predrisk) %>% 
      dplyr::summarize(n=sum(U19_weight),
                       observed_risk=sum(get(LC_event_time), na.rm = TRUE)/sum(U19_weight),
                       lowerCI=observed_risk-1.96*sqrt((observed_risk*(1-observed_risk))/n),
                       upperCI=observed_risk+1.96*sqrt((observed_risk*(1-observed_risk))/n),
                       predicted_risk=sum(get(risk_predtime)*U19_weight)/sum(U19_weight))%>% # i am using mean, should i be using medians?
      mutate(lowerCI=ifelse(lowerCI<0, 0, lowerCI),
             upperCI=ifelse(upperCI>1,1, upperCI),#cap confidence intervals for risk
             E_O_overall=summary_EO_data$cal_CI, 
             time=t) 
    data_plot <- bind_rows(data_plot, data_plot_t)
  }
  
  return(data_plot)
}


calplot_function <- function(datacalplot_t, line_width=0.75,
                             annotsize=4,axissize=12,axistitlesize=15){
  observed_risk_v <- max(datacalplot_t$observed_risk)#+0.0015
  max_y <- max(datacalplot_t$upperCI)#+0.001
  max_x <- max(datacalplot_t$predicted_risk)#+0.0007
  max_eq <- max(max_y, max_x)
  
  plot <- datacalplot_t %>%  
    ggplot(aes(y=observed_risk, x=predicted_risk))+geom_point(color="#BC3C29FF")+ #color=model
    geom_pointrange(aes(ymin=lowerCI, ymax=upperCI), color="#BC3C29FF")+
    geom_smooth(method="lm", se=FALSE, color="#BC3C29FF", linewidth = line_width)+
    geom_abline(slope=1, col="grey")+
    labs(y="Observed risk", x="Predicted risk", title=paste0("t=",t), color="")+
    coord_equal()+
    theme_bw()+
    
    theme(legend.position="none",
          axis.text = element_text(size=axissize),
          axis.title=element_text(size=axistitlesize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    annotate("text", x=0, y=max_eq ,
             label=unique(datacalplot_t$E_O_overall),size=annotsize,
             hjust=0, vjust=1.5)+
    scale_y_continuous(labels = scales::percent, expand = c(0, 0),  limits = c(0, max_eq),
                       breaks = function(x) scales::breaks_extended(n = 4)(x)[scales::breaks_extended(n = 4)(x) != 0])+
    scale_x_continuous(labels = scales::percent, expand = c(0, 0),  limits = c(0, max_eq),
                       breaks = scales::breaks_extended(n = 4))
  return(plot)
}

## Performance USPSTF ###
perf_uspstf <- function(datax, uspstfy,LC_event_t){
  TP <- datax %>% filter(!!sym(LC_event_t)==1 & !!sym(uspstfy)==1) %>% nrow()
  FN <- datax %>% filter(!!sym(LC_event_t)==1 & !!sym(uspstfy)==0) %>% nrow() 
  FP <- datax %>% filter(!!sym(LC_event_t)==0 & !!sym(uspstfy)==1) %>% pull(U19_weight) %>% sum() 
  TN <- datax %>% filter(!!sym(LC_event_t)==0 & !!sym(uspstfy)==0) %>% pull(U19_weight) %>% sum()  
  
  TPR <- TP/(TP+FN)
  FPR <- FP/(FP+TN)
  return(data.frame(TPR=TPR, FPR=FPR))  
}


## venn diagrams ###
#library(VennDiagram)         
venn_cases_time <- function(datax, t, protein_thresh, plco_thresh, USPSTFY_t_yn,  pop_text, print.mode){ # make sure to only include cases
  
  if(USPSTFY_t_yn=="n"){
    
    ggvenn_list_cases_uw <- list(`PLCO cases`=rep(pull(filter(datax, plcom2012>=plco_thresh), U19_id),
                                                  round(pull(filter(datax, plcom2012>=plco_thresh), U19_weight),0)), 
                                 `Protein cases`=rep(pull(filter(datax, !!sym(paste0("proteinmodelrisk_",t,"y"))>=protein_thresh), U19_id),
                                                     round(pull(filter(datax, !!sym(paste0("proteinmodelrisk_",t,"y"))>=protein_thresh), U19_weight),0)))
    print("ok1")
    ggvenn_list_cases <- list(PLCO=paste0(ggvenn_list_cases_uw$`PLCO cases`, "_", sequence(rle(ggvenn_list_cases_uw$`PLCO cases`)$lengths)),
                              protein=paste0(ggvenn_list_cases_uw$`Protein cases`, "_", sequence(rle(ggvenn_list_cases_uw$`Protein cases`)$lengths))) #sequence rle adds id idenx so that the same id is counted each time 
    print("ok2")
    VennDiagram_list<- venn.diagram(x = ggvenn_list_cases,filename = NULL, fill=c("#B3BBDF", "#E6ACA4"), 
                                    euler.d=TRUE, cat.col=c("#7D829C", "#B88983"),
                                    category.names = c(eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nPLCOm2012'))"))),
                                                       eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nINTEGRAL Risk\nmodel'))")))),
                                    cat.cex=1.25, cat.default.pos="outer",
                                    cat.pos=c(-45,45), cat.dist=0.095,
                                    print.mode=print.mode, cex=1.5,sigdigs=4, 
                                    ext.text=TRUE,
                                    margin=0.08)
    
    print(cowplot::plot_grid(VennDiagram_list))
  }
  else if(USPSTFY_t_yn=="USPSTF2021"){
    ggvenn_list_cases_uw <- list(`PLCO cases`=rep(pull(filter(datax, plcom2012>=plco_thresh), U19_id),
                                                  round(pull(filter(datax, plcom2012>=plco_thresh), U19_weight),0)), 
                                 `Protein cases`=rep(pull(filter(datax, !!sym(paste0("proteinmodelrisk_",t,"y"))>=protein_thresh), U19_id),
                                                     round(pull(filter(datax, !!sym(paste0("proteinmodelrisk_",t,"y"))>=protein_thresh), U19_weight),0)),
                                 `USPSTF2021`=rep(pull(filter(datax, USPSTF2021==1), U19_id),
                                                  round(pull(filter(datax, USPSTF2021==1), U19_weight),0)))
    ggvenn_list_cases <- list(PLCO=paste0(ggvenn_list_cases_uw$`PLCO cases`, "_", sequence(rle(ggvenn_list_cases_uw$`PLCO cases`)$lengths)),
                              protein=paste0(ggvenn_list_cases_uw$`Protein cases`, "_", sequence(rle(ggvenn_list_cases_uw$`Protein cases`)$lengths)),
                              USPSTF2021= paste0(ggvenn_list_cases_uw$USPSTF2021, "_", sequence(rle(ggvenn_list_cases_uw$USPSTF2021)$lengths))) #sequence rle adds id idenx so that the same id is counted each time 
    VennDiagram_list<-  venn.diagram(x = ggvenn_list_cases,filename = NULL, fill=c("#B3BBDF", "#E6ACA4", "#A4D7C6"), 
                                     euler.d=TRUE, cat.col=c("#7D829C", "#B88983","#a9c2ba"),
                                     category.names = c(eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nPLCOm2012'))"))),
                                                        eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nINTEGRAL Risk\nmodel'))"))),
                                                        eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nUSPSTF2021'))")))),
                                     cat.cex=1.25, cat.default.pos="outer",
                                     cat.pos=c(-80,80,180), cat.dist=0.15,
                                     print.mode=print.mode, cex=1.5,sigdigs=4, 
                                     ext.text=TRUE,
                                     margin=0.08)
    print(cowplot::plot_grid(VennDiagram_list))
    
  }
  
  else if(USPSTFY_t_yn=="USPSTF2013"){
    ggvenn_list_cases_uw <- list(`PLCO cases`=rep(pull(filter(datax, plcom2012>=plco_thresh), U19_id),
                                                  round(pull(filter(datax, plcom2012>=plco_thresh), U19_weight),0)), 
                                 `Protein cases`=rep(pull(filter(datax, !!sym(paste0("proteinmodelrisk_",t,"y"))>=protein_thresh), U19_id),
                                                     round(pull(filter(datax, !!sym(paste0("proteinmodelrisk_",t,"y"))>=protein_thresh), U19_weight),0)),
                                 `USPSTF2013`=rep(pull(filter(datax, USPSTF2013==1), U19_id),
                                                  round(pull(filter(datax, USPSTF2013==1), U19_weight),0)))
    ggvenn_list_cases <- list(PLCO=paste0(ggvenn_list_cases_uw$`PLCO cases`, "_", sequence(rle(ggvenn_list_cases_uw$`PLCO cases`)$lengths)),
                              protein=paste0(ggvenn_list_cases_uw$`Protein cases`, "_", sequence(rle(ggvenn_list_cases_uw$`Protein cases`)$lengths)),
                              USPSTF2013= paste0(ggvenn_list_cases_uw$USPSTF2013, "_", sequence(rle(ggvenn_list_cases_uw$USPSTF2021)$lengths)))
    
    VennDiagram_list<-  venn.diagram(x = ggvenn_list_cases,filename = NULL, fill=c("#B3BBDF", "#E6ACA4", "#A4D7C6"), 
                                     euler.d=TRUE, cat.col=c("#7D829C", "#B88983","#a9c2ba"),
                                     category.names = c(eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nPLCOm2012'))"))),
                                                        eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nINTEGRAL Risk\nmodel'))"))),
                                                        eval(parse(text=paste0("expression(bold('",pop_text,"\nscreened by\nnUSPSTF2013'))")))),
                                     cat.cex=1.25, cat.default.pos="outer",
                                     cat.pos=c(-80,80,180), cat.dist=0.15,
                                     print.mode=print.mode, cex=1.5,sigdigs=4, 
                                     ext.text=TRUE,
                                     margin=0.08)
    print(cowplot::plot_grid(VennDiagram_list))
  }
  # sigdig only does significant digits in frequencies-> need to adjust manually the labels 
  x <- !sapply( sapply(VennDiagram_list,'[[',"label") , is.null )
  ind <- seq.int(length(x))
  ind <- ind[x]
  vec_round <- c(NULL)
  vec_raw <- c(NULL)
  vec_raw_perc <- c(NULL)
  # length labels group 
  labels <- sapply(VennDiagram_list, function(x) x$label)
  labels <- labels[!is.na(labels)]
  labels <- labels[grepl("%", labels)]
  for (i in 1:length(labels)) { #length labels
    print(labels[i])
    full_label <- `[[`(`[[`( VennDiagram_list , ind[i] ), "label") 
    print("full_label")
    print(full_label)
    raw_number <- str_remove(full_label, "\\n.*")
    print("raw_numner")
    print(raw_number)
    perc <- as.numeric(str_match(full_label, "\\(*(.*?)%")[,2])
    round_per <- round(perc, 0)
    print("roundper")
    print(round_per)
    if(any(grepl("raw",print.mode))){
      full_label_updated <- paste0(raw_number,"\n(",round_per,"%)")
    }
    else{
      full_label_updated <- paste0(round_per,"%")
    }
    # update labels in venndiag
    print("label updated")
    print(full_label_updated)
    `[[`(`[[`( VennDiagram_list , ind[i] ), "label") <- full_label_updated
    # to test if sum is 100 
    vec_raw_perc <- c(vec_raw_perc, perc)
    vec_round <- c(vec_round, round_per)
    vec_raw <- c(vec_raw, raw_number)
  }
  print("vec_raw_perc")
  print(vec_raw_perc)
  print("vec round")
  print(vec_round)
  print("vec raw")
  print(vec_raw)
  sum_round_per <- sum(vec_round)
  print(sum_round_per)
  if(sum_round_per>100){
    dec <- str_replace(vec_raw_perc, ".*\\.","0.") # precentages not rounded -> 
    ind_change <- which.min(dec)
    print("ind change")
    print(ind_change)
    # change the indexed label to round down if upper than 100 
    perc <- vec_raw_perc[ind_change]
    round_per <- floor(perc)
    raw_number_up <- vec_raw[ind_change]
    if(any(grepl("raw",print.mode))){
      full_label_updated <- paste0(raw_number_up,"\n(",round_per,"%)")
    }
    else{
      full_label_updated <- paste0(round_per,"%")
    }
    # update labels in venndiag
    print("label updated")
    print(full_label_updated)
    `[[`(`[[`( VennDiagram_list , ind[ind_change] ), "label") <- full_label_updated
  }
  if(sum_round_per<100){
    dec <- str_remove(vec_raw_perc, ".*\\.")
    ind_change <- which.max(dec)
    # change the indexed label to round up if lower than 100 
    perc <- vec_raw_perc[ind_change]
    print("if less 100")
    print(perc)
    round_per <- ceiling(perc)
    print(round_per)
    raw_number_up <- vec_raw[ind_change]
    
    if(any(grepl("raw",print.mode))){
      full_label_updated <- paste0(raw_number_up,"\n(",round_per,"%)") #raw_number
    }
    else{
      full_label_updated <- paste0(round_per,"%")
    }
    # update labels in venndiag
    print("label updated")
    print(full_label_updated)
    `[[`(`[[`( VennDiagram_list , ind[ind_change] ), "label") <- full_label_updated
  }
  
  venn_plot <- cowplot::plot_grid(VennDiagram_list)
  return(venn_plot)
}


getdata_lorezplot <- function(datax,riskmodel_time, LC_event_time,t){
  # REF USPSTF pop elig & cases elig 
  pop_eligUSPSTF <- datax %>% filter(LC_event==0, !is.na(!!sym(LC_event_time))) %>% 
    dplyr::select(c("USPSTF2013", "USPSTF2021", LC_event_time, "U19_weight")) %>% 
    pivot_longer(c("USPSTF2013", "USPSTF2021"), names_to="USPSTF_Y", values_to="elig_yn") %>% 
    group_by(USPSTF_Y) %>% 
    dplyr::summarise(n=sum(U19_weight), pop_elig=sum(elig_yn*U19_weight),
                     prop_pop_elig=pop_elig/n)
  
  toteligUSPSTF_cases <-datax %>% filter(!is.na(!!sym(LC_event_time)), !!sym(LC_event_time)==1) %>% 
    dplyr::select(c("USPSTF2013", "USPSTF2021", !!sym(LC_event_time))) %>% 
    pivot_longer(c("USPSTF2013", "USPSTF2021"), names_to="USPSTF_Y", values_to="elig_yn") %>% 
    group_by(USPSTF_Y) %>% 
    dplyr::summarise(ncases=sum(!!sym(LC_event_time)), 
                     cases_elig=sum(elig_yn*!!sym(LC_event_time)),
                     prop_cases_elig=cases_elig/ncases)
  
  pop_eligUSPSTF <- merge(pop_eligUSPSTF,toteligUSPSTF_cases, by="USPSTF_Y")
  
  # Annot thresholds based on USPSTFS for both models
  lorenz_data_baseline <- datax %>%filter(!is.na(!!sym(LC_event_time)), LC_event==0)  %>% 
    pivot_longer(c("plcom2012",riskmodel_time), names_to="model", values_to="model risk") %>% 
    dplyr::select(model, `model risk`, U19_weight) %>% 
    group_by(model) %>% 
    arrange(desc(`model risk`)) %>% mutate(cumsum_elig=cumsum(U19_weight),
                                           total_elig=sum(U19_weight), 
                                           elig_prop=cumsum_elig/total_elig) %>% ungroup()  %>% 
    mutate(diff_USPSTF21=cumsum_elig -pop_eligUSPSTF$pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"],#abs(elig_prop -pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"]),
           diff_USPSTF13=cumsum_elig -pop_eligUSPSTF$pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"]) %>% 
    pivot_longer(c("diff_USPSTF21", "diff_USPSTF13"), names_to="diff set", values_to = "valuediff") 
  
  risks_prot_annot_lower <- lorenz_data_baseline%>% 
    group_by(model,`diff set`) %>% filter(valuediff<=0) %>% slice_min(abs(valuediff),n=1)%>%
    ungroup() %>% 
    dplyr::select(c("model", "model risk", "diff set", "cumsum_elig" )) 
  
  risks_prot_annot_upper <- lorenz_data_baseline%>% 
    group_by(model,`diff set`) %>% filter(valuediff>0) %>% slice_min(valuediff,n=1)%>%
    ungroup() %>% 
    dplyr::select(c("model", "model risk", "diff set", "cumsum_elig" )) 
  
  risks_prot_annot <- bind_rows(risks_prot_annot_lower,risks_prot_annot_upper)
  risks_prot_annot <- risks_prot_annot %>% 
    mutate(USPSTF_elig=case_when(`diff set`=="diff_USPSTF13"~pop_eligUSPSTF$pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"],
                                 `diff set`=="diff_USPSTF21"~pop_eligUSPSTF$pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"])) %>% 
    group_by(model, `diff set`) %>% 
    mutate(riskthreshold_pooled=min(`model risk`)+
             (max(`model risk`)- min(`model risk`))*((USPSTF_elig-max(cumsum_elig))/(min(cumsum_elig)-max(cumsum_elig))))
  #R= r1+(r2-r1)*((ctarget-c1)/(c2-c1))
  risks_prot_annot <- risks_prot_annot %>% distinct(riskthreshold_pooled, .keep_all = T) %>% ungroup()
  
  #add how many cases would be eligible
  total_ycases <- sum(datax[,LC_event_time], na.rm=TRUE)
  cases_thresh_USPSTF13 <- nrow(datax %>% 
                                  filter(!!sym(riskmodel_time)>=risks_prot_annot$riskthreshold_pooled[risks_prot_annot$model==riskmodel_time & risks_prot_annot$`diff set`=="diff_USPSTF13"],
                                         !!sym(LC_event_time)==1))
  
  cases_thresh_USPSTF21 <- nrow(datax %>% 
                                  filter(!!sym(riskmodel_time)>=risks_prot_annot$riskthreshold_pooled[risks_prot_annot$model==riskmodel_time &risks_prot_annot$`diff set`=="diff_USPSTF21"],
                                         !!sym(LC_event_time)==1))
  
  cases_plco_13 <- nrow(datax %>% filter(plcom2012>=risks_prot_annot$riskthreshold_pooled[risks_prot_annot$model=="plcom2012" & risks_prot_annot$`diff set`=="diff_USPSTF13"], 
                                         !!sym(LC_event_time)==1))
  cases_plco_21 <- nrow(datax %>% filter(plcom2012>=risks_prot_annot$riskthreshold_pooled[risks_prot_annot$model=="plcom2012" & risks_prot_annot$`diff set`=="diff_USPSTF21"], 
                                         !!sym(LC_event_time)==1))
  
  risks_prot_annot <- risks_prot_annot %>% 
    mutate(cases_elig=case_when(model=="plcom2012" &`diff set`=="diff_USPSTF13"~cases_plco_13, 
                                model=="plcom2012" &`diff set`=="diff_USPSTF21"~cases_plco_21,
                                model==riskmodel_time &`diff set`=="diff_USPSTF13"~cases_thresh_USPSTF13, 
                                model==riskmodel_time &`diff set`=="diff_USPSTF21"~cases_thresh_USPSTF21),
           prop_caseselig=cases_elig/total_ycases)
  
  prot_thresh_f <- risks_prot_annot$riskthreshold_pooled[risks_prot_annot$model==riskmodel_time &risks_prot_annot$`diff set`=="diff_USPSTF21"]
  thresh_plco <- risks_prot_annot$riskthreshold_pooled[risks_prot_annot$model=="plcom2012" &risks_prot_annot$`diff set`=="diff_USPSTF21"]
  
  lorenz_data_plot <-  datax %>%filter(!is.na(!!sym(LC_event_time))) %>%  
    pivot_longer(c("plcom2012", riskmodel_time), names_to="model", values_to="model risk") %>% 
    dplyr::select(LC_event_time, model, `model risk`, U19_weight) %>% 
    group_by(model) %>% 
    arrange(desc(`model risk`)) %>% mutate(cumsum_elig=cumsum(U19_weight),
                                           cumsum_cases=cumsum(!!sym(LC_event_time)),
                                           total_elig=sum(U19_weight), 
                                           total_cases=sum(!!sym(LC_event_time), na.rm = TRUE), 
                                           elig_prop=cumsum_elig/total_elig, 
                                           cases_elig_prop=cumsum_cases/total_cases) %>% ungroup() 
  
  lorenz_data_plot  <- risks_prot_annot %>% 
    merge(lorenz_data_plot, by=c("model", "model risk"), all.y=TRUE) %>% 
    mutate(annot_risk=ifelse(!is.na(`diff set`), "yes", "no"), 
           annot_risk_full=ifelse(annot_risk=="yes",
                                  paste0(ifelse(model=="plcom2012","PLCOm2012: ", "INTEGRAL Risk model: "), 
                                         round(100*cases_elig_prop,0),"% of cancers",
                                         ifelse(model=="plcom2012"," (6y risk\u2265", paste0(" (",t,"y risk\u2265")),
                                         round(100*riskthreshold_pooled,3), "%) "),
                                  NA)) %>% 
    group_by(`diff set`) %>% 
    mutate(annot_risk_comb=paste(annot_risk_full, collapse = "\n"))
  
  return(list(lorenz_data_plot,pop_eligUSPSTF,prot_thresh_f,thresh_plco))
}
#

lorenz_plot_function <- function(datax, riskmodel_time, LC_event_time,t,size_txt){
  full_plotdata <- getdata_lorezplot(datax,riskmodel_time, LC_event_time,t)
  lorenz_data_plot <- full_plotdata[[1]]
  pop_eligUSPSTF <- full_plotdata[[2]]
  v_lines <- round(pop_eligUSPSTF$prop_pop_elig, digits=3)
  annot_x_US13 <- round(pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"], digits=2)                         
  annot_y_US13 <- round(pop_eligUSPSTF$prop_cases_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"], digits=2) 
  annot_x_US21 <- round(pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"], digits=2) 
  annot_y_US21 <- round(pop_eligUSPSTF$prop_cases_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"], digits=2) 
  
  annot_text_US13_raw <- paste0(unique(lorenz_data_plot[lorenz_data_plot$`diff set`=="diff_USPSTF13" & !is.na(lorenz_data_plot$`diff set`),
                                                        "annot_risk_comb", drop=T]),
                                "\nUSPSTF2013: ",round(pop_eligUSPSTF$prop_cases_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"]*100,0), "% of cancers ")
  annot_text_US13_split <- stringr::str_split(annot_text_US13_raw, "\n")[[1]]
  annot_text_US13 <- eval(parse(text=paste0("'<span style=\"color:#bc3c29ff\">",
                                            annot_text_US13_split[grepl("INTEGRAL",annot_text_US13_split)],
                                            "</span><br><span style=\"color:#0072b5ff\">",
                                            annot_text_US13_split[grepl("PLCOm2012",annot_text_US13_split)],
                                            "</span><br><span style=\"color:#E18727FF\">",
                                            annot_text_US13_split[grepl("USPSTF",annot_text_US13_split)],
                                            "</span><br><span style=\"color:gray21\">",
                                            round(pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"]*100,0),"% population screened '")))
  
  
  annot_text_US21_raw <- paste0(unique(lorenz_data_plot[lorenz_data_plot$`diff set`=="diff_USPSTF21" & !is.na(lorenz_data_plot$`diff set`) ,
                                                        "annot_risk_comb", drop=T]),
                                "\nUSPSTF2021: ",round(pop_eligUSPSTF$prop_cases_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"]*100,0), "% of cancers")
  annot_text_US21_split <- stringr::str_split(annot_text_US21_raw, "\n")[[1]]
  annot_text_US21 <- eval(parse(text=paste0("'<span style=\"color:#bc3c29ff\">",
                                            annot_text_US21_split[grepl("INTEGRAL",annot_text_US21_split)],
                                            "</span><br><span style=\"color:#0072b5ff\">",
                                            annot_text_US21_split[grepl("PLCOm2012",annot_text_US21_split)],
                                            "</span><br><span style=\"color:#E18727FF\">",
                                            annot_text_US21_split[grepl("USPSTF",annot_text_US21_split)],
                                            "</span><br><span style=\"color:gray21\">",
                                            round(pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"]*100,0),"% population screened'")))
  
  plot <- lorenz_data_plot %>% 
    ggplot(aes(x=elig_prop, y=cases_elig_prop, color=model))+
    geom_point(data=. %>% filter(annot_risk=="yes"))+
    geom_path()+
    geom_vline(xintercept =v_lines, linetype = "dashed", color="gray21" )+
    #segments plco to link to annotation to USPSTF 21
    geom_segment(x = lorenz_data_plot$elig_prop[lorenz_data_plot$model=="plcom2012"& lorenz_data_plot$`diff set`=="diff_USPSTF21" & !is.na(lorenz_data_plot$`diff set`)],
                 y =  lorenz_data_plot$cases_elig_prop[lorenz_data_plot$model=="plcom2012"& lorenz_data_plot$`diff set`=="diff_USPSTF21" & !is.na(lorenz_data_plot$`diff set`)],
                 xend =pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"]+0.02, yend = annot_y_US21, color="#0072b5ff")+
    # # plco to USPSYF 13
    geom_segment(x = lorenz_data_plot$elig_prop[lorenz_data_plot$model=="plcom2012"& lorenz_data_plot$`diff set`=="diff_USPSTF13" & !is.na(lorenz_data_plot$`diff set`)],
                 y =  lorenz_data_plot$cases_elig_prop[lorenz_data_plot$model=="plcom2012"& lorenz_data_plot$`diff set`=="diff_USPSTF13" & !is.na(lorenz_data_plot$`diff set`)],
                 xend = pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"]+0.02, yend = annot_y_US13, color="#0072b5ff")+
    # segments protein model to link to annotation
    # to USPSTF 21
    geom_segment(x = lorenz_data_plot$elig_prop[lorenz_data_plot$model==riskmodel_time& lorenz_data_plot$`diff set`=="diff_USPSTF21" & !is.na(lorenz_data_plot$`diff set`)],
                 y =  lorenz_data_plot$cases_elig_prop[lorenz_data_plot$model==riskmodel_time& lorenz_data_plot$`diff set`=="diff_USPSTF21" & !is.na(lorenz_data_plot$`diff set`)],
                 xend = pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"]+0.02, yend = annot_y_US21, color="#bc3c29ff")+
    # # to uspsyf 13
    geom_segment(x = lorenz_data_plot$elig_prop[lorenz_data_plot$model==riskmodel_time& lorenz_data_plot$`diff set`=="diff_USPSTF13" & !is.na(lorenz_data_plot$`diff set`)],
                 y =  lorenz_data_plot$cases_elig_prop[lorenz_data_plot$model==riskmodel_time& lorenz_data_plot$`diff set`=="diff_USPSTF13" & !is.na(lorenz_data_plot$`diff set`)],
                 xend = pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"]+0.02, yend = annot_y_US13, color="#bc3c29ff")+
    # # segments USPSTF 13 to annot
    geom_segment(x=pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"],
                 y=pop_eligUSPSTF$prop_cases_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"],
                 xend=pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2013"]+0.02,
                 yend=annot_y_US13, color="#E18727FF")+
    # # segments USPSTF 21 to annot
    geom_segment(x=pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"],
                 y=pop_eligUSPSTF$prop_cases_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"],
                 xend=pop_eligUSPSTF$prop_pop_elig[pop_eligUSPSTF$USPSTF_Y=="USPSTF2021"]+0.02,
                 yend=annot_y_US21, color="#E18727FF")+
    
    # annotation text 
    annotate(geom="richtext", x=c(annot_x_US21+0.02, annot_x_US13+0.02), 
             y=c(annot_y_US21, annot_y_US13), 
             label=c(annot_text_US21,annot_text_US13), 
             hjust=0, size=4.5)+
    annotate("point",x = pop_eligUSPSTF$prop_pop_elig, y = pop_eligUSPSTF$prop_cases_elig,shape = 23,color="#E18727FF",fill="#E18727FF", size=5)+
    #annotate("label", x=USPSTF_perelig+0.02, y=USPSTF_caseelig, label=USPSTF_annot, color="#E18727FF", hjust=0)+
    theme_bw()+
    labs(x="Eligible population", y="Cases eligible", color="Model used")+
    scale_x_continuous(labels = scales::percent, expand = c(0, 0),limits = c(0,NA),
                       breaks=c(0,0.25,0.50,0.75,1))+
    scale_y_continuous(labels = scales::percent, expand = c(0, 0),limits = c(0,NA),
                       breaks=c(0.25,0.50,0.75,1))+
    scale_color_manual(values=c("#0072b5ff","#bc3c29ff"), 
                       labels=c("PLCOm2012", "INTEGRAL Risk Model"))
  
  return(plot)
}


# threshold identification 
get_prot_threshold <- function(datax,nbscreenedref,model_time ){
  prot_tresh <- datax %>%
    arrange(desc(!!sym(model_time))) %>%
    dplyr::select(all_of(c("U19_id","plcom2012", model_time, "U19_weight"))) %>%
    mutate(csum=cumsum(U19_weight)) %>% 
    mutate(diff_toref=abs(nbscreenedref-csum)) %>%
    slice_min(diff_toref, n=1) %>% pull(!!sym(model_time))
  return(prot_tresh)
}


# characteristics population selected 
wtd_distributions <- function(var, weight, rounding=NULL, brk=F){
  if(is.null(rounding)){
    rounding=2
  }
  mean <- formatC(Hmisc::wtd.mean(var,weight, na.rm = TRUE), digits=rounding, format="f")
  Q1 <- formatC(Hmisc::wtd.quantile(var,weight,probs=0.25, na.rm = TRUE)[[1]], digits=rounding, format = "f")
  Q3 <- formatC(Hmisc::wtd.quantile(var,weight,probs=0.75, na.rm = TRUE)[[1]], digits=rounding, format = "f")
  
  if(isTRUE(brk)){mean_Q <- paste0(mean, "\n(", Q1, "-", Q3, ")")}
  else {mean_Q <- paste0(mean, " (", Q1, "-", Q3, ")")}
    
  return(mean_Q)
}

summary_eligibility_charac <- function(data_charac_x,t,case_status, rounding, brk=F){
  LC_event_t <- paste0("LC_event",t,"y")
  model_prot_y <- paste0("elig_prot",t,"y")
  tot_cases <- sum(data_charac_x[,LC_event_t], na.rm=TRUE)
  data_charac_x <- data_charac_x %>% mutate(quit_years_tot=quit_years,
                                            quit_years=ifelse(smoke_status==3, NA, quit_years)) %>% 
    filter(!is.na(!!sym(LC_event_t))) %>% 
    pivot_longer(cols=all_of(c("PLCO_08", "elig_USPSTF2021", model_prot_y)), 
                 names_to = "elig_by", values_to = "elig_stat") %>% 
    filter(elig_stat!="Non elig") 
  if(case_status=="cases"){
    data_charac_x <- data_charac_x %>% filter(!!sym(LC_event_t)==1)
  }
  if (case_status=="baseline"){
    data_charac_x <- data_charac_x %>% mutate(U19_weight=ifelse(LC_event==1,0,U19_weight))
  }
  
  data_charc <- data_charac_x%>% 
    group_by(elig_by) %>% 
    dplyr::summarise(N_screened=round(sum(U19_weight), rounding), 
                     N_cases_screened=sum(!!sym(LC_event_t)), 
                     Age=wtd_distributions(age, U19_weight, rounding,brk), 
                     CPD=wtd_distributions(intensity, U19_weight,rounding,brk),
                     `Years smoked`=wtd_distributions(years_smoked, U19_weight,rounding,brk),
                     quit_years_former=wtd_distributions(quit_years, U19_weight, rounding,brk), 
                     quit_years_all=wtd_distributions(quit_years_tot, U19_weight, rounding,brk),
                     per_former=100*sum((smoke_status==2)*U19_weight)/N_screened) %>% 
    mutate(time=t,per_total_cases=formatC(100*N_cases_screened/tot_cases, digits=0, format="f"),
           per_former=paste0(round(per_former,0), "%"),
           N_cases_screened_per=paste0(N_cases_screened, " (", per_total_cases, "%)")) %>% 
    dplyr::select(elig_by,time, N_screened, N_cases_screened_per,Age, CPD, `Years smoked`, 
                  quit_years_former, per_former,quit_years_all)
  print(data_charc)
  return(data_charc)
}
