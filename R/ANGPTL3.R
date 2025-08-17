#prevalent和incident
ukbb_angptl3_phenotype_time5 <- read.delim("/processed_data/medenlian_all_gene/two_sample_mr/hepatokine/angptl3_ukbb/ukbb_angptl3_phenotype_time5.txt")
#date formating
ukbb_angptl3_phenotype_time5$Time_assess <- as.Date(ukbb_angptl3_phenotype_time5$Time_assess)
ukbb_angptl3_phenotype_time5$Time_blood_sample_collected0 <- as.Date(ukbb_angptl3_phenotype_time5$Time_blood_sample_collected0)
ukbb_angptl3_phenotype_time5$date_last_personal_contact <- as.Date(ukbb_angptl3_phenotype_time5$date_last_personal_contact)

ukbb_angptl3_phenotype_time5$dialysis_time <- as.Date(ukbb_angptl3_phenotype_time5$dialysis_time)

##distingulish incident,prevalent disease
create_ckd_status <- function(data, 
                              ckd_time_var, 
                              assess_time_var = "Time_assess") {
  # Convert dates and calculate gap
  data[[ckd_time_var]] <- as.Date(data[[ckd_time_var]])
  data[[assess_time_var]] <- as.Date(data[[assess_time_var]])
  gap_var <- paste0(ckd_time_var, "_gap")
  
  # generate incident 和 prevalent colum
  incident_var <- paste0(ckd_time_var, "_incident")
  prevalent_var <- paste0(ckd_time_var, "_prevalent")
  
  data <- data %>%
    mutate(
      !!gap_var := .data[[assess_time_var]] - .data[[ckd_time_var]],
      
      # Prevalent cases (diagnosis before/at assessment)
      !!prevalent_var := case_when(
        .data[[gap_var]] >= 0 ~ 1L,
        .data[[gap_var]] < 0 ~ 0L,
        is.na(.data[[gap_var]]) ~ 0L
      ),
      
      # Incident cases (diagnosis after assessment)
      !!incident_var := case_when(
        .data[[gap_var]] >= 0 ~ 0L,
        .data[[gap_var]] < 0 ~ 1L, 
        is.na(.data[[gap_var]]) ~ 0L
      )
    )
  
  return(data)
}

# apply
ukbb_angptl3_phenotype_time5 <- create_ckd_status(ukbb_angptl3_phenotype_time5, "Angina_time")

#logistic model
model <- glm(dialysis_time_prevalent ~ Fasting_time + Ethnic + smoke + Age + Sex + deprivation_index + angptl3_olink, 
             family = binomial, 
             data = ukbb_angptl3_phenotype_time5)

model_summary <- summary(model)$coefficients

# calculate OR and 95%CI
or_results <- exp(cbind(OR = coef(model), confint(model)))
or_results


#generate incident disease duration
incident_duration <- function(data, 
                              time_var,
                              assess_var = "Time_assess",
                              last_contact_var = "date_last_personal_contact",
                              prevalent_var = NULL,
                              incident_var = NULL) {
  
  if (is.null(prevalent_var)) prevalent_var <- paste0(time_var, "_prevalent")
  if (is.null(incident_var)) incident_var <- paste0(time_var, "_incident")
  
  # 检查必要列是否存在
  required_cols <- c(time_var, assess_var, last_contact_var, prevalent_var, incident_var)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("缺失必要列:", paste(missing_cols, collapse = ", ")))
  }
  
  # 1. delete prevelant cases（prevalent=1）
  data_filtered <- data[data[[prevalent_var]] != 1, ]
  
  # 2. calculate duration 
  duration_var <- paste0(time_var, "_duration")
  data_filtered[[duration_var]] <- ifelse(
    data_filtered[[incident_var]] == 1,
    as.numeric(data_filtered[[time_var]] - data_filtered[[assess_var]]),
    as.numeric(data_filtered[[last_contact_var]] - data_filtered[[assess_var]])
  )
  
  return(data_filtered)
}

##apply
AKI_time_cox <- incident_duration(
  data = ukbb_angptl3_phenotype_time5,
  time_var = "AKI_time"  # 自动匹配dialysis_time_prevalent/incident
)

##Cox analysis
cox_model <- coxph(
  Surv(AKI_time_duration, AKI_time_incident) ~ Fasting_time + Ethnic + smoke + Age + Sex + deprivation_index + angptl3_olink , 
  data = AKI_time_cox
)

summary(cox_model)
coef_summary <- summary(cox_model)$coefficients
coef_summary[, "Pr(>|z|)"] <- format(coef_summary[, "Pr(>|z|)"], scientific = FALSE, digits = 20)
print(coef_summary)

# calculate OR and 95%CI
or_results <- exp(cbind(OR = coef(cox_model), confint(cox_model)))
or_results










