pkgbuild::has_build_tools(debug = TRUE)
install.packages("glmmTMB")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("lme4")
install.packages("emmeans")
install.packages("DHARMa")
install.packages("performance")
install.packages("ggeffects")
install.packages("statmod")


library(glmmTMB)
library(readxl)
library(dplyr)
library(emmeans)
library(ggplot2)
library(DHARMa)
library(performance)
library(statmod)
library(tidyr)
library(ggeffects)


#SCETION 1: PREPARING THE DATA
file_1 <- read_excel("C:/Users/trbag//Downloads/BioMaths/Hopefully final/TIR every day.xlsx") #loads data for our main analysis
file_1 <- filter(file_1, TrtGroup != "BPFiasp") 
#just make sure to remove all patients with BPFiasp as a treatment group otherwise the summary will say there are 60 patients without an age at enrollment 
summary(file_1)

Daily_vals <- data.frame(file_1) #create a data frame that contains all columns of data from file_1
file_2 <- read_excel("C:/Users/trbag/Downloads/BioMaths/Hopefully final/BPFIASP_TIR.xlsx")
file_2 <- filter(file_2, TrtGroup != "Control") #exclude all control patuents from our file
summary(file_2)
BP_Data <- data.frame(file_2)


#We use factor to specify that the strings in certain columns (diabetes treatment + study group) are actually elements of different groups not just strings of text 
#factor tells R that something is a categorical variable
#rlevel is used when we we want to compare elements in a certain column/factor 
#we use rlevel to manually specify what we want the reference level of a specific factor to be  
Daily_vals$StudyGroup <- factor(Daily_vals$StudyGroup)
Daily_vals$StudyGroup <- relevel(Daily_vals$StudyGroup, ref = "Continued non-closed-loop treatment")
#We use the Control group as the thing with which other elements in that factor (i.e the BP group) will be compared 
Daily_vals$DiabTrt <- factor(Daily_vals$DiabTrt)
Daily_vals$AgeGroup <- factor(Daily_vals$AgeGroup)
Daily_vals$StudyDayFactor <- factor(
  Daily_vals$StudyDay,
  levels = seq( #need to define levels of the Study Day factor as different, equally-spaced days of the study
    #seq fills in intermediate values between the first (min) and last (max) study days 
    min(Daily_vals$StudyDay, na.rm = TRUE), #na.rm = TRUE ignores any missing days 
    max(Daily_vals$StudyDay, na.rm = TRUE)
  ) #missing days still included in the sequence because AR(1) needs to calculate the correlation between consecutive days
)
#----------------------------------------------------------------------------------

#SECTION 2: FINDING THE RIGHT MODEL

#1. CHECK FOR OVERDISPERSION IN DATA
model_BINOM<- glmmTMB(
  Ratio ~ StudyGroup + DiabTrt + Age + (1 | PtID), #Ratio is specified as the response variable 
  # We define our fixed effects above and (1|PtID) is our patient-specific random effect
  data   = Daily_vals, #specify the dataframe to be used 
  weights = TotalReadings, ##Specify weights to define the total readings as the total number of trials
  family = binomial(link = "logit") #specify the distribution and the link function 
)
check_overdispersion(simulate_residuals(model_BINOM)) #Test for overedispersion in initial model


#2. FIRST MODEL ATTEMPT: standard BB distribution

model_ONE <- glmmTMB(
  Ratio ~ StudyGroup + DiabTrt + Age + (1 | PtID),
  data = Daily_vals,
  weights = TotalReadings,
  family = betabinomial(link = "logit") #change family to betabinomial 
)


#3. SECOND MODEL ATTEMPT: ADDING A TEMPORAL AUTOCORRELATION COEFFICIENT TO ACCOUNT FOR CORRELATION BETWEEN DIFFERENT DAYS OF READINGS FOR A PATIENT
model_TWO <- glmmTMB(
  Ratio ~ StudyGroup + DiabTrt + Age +  ar1(StudyDayFactor + 0 | PtID) + (1 | PtID), #AR(1)/Temporal autocorrelation coefficient using StudyDay as the time factor and PtID as the group factor
  data = Daily_vals,
  weights = TotalReadings,
  family = betabinomial(link = "logit")
)


#4. FINAL MODEL ATTEMPT: ALLOWING THE DISPERSION PARAMETER FROM THE BB TO VARY BETWEEN DIFFERENT STUDY GROUPS
model_FINAL <- glmmTMB(
  Ratio ~ StudyGroup + DiabTrt + Age +
    (1 | PtID) +
    ar1(StudyDayFactor + 0 | PtID), 
  data = Daily_vals,
  weights = TotalReadings,
  dispformula =  ~ StudyGroup, #Allow dispersion to vary between different study groups 
  family = betabinomial(link = "logit")
)


#5. COMPARE AIC VALUES TO DETERMINE THE BEST MODEL
AIC(model_BINOM,model_ONE,model_TWO,model_FINAL)

#-----------------------------------------------------------------------------------

#SECTION 3: HYPOTHESIS 1 TESTING

#6. OBTAIN RESULTS FROM FINAL MODEL AND CHECK COLINEARITY
check_collinearity(model_FINAL)
summary(model_FINAL)
emm <- emmeans(model_FINAL, ~ StudyGroup, type = "response") #Calculates the estimated mean TIR for each study group (BP vs Control) as a proportion by setting type="response"
summary(emm)

#7. Use Contrast to decide whether to accept our reject our H0 
contrast(emm, method = "pairwise") #Compares the estimated BP vs control groups and uses Wald test to determine whether the difference is statistically significant
confint(contrast(emm, method = "pairwise")) #Calculate the 95% confidence interval from the previous line (95% is the default for this function)


# -------------------------------------------------------------------------
#SECTION 4:TESTING HYPOTHESIS 2: PATIENTS OVER AGE OF 65


#NOW WE WANT TO LOOK AT TIR FOR PATIENTS UNDER + OVER 65
#We now want to compare the groups Undr65 and Ovr65 and we want to look at the proportion of out of range readings which exhibit "extreme hypoglycemia"
#Hypothesis: patients over the age of 65 will experience hypoglycemia more than patients under the age of 65 

Daily_vals$AgeGroup <- relevel(Daily_vals$AgeGroup, ref = "Undr65")
model_AGE <- glmmTMB(
  RatioLow ~ AgeGroup + StudyGroup + DiabTrt + ar1(StudyDayFactor + 0 | PtID) + (1 | PtID), #Change response variable to RatioLow (Readingsbelow54/totalreadings)
  #Use AgeGroup as the fixed effect instead of Age
  data   = Daily_vals,
  weights = TotalReadings,
  dispformula =  ~ StudyGroup,
  family = betabinomial(link = "logit")
)
summary(model_AGE)
emm_AGE <- emmeans(model_AGE, ~ AgeGroup, type = "response") #Calculates mean estimated proportion of readings below 54mg.dl for the Undr65 and Ovr65 age groups
summary(emm_AGE)

## ── Difference between groups with p-value ───────────────────
contrast(emm_AGE, method = "pairwise")


#-------------------------------------------------------------------------
#HYPOTHESIS: BPFIASP
BP_Data$TrtGroup <- factor(BP_Data$TrtGroup) 
BP_Data$TrtGroup <- relevel(BP_Data$TrtGroup, ref = "BP") #Manually set BP as our reference group for the TrtGroup factor 
#We use the Control group as the thing with which other elements in that factor (i.e the BP group) will be compared 
BP_Data$DiabTrt <- factor(BP_Data$DiabTrt)
BP_Data$DiabTrt <- factor(BP_Data$DiabTrt)
BP_Data$StudyDayFactor <- factor(
  BP_Data$StudyDay,
  levels = seq(
    min(BP_Data$StudyDay, na.rm = TRUE),
    max(BP_Data$StudyDay, na.rm = TRUE)
  )
)



model_FIASP <- glmmTMB(
  Ratio ~ TrtGroup + DiabTrt  + Age + ar1(StudyDayFactor + 0 | PtID) + (1 | PtID), #keep the response variable as the same Ratio from hypothesis 1
  #Use TrtGroup here instead of Study group because the study group for BP and BPFiasp patients is "Bionic Pancreas" but TrtGroup specifies BP and BPFiasp
  data   = BP_Data,
  weights = TotalReadings,
  dispformula =  ~ TrtGroup, #Allow dispersion to vary between treatment groups (equivalent to allowing it to vary between study groups in the other hypotheses)
  family = betabinomial(link = "logit")
)
summary(model_FIASP)

## ── Estimated TIR % per group ─────────────────────────────────
emm <- emmeans(model_FIASP, ~ TrtGroup, type = "response") #Calculates estimated mean TIR for each treatment group
summary(emm)

## ── Difference between groups with p-value ───────────────────
contrast(emm, method = "pairwise")


#-------------------------------------------------------------


#GRAPHS SECTION



#GRAPH 1: ESTIMATED MEAN TIR for Hypothesis 1 
    
      
      # Repeats the emmeans calculation from hypothesis 1 but keeps it on the log-odds scale
      emm <- emmeans(model_FINAL, ~ StudyGroup)
      
      #Takes the estimated TIR and the 95% confidence intervals from the emm summary
      emm_results <- as.data.frame(
        summary(
          emm,
          type = "response", #converts log-odds values to proportions 
          infer = c(TRUE, FALSE), #displays the confidence intervals but not the standard errors or degrees of freedom
          level = 0.95 #sets confidence level to 95%
        )
      )
      
      #Performs the same contrast calculation from hypothesis 1 BUT places results in a data frame 
      contrast_results <- as.data.frame(
        summary(
          contrast(emm, method = "pairwise"),
          type = "response" #results converted from log-odds to proportions
        )
      )
      
      #Sets positions where bars corresponding to each study group will be and converts mean TIR and upper + lower confidence intervals from proportions to percentages
      emm_results <- emm_results %>% 
        mutate( #Adds 4 columns to emm_results
          x_pos = seq_len(n()), #counts number of rows in emm_results data frame (2) so x_pos is the sequence (1,2)
          mean_TIR = prob * 100, 
          lower = asymp.LCL * 100,
          upper = asymp.UCL * 100
        )
      
      # Store and extract the name of each study group
      group_names <- as.character(emm_results$StudyGroup)
      
      control_group <- group_names[ #function returns the full name for the control group found in group_names
        grepl("Continued", group_names, ignore.case = TRUE) #Identifies which names in group_names contain "Continued"
      ]
      
      bp_group <- group_names[ #function returns the full name for the BP group found in group_names
        grepl("Bionic", group_names, ignore.case = TRUE)#Identifies which names in group_names contain "Bionic"
      ]
      
      bar_colours <- setNames( #Assigns a specific colour to each study group 
        c("#E89088", "#59C3C9"),
        c(control_group, bp_group)
      )
      
      ci_colours <- setNames( #Assigns a specific colour to the separate confidence intervals for each study group
        c("red3", "dodgerblue3"),
        c(control_group, bp_group)
      )
      
      #Sort confidence interval values across study groups from lowest to highest to put on y-axis
      ci_values <- sort(c(emm_results$lower, emm_results$upper))
      
      # Plot
      ggplot(
        emm_results,
        aes(
          x = x_pos,
          y = mean_TIR,
          fill = StudyGroup 
          
        )
      ) +
        
        # Specify size of bars for each study group 
        geom_col(
          width = 0.55,
          alpha = 0.9 
        ) +
        
        # Draw lines representing confidence intervals
        geom_errorbar(
          aes(
            ymin = lower,
            ymax = upper,
            colour = StudyGroup
          ),
          width = 0.15,
          linewidth = 0.9
        ) +
        
        #Draw dotted line from the lower confidence interval endpoints to the y-axis (in same colour as confidence interval lines)
        geom_segment(
          aes(
            x = 0.5,
            xend = x_pos,
            y = lower,
            yend = lower,
            colour = StudyGroup
          ),
          linetype = "dotted",
          linewidth = 0.8
        ) +
        
        # Draw dotted line from upper confidence interval endpoints to the y-axis (in same colour as confidence interval lines)
        geom_segment(
          aes(
            x = 0.5,
            xend = x_pos,
            y = upper,
            yend = upper,
            colour = StudyGroup
          ),
          linetype = "dotted",
          linewidth = 0.8
        ) +
        
        # Write the estimated mean TIR (from emmeans as a percentage) in the middle of each bar
        geom_text(
          aes(
            y = mean_TIR / 2,
            label = paste0(round(mean_TIR, 1), "%")
          ),
          colour = "white",
          size = 5,
          fontface = "bold"
        ) +
        
        # Add text showing the odds ratio and p-value from emmeans 
        annotate(
          "text",
          x = 1.5,
          y = 91,
          label = paste0(
            "Odds ratio = ",
            round(contrast_results$odds.ratio[1], 3),
            "\np < 0.0001"
          ),
          size = 4
        ) +
        
        # Specify the bar colours
        scale_fill_manual(
          values = bar_colours
        ) +
        
        # Specify the colours for the confidence interval solid and dotted lines
        scale_colour_manual(
          values = ci_colours
        ) +
        
        # Format x-axis: Specify length of x-axis (without showing values), label position 1 and 2 for each study group, set limits for x-axis and remove excess space 
        scale_x_continuous( #allows us to assign text to numerical positions on the x-axis
          breaks = emm_results$x_pos,
          labels = emm_results$StudyGroup,
          limits = c(0.5, 2.5),
          expand = c(0, 0) #Don't exclude space before or after limit values 
        ) +
        
        # Label y-axis with possible TIR percentages and add all confidence interval endpoints
        scale_y_continuous(
          breaks = sort(unique(c(
            seq(0, 100, 25),
            ci_values
          ))),
          labels = function(x) {
            paste0(round(x, 1), "%")
          },
          limits = c(0, 100)
        ) +
        
        labs( #Label the graph and all axes 
          title = "Estimated Mean TIR by Treatment Group",
          subtitle = "Model-adjusted estimates with 95% confidence intervals",
          x = "Treatment Group",
          y = "Estimated Mean TIR (%)"
        ) +
        
        guides( #Hides colour legends)
          fill = "none",
          colour = "none"
        ) +
        
        theme_minimal()

#----------------------------------------------------------------------------------
#GRAPH 2: PREDICTED TIR FOR EACH AGE IN EACH TREATMENT GROUP 
      
      
      # Predict TIR across age for both treatment groups
      pred_age <- ggpredict(
        model_FINAL,
        terms = c("Age [all]", "StudyGroup") #Based on model_FINAL, predicts daily TIR for each unique age
      )
      
      # Plot
      plot(pred_age) +
        scale_y_continuous(
          labels = scales::percent,
          limits = c(0, 1)
        ) +
        labs(
          title = "Predicted TIR by Age and Treatment Group",
          subtitle = "Fixed-effect predictions with 95% confidence intervals",
          x = "Age (years)",
          y = "Predicted TIR (%)",
          colour = "Treatment Group",
          fill = "Treatment Group"
        ) +
        theme_minimal()

#-------------------------------------------------------
#POSSIBLE GRAPH 3: OVERALL TIR PER PATIENT PER STUDY
      
      
      # Each patient's overall time in range, weighted by readings per day
      pt_tir <- Daily_vals %>%
        group_by(PtID, StudyGroup) %>%
        summarise(TIR = sum(Ratio * TotalReadings) / sum(TotalReadings),
                  .groups = "drop")
      
      group_means <- pt_tir %>%
        group_by(StudyGroup) %>%
        summarise(mean_TIR = mean(TIR))
      
      # Rank patients within each group so the dots form a sorted curve
      pt_tir <- pt_tir %>%
        group_by(StudyGroup) %>%
        arrange(TIR, .by_group = TRUE) %>%
        mutate(rank = row_number()) %>%
        ungroup()
      
      ggplot(pt_tir, aes(x = rank, y = TIR, colour = StudyGroup)) +
        geom_point(size = 2, alpha = 0.8) +
        geom_hline(data = group_means,
                   aes(yintercept = mean_TIR, colour = StudyGroup),
                   linetype = "dashed", linewidth = 0.8) +
        facet_wrap(~ StudyGroup, scales = "free_x") +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        labs(x = "Patients (ordered from lowest to highest time in range)",
             y = "Time in range",
             title = "Time in range for each patient by study group",
             caption = "Each dot is one patient; dashed line is the group mean") +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "none")
      
      ggsave("fig_patient_tir_ranked.png", width = 9, height = 5, dpi = 300)


#-----------------------------------------------------------------------------------

#GRAPH 4: SHOW OBSERVED VS PREDICTED OVERALL TIR THROUGHOUT STUDY
      library(ggplot2)
      library(dplyr)
      
      # 1. Prepare the data
      # Assumes model_FINAL was fitted using all rows of Daily_vals
      # in the same order
      
      plot_data <- Daily_vals
      
      # Model-fitted daily TIR (including random effects)
      plot_data$FittedTIR <- predict(
        model_FINAL,
        type = "response"
      ) * 100
      
      # Observed daily TIR
      plot_data$ActualTIR <- plot_data$Ratio * 100
      
      
      # 2. Calculate observed and fitted TIR for each patient
      
      patient_comparison <- plot_data %>%
        group_by(PtID, StudyGroup) %>%
        summarise(
          Observed = weighted.mean(
            ActualTIR,
            TotalReadings,
            na.rm = TRUE
          ),
          Fitted = weighted.mean(
            FittedTIR,
            TotalReadings,
            na.rm = TRUE
          ),
          .groups = "drop"
        )
      
      
      # 3. Create data for the perfect-agreement line
      
      agreement_line <- data.frame(
        Observed = c(0, 100),
        Fitted = c(0, 100),
        Line = "Actual TIR = Predicted TIR"
      )
      
      
      # 4. Create the graph
      
      ggplot(
        patient_comparison,
        aes(
          x = Observed,
          y = Fitted,
          colour = StudyGroup
        )
      ) +
        
        # Each point represents one patient
        geom_point(
          size = 2.5,
          alpha = 0.8
        ) +
        
        # Dashed line: Actual TIR = Predicted TIR
        geom_line(
          data = agreement_line,
          aes(
            x = Observed,
            y = Fitted,
            linetype = Line
          ),
          inherit.aes = FALSE,
          colour = "black",
          linewidth = 0.7
        ) +
        
        # Treatment group colours
        scale_colour_manual(
          name = "Treatment Group",
          values = c(
            "Continued non-closed-loop treatment" = "#E89088",
            "Bionic Pancreas" = "#59C3C9"
          )
        ) +
        
        # Add dashed line to legend
        scale_linetype_manual(
          name = NULL,
          values = c(
            "Actual TIR = Predicted TIR" = "dashed"
          )
        ) +
        
        # Set equal axes
        coord_equal(
          xlim = c(0, 100),
          ylim = c(0, 100)
        ) +
        
        # Titles and labels
        labs(
          title = "Observed vs Model-Fitted TIR by Patient",
          subtitle = "Each point represents one patient",
          x = "Observed TIR (%)",
          y = "Model-Fitted TIR (%)"
        ) +
        
        # Appearance
        theme_minimal(base_size = 13) +
        
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_line(
            colour = "grey95"
          )
        )
#-----------------------------------------------------------

#POSSIBLE GRAPH 5: ESTIMATED MEAN PROPORTION OF READINGS BELOW 54mg.dl
      library(ggplot2)
      library(dplyr)
      library(emmeans)
      
      # 1. Obtain estimated proportions
      
      emm_age <- emmeans(model_AGE, ~ AgeGroup)
      
      age_results <- as.data.frame(
        summary(emm_age, type = "response")
      ) %>%
        mutate(
          AgeGroup = factor(
            as.character(AgeGroup),
            levels = c("Undr65", "Ovr65")
          ),
          Mean = prob * 100,
          Lower = asymp.LCL * 100,
          Upper = asymp.UCL * 100
        ) %>%
        arrange(AgeGroup) %>%
        mutate(x_pos = row_number())
      
      print(age_results)
      
      # 2. Odds ratio and p-value
      
      age_contrast <- as.data.frame(
        summary(
          contrast(emm_age, method = "pairwise"),
          type = "response"
        )
      )
      
      # 3. Y-axis labels
      
      ci_ticks_age <- sort(c(
        age_results$Lower,
        age_results$Upper
      ))
      
      regular_ticks_age <- seq(0, 0.6, 0.1)
      
      # Remove regular ticks within 0.04 percentage points
      # of any confidence interval endpoint
      regular_ticks_age <- regular_ticks_age[
        sapply(regular_ticks_age, function(x) {
          all(abs(x - ci_ticks_age) >= 0.04)
        })
      ]
      
      age_ticks <- sort(unique(c(
        regular_ticks_age,
        ci_ticks_age
      )))
      
      # 4. Graph
      
      graph_age <- ggplot(
        age_results,
        aes(x = x_pos, y = Mean, fill = AgeGroup)
      ) +
        
        geom_col(
          width = 0.55,
          alpha = 0.9
        ) +
        
        # Coloured confidence interval error bars
        geom_errorbar(
          aes(
            ymin = Lower,
            ymax = Upper,
            colour = AgeGroup
          ),
          width = 0.15,
          linewidth = 0.9
        ) +
        
        # Lower CI dotted lines
        geom_segment(
          aes(
            x = 0.5,
            xend = x_pos,
            y = Lower,
            yend = Lower,
            colour = AgeGroup
          ),
          linetype = "dotted",
          linewidth = 0.7
        ) +
        
        # Upper CI dotted lines
        geom_segment(
          aes(
            x = 0.5,
            xend = x_pos,
            y = Upper,
            yend = Upper,
            colour = AgeGroup
          ),
          linetype = "dotted",
          linewidth = 0.7
        ) +
        
        # Percentages inside bars
        geom_text(
          aes(
            y = Mean / 2,
            label = sprintf("%.3f%%", Mean)
          ),
          colour = "white",
          fontface = "bold",
          size = 5
        ) +
        
        # Odds ratio and p-value
        annotate(
          "text",
          x = 1.5,
          y = 0.53,
          label = paste0(
            "Odds ratio = ",
            round(age_contrast$odds.ratio[1], 3),
            "\np = ",
            signif(age_contrast$p.value[1], 3)
          ),
          size = 4
        ) +
        
        scale_fill_manual(
          values = c(
            "Undr65" = "#E89088",
            "Ovr65" = "#59C3C9"
          )
        ) +
        
        scale_colour_manual(
          values = c(
            "Undr65" = "red3",
            "Ovr65" = "dodgerblue3"
          )
        ) +
        
        scale_x_continuous(
          breaks = age_results$x_pos,
          labels = as.character(age_results$AgeGroup),
          limits = c(0.5, 2.5)
        ) +
        
        scale_y_continuous(
          breaks = age_ticks,
          labels = function(x) sprintf("%.3f%%", x),
          minor_breaks = NULL,
          expand = expansion(mult = c(0, 0.03))
        ) +
        
        coord_cartesian(
          ylim = c(0, 0.6)
        ) +
        
        labs(
          title = "Estimated Severe Hypoglycaemia by Age Group",
          subtitle = "Model-adjusted estimates with 95% confidence intervals",
          x = "Age Group",
          y = "Estimated Readings Below 54 mg/dL (%)"
        ) +
        
        guides(fill = "none", colour = "none") +
        
        theme_minimal(base_size = 13) +
        
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 10),
          panel.grid.minor = element_blank()
        )
      
      print(graph_age)

#-----------------------------------------------------------


#GRAPH 6: BP vs BPFiasp
      library(ggplot2)
      library(dplyr)
      library(emmeans)
      
      # 1. Obtain model estimates
      
      emm_fiasp <- emmeans(model_FIASP, ~ TrtGroup)
      
      fiasp_results <- as.data.frame(
        summary(emm_fiasp, type = "response")
      ) %>%
        mutate(
          TrtGroup = factor(
            as.character(TrtGroup),
            levels = c("BP", "BPFiasp")
          ),
          Mean = prob * 100,
          Lower = asymp.LCL * 100,
          Upper = asymp.UCL * 100
        ) %>%
        arrange(TrtGroup) %>%
        mutate(x_pos = row_number())
      
      # 2. Odds ratio and p-value
      
      fiasp_contrast <- as.data.frame(
        summary(
          contrast(emm_fiasp, method = "pairwise"),
          type = "response"
        )
      )
      
      # 3. Create graph
      
      graph_fiasp <- ggplot(
        fiasp_results,
        aes(x = x_pos, y = Mean, fill = TrtGroup)
      ) +
        
        # Bars
        geom_col(
          width = 0.55,
          alpha = 0.9
        ) +
        
        # Coloured error bars
        geom_errorbar(
          aes(
            ymin = Lower,
            ymax = Upper,
            colour = TrtGroup
          ),
          width = 0.15,
          linewidth = 0.9
        ) +
        
        # Dotted lower CI lines
        geom_segment(
          aes(
            x = 0.5,
            xend = x_pos,
            y = Lower,
            yend = Lower,
            colour = TrtGroup
          ),
          linetype = "dotted",
          linewidth = 0.7
        ) +
        
        # Dotted upper CI lines
        geom_segment(
          aes(
            x = 0.5,
            xend = x_pos,
            y = Upper,
            yend = Upper,
            colour = TrtGroup
          ),
          linetype = "dotted",
          linewidth = 0.7
        ) +
        
        # Lower CI labels: below the top of each bar
        geom_text(
          aes(
            y = Lower - 3,
            label = paste0(
              "Lower CI: ",
              sprintf("%.1f%%", Lower)
            )
          ),
          colour = "black",
          size = 3.5
        ) +
        
        # Upper CI labels: above each bar
        geom_text(
          aes(
            y = Upper + 3,
            label = paste0(
              "Upper CI: ",
              sprintf("%.1f%%", Upper)
            )
          ),
          colour = "black",
          size = 3.5
        ) +
        
        # Estimated means inside bars
        geom_text(
          aes(
            y = Mean / 2,
            label = sprintf("%.1f%%", Mean)
          ),
          colour = "white",
          fontface = "bold",
          size = 5
        ) +
        
        # Odds ratio and p-value
        annotate(
          "text",
          x = 1.5,
          y = 92,
          label = paste0(
            "Odds ratio = ",
            round(fiasp_contrast$odds.ratio[1], 3),
            "\np = ",
            signif(fiasp_contrast$p.value[1], 3)
          ),
          size = 4
        ) +
        
        # Bar colours
        scale_fill_manual(
          values = c(
            "BP" = "#E89088",
            "BPFiasp" = "#59C3C9"
          )
        ) +
        
        # Error bar and dotted line colours
        scale_colour_manual(
          values = c(
            "BP" = "red3",
            "BPFiasp" = "dodgerblue3"
          )
        ) +
        
        # X-axis labels
        scale_x_continuous(
          breaks = fiasp_results$x_pos,
          labels = as.character(fiasp_results$TrtGroup),
          limits = c(0.5, 2.5)
        ) +
        
        # Regular y-axis ticks only
        scale_y_continuous(
          breaks = seq(0, 100, 25),
          labels = function(x) paste0(x, "%")
        ) +
        
        coord_cartesian(ylim = c(0, 100)) +
        
        labs(
          title = "Estimated Mean TIR: BP vs BPFiasp",
          subtitle = "Model-adjusted estimates with 95% confidence intervals",
          x = "Treatment Group",
          y = "Estimated Mean TIR (%)"
        ) +
        
        guides(
          fill = "none",
          colour = "none"
        ) +
        
        theme_minimal(base_size = 13) +
        
        theme(
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
      
      print(graph_fiasp)
      
      

#-------------------------------------

#Graph 7: Compare patients actual vs fitted daily TIR
      library(ggplot2)
      library(dplyr)
      library(tidyr)
      
      # 1. Prepare data
      # Assumes model_FINAL was fitted using all rows
      # of Daily_vals in the same order
      
      plot_data <- Daily_vals
      
      # Model-fitted daily TIR
      plot_data$FittedTIR <- predict(
        model_FINAL,
        type = "response"
      ) * 100
      
      # Observed daily TIR
      plot_data$ActualTIR <- plot_data$Ratio * 100
      
      # Convert study day to numeric
      plot_data$StudyDay <- as.numeric(
        as.character(plot_data$StudyDayFactor)
      )
      
      
      # 2. Randomly select 2 patients per treatment group
      
      set.seed(123)
      
      selected_patients <- plot_data %>%
        distinct(StudyGroup, PtID) %>%
        group_by(StudyGroup) %>%
        slice_sample(n = 2) %>%
        ungroup()
      
      print(selected_patients)
      
      
      # 3. Keep selected patients
      
      selected_data <- plot_data %>%
        semi_join(
          selected_patients,
          by = c("StudyGroup", "PtID")
        )
      
      
      # 4. Convert to long format
      
      trajectory_data <- selected_data %>%
        select(
          PtID,
          StudyGroup,
          StudyDay,
          ActualTIR,
          FittedTIR
        ) %>%
        pivot_longer(
          cols = c(ActualTIR, FittedTIR),
          names_to = "Type",
          values_to = "TIR"
        ) %>%
        mutate(
          Type = factor(
            Type,
            levels = c("ActualTIR", "FittedTIR"),
            labels = c("Observed", "Model-fitted")
          ),
          
          PatientLabel = paste0(
            ifelse(
              StudyGroup == "Bionic Pancreas",
              "BP",
              "Control"
            ),
            " - Patient ", PtID
          ),
          
          LineColour = case_when(
            Type == "Model-fitted" ~ "Model-fitted",
            StudyGroup == "Bionic Pancreas" ~ "BP observed",
            TRUE ~ "Control observed"
          )
        ) %>%
        arrange(PtID, Type, StudyDay)
      
      
      # 5. Create four-panel graph
      
      graph_trajectories <- ggplot(
        trajectory_data,
        aes(
          x = StudyDay,
          y = TIR,
          colour = LineColour,
          linetype = Type,
          group = interaction(PtID, Type)
        )
      ) +
        
        # Patient trajectories
        geom_line(
          linewidth = 0.7,
          alpha = 0.85
        ) +
        
        # Four separate patient panels
        facet_wrap(
          ~ PatientLabel,
          ncol = 2
        ) +
        
        # Line colours
        scale_colour_manual(
          name = NULL,
          breaks = c(
            "BP observed",
            "Control observed",
            "Model-fitted"
          ),
          values = c(
            "BP observed" = "#59C3C9",
            "Control observed" = "#E89088",
            "Model-fitted" = "black"
          )
        ) +
        
        # Line types
        scale_linetype_manual(
          values = c(
            "Observed" = "solid",
            "Model-fitted" = "dashed"
          ),
          guide = "none"
        ) +
        
        # One combined legend
        guides(
          colour = guide_legend(
            override.aes = list(
              linetype = c(
                "solid",
                "solid",
                "dashed"
              ),
              alpha = 1,
              linewidth = 0.9
            )
          )
        ) +
        
        # Y-axis
        scale_y_continuous(
          breaks = seq(0, 100, 20),
          labels = function(x) paste0(x, "%")
        ) +
        
        coord_cartesian(
          ylim = c(0, 100)
        ) +
        
        # Labels
        labs(
          title = "Observed vs Model-Fitted Daily TIR",
          subtitle = "Two randomly selected patients per treatment group",
          x = "Study Day",
          y = "Daily TIR (%)"
        ) +
        
        # Appearance
        theme_minimal(base_size = 12) +
        
        theme(
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()
        )
      
      print(graph_trajectories)
      
      
      # 6. Save high-resolution figure for Word
      
      ggsave(
        filename = "Figure4_Patient_TIR_Trajectories.png",
        plot = graph_trajectories,
        width = 11,
        height = 7,
        dpi = 300
      )
      



###-------------------------------------
#GRAPH 8: DHARMa
      library(DHARMa)
      
      set.seed(123)
      
      # Simulate 1000 datasets from your fitted model
      res_FINAL <- simulateResiduals(
        fittedModel = model_FINAL,
        n = 1000,
        refit = FALSE
      )
      
      # Standard diagnostic plots
      plot(res_FINAL)
      
      # Formal diagnostic tests
      testUniformity(res_FINAL)
      testDispersion(res_FINAL)
      testOutliers(res_FINAL)
