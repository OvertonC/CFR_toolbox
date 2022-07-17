# CFR risk ratio script

library(epitools)
library(dplyr)
library(ggplot2)
setwd("C:/Users/overt/OneDrive/Documents/Work/Manchester_postdoc/R/CareHomeCFR/output_files/EN_data")


#################### AGE AGE AGE #############################################################################
df_age_4 <- read.csv("output_over_65_28day_allbackwardsweekly4variable_lag.csv") %>% mutate(Age = "85+")
df_age_3 <- read.csv("output_over_65_28day_allbackwardsweekly3variable_lag.csv") %>% mutate(Age = "75 to 85")
df_age_2 <- read.csv("output_over_65_28day_allbackwardsweekly2variable_lag.csv") %>% mutate(Age = "65 to 75")
df_ages <- rbind(df_age_2,df_age_3,df_age_4)


output <- data.frame()
age_list <- unique(df_ages$Age)
for (age_1 in age_list) {
  for (age_2 in age_list){
    df1 <- df_ages %>% filter(Age == age_1)
    df2 <- df_ages %>% filter(Age == age_2)
    rr <- matrix(0,nrow(df_age_2),1)
    rrl <- matrix(0,nrow(df_age_2),1)
    rru <- matrix(0,nrow(df_age_2),1)
    for (i in 1:nrow(df_age_2)){
      if (round(df1$X7[i]) == 0 | round(df2$X7[i]) == 0){
        print("here")
      } else {
      rr[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[2]
      rrl[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[4]
      rru[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[6]
      }
    }
    
    output_temp <- as.data.frame(matrix(c(rr,rrl,rru),ncol = 3))
    output_temp$age_1 <- age_1
    output_temp$age_2 <- age_2
    output <- rbind(output,output_temp)
  }
}

output$age_comb <- paste0(output$age_2," relative to ",output$age_1)
output <- output %>% group_by(age_comb) %>% mutate(time = as.Date("2020-03-01")-1 + 7*seq(1,92,1)) %>% ungroup() %>% 
  filter(age_comb %in% c("65 to 75 relative to 75 to 85","65 to 75 relative to 85+","75 to 85 relative to 85+")) %>% filter(
    time < as.Date("2021-09-22")
  )

plt <- ggplot(data = output,aes(x = time)) + 
  geom_line(aes(y = V1),color = "red") +
  geom_hline(yintercept = 1) + 
  geom_ribbon(aes(ymin = V2, ymax = V3),alpha =0.1,fill = "red") +
  coord_cartesian(ylim = c(0,2)) + 
  facet_wrap(~age_comb) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  labs(x = "Time",y= "Risk ratio")
plt
ggsave(filename = "ages_risk_ratio_CFR.png",plot = plt, width = 12, height = 7)
################################################################################################





#################### CARE TYPE #############################################################################
df1 <- read.csv("output_over_65_28day_withoutbackwardsoverall1variable_lag.csv") %>% mutate(care = "Without",method = "overall")
df2 <- read.csv("output_over_65_28day_withbackwardsoverall1variable_lag.csv") %>% mutate(care = "With",method = "overall")
df3 <- read.csv("output_over_65_28day_withoutbackwardsweekly1variable_lag.csv") %>% mutate(care = "Without",method = "weekly")
df4 <- read.csv("output_over_65_28day_withbackwardsweekly1variable_lag.csv") %>% mutate(care = "With",method = "weekly")
df5 <- read.csv("output_over_65_28day_withoutbackwardsdaily1variable_lag.csv") %>% mutate(care = "Without",method = "daily")
df6 <- read.csv("output_over_65_28day_withbackwardsdaily1variable_lag.csv") %>% mutate(care = "With",method = "daily")
df_care_type <- rbind(df1,df2,df3,df4,df5,df6)
df_care_type$group <- paste0(df_care_type$care,"_",df_care_type$method)


output <- data.frame()
care_list <- unique(df_care_type$care)
method_list <- unique(df_care_type$method)
for (method_ in method_list){
  for (care_1 in care_list) {
    for (care_2 in care_list){
      df1 <- df_care_type %>% filter(care == care_1,method == method_)
      df2 <- df_care_type %>% filter(care == care_2,method == method_)
      rr <- matrix(0,nrow(df1),1)
      rrl <- matrix(0,nrow(df1),1)
      rru <- matrix(0,nrow(df1),1)
      for (i in 1:nrow(df1)){
        if (round(df1$X7[i]) == 0 | round(df2$X7[i]) == 0){
          print("here")
        } else {
          rr[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[2]
          rrl[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[4]
          rru[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[6]
        }
      }
      
      output_temp <- as.data.frame(matrix(c(rr,rrl,rru),ncol = 3))
      output_temp$care_1 <- care_1
      output_temp$care_2 <- care_2
      output_temp$method <- method_
      output <- rbind(output,output_temp)
    }
  }
}

output$care_comb <- paste0(output$care_2," relative to ",output$care_1, " - ",output$method)
output <- output %>% mutate(scale = case_when(
  method == "daily" ~ 1,
  method == "overall" ~ 1,
  TRUE ~ 7
), length_ = case_when(
  method == "daily" ~ 650,
  method == "overall" ~ 650,
  TRUE ~ 92
))
output <- output %>% group_by(care_comb) %>% mutate(time = as.Date("2020-03-01")-scale + scale*seq(1,n(),1)) %>% ungroup() %>% 
  filter(care_comb %in% c("With relative to Without - overall", "With relative to Without - weekly","With relative to Without - daily")) %>% filter(
    time < as.Date("2021-09-22")
  )

plt <- ggplot(data = output,aes(x = time)) + 
  geom_line(aes(y = V1),color = "red") +
  geom_hline(yintercept = 1) + 
  geom_ribbon(aes(ymin = V2, ymax = V3),alpha =0.1,fill = "red") +
  coord_cartesian(ylim = c(0,2)) + 
  facet_wrap(~care_comb) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  labs(x = "Time",y= "Risk ratio")
plt
ggsave(filename = "care_type_risk_ratio_CFR.png",plot = plt, width = 12, height = 7)
################################################################################################




#################### Regions #############################################################################
setwd("C:/Users/overt/OneDrive/Documents/Work/Manchester_postdoc/R/CareHomeCFR/output_files/region_data")
df1 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagWest Midlands_Regions.csv") %>% mutate(region = "West Midlands")
df2 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagYorkshire and Humber_Regions.csv") %>% mutate(region = "Yorkshire and Humber")
df3 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagEast Midlands_Regions.csv") %>% mutate(region = "East Midlands")
df4 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagEast of England_Regions.csv") %>% mutate(region = "East of England")
df5 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagLondon_Regions.csv") %>% mutate(region = "London")
df6 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagNorth East_Regions.csv") %>% mutate(region = "North East")
df7 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagNorth West_Regions.csv") %>% mutate(region = "North West")
df8 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagSouth East_Regions.csv") %>% mutate(region = "South East")
df9 <- read.csv("output_over_65_28day_allbackwardsweekly1variable_lagSouth West_Regions.csv") %>% mutate(region = "South West")
df10 <- read.csv("C:/Users/overt/OneDrive/Documents/Work/Manchester_postdoc/R/CareHomeCFR/output_files/EN_data/output_over_65_28day_allbackwardsweekly1variable_lag.csv") %>% mutate(region = "England")
df_region <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)


output <- data.frame()
region_list <- unique(df_region$region)
for (care_1 in region_list) {
  df1 <- df_region %>% filter(region == care_1)
  df2 <- df_region %>% filter(region == "England")
  rr <- matrix(0,nrow(df1),1)
  rrl <- matrix(0,nrow(df1),1)
  rru <- matrix(0,nrow(df1),1)
  for (i in 1:nrow(df1)){
    if (round(df1$X7[i]) == 0 | round(df2$X7[i]) == 0){
      print("here")
    } else {
      rr[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[2]
      rrl[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[4]
      rru[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[6]
    }
  }
  
  output_temp <- as.data.frame(matrix(c(rr,rrl,rru),ncol = 3))
  output_temp$region <- care_1
  output <- rbind(output,output_temp)
}


output$region_comb <- paste0("England"," relative to ",output$region)
# output <- output %>% mutate(scale = case_when(
#   method == "daily" ~ 1,
#   method == "overall" ~ 1,
#   TRUE ~ 7
# ), length_ = case_when(
#   method == "daily" ~ 650,
#   method == "overall" ~ 650,
#   TRUE ~ 92
# ))
output <- output %>% group_by(region_comb) %>% mutate(time = as.Date("2020-03-01")-7 + 7*seq(1,n(),1)) %>% ungroup() %>%
  filter(region != "England")  %>% filter(
    time < as.Date("2021-09-22")
  )
output$V1[is.infinite(output$V1)] <- NA

plt <- ggplot(data = output,aes(x = time)) + 
  geom_line(aes(y = V1),color = "red") +
  geom_hline(yintercept = 1) + 
  geom_ribbon(aes(ymin = V2, ymax = V3),alpha =0.1,fill = "red") +
  coord_cartesian(ylim = c(0,2)) + 
  facet_wrap(~region_comb) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  labs(x = "Time",y= "Risk ratio")
plt
ggsave(filename = "region_PHE_risk_ratio_CFR.png",plot = plt, width = 12, height = 7)
################################################################################################




#################### Regions #############################################################################
setwd("C:/Users/overt/OneDrive/Documents/Work/Manchester_postdoc/R/CareHomeCFR/output_files/region_data")
df1 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagWest Midlands_Regions.csv") %>% mutate(region = "West Midlands")
df2 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagYorkshire and Humber_Regions.csv") %>% mutate(region = "Yorkshire and Humber")
df3 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagEast Midlands_Regions.csv") %>% mutate(region = "East Midlands")
df4 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagEast of England_Regions.csv") %>% mutate(region = "East of England")
df5 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagLondon_Regions.csv") %>% mutate(region = "London")
df6 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagNorth East_Regions.csv") %>% mutate(region = "North East")
df7 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagNorth West_Regions.csv") %>% mutate(region = "North West")
df8 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagSouth East_Regions.csv") %>% mutate(region = "South East")
df9 <- read.csv("output_over_65_Confirmedday_allbackwardsweekly1variable_lagSouth West_Regions.csv") %>% mutate(region = "South West")
df10 <- read.csv("C:/Users/overt/OneDrive/Documents/Work/Manchester_postdoc/R/CareHomeCFR/output_files/EN_data/output_over_65_28day_allbackwardsweekly1variable_lag.csv") %>% mutate(region = "England")
df_region <- rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)


output <- data.frame()
region_list <- unique(df_region$region)
for (care_1 in region_list) {
  df1 <- df_region %>% filter(region == care_1)
  df2 <- df_region %>% filter(region == "England")
  rr <- matrix(0,nrow(df1),1)
  rrl <- matrix(0,nrow(df1),1)
  rru <- matrix(0,nrow(df1),1)
  for (i in 1:nrow(df1)){
    if (round(df1$X7[i]) == 0 | round(df2$X7[i]) == 0){
      print("here")
    } else {
      rr[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[2]
      rrl[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[4]
      rru[i] <- riskratio(matrix(c(df1$X6[i],df2$X6[i],round(df1$X7[i]),round(df2$X7[i])),nrow=2),conf.level = 0.95,rev= "c")$measure[6]
    }
  }
  
  output_temp <- as.data.frame(matrix(c(rr,rrl,rru),ncol = 3))
  output_temp$region <- care_1
  output <- rbind(output,output_temp)
}


output$region_comb <- paste0("England"," relative to ",output$region)
# output <- output %>% mutate(scale = case_when(
#   method == "daily" ~ 1,
#   method == "overall" ~ 1,
#   TRUE ~ 7
# ), length_ = case_when(
#   method == "daily" ~ 650,
#   method == "overall" ~ 650,
#   TRUE ~ 92
# ))
output <- output %>% group_by(region_comb) %>% mutate(time = as.Date("2020-03-01")-7 + 7*seq(1,n(),1)) %>% ungroup() %>%
  filter(region != "England")  %>% filter(
    time < as.Date("2021-09-22")
  )
output$V1[is.infinite(output$V1)] <- NA

plt <- ggplot(data = output,aes(x = time)) + 
  geom_line(aes(y = V1),color = "red") +
  geom_hline(yintercept = 1) + 
  geom_ribbon(aes(ymin = V2, ymax = V3),alpha =0.1,fill = "red") +
  coord_cartesian(ylim = c(0,2)) + 
  facet_wrap(~region_comb) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  labs(x = "Time",y= "Risk ratio")
plt
ggsave(filename = "region_CQC_risk_ratio_CFR.png",plot = plt, width = 12, height = 7)
################################################################################################



############# Method scoring ###################################################################################
setwd("C:/Users/overt/OneDrive/Documents/Work/Manchester_postdoc/R/CareHomeCFR/output_files/EN_data")
df1 <- (read.csv("output_over_65_28day_allforwardsdaily1variable_lag.csv") %>% mutate(method = "Forward"))[,c(1,4,9)]
names(df1) <- c("time","estimate","method")
df2 <- (read.csv("output_over_65_28day_allbackwardsdaily1variable_lag.csv") %>% mutate(method = "Backward"))[,c(1,4,9)]
names(df2) <- c("time","estimate","method")
df3 <- read.csv("cohort_CFR_current.csv") %>% mutate(method = "Cohort")
df3 <- df3 %>% mutate(time = seq(0,nrow(df3)-1,1)) %>% rename("estimate"="X0")
df_method <- rbind(df1,df2,df3)

df_method <- df_method %>% mutate(time = as.Date("2020-03-01") + time) %>% filter(time < as.Date("2021-09-22"))
df_cohort <- df_method %>% filter(method == "Cohort")
df_forward <- df_method %>% 
  left_join(df_cohort, by = "time") %>%
  filter(method.x == "Forward") %>% mutate(score = (estimate.y-estimate.x)/estimate.y) %>% filter(time <= max(time)-10)

df_backward <- df_method %>% 
  left_join(df_cohort, by = "time") %>%
  filter(method.x == "Backward") %>% mutate(score = (estimate.y-estimate.x)/estimate.y) %>% filter(time <= max(time)-10)

df_backward_shift <- df_method %>% mutate(time = time-10) %>% 
  left_join(df_cohort, by = "time") %>%
filter(method.x == "Backward") %>% mutate(score = (estimate.y-estimate.x)/estimate.y)

df_forward_mean <- df_forward %>% filter(!is.na(score),!is.infinite(score))
df_backward_mean <- df_backward %>% filter(!is.na(score),!is.infinite(score))
df_backward_shift_mean <- df_backward_shift %>% filter(!is.na(score),!is.infinite(score)) %>% mutate(method.x = "Backward (shifted)")

mean(df_forward_mean$score,na.rm=TRUE)
mean(df_backward_mean$score,na.rm=TRUE)
mean(df_backward_shift_mean$score,na.rm=TRUE)

df_plotting <- rbind(df_forward_mean,df_backward_mean,df_backward_shift_mean) %>% rename("Method" = "method.x")

plt <- ggplot(df_plotting,aes(x = time)) +
  geom_point(aes(y=score,color = Method)) + 
  geom_hline(aes(yintercept = 0)) + 
  labs(x = "Time",y= "Relative error")
plt
ggsave(filename = "relative_error_CFR.png",plt, width = 12, height = 7)
