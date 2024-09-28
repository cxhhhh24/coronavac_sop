rm(list = ls())

library(haven)
library(mgcv)
library(DescTools)
library(ggplot2)
library(patchwork)
library(reshape2)
library(readxl)
library(dplyr)
library(tidyverse)


dat <- read.csv("output/dat_cleaned.csv")
dat2 <- read.csv("output/dat_analysis.csv",fileEncoding="GBK")

#- formating date
dat <- dat[, -1]

# titer
for (i in c('V1_nt','SI1_nt', 'SI2_nt','I1_nt', 'I2_nt', 'I3_nt', 'I4_nt')) {
  dat[, i] <- str_replace_all(dat[, i], "<4", "2")
  dat[, i] <- str_replace_all(dat[, i], ">8192", "8192")
  dat[, i] <- as.numeric(dat[, i])
}
#date
for (i in c(names(dat)[str_detect(names(dat), "date")])) {
  dat[, i] <- as.Date(dat[, i], "%Y-%m-%d")
}


#- ajusting variables
dat$dose <- as.numeric(dat$dose)
dat$dose_interval <- as.numeric(dat$dose2_date-dat$dose1_date)
dat$pb14_dose2_interval <- dat$SI1_date-dat$dose2_date
dat$pb28_dose2_interval <- dat$SI2_date-dat$dose2_date
dat$age_group <- factor(dat$age_group, levels = c("young", "old"))
dat$sex <- factor(dat$sex, levels = c("f", "m"))

dat$race <- as.character(dat$race)
dat$race[dat$race == "白种人"] <- "white"
dat$race[dat$race == "多种族"] <- "multi_race"
dat$race[dat$race == "黑种人或非裔美洲人"] <- "black_or_african_american"
dat$race[dat$race == "美洲印第安人或阿拉斯加本地人"] <- "indian_or_native_of_Alaska"
dat$race[dat$race == "亚裔"] <- "asian"
dat$race <- factor(dat$race, levels = c("white", "black_or_african_american", 
                                        "indian_or_native_of_Alaska", "asian"  ,"multi_race"))
table(dat$race)

dat$bmi_binary[is.na(dat$bmi_binary)] <- "unk"
dat$bmi_binary <- factor(dat$bmi_binary, levels = c("0", "1","unk"))
table(dat$bmi_binary)


for (i in c("obesity", "malignant", "cardiovascular", 
            "chroniclung", "diabetes")) {
  str_replace_all(dat[, i], "Y", "1") -> dat[, i]
  str_replace_all(dat[, i], "N", "0") -> dat[, i]
  dat[, i] <- as.numeric(dat[, i])
}
dat$type_comorbidy[dat$obesity == 1] <- "obesity"
dat$type_comorbidy[dat$malignant == 1] <- "malignant"
dat$type_comorbidy[dat$cardiovascular == 1] <- "cardiovascular"
dat$type_comorbidy[dat$chroniclung == 1] <- "chroniclung"
dat$type_comorbidy[dat$diabetes == 1] <- "diabetes"
dat$type_comorbidy[is.na(dat$type_comorbidy) & dat$comorbidities == "Y"] <- "other"
dat$type_comorbidy <- factor(dat$type_comorbidy, levels = c("obesity", "malignant", "cardiovascular", "chroniclung", "diabetes", "other"))
dat$whe_comorbity <- factor(dat$whe_comorbity, levels = c(0, 1))

dat$whe_prior_infection <- as.character(dat$whe_prior_infection)
dat$whe_prior_infection[is.na(dat$whe_prior_infection)] <- "unk"
dat$whe_prior_infection <- factor(dat$whe_prior_infection, levels = c("Y", "N", "unk"))
table(dat$whe_prior_infection)

dat$whe_nt_available[!is.na(dat$V1_nt) | !is.na(dat$SI1_nt) | !is.na(dat$SI2_nt) | !is.na(dat$I1_nt) | 
                       !is.na(dat$I2_nt) | !is.na(dat$I3_nt) | !is.na(dat$I4_nt)] <- 1
dat$whe_nt_available[is.na(dat$whe_nt_available)] <- 0
table(dat$whe_nt_available)

dat$no_onset <- apply(!is.na(dat[, c("onset_date_1", "onset_date_2", "onset_date_3")]), 1, sum)


#---- selecting included participants
dat <- dat[dat$id %in% dat2$id, ]


########## PB 14 #####################
#- define cases
#- pb14 cases
follow_up_period <- 240
dat$pb14_case[dat$whe_CN_std_1 == "Y" & (dat$onset_date_1-dat$SI1_date >= 7) & (dat$onset_date_1-dat$SI1_date <= follow_up_period) & dat$treat_group == "Vaccine"] <- 1 #case
dat$pb14_case[dat$whe_CN_std_1 == "Y" & (dat$onset_date_1-dat$SI1_date >= 7) & (dat$onset_date_1-dat$SI1_date <= follow_up_period) & dat$treat_group == "Placebo" & dat$onset_date_1 < dat$sup_dose1_date] <- 1
# dat$pb14_case[dat$whe_CN_std_1 == "Y" & (dat$onset_date_1-dat$SI1_date < 7)  & (dat$onset_date_1 < dat$sup_dose1_date)] <- 2 #intercurrent-case
dat$pb14_case[dat$whe_CN_std_1 == "N" & is.na(dat$onset_date_1) == T] <- 2 # asym
dat$pb14_case[dat$whe_CN_std_1 == "N" & is.na(dat$onset_date_1) == F] <- 0 #noncase
# dat$pb14_case[dat$whe_CN_std_1 == "N"] <- 0 #noncase
dat$pb14_case[is.na(dat$pb14_case)] <- 0 #noncase
table(dat$pb14_case)


#-define asym cases 
tmp <- dat[, c("id", "SI2_nt", "I1_nt", "I2_nt","I3_nt", "I4_nt", "dose2_date",  "SI2_date", "I1_date", "I2_date", "I3_date", "I4_date")]
tmp$due_date <- tmp$dose2_date + follow_up_period
col <- c()
for (i in c(seq(8,12,1))) {
  new <- tmp[,i] - tmp[,ncol(tmp)]
  col <- cbind(col, new)
}

col <- as.data.frame(1* (col < 0))
col$time <- apply(col, 1, sum, na.rm =T)
max(col$time)
tmp <- tmp[, 1:(max(col$time)+1)]


asym_id1 <- c()
for (i in 2:max(col$time)) {
  for (j in (i+1):(max(col$time)+1)) {
    tmp1 <- tmp[, c(1,i,j)]
    tmp1$r <- tmp1[,3]/tmp1[,2]
    tmp2 <- unique(tmp1$id[tmp1$r >= 4])
    asym_id1 <- c(asym_id1, tmp2)
  }
}

asym_id2 <- c()
for (i in 2:max(col$time)) {
  for (j in (i+1):(max(col$time)+1)) {
    tmp1 <- tmp[, c(1,i,j)]
    tmp2 <- unique(tmp1$id[tmp1[,2] < 8 & tmp1[,3] >= 8])
    asym_id2 <- c(asym_id2, tmp2)
  }
}

asym_id <- c(asym_id1, asym_id2)
asym_id <- asym_id[!duplicated(asym_id)]
dat$pb14_case[dat$id %in% asym_id & dat$pb14_case != 1] <- 2
table(dat$pb14_case)



#--- estimate placebo average infection risk
dat_pb <- dat %>% filter(treat_group == "Placebo")
dat_pb <- dat_pb %>% filter(dose == 2) # 5838    # criteria 1
dat_pb <- dat_pb %>% filter(V1_nt == 2 & V1_igg < 0.8 & V1_igt < 1) # 1216   # criteria 2
dat_pb <- dat_pb %>% filter((dat_pb$dose_interval<=28  & dat_pb$dose_interval>=14)) # 1202   # criteria 3

dat_pb <- dat_pb[!is.na(dat_pb$SI1_nt) & !is.na(dat_pb$pb14_case), ] #1041 criteria 4
dat_pb <- dat_pb %>% filter(pb14_dose2_interval<=28 & pb14_dose2_interval>=14) #435 criteria 5
table(dat_pb$pb14_case)



#overall infection risk
placebo_risk_inf <- length(unique(dat_pb$id[dat_pb$pb14_case %in% c(1,2)]))/length(unique(dat_pb$id)) #0.4666667
#overall sym infection risk
placebo_risk_sym <- length(unique(dat_pb$id[dat_pb$pb14_case %in% c(1)]))/length(unique(dat_pb$id))  #0.4666667
#overall asym infection risk
placebo_risk_asym <- length(unique(dat_pb$id[dat_pb$pb14_case %in% c(2)]))/length(unique(dat_pb$id))  #0.4666667

#---estimate overall baseline exposure risk --- 
table(dat$pb14_case)
dat$pb14_case_overall <- dat$pb14_case
dat$pb14_case_overall[dat$pb14_case == 2] <- 1
base_risk_placebo <- glm(pb14_case_overall ~ age + sex + whe_white + bmi_binary + whe_comorbity, 
                         data = dat %>% filter(treat_group == "Placebo"), family = binomial())

summary(base_risk_placebo)

dat$pb14_case_sym <- dat$pb14_case
dat$pb14_case_sym[dat$pb14_case_sym == 2] <- NA
base_risk_placebo_sym <- glm(pb14_case_sym ~ age + sex + whe_white + bmi_binary + whe_comorbity, 
                             data = dat %>% filter(treat_group == "Placebo"), family = binomial())

summary(base_risk_placebo_sym)

dat$pb14_case_asym <- dat$pb14_case
dat$pb14_case_asym[dat$pb14_case_sym == 1] <- NA
dat$pb14_case_asym[dat$pb14_case_asym == 2] <- 1
base_risk_placebo_asym <- glm(pb14_case_asym ~ age + sex + whe_white + bmi_binary + whe_comorbity, 
                              data = dat %>% filter(treat_group == "Placebo"), family = binomial())

summary(base_risk_placebo_asym)


dat_vac <- dat %>% filter(treat_group == "Vaccine") # 6063 
dat_vac$base_risk_overall <- predict(base_risk_placebo, newdata = dat_vac, type = "response")
dat_vac$base_risk_sym <- predict(base_risk_placebo_sym, newdata = dat_vac, type = "response")
dat_vac$base_risk_asym <- predict(base_risk_placebo_asym, newdata = dat_vac, type = "response")




#--- define correlation populaiton  --- 
# criteria 1: receive 2 dose of coronaVac
# criteria 2: baseline seronegative
# criteria 3: dose 2 were not overwindow
dat_vac <- dat_vac %>% filter(dose == 2) # 6063    # criteria 1
table(dat_vac$V1_nt)
dat_vac <- dat_vac %>% filter(V1_nt == 2) # 1775   # criteria 2
dat_vac <- dat_vac %>% filter((dat_vac$dose_interval<=28  & dat_vac$dose_interval>=14)) # 1742   # criteria 3

#--- inverse probability sampling
dat_vac$interval_dose1_dose2 <-  dat_vac$dose2_date - dat_vac$dose1_date
dat_vac$whe_pb14titer <- 1*(!is.na(dat_vac$SI1_nt))
ips <- glm(whe_pb14titer ~ age_group + interval_dose1_dose2 + pb14_case,family = 'binomial', data = dat_vac)
dat_vac$ips <- predict(ips, newdata = dat_vac, type = "response")

#--- define correlation cohort  --- 
# criteria 4: have pb14 nt result
# criteria 5: pb14 were not overwindow (2 week window)
# criteria 6: remove intercurrent cases
dat_vac <- dat_vac[!is.na(dat_vac$SI1_nt) & !is.na(dat_vac$pb14_case), ] #1375 criteria 4
dat_vac <- dat_vac %>% filter(pb14_dose2_interval<=28 & pb14_dose2_interval>=14) #871 criteria 5
# dat_vac <- dat_vac[dat_vac$pb14_case != 2, ] #836 criteria 6
sum(is.na(dat_vac$SI1_nt))
table(dat_vac$pb14_case)



# distribution of NT  
as.data.frame.matrix(table(dat_vac$SI1_nt, dat_vac$pb14_case)) -> a
a$prop_inf <- round((a$`1`+a$`2`)/(a$`0`+a$`1`+a$`2`) * 100, 2)
a$prop_sym <- round(a$`1`/(a$`0`+a$`1`+a$`2`) * 100, 2)
a$prop_asym <- round(a$`2`/(a$`0`+a$`1`+a$`2`) * 100, 2)
a$titer <- row.names(a)
colnames(a)[1:3] <- c("noncase", "sym_case", "asym_case")
a <- a[, c("titer", "noncase", "sym_case", "asym_case", "prop_inf", "prop_sym", "prop_asym")]

# write.csv(a, "output/pb14_nt_dis.csv")



#--- grouping dataset based on clinicial endpoints
dat_sym <- dat_vac[dat_vac$pb14_case != 2, ]
table(dat_sym$pb14_case)

dat_asym <- dat_vac[dat_vac$pb14_case != 1, ]
dat_asym$pb14_case[dat_asym$pb14_case == 2] <- 1
table(dat_asym$pb14_case)

dat_inf <- dat_vac
dat_inf$pb14_case[dat_inf$pb14_case == 2] <- 1
table(dat_inf$pb14_case)


##### Sym ######
#--- estimating abosulte risk
dat_sym <- dat_sym[order(dat_sym$SI1_nt), ] 
gam_model <- gam(pb14_case ~ s(log(SI1_nt,2), bs="cr", k = 3) + base_risk_sym, data = dat_sym , family = "binomial", weights = 1/ips)
dat_sym$abs_risk <- predict(gam_model, newdata = dat_sym, type = "response")

#-- boostrap
set.seed(1234)
n=100
pred_sym <- c()
for (i in 1:n) {
  sample <- dat_sym[sample(1:nrow(dat_sym), nrow(dat_sym), replace = T), ]
  gam_model <- gam(pb14_case ~ s(log(SI1_nt,2), bs="cr", k = 3) + base_risk_sym, data = sample, family = "binomial", weights = 1/ips)
  tmp <- data.frame(no = i, SI1_nt = 2:8192, base_risk_sym = median(dat_sym$base_risk_sym, na.rm = T))
  tmp$abs_risk <- predict(gam_model, newdata = tmp, type = "response")
  pred_sym <- rbind(pred_sym, tmp)
}

pred_sym <- dcast(pred_sym, SI1_nt ~ no, value.var = "abs_risk")
pred_sym$abs_risk <- apply(pred_sym[, 2:n+1], 1, quantile, prob = 0.5, na.rm = T)
pred_sym$abs_risk_lci <- apply(pred_sym[, 2:n+1], 1, quantile, prob = 0.025, na.rm = T)
pred_sym$abs_risk_uci <- apply(pred_sym[, 2:n+1], 1, quantile, prob = 0.975, na.rm = T)
pred_sym <- pred_sym[, c("SI1_nt", "abs_risk", "abs_risk_lci", "abs_risk_uci")]


pred_sym$rel_risk <- pred_sym$abs_risk/placebo_risk_sym
pred_sym$rel_risk_lci <- pred_sym$abs_risk_lci/placebo_risk_sym
pred_sym$rel_risk_uci <- pred_sym$abs_risk_uci/placebo_risk_sym

pred_sym$ve <- 1-pred_sym$rel_risk
pred_sym$ve_lci <- 1-pred_sym$rel_risk_uci
pred_sym$ve_uci <- 1-pred_sym$rel_risk_lci
pred_sym$ve_lci[pred_sym$ve_lci < 0] <- 0


##### Asym ######
#--- estimating abosulte risk
dat_asym <- dat_asym[order(dat_asym$SI1_nt), ] 
gam_model <- gam(pb14_case ~ s(log(SI1_nt,2), bs="cr", k = 3) + base_risk_asym, data = dat_asym , family = "binomial", weights = 1/ips)
dat_asym$abs_risk <- predict(gam_model, newdata = dat_asym, type = "response")

#-- boostrap
set.seed(1234)
n=100
pred_asym <- c()
for (i in 1:n) {
  sample <- dat_asym[sample(1:nrow(dat_asym), nrow(dat_asym), replace = T), ]
  gam_model <- gam(pb14_case ~ s(log(SI1_nt,2), bs="cr", k = 3) + base_risk_asym, data = sample, family = "binomial", weights = 1/ips)
  tmp <- data.frame(no = i, SI1_nt = 2:8192, base_risk_asym = median(dat_asym$base_risk_asym, na.rm = T))
  tmp$abs_risk <- predict(gam_model, newdata = tmp, type = "response")
  pred_asym <- rbind(pred_asym, tmp)
}

pred_asym <- dcast(pred_asym, SI1_nt ~ no, value.var = "abs_risk")
pred_asym$abs_risk <- apply(pred_asym[, 2:n+1], 1, quantile, prob = 0.5, na.rm = T)
pred_asym$abs_risk_lci <- apply(pred_asym[, 2:n+1], 1, quantile, prob = 0.025, na.rm = T)
pred_asym$abs_risk_uci <- apply(pred_asym[, 2:n+1], 1, quantile, prob = 0.975, na.rm = T)
pred_asym <- pred_asym[, c("SI1_nt", "abs_risk", "abs_risk_lci", "abs_risk_uci")]


pred_asym$rel_risk <- pred_asym$abs_risk/placebo_risk_asym
pred_asym$rel_risk_lci <- pred_asym$abs_risk_lci/placebo_risk_asym
pred_asym$rel_risk_uci <- pred_asym$abs_risk_uci/placebo_risk_asym

pred_asym$ve <- 1-pred_asym$rel_risk
pred_asym$ve_lci <- 1-pred_asym$rel_risk_uci
pred_asym$ve_uci <- 1-pred_asym$rel_risk_lci
pred_asym$ve_lci[pred_asym$ve_lci < 0] <- 0


##### Inf ######
#--- estimating abosulte risk
dat_inf <- dat_inf[order(dat_inf$SI1_nt), ] 
gam_model <- gam(pb14_case ~ s(log(SI1_nt,2), bs="cr", k = 3) + base_risk_overall, data = dat_inf , family = "binomial", weights = 1/ips)
dat_inf$abs_risk <- predict(gam_model, newdata = dat_inf, type = "response")

#-- boostrap
set.seed(1234)
n=100
pred_inf <- c()
for (i in 1:n) {
  sample <- dat_inf[sample(1:nrow(dat_inf), nrow(dat_inf), replace = T), ]
  gam_model <- gam(pb14_case ~ s(log(SI1_nt,2), bs="cr", k = 3) + base_risk_overall, data = sample, family = "binomial", weights = 1/ips)
  tmp <- data.frame(no = i, SI1_nt = 2:8192, base_risk_overall = median(dat_inf$base_risk_overall, na.rm = T))
  tmp$abs_risk <- predict(gam_model, newdata = tmp, type = "response")
  pred_inf <- rbind(pred_inf, tmp)
}

pred_inf <- dcast(pred_inf, SI1_nt ~ no, value.var = "abs_risk")
pred_inf$abs_risk <- apply(pred_inf[, 2:n+1], 1, quantile, prob = 0.5, na.rm = T)
pred_inf$abs_risk_lci <- apply(pred_inf[, 2:n+1], 1, quantile, prob = 0.025, na.rm = T)
pred_inf$abs_risk_uci <- apply(pred_inf[, 2:n+1], 1, quantile, prob = 0.975, na.rm = T)
pred_inf <- pred_inf[, c("SI1_nt", "abs_risk", "abs_risk_lci", "abs_risk_uci")]


pred_inf$rel_risk <- pred_inf$abs_risk/placebo_risk_inf
pred_inf$rel_risk_lci <- pred_inf$abs_risk_lci/placebo_risk_inf
pred_inf$rel_risk_uci <- pred_inf$abs_risk_uci/placebo_risk_inf

pred_inf$ve <- 1-pred_inf$rel_risk
pred_inf$ve_lci <- 1-pred_inf$rel_risk_uci
pred_inf$ve_uci <- 1-pred_inf$rel_risk_lci
pred_inf$ve_lci[pred_inf$ve_lci < 0] <- 0

# combine data
pred_asym$type <- 1
pred_sym$type <- 2
pred_inf$type <- 3
pred <- rbind(pred_asym, pred_sym,pred_inf)
pred$type <- factor(pred$type)

dat_asym$type <- 1
dat_sym$type <- 2
dat_inf$type <- 3
dat_vac <- rbind(dat_asym, dat_sym, dat_inf)
dat_vac$type <- factor(dat_vac$type)




#--- plot 
theme(legend.position = "none",
      legend.justification = c(0,0),
      plot.margin = margin(1,1,0.5,0.5,unit="cm"),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.ticks = element_line(colour="black",size=0.4),
      axis.line = element_line(colour="black",size=1 ),
      axis.title.x = element_text(color="black",size=16,margin = margin(0, 0, 0, 0),face="bold", vjust = -3),
      axis.title.y = element_text(color="black",size=16,margin = margin(0, 0, 0, 0),face="bold", vjust = 3),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"), # adjust spacing between two faceted plots
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.background = element_rect(color = "white"),
      axis.text.x  = element_text(color="black",size=16,margin = margin(4, 0, 0, 0),face="bold"),
      axis.text.y  = element_text(color="black",size=16,margin = margin(0, 4, 0, 0),face="bold"),
      panel.grid = element_blank(),
      panel.border=element_blank(),
      strip.background = element_rect(colour = "black",fill = "#6996AA"),
      strip.text =  element_text(colour = "white",size=10,face="bold")) -> theme1



ggplot()+
  geom_point(data = dat_vac %>% filter(type != 1), aes(x = log(SI1_nt,2), y = abs_risk, color = type),alpha=0.25) + 
  geom_line(data = pred %>% filter(type != 1), aes(x = log(SI1_nt,2), y = abs_risk, color = type),size = 1) +
  geom_ribbon(data = pred %>% filter(type != 1), aes(x = log(SI1_nt,2), ymin = abs_risk_lci, ymax = abs_risk_uci, fill = type), alpha=0.15, show.legend = F) +
  scale_x_continuous(name = "Neutralization titer", expand = c(0,0),breaks = seq(1, 13, 2), 
                     labels = 2^seq(1, 13, 2), limits = c(-1, max(seq(1, 13, 2)))) +
  scale_y_continuous(name = "Absolute risk", limits = c(0,0.8), expand = c(0,0.001)) +
  scale_color_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"), labels = c("Symptomatic infections", "All infections"))+
  scale_fill_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"))+
  coord_cartesian(xlim = c(0,max(seq(1, 13, 2))), ylim = c(0,0.8), clip = "off") +
  geom_hline(yintercept =placebo_risk_sym, linetype = "dashed",  color = "black", size = 1) + 
  geom_hline(yintercept =placebo_risk_inf, linetype = "dashed",  color = "black", size = 1) + 
  annotate("text", x =9.5, y = 0.5, label = "Placebo overall infection risk") +
  annotate("text", x = 9, y = 0.43, label = "Placebo overall symptomatic risk") +
  theme(legend.position = "top",
        legend.justification = c(0,0),
        plot.margin = margin(1,1,0.5,0.5,unit="cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=1 ),
        axis.title.x = element_text(color="black",size=16,margin = margin(0, 0, 0, 0),face="bold",  vjust = -3, hjust = 0.5),
        axis.title.y = element_text(color="black",size=16,margin = margin(0, 0, 0, 0),face="bold", vjust = 3, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.8, "lines"), # adjust spacing between two faceted plots
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.background = element_rect(color = "white"),
        axis.text.x  = element_text(color="black",size=16,margin = margin(4, 0, 0, 0),face="bold"),
        axis.text.y  = element_text(color="black",size=16,margin = margin(0, 4, 0, 0),face="bold"),
        panel.grid = element_blank(),
        panel.border=element_blank(),
        strip.background = element_rect(colour = "black",fill = "#6996AA"),
        strip.text =  element_text(colour = "white",size=10,face="bold"),
        title = element_text(hjust = 0, size = 14, vjust = 0,face = "bold"))+
   ggtitle(label = "A. 14 days post the second dose")-> p1


ggplot()+
  # geom_density(data = dat_vac %>% filter(type == 1), aes(x = log(SI1_nt,2)), alpha=0.15, color = "#66CC99", fill = "#66CC99") +
  geom_density(data = dat_vac %>% filter(type == 2), aes(x = log(SI1_nt,2)), alpha=0.15,color = "#E1ADA7", fill = "#E1ADA7") +
  geom_density(data = dat_vac %>% filter(type == 3), aes(x = log(SI1_nt,2)), alpha=0.15,color = "#1f78b4", fill = "#1f78b4") +
  geom_line(data = pred %>% filter(type != 1), aes(x = log(SI1_nt,2), y = rel_risk/4, color = type),size = 1) +
  geom_ribbon(data =  pred %>% filter(type != 1), aes(x = log(SI1_nt,2), ymin = rel_risk_lci/4, ymax = rel_risk_uci/4, fill = type), alpha=0.15) +
  scale_x_continuous(name = "Neutralization titer", expand = c(0,0),breaks = seq(1, 13, 2), 
                     labels = 2^seq(1, 13, 2), limits = c(0, max(seq(1, 13, 2)))) +
  scale_y_continuous(name = "Relative risk", expand = c(0,0), breaks = c(0, 0.4,  0.8,  1.2, 1.6)/4,limits = c(0,1.6)/4,
                     labels =c("0", "0.4","0.8" ,"1.2", "1.6"),
                     sec.axis = sec_axis(trans = ~.*1, name = "Density\n ", breaks = c(0, 0.1, 0.2, 0.3, 0.4))) +
  scale_color_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"), labels = c("Symptomatic infections", "All infections"))+
  scale_fill_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"))+
  coord_cartesian(xlim = c(0,max(seq(1, 13, 2))), clip = "off") +
  theme1-> p2


ggplot()+
  # geom_density(data = dat_vac %>% filter(type == 1), aes(x = log(SI1_nt,2)), alpha=0.15, color = "#66CC99", fill = "#66CC99") +
  geom_density(data = dat_vac %>% filter(type == 2), aes(x = log(SI1_nt,2)), alpha=0.15,color = "#E1ADA7", fill = "#E1ADA7") + 
  geom_density(data = dat_vac %>% filter(type == 3), aes(x = log(SI1_nt,2)), alpha=0.15,color = "#1f78b4", fill = "#1f78b4") + 
  geom_line(data = pred %>% filter(type != 1), aes(x = log(SI1_nt,2), y = ve/3.3, color = type),size = 1) +
  geom_ribbon(data = pred %>% filter(type != 1), aes(x = log(SI1_nt,2), ymin = ve_lci/3.3, ymax = ve_uci/3.3, fill = type), alpha=0.15) +
  scale_x_continuous(expand = c(0,-0.5), name = "Neutralization titer",breaks = seq(1, 13, 2),
                     labels = 2^seq(1, 13, 2), limits = c(1, max(seq(1, 13, 2)))) +
  scale_y_continuous(name = "Vaccine efficacy (%)", expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75, 1.01)/3.3,limits = c(0,1.01)/3.3,
                     labels =c("0", "25", "50", "75","100"),
                     sec.axis = sec_axis(trans = ~.*1, name = "Density\n ", breaks = c(0, 0.1, 0.2, 0.3))) +
  scale_color_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"), labels = c("Symptomatic infections", "All infections"))+
  scale_fill_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"))+
  coord_cartesian(xlim = c(0,max(seq(1, 13, 2))+1), clip = "off") +
  annotate("segment", x = log(0,2), xend = log(76,2), y = 0.5/3.3, yend = 0.5/3.3, size=1, color="grey", linetype = 2)+
  annotate("segment", x = log(28,2), xend = log(28,2), y = 0, yend = 0.5/3.3, size=1, color="grey", linetype = 2)+
  annotate("segment", x = log(76,2), xend = log(76,2), y = 0, yend = 0.5/3.3, size=1, color="grey", linetype = 2)+
  annotate("text", x = 4.3, y = 0.015, label = "28", size = 5) +
  annotate("text", x = 6.8, y = 0.015, label = "76", size = 5) +
  theme1 -> p3




########## PB 28 #####################
#- define cases
#- pb28 cases
# follow_up_period <- 365
dat$pb28_case[dat$whe_CN_std_1 == "Y" & (dat$onset_date_1-dat$SI2_date >= 7) & (dat$onset_date_1-dat$SI2_date <= follow_up_period) & dat$treat_group == "Vaccine"] <- 1 #case
dat$pb28_case[dat$whe_CN_std_1 == "Y" & (dat$onset_date_1-dat$SI2_date >= 7) & (dat$onset_date_1-dat$SI2_date <= follow_up_period) & dat$treat_group == "Placebo" & dat$onset_date_1 < dat$sup_dose1_date] <- 1
# dat$pb28_case[dat$whe_CN_std_1 == "Y" & (dat$onset_date_1-dat$SI1_date < 7)  & (dat$onset_date_1 < dat$sup_dose1_date)] <- 2 #intercurrent-case
dat$pb28_case[dat$whe_CN_std_1 == "N" & is.na(dat$onset_date_1) == T] <- 2 # asym
dat$pb28_case[dat$whe_CN_std_1 == "N" & is.na(dat$onset_date_1) == F] <- 0 #noncase
# dat$pb28_case[dat$whe_CN_std_1 == "N"] <- 0 #noncase
dat$pb28_case[is.na(dat$pb28_case)] <- 0 #noncase
table(dat$pb28_case)


#-define asym cases 
tmp <- dat[, c("id", "SI2_nt", "I1_nt", "I2_nt","I3_nt", "I4_nt", "dose2_date",  "SI2_date", "I1_date", "I2_date", "I3_date", "I4_date")]
tmp$due_date <- tmp$dose2_date + follow_up_period
col <- c()
for (i in c(seq(8,12,1))) {
  new <- tmp[,i] - tmp[,ncol(tmp)]
  col <- cbind(col, new)
}


col <- as.data.frame(1* (col < 0))
col$time <- apply(col, 1, sum, na.rm =T)
max(col$time)
tmp <- tmp[, 1:(max(col$time)+1)]

asym_id1 <- c()
for (i in 2:max(col$time)) {
  for (j in (i+1):(max(col$time)+1)) {
    tmp1 <- tmp[, c(1,i,j)]
    tmp1$r <- tmp1[,3]/tmp1[,2]
    tmp2 <- unique(tmp1$id[tmp1$r >= 4])
    asym_id1 <- c(asym_id1, tmp2)
  }
}

asym_id2 <- c()
for (i in 2:max(col$time)) {
  for (j in (i+1):(max(col$time)+1)) {
    tmp1 <- tmp[, c(1,i,j)]
    tmp2 <- unique(tmp1$id[tmp1[,2] < 8 & tmp1[,3] >= 8])
    asym_id2 <- c(asym_id2, tmp2)
  }
}

asym_id <- c(asym_id1, asym_id2)
asym_id <- asym_id[!duplicated(asym_id)]
dat$pb28_case[dat$id %in% asym_id & dat$pb28_case != 1] <- 2
table(dat$pb28_case)


#--- estimate placebo average risk
dat_pb <- dat %>% filter(treat_group == "Placebo")
dat_pb <- dat_pb %>% filter(dose == 2) # 5836    # criteria 1
dat_pb <- dat_pb %>% filter(V1_nt == 2 & V1_igg < 0.8 & V1_igt < 1) # 1293   # criteria 2
dat_pb <- dat_pb %>% filter((dat_pb$dose_interval<=28  & dat_pb$dose_interval>=14)) # 1278   # criteria 3

dat_pb <- dat_pb[!is.na(dat_pb$SI2_nt) & !is.na(dat_pb$pb28_case), ] #1011 criteria 4
dat_pb <- dat_pb %>% filter(pb28_dose2_interval<=42 & pb28_dose2_interval>=28) #436 criteria 5
table(dat_pb$pb28_case)

#overall infection risk
placebo_risk_inf <- length(unique(dat_pb$id[dat_pb$pb28_case %in% c(1,2)]))/length(unique(dat_pb$id)) #0.4628297
#overall sym infection risk
placebo_risk_sym <- length(unique(dat_pb$id[dat_pb$pb28_case %in% c(1)]))/length(unique(dat_pb$id)) #0.4628297


#--- estimate overall baseline exposure risk --- 
table(dat$pb28_case)
dat$pb28_case_overall <- dat$pb28_case
dat$pb28_case_overall[dat$pb28_case == 2] <- 1
base_risk_placebo <- glm(pb28_case_overall ~ age + sex + whe_white + bmi_binary + whe_comorbity, 
                         data = dat %>% filter(treat_group == "Placebo"), family = binomial())


dat$pb28_case_sym <- dat$pb28_case
dat$pb28_case_sym[dat$pb28_case_sym == 2] <- NA
base_risk_placebo_sym <- glm(pb28_case_sym ~ age + sex + whe_white + bmi_binary + whe_comorbity, 
                             data = dat %>% filter(treat_group == "Placebo"), family = binomial())

dat$pb28_case_asym <- dat$pb28_case
dat$pb28_case_asym[dat$pb28_case_sym == 1] <- NA
dat$pb28_case_asym[dat$pb28_case_asym == 2] <- 1
base_risk_placebo_asym <- glm(pb28_case_asym ~ age + sex + whe_white + bmi_binary + whe_comorbity, 
                              data = dat %>% filter(treat_group == "Placebo"), family = binomial())

dat_vac <- dat %>% filter(treat_group == "Vaccine") # 6063 
dat_vac$base_risk_overall <- predict(base_risk_placebo, newdata = dat_vac, type = "response")
dat_vac$base_risk_sym <- predict(base_risk_placebo_sym, newdata = dat_vac, type = "response")
dat_vac$base_risk_asym <- predict(base_risk_placebo_asym, newdata = dat_vac, type = "response")






#--- define correlation populaiton  --- 
# criteria 1: receive 2 dose of coronaVac
# criteria 2: baseline seronegative
# criteria 3: dose 2 were not overwindow
dat_vac <- dat_vac %>% filter(dose == 2) # 6021    # criteria 1
dat_vac <- dat_vac %>% filter(V1_nt == 2) # 1766   # criteria 2
dat_vac <- dat_vac %>% filter((dat_vac$dose_interval<=28  & dat_vac$dose_interval>=14)) # 1743  # criteria 3
dat_vac$interval_dose1_dose2 <-  dat_vac$dose2_date - dat_vac$dose1_date


#--- inverse probability sampling
dat_vac$interval_dose1_dose2 <-  dat_vac$dose2_date - dat_vac$dose1_date
dat_vac$whe_pb28titer <- 1*(!is.na(dat_vac$SI2_nt))
ips <- glm(whe_pb28titer ~ age_group + interval_dose1_dose2 + pb28_case,family = 'binomial', data = dat_vac)
dat_vac$ips <- predict(ips, newdata = dat_vac, type = "response")


#--- define correlation cohort  --- 
# criteria 4: have pb28 nt result
# criteria 5: pb28 were not overwindow (2 week window)
# criteria 6: remove intercurrent cases
dat_vac <- dat_vac[!is.na(dat_vac$SI2_nt) & !is.na(dat_vac$pb28_case), ] #1279 criteria 4
dat_vac <- dat_vac %>% filter(pb28_dose2_interval<=42  & pb28_dose2_interval>=28) #799 criteria 5
# dat_vac <- dat_vac[dat_vac$pb28_case != 2, ] #714 criteria 6
sum(is.na(dat_vac$SI2_nt))
table(dat_vac$pb28_case)


#--- grouping dataset based on clinicial endpoints
dat_sym <- dat_vac[dat_vac$pb28_case != 2, ]

#sensitivity
# dat_sym <- dat_vac
# dat_sym$pb28_case[dat_sym$pb28_case == 2] <- 0
# table(dat_sym$pb28_case)


dat_asym <- dat_vac[dat_vac$pb28_case != 1, ]
dat_asym$pb28_case[dat_asym$pb28_case == 2] <- 1
table(dat_asym$pb28_case)

dat_inf <- dat_vac
dat_inf$pb28_case[dat_inf$pb28_case == 2] <- 1
table(dat_inf$pb28_case)


##### Sym ######
#--- estimating abosulte risk
dat_sym <- dat_sym[order(dat_sym$SI2_nt), ] 
gam_model <- gam(pb28_case ~ s(log(SI2_nt,2), bs="cr", k = 3) + base_risk_sym, data = dat_sym, family = "binomial", weights = 1/ips)
dat_sym$abs_risk <- predict(gam_model, newdata = dat_sym, type = "response")

#-- boostrap
set.seed(1234)
n=100
pred_sym <- c()
for (i in 1:n) {
  sample <- dat_sym[sample(1:nrow(dat_sym), nrow(dat_sym), replace = T), ]
  gam_model <- gam(pb28_case ~ s(log(SI2_nt,2), bs="cr", k = 3) + base_risk_sym, data = sample, family = "binomial", weights = 1/ips)
  tmp <- data.frame(no = i, SI2_nt = 2:8192, base_risk_sym = median(dat_sym$base_risk_sym, na.rm = T))
  tmp$abs_risk <- predict(gam_model, newdata = tmp, type = "response")
  pred_sym <- rbind(pred_sym, tmp)
}

pred_sym <- dcast(pred_sym, SI2_nt ~ no, value.var = "abs_risk")
pred_sym$abs_risk <- apply(pred_sym[, 2:n+1], 1, quantile, prob = 0.5, na.rm = T)
pred_sym$abs_risk_lci <- apply(pred_sym[, 2:n+1], 1, quantile, prob = 0.025, na.rm = T)
pred_sym$abs_risk_uci <- apply(pred_sym[, 2:n+1], 1, quantile, prob = 0.975, na.rm = T)
pred_sym <- pred_sym[, c("SI2_nt", "abs_risk", "abs_risk_lci", "abs_risk_uci")]


pred_sym$rel_risk <- pred_sym$abs_risk/placebo_risk_sym
pred_sym$rel_risk_lci <- pred_sym$abs_risk_lci/placebo_risk_sym
pred_sym$rel_risk_uci <- pred_sym$abs_risk_uci/placebo_risk_sym

pred_sym$ve <- 1-pred_sym$rel_risk
pred_sym$ve_lci <- 1-pred_sym$rel_risk_uci
pred_sym$ve_uci <- 1-pred_sym$rel_risk_lci

pred_sym$ve_lci[pred_sym$ve_lci < 0] <- 0
# write.csv(pred_sym, "pred_sym.csv")



##### Asym ######
#--- estimating abosulte risk
dat_asym <- dat_asym[order(dat_asym$SI2_nt), ] 
gam_model <- gam(pb28_case ~ s(log(SI2_nt,2), bs="cr", k = 3) + base_risk_asym, data = dat_asym, family = "binomial", weights = 1/ips)
dat_asym$abs_risk <- predict(gam_model, newdata = dat_asym, type = "response")

#-- boostrap
set.seed(1234)
n=100
pred_asym <- c()
for (i in 1:n) {
  sample <- dat_asym[sample(1:nrow(dat_asym), nrow(dat_asym), replace = T), ]
  gam_model <- gam(pb28_case ~ s(log(SI2_nt,2), bs="cr", k = 3) + base_risk_asym, data = sample, family = "binomial", weights = 1/ips)
  tmp <- data.frame(no = i, SI2_nt = 2:8192, base_risk_asym = median(dat_asym$base_risk_asym, na.rm = T))
  tmp$abs_risk <- predict(gam_model, newdata = tmp, type = "response")
  pred_asym <- rbind(pred_asym, tmp)
}

pred_asym <- dcast(pred_asym, SI2_nt ~ no, value.var = "abs_risk")
pred_asym$abs_risk <- apply(pred_asym[, 2:n+1], 1, quantile, prob = 0.5, na.rm = T)
pred_asym$abs_risk_lci <- apply(pred_asym[, 2:n+1], 1, quantile, prob = 0.025, na.rm = T)
pred_asym$abs_risk_uci <- apply(pred_asym[, 2:n+1], 1, quantile, prob = 0.975, na.rm = T)
pred_asym <- pred_asym[, c("SI2_nt", "abs_risk", "abs_risk_lci", "abs_risk_uci")]



pred_asym$rel_risk <- pred_asym$abs_risk/placebo_risk_asym
pred_asym$rel_risk_lci <- pred_asym$abs_risk_lci/placebo_risk_asym
pred_asym$rel_risk_uci <- pred_asym$abs_risk_uci/placebo_risk_asym

pred_asym$ve <- 1-pred_asym$rel_risk
pred_asym$ve_lci <- 1-pred_asym$rel_risk_uci
pred_asym$ve_uci <- 1-pred_asym$rel_risk_lci

pred_asym$ve_lci[pred_asym$ve_lci < 0] <- 0
# write.csv(pred_asym, "pred_asym.csv")


##### Inf ######
#--- estimating abosulte risk
dat_inf <- dat_inf[order(dat_inf$SI2_nt), ] 
gam_model <- gam(pb28_case ~ s(log(SI2_nt,2), bs="cr", k = 3) + base_risk_overall, data = dat_inf, family = "binomial", weights = 1/ips)
dat_inf$abs_risk <- predict(gam_model, newdata = dat_inf, type = "response")

#-- boostrap
set.seed(1234)
n=100
pred_inf <- c()
for (i in 1:n) {
  sample <- dat_inf[sample(1:nrow(dat_inf), nrow(dat_inf), replace = T), ]
  gam_model <- gam(pb28_case ~ s(log(SI2_nt,2), bs="cr", k = 3) + base_risk_overall, data = sample, family = "binomial", weights = 1/ips)
  tmp <- data.frame(no = i, SI2_nt = 2:8192, base_risk_overall = median(dat_inf$base_risk_overall, na.rm = T))
  tmp$abs_risk <- predict(gam_model, newdata = tmp, type = "response")
  pred_inf <- rbind(pred_inf, tmp)
}

pred_inf <- dcast(pred_inf, SI2_nt ~ no, value.var = "abs_risk")
pred_inf$abs_risk <- apply(pred_inf[, 2:n+1], 1, quantile, prob = 0.5, na.rm = T)
pred_inf$abs_risk_lci <- apply(pred_inf[, 2:n+1], 1, quantile, prob = 0.025, na.rm = T)
pred_inf$abs_risk_uci <- apply(pred_inf[, 2:n+1], 1, quantile, prob = 0.975, na.rm = T)
pred_inf <- pred_inf[, c("SI2_nt", "abs_risk", "abs_risk_lci", "abs_risk_uci")]



pred_inf$rel_risk <- pred_inf$abs_risk/placebo_risk_inf
pred_inf$rel_risk_lci <- pred_inf$abs_risk_lci/placebo_risk_inf
pred_inf$rel_risk_uci <- pred_inf$abs_risk_uci/placebo_risk_inf

pred_inf$ve <- 1-pred_inf$rel_risk
pred_inf$ve_lci <- 1-pred_inf$rel_risk_uci
pred_inf$ve_uci <- 1-pred_inf$rel_risk_lci

pred_inf$ve_lci[pred_inf$ve_lci < 0] <- 0
# write.csv(pred_inf, "pred_inf.csv")



# combine data
pred_asym$type <- 1
pred_sym$type <- 2
pred_inf$type <- 3
pred <- rbind(pred_asym, pred_sym,pred_inf)
pred$type <- factor(pred$type)

dat_asym$type <- 1
dat_sym$type <- 2
dat_inf$type <- 3
dat_vac <- rbind(dat_asym, dat_sym, dat_inf)
dat_vac$type <- factor(dat_vac$type)



#--- plot 
theme(legend.position = "none",
      legend.justification = c(0,0),
      plot.margin = margin(1,1,0.5,0.5,unit="cm"),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.ticks = element_line(colour="black",size=0.4),
      axis.line = element_line(colour="black",size=1 ),
      axis.title.x = element_text(color="black",size=16,face="bold", vjust = -3, hjust = 0.5),
      axis.title.y = element_text(color="black",size=16,face="bold", vjust = 3, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"), # adjust spacing between two faceted plots
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.background = element_rect(color = "white"),
      axis.text.x  = element_text(color="black",size=16,margin = margin(4, 0, 0, 0),face="bold"),
      axis.text.y  = element_text(color="black",size=16,margin = margin(0, 4, 0, 0),face="bold"),
      panel.grid = element_blank(),
      panel.border=element_blank(),
      strip.background = element_rect(colour = "black",fill = "#6996AA"),
      strip.text =  element_text(colour = "white",size=10,face="bold"),
      title = element_text(hjust = 0, size = 14, vjust = 0,face = "bold")) -> theme1

ggplot()+
  geom_point(data = dat_vac %>% filter(type != 1), aes(x = log(SI2_nt,2), y = abs_risk, color = type),alpha=0.25) + 
  geom_line(data = pred %>% filter(type != 1), aes(x = log(SI2_nt,2), y = abs_risk, color = type),size = 1) +
  geom_ribbon(data = pred %>% filter(type != 1), aes(x = log(SI2_nt,2), ymin = abs_risk_lci, ymax = abs_risk_uci, fill = type), 
              alpha=0.15, show.legend = F) +
  scale_x_continuous(name = "Neutralization titer", expand = c(0,0),breaks = seq(1, 13, 2), 
                     labels = 2^seq(1, 13, 2), limits = c(-15, max(seq(1, 13, 2))+5)) +
  scale_y_continuous(name = "Absolute risk", expand = c(0,0.001)) +
  scale_color_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"), labels = c("Symptomatic infections", "All infections"))+
  scale_fill_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"))+
  coord_cartesian(xlim = c(0,max(seq(1, 13, 2))), ylim = c(0,0.8), clip = "off") +
  geom_hline(yintercept =placebo_risk_sym, linetype = "dashed",  color = "black", size = 1) + 
  geom_hline(yintercept =placebo_risk_inf, linetype = "dashed",  color = "black", size = 1) + 
  annotate("text", x =9.5, y = 0.5, label = "Placebo overall infection risk") +
  annotate("text", x = 9, y = 0.43, label = "Placebo overall symptomatic risk") +
  annotate("rect", xmin = 1.2, xmax = 2, ymin = 0.74, ymax = 0.8, fill = "#E1ADA7", color = "black") + 
  annotate("rect", xmin = 1.2, xmax = 2, ymin = 0.67, ymax = 0.73, fill = "#1f78b4", color = "black") + 
  annotate("text", x =5.2, y = 0.77, label = "Symptomatic infections") +
  annotate("text", x =4, y = 0.7, label = "All infectinos") +
  theme1 +
  ggtitle(label = "A") -> p4



ggplot()+
  # geom_density(data = dat_vac %>% filter(type == 1), aes(x = log(SI2_nt,2)), alpha=0.15, color = "#1f78b4", fill = "#1f78b4") +
  geom_density(data = dat_vac %>% filter(type == 2), aes(x = log(SI2_nt,2)), alpha=0.15,color = "#E1ADA7", fill = "#E1ADA7") +
  geom_density(data = dat_vac %>% filter(type == 3), aes(x = log(SI2_nt,2)), alpha=0.15,color = "#1f78b4", fill = "#1f78b4") +
  geom_line(data = pred %>% filter(type != 1), aes(x = log(SI2_nt,2), y = rel_risk/4, color = type),size = 1) +
  geom_ribbon(data = pred %>% filter(type != 1), aes(x = log(SI2_nt,2), ymin = rel_risk_lci/4, ymax = rel_risk_uci/4, fill = type), alpha=0.15) +
  scale_x_continuous(name = "Neutralization titer", expand = c(0,0),breaks = seq(1, 13, 2), 
                     labels = 2^seq(1, 13, 2), limits = c(0, max(seq(1, 13, 2)))) +
  scale_y_continuous(name = "Relative risk", expand = c(0,0), breaks = c(0, 0.4, 0.8, 1.2, 1.6)/4,limits = c(0,1.6)/4,
                     labels =c("0", "0.4", "0.8", "1.2","1.6"),
                     sec.axis = sec_axis(trans = ~.*1, name = "Density\n ", breaks = c(0, 0.1, 0.2, 0.3, 0.4))) +
  scale_color_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"), labels = c("Symptomatic infections", "All infections"))+
  scale_fill_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"))+
  coord_cartesian(xlim = c(0,max(seq(1, 13, 2))), clip = "off") +
  theme1 + ggtitle(label = "B") -> p5


  

ggplot()+
  # geom_density(data = dat_vac %>% filter(type == 1), aes(x = log(SI2_nt,2)), alpha=0.15, color = "#1f78b4", fill = "#1f78b4") +
  geom_density(data = dat_vac %>% filter(type == 2), aes(x = log(SI2_nt,2)), alpha=0.15,color = "#E1ADA7", fill = "#E1ADA7") + 
  geom_density(data = dat_vac %>% filter(type == 3), aes(x = log(SI2_nt,2)), alpha=0.15,color = "#1f78b4", fill = "#1f78b4") + 
  geom_line(data = pred %>% filter(type != 1), aes(x = log(SI2_nt,2), y = ve/3.3, color = type),size = 1) +
  geom_ribbon(data = pred %>% filter(type != 1), aes(x = log(SI2_nt,2), ymin = ve_lci/3.3, ymax = ve_uci/3.3, fill = type), alpha=0.15) +
  scale_x_continuous(expand = c(0,-0.5), name = "Neutralization titer",breaks = seq(1, 13, 2),
                     labels = 2^seq(1, 13, 2), limits = c(1, max(seq(1, 13, 2)))) +
  scale_y_continuous(name = "Vaccine efficacy (%)", expand = c(0,0), breaks = c(0, 0.25, 0.5, 0.75, 1.01)/3.3,limits = c(0,1.01)/3.3,
                     labels =c("0", "25", "50", "75","100"),
                     sec.axis = sec_axis(trans = ~.*1, name = "Density\n ", breaks = c(0, 0.1, 0.2, 0.3))) +
  scale_color_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"), labels = c("Symptomatic infections", "All infections"))+
  scale_fill_manual(values = c("#E1ADA7","#1f78b4", "#66CC99"))+
  coord_cartesian(xlim = c(0,max(seq(1, 13, 2))+1), clip = "off") +
  annotate("segment", x = log(0,2), xend = log(42,2), y = 0.5/3.3, yend = 0.5/3.3, size=1, color="grey", linetype = 2)+
  annotate("segment", x = log(30,2), xend = log(30,2), y = 0, yend = 0.5/3.3, size=1, color="grey", linetype = 2)+
  annotate("segment", x = log(42,2), xend = log(42,2), y = 0, yend = 0.5/3.3, size=1, color="grey", linetype = 2)+
  annotate("text", x = 4.2, y = 0.015, label = "30", size = 5) +
  annotate("text", x = 6, y = 0.015, label = "42", size = 5) +  theme1 + ggtitle(label = "C")  ->  p6



# tiff("NT_COP.tiff", width = 42, height = 22, units = "cm", res = 400)
# (p1+p2+p3) / (p4+p5+p6)
# dev.off()


tiff("NT_COP_sen1.tiff", width = 42, height = 12, units = "cm", res = 400)
(p4+p5+p6)
dev.off()

