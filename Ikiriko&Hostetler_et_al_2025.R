#Data Analysis and Figure Generation for publication associated with Ikiriko and Hostetler, 2025#
library(dplyr)
citation("dplyr")
packageVersion("dplyr")
library(ggplot2)
citation("ggplot2")
packageVersion("ggplot2")
library(tidyr)
citation("tidyr")
packageVersion("tidyr")
library(cowplot)
citation("cowplot")
packageVersion("cowplot")
library(agricolae)
citation("agricolae")
packageVersion("agricolae")
library(rcompanion)
citation("rcompanion")
packageVersion("rcompanion")
library(lubridate) #monthday
R.version

cat("\014")
rm(list=ls()) 
ls() 
setwd(dir = "/Users/ashley/Desktop/BMMB/Data/")
getwd()
#Figure 1####
data0 = read.csv("PusherDatabase_BMMB.csv", header = TRUE, na.strings = "NA")
data=data0
unique(data$Additional.Notes)
data = subset(data, Experiment == "BMMB")
data = subset(data, Accession == "B73" | Accession == "Mo17" | Accession == "CML258")
data = subset(data, Year == "2021" | Year == "2023" | Year == "2024")
data = subset(data, Plot.ID !="202" & Plot.ID !="203" & Plot.ID != "204"  & Plot.ID !="206"  & Plot.ID !="245"  & Plot.ID !="345"  & Plot.ID !="352"  & Plot.ID !="354"  & Plot.ID !="356") #These plots were removed because they were part of the "Tassling objective
data$Plant.Number = replace(data$Plant.Number, data$Plant.Number == "1.1", "7")
data$TempID = paste(data$Plot.ID, data$Plant.Number, data$Date, data$TimeOfTesting, sep = "_")
data = subset(data, TempID != "964_5_9/27/21_12:22 PM") #This was a repeat save of Plant 4 without repeat testing
data$TempID = NULL
data$Condition = paste(data$Year, data$Accession, sep = "_")
data$ID = paste(data$Plot.ID, data$Plant.Number, data$Year, sep = "_")
data = as.data.frame(data[,c(20, 21, 2, 5, 3, 7, 8, 15, 19)])
colnames(data)[8] = "EI"
data$PlantingDate = data$Year
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2021", "5/20/2021")
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2023", "5/17/2023")
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2024", "5/20/2024")
data$PlantingDate <- as.Date(data$PlantingDate, format="%m/%d/%Y")
data$Date = as.Date(data$Date, format="%m/%d/%y")
data$DAP = difftime(data$Date, data$PlantingDate, units = "days")
data$DAP = as.numeric(data$DAP)
head(data)
X = data %>%
  group_by(Condition) %>%
  summarise(unique_values = list(unique(DAP)))
X = as.data.frame(X)
X
FullPusherDataFile = data
counts = data %>%
  group_by(ID) %>%
  summarise(count = n())
summary_counts = counts %>%
  group_by(count) %>%
  summarise(rows_with_count = n())
summary_counts
data = data %>%
  group_by(ID) %>%
  filter(n_distinct(EI) > 5) %>%
  ungroup()
data$Accession <- factor(data$Accession, levels = c("B73","Mo17", "CML258"))
data %>%
  group_by(Year) %>%
  summarize(
    First_Date = min(Date), 
    Last_Date = max(Date)
  )
result = data %>%
  distinct(Condition, Date, .keep_all = TRUE) %>%
  group_by(Condition) %>%
  arrange(DAP, .by_group = TRUE) %>%
  mutate(Diff_Days = DAP - lag(DAP)) %>%
  ungroup()
result2 = result %>%
  group_by(Condition) %>%
  summarise(
    average_DAP = mean(Diff_Days, na.rm = TRUE),
    range_DAP = max(Diff_Days, na.rm = TRUE) - min(Diff_Days, na.rm = TRUE),
    min_DAP = min(Diff_Days, na.rm = TRUE),
    max_DAP = max(Diff_Days, na.rm = TRUE)
  )
print(result2)
PlotA = ggplot(data, aes(x=DAP, y=EI, color = Accession)) +
  geom_line(aes(x=DAP, y=EI, group = ID, color = Accession), linewidth=0.25)+
  geom_point(aes(x=DAP, y=EI), color="black", size=0.1)+
  ylab("Stalk Flexural Stiffness Nm2")+
  xlab("Days After Planting")+
  scale_color_manual(values = c("pink4", "slategray3", "mistyrose3"))+
  scale_y_continuous(limits=c(0,150), breaks=seq(0,150,25))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(Year~Accession)
PlotA
results_df = data.frame() 
for (i in unique(data$Condition)) { 
  subset_data = data[which(data$Condition == i),]
  model = nls(EI ~ SSlogis(DAP, Asym, xmid, scal), subset_data)
  subset_data$pred_EI = predict(model)
  coeffs = coef(model) 
  Asym = coeffs['Asym'] 
  xmid = coeffs['xmid'] 
  scal = coeffs['scal'] 
  temp_df <- data.frame(
    Condition = i,
    Year = subset_data$Year,
    Accession = subset_data$Accession,
    DAP = subset_data$DAP,
    Pred_EI = subset_data$pred_EI,
    Asymp = Asym,
    Xmid = xmid,
    scale = scal
  )
  temp_df = temp_df[!duplicated(temp_df$DAP), ]
  results_df <- rbind(results_df, temp_df)
}
results_df$ID = results_df$Condition
results_df=results_df %>% separate(ID, c("Year", "Genotype"), "_")
PlotB = ggplot(results_df, aes(x=DAP, color = Accession)) +
  geom_point(aes(y = Pred_EI, group = Condition), color="black") +  # Regular points
  geom_line(aes(y = Pred_EI, group = Condition, color = Accession)) +
  geom_hline(data = results_df, aes(yintercept = Asymp, group = Condition), color = "black", linewidth=0.25) +  # Asymptote
  scale_color_manual(values = c("pink4", "slategray3", "mistyrose3"))+
  xlab("Days After Planting")+
  ylab("")+
  scale_y_continuous(limits=c(0,150), breaks=seq(0,150,25))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(Year~Accession)
PlotB
result3 = results_df %>%
  group_by(Condition) %>%
  filter(!duplicated(Asymp)) %>%
  ungroup()
result3 = result3[,c(1,2,5:9)]
result3$change = (result3$scale/4)*result3$Asymp
result3
PusherData = data
ModelData = results_df
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PlotA","PlotB")))
#Pythium#
data0 = read.csv("PusherDatabase_BMMB.csv", header = TRUE, na.strings = "NA")
data=data0
unique(data$Experiment)
data = subset(data, Experiment == "Pythium")
data$filename = data$Raw.Data.Label
colnames(data)
data = data[,c(1:17)]
data$Additional.Notes = replace(data$Additional.Notes, data$Additional.Notes == "Pythium", "P")
data$Additional.Notes = replace(data$Additional.Notes, data$Additional.Notes == "Control", "C")
data$Additional.Notes = replace(data$Additional.Notes, data$Additional.Notes == "c", "C")
data$Additional.Notes = replace(data$Additional.Notes, data$Additional.Notes == "p", "P")
data$Additional.Notes = replace(data$Additional.Notes, data$Additional.Notes == "p'", "P")
data$Treatment = data$Additional.Notes
df = subset(data, Year == "2023")
colnames(df)[7] = "Field"
colnames(df)[8] = "Pair"
df$PlantingDate = df$Field
df$PlantingDate = replace(df$PlantingDate, df$PlantingDate == "1", "4/20/2023")
df$PlantingDate = replace(df$PlantingDate, df$PlantingDate == "2", "4/19/2023")
df$PlantingDate = replace(df$PlantingDate, df$PlantingDate == "3", "4/18/2023")
df$PlantingDate <- as.Date(df$PlantingDate, format="%m/%d/%Y")
df$Date = as.Date(df$Date, format="%m/%d/%y")
df$DAP = difftime(df$Date, df$PlantingDate, units = "days")
df$DAP = as.numeric(df$DAP)
X = df %>%
  group_by(Field) %>%
  summarise(unique_values = list(unique(DAP)))
X = as.data.frame(X)
X
data = df
data$ID = paste(data$Field, data$Pair, data$Treatment, sep = "_")
counts = data %>%
  group_by(ID) %>%
  summarise(count = n())
summary_counts = counts %>%
  group_by(count) %>%
  summarise(rows_with_count = n())
summary_counts
X = df %>%
  group_by(Field, Pair) %>%
  summarise(unique_values = list(unique(DAP)))
colnames(data)[15] = "EI"
data = data %>%
  group_by(ID) %>%
  filter(n_distinct(EI) > 2) %>%
  ungroup()
data = as.data.frame(data)
PlotC= ggplot(data, aes(x=DAP, y=EI, color = Treatment)) +
  geom_line(aes(x=DAP, y=EI, group = ID, color = Treatment), linewidth=0.25)+
  geom_point(aes(x=DAP, y=EI), color="black", size=0.1)+
  ylab("Stalk Flexural Stiffness Nm2")+
  xlab("Days After Planting")+
  scale_color_manual(values = c("indianred4", "burlywood"))+
  scale_y_continuous(limits=c(0,55), breaks=seq(0,55,10))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(.~Field)
PlotC
results_df = data.frame() 
head(data)
data$Condition = paste(data$Field, data$Treatment, sep = "_")
for (i in unique(data$Condition)) { 
  subset_data = data[which(data$Condition == i),]
  model = nls(EI ~ SSlogis(DAP, Asym, xmid, scal), subset_data)
  subset_data$pred_EI = predict(model)
  coeffs = coef(model) 
  Asym = coeffs['Asym'] 
  xmid = coeffs['xmid'] 
  scal = coeffs['scal'] 
  temp_df <- data.frame(
    Condition = i,
    Treatment = subset_data$Treatment, 
    Field = subset_data$Field,
    DAP = subset_data$DAP,
    Pred_EI = subset_data$pred_EI,
    Asymp = Asym,
    Xmid = xmid,
    scale = scal
  )
  temp_df = temp_df[!duplicated(temp_df$DAP), ]
  results_df <- rbind(results_df, temp_df)
}
results_df$ID = results_df$Condition
PlotD = ggplot(results_df, aes(x=DAP, color = Treatment)) +
  geom_point(aes(y = Pred_EI, group = Treatment), color="black") +  # Regular points
  geom_line(aes(y = Pred_EI, group = Treatment, color = Treatment)) +
  geom_hline(data = results_df, aes(yintercept = Asymp, group = Treatment), color = "black", linewidth=0.25) +  # Asymptote
  scale_color_manual(values = c("indianred4", "burlywood"))+
  xlab("Days After Planting")+
  ylab("")+
  scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(.~Field)
PlotD
result3 = results_df %>%
  group_by(Condition) %>%
  filter(!duplicated(Asymp)) %>%
  ungroup()
result3
result3$change = (result3$scale/4)*result3$Asymp
result3
ggdraw() +
  draw_plot(PlotA, x = 0, y = 0.3, width = 0.5, height = 0.70) +
  draw_plot(PlotB, x = 0.5, y = 0.3, width = 0.5, height = 0.70) +
  draw_plot(PlotC, x = 0, y = 0, width = 0.5, height = 0.30) +
  draw_plot(PlotD, x = 0.5, y = 0, width = 0.5, height = 0.30) +
  draw_plot_label(label = c("A", "C", "B", "D"), 
                  size = 12,
                  x = c(0,0.5,0,0.5), 
                  y = c(1,1,0.3,0.3))
PythiumPusherData = data
PythiumModelData = results_df
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData")))
#Supplemental Figure 1####
df1 = read.csv("GDD_2021.csv", header = TRUE, na.strings = "NA")
df2 = read.csv("GDD_2023.csv", header = TRUE, na.strings = "NA")
df3 = read.csv("GDD_2024.csv", header = TRUE, na.strings = "NA")
data = rbind(df1, df2, df3)
data$Date = as.Date(data$Date, format="%m/%d/%Y")
data$month_day = format(data$Date, "%m/%d")
data$Year = format(data$Date, "%Y")
DateData = PusherData[,c(3,4,5,11)]
data = merge(data, DateData, by = c("Year", "Date"), all.x = TRUE)
data = data %>%
  distinct(Date, Accession, .keep_all = TRUE)
data = data %>%
  mutate(FTAccession = case_when(
    Year == 2021 & month_day == "08/11" ~ "CML258",
    Year== 2023 & month_day == "08/31" ~ "CML258",
    Year == 2023 & month_day == "08/02" ~ "B73",
    Year == 2024 & month_day == "08/20" ~ "CML258",
    Year == 2024 & month_day == "07/22" ~ "B73",
    TRUE ~ NA_character_  # If none of the conditions match, assign NA
  ))
data = data %>%
  mutate(FTAccession = case_when(
    Year== 2023 & month_day == "08/02"& Accession == "Mo17" ~ "Mo17",
    Year== 2023 & month_day == "08/02"& Accession == "CML258" ~ NA_character_,
    Year== 2023 & month_day == "08/31" & Accession == "CML258" ~ "CML258",
    Year== 2023 & month_day == "08/31" & Accession == "B73" ~ NA_character_,
    Year== 2023 & month_day == "08/31" & Accession == "Mo17" ~ NA_character_,
    Year== 2024 & month_day == "07/22" & Accession == "Mo17" ~ "Mo17",
    TRUE ~ FTAccession  # Keep existing value in column G if no condition is met
  ))
df = subset(data, FTAccession != "NA")
df = subset(df, Year == "2023" | Year == "2024")
df1 = df %>%
  group_by(Accession) %>%
  summarize(Cumulative = mean(Cumulative, na.rm = TRUE))
df1 = data.frame(
  Condition = c("2021_CML258","2023_Mo17","2023_CML258","2023_B73","2024_Mo17","2024_CML258","2024_B73"),
  Cumulative = c(2312, 1588, 2312, 1588, 1588, 2312, 1588),
  Year = c("2021", "2023", "2023", "2023", "2024", "2024", "2024"),
  Accession = c("CML258", "Mo17", "CML258", "B73", "Mo17", "CML258", "B73")
)
PlotA = ggplot(data, aes(x=month_day, y=Cumulative, color = Year)) +
  geom_line(aes(x=month_day, y=Cumulative, group = Year, color = Year))+
  geom_point(aes(x=month_day, y=Cumulative, group = Year, color = Year, shape = Accession), color="black", size=1)+
  geom_hline(data = df1, aes(yintercept = Cumulative, group = Year), color = "black", linewidth=0.25) +  # Asymptote
  ylab("Cumulative Growing Degree Days")+
  xlab("Date (Month/Day)")+
  ggtitle("A")+
  scale_color_manual(values = c("mistyrose3", "palegreen4", "indianred4"))+
  scale_shape_manual(values = c(19,19,19,19))+  # Manually setting shapes (16 = circle, 17 = triangle)
  theme_bw()+
  theme(
    axis.text.x = element_text(size=6, angle = 90, vjust=0.5),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
PlotA
result = data %>%
  group_by(Year) %>%
  mutate(
    first_date = min(Date),
    last_date = max(Date)
  ) %>%
  filter(Date >= first_date & Date <= last_date) %>%
  summarize(
    avg_DD = mean(Daily.DD, na.rm = TRUE)
  )
head(data)
GDD = data
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData","GDD","PlotA","df1")))
results_df = ModelData
results_df$PlantingDate = results_df$Year
results_df$PlantingDate = replace(results_df$PlantingDate, results_df$PlantingDate == "2021", "5/20/2021")
results_df$PlantingDate = replace(results_df$PlantingDate, results_df$PlantingDate == "2023", "5/17/2023")
results_df$PlantingDate = replace(results_df$PlantingDate, results_df$PlantingDate == "2024", "5/20/2024")
results_df$PlantingDate <- as.Date(results_df$PlantingDate, format="%m/%d/%Y")
results_df$Date = results_df$PlantingDate + results_df$DAP
GDD = GDD[,c(1,2,4)]
results_df = merge(results_df, GDD, by = c("Date","Year"))
results_df$Accession = factor(results_df$Accession, levels = c("B73","Mo17", "CML258"))
df1$Accession = factor(df1$Accession, levels = c("B73","Mo17", "CML258"))
PlotB = ggplot(results_df) +
  geom_point(aes(x=Cumulative, y = Pred_EI, group = Condition), color="black") +  # Regular points
  geom_line(aes(x=Cumulative, y = Pred_EI, group = Condition, color=Accession)) +
  geom_hline(data = results_df, aes(yintercept = Asymp, group = Condition), color = "black", linewidth=0.25) +  # Asymptote
  geom_vline(data = df1, aes(xintercept = Cumulative, group = Condition), color = "black", linewidth=0.25) +  # Asymptote
  scale_color_manual(values = c("pink4", "slategray3", "burlywood", "black"))+
  xlab("Cumulative Growing Degree Days")+
  ylab("")+
  ggtitle("B")+
  scale_y_continuous(limits=c(0,150), breaks=seq(0,150,25))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(Year~Accession, drop = TRUE)
PlotB
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData")))

#Supplemental Figure 3####
data = FullPusherDataFile
Control = data %>%
  group_by(ID) %>%
  filter(n_distinct(EI) < 2) %>%
  ungroup()
Control = as.data.frame(Control)
Control$Treatment = "Control"
Control = subset(Control, Year != "2024")
Control_sorted = Control[order(Control$Year, Control$Plot.ID, Control$Plant.Number), ]
unique(Control_sorted$Plot.ID)
Weekly = data %>%
  group_by(ID) %>%
  filter(n_distinct(EI) > 5) %>%
  ungroup()
Weekly = as.data.frame(Weekly)
Weekly$Treatment = "Weekly"
Weekly = subset(Weekly, Condition != "2024_CML258")
Weekly = subset(Weekly, Year != "2024")
Weekly_sorted = Weekly[order(Weekly$Year, Weekly$Plot.ID, Weekly$Plant.Number), ]
unique(Weekly_sorted$Plot.ID)
result = Control %>%
  group_by(Plot.ID, Plant.Number, Accession, Year) %>%
  summarise(RowCount = n(), .groups = "drop")
result = Weekly %>%
  group_by(Plot.ID, Plant.Number, Accession, Year) %>%
  summarise(RowCount = n(), .groups = "drop")
result$ID = paste(result$Plot.ID, result$Plant.Number, result$Year, sep="_")
Weekly = Weekly %>%
  group_by(ID) %>%
  filter(DAP == max(DAP, na.rm = TRUE)) %>%
  ungroup()
Weekly = as.data.frame(Weekly)
Weekly = semi_join(Weekly, Control, by = "Plot.ID")
Weekly_sorted = Weekly[order(Weekly$Year, Weekly$Plot.ID, Weekly$Plant.Number), ]
Control = semi_join(Control, Weekly, by = "Plot.ID")
Control_sorted = Control[order(Control$Year, Control$Plot.ID, Control$Plant.Number), ]
data = rbind(Weekly, Control)
data$Accession = factor(data$Accession, levels = c("B73","Mo17", "CML258"))
ggplot(data, aes(x=Treatment, y=EI, fill = Treatment)) +
  geom_boxplot() +
  ylab("Stalk Flexural Stiffness Nm2")+
  xlab("Testing")+
  scale_fill_manual(values = c("mistyrose3", "indianred4"))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~ Condition, nrow = 2, ncol = 2)
#Set up data frame to fill in
SW <- matrix(NA,nrow=2,ncol=2)
rownames(SW) <- c("W","pvalue")
SW[1,1] = "W"
SW[2,1] = "pvalue"
Z = "" #placeholder
condition=unique(data$Condition)
for (i in condition){
  a = i
  subset_data = subset(data, Condition == a)
  aov = lm(EI ~ Treatment, data = subset_data)
  hsd = aov(aov)              
  resid = residuals(object = hsd)
  shap = shapiro.test(x=resid) #Run a shapiro test on residuals 
  SW[1,2] = shap$statistic #extract information
  SW[2,2] = shap$p.value #extract information
  colnames(SW)[2] = a #extract information
  SW = as.data.frame(SW)
  write.table(a, file = "ShapiroResults.txt", sep = "\t",
              row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(SW, file = "ShapiroResults.txt", sep = "\t",
              row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(Z, file = "ShapiroResults.txt", sep = "\t",
              row.names = FALSE, col.names = FALSE, append=TRUE)
  if (shap$p.value > 0.05){  
    print(shapiro.test(x=resid)$p.value)
    aov = lm(EI ~ Treatment, data = subset_data)
    write.table(a, file = "ANOVA.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova(aov) 
    anova = anova(aov)
    anova = as.data.frame(anova)
    Z = ""
    write.table(anova, file = "ANOVA.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "ANOVA.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("Treatment"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "means.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "means.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "groups.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "groups.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  #for data that needs normalization
  if (shap$p.value < 0.05){
    print(shapiro.test(x=resid)$p.value)
    norm_col_name <- paste0("norm_", a)
    par(mfrow=c(3,1))
    subset_data$norm_col_name = transformTukey(subset_data$EI,
                                               start = -10,
                                               end = 10,
                                               int = 0.025,
                                               plotit = TRUE, #can be false 
                                               verbose = FALSE,
                                               quiet = FALSE,
                                               statistic = 1,
                                               returnLambda = FALSE  )
    aov = lm(norm_col_name ~ Treatment, data = subset_data)
    col = norm_col_name
    write.table(col, file = "ANOVA_TRANS.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    write.table(anova, file = "ANOVA_TRANS.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "ANOVA_TRANS.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)              
    hsd_test = HSD.test(hsd, trt = c("Treatment"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = col 
    write.table(means1, file = "means_trans.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "means_trans.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "groups_trans.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "groups_trans.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
}
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData")))
#Figure 2####
data0 = read.csv("PhenotypesDatabase_BMMB.csv", header = TRUE, na.strings = "NA")
data=data0
data = subset(data, Accession == "B73" | Accession == "Mo17" | Accession == "CML258")
data = subset(data, Year == "2021" | Year == "2023" | Year == "2024")
data = subset(data, Treatment !="PT" & Treatment !="ER" & Treatment !="DT" )
data$AvgDiameter = (data$Skinny.Stem.Diameter..mm. + data$Fat.Stem.Diameter..mm.) /2
data$Diameter..mm. = as.numeric(data$Diameter..mm.)
data = data %>%
  mutate(StalkDiameter = ifelse(!is.na(Diameter..mm.), Diameter..mm., AvgDiameter))
data = data[,c(1,3,7:12,32)]
data$Condition = paste(data$Year, data$Accession, sep = "_")
data$ID = paste(data$Plot, data$Plant, data$Year, sep = "_")
data$PlantingDate = data$Year
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2021", "5/20/2021")
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2023", "5/17/2023")
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2024", "5/20/2024")
data$PlantingDate <- as.Date(data$PlantingDate, format="%m/%d/%Y")
data$CollectionDate = as.Date(data$CollectionDate, format="%m/%d/%Y")
data$CollectionDate = replace(data$CollectionDate, data$CollectionDate == "2024-01-10", "2024-10-01")
data$DAP = difftime(data$CollectionDate, data$PlantingDate, units = "days")
data$DAP = as.numeric(data$DAP)
counts = data %>%
  group_by(ID) %>%
  summarise(count = n())
summary_counts = counts %>%
  group_by(count) %>%
  summarise(rows_with_count = n())
summary_counts
counts = subset(counts, count > 8)
counts = as.data.frame(counts)
data = merge(data, counts, by = "ID")
data = data[,c(1:13)]
colnames(data)[9] = "PH"
data$Accession = factor(data$Accession, levels = c("B73","Mo17", "CML258"))
PlotA = ggplot(data, aes(x=DAP, y=StalkDiameter, color = Accession)) +
  geom_line(aes(x=DAP, y=StalkDiameter, group = ID, color = Accession), linewidth=0.25, na.rm = TRUE )+
  geom_point(aes(x=DAP, y=StalkDiameter), color="black", size=0.1)+
  ylab("Stalk diameter (mm)")+
  xlab("Days After Planting")+
  scale_color_manual(values = c("pink4", "slategray3", "mistyrose3"))+
  scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(Year~Accession)
PlotA
df = data
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData","df", "PlotA")))
df = as.data.frame(df)
df = as.data.frame(df)
colnames(df)[3] = "Date"
data = PusherData
df$Date = replace(df$Date, df$Date == "2023-07-18", "2023-07-20") #Plant height data was taken on 07/18/2023 but Pusher was run on 07/20/2023 so I changed this date to make them merge-able
df$Date = replace(df$Date, df$Date == "2021-09-08", "2021-09-10") #Final Plant height data was taken on 09/08/2021 but Pusher was run on 09/10/2021 so I changed this date to make them merge-able
df$Date = replace(df$Date, df$Date == "2023-09-14", "2023-09-15") #Plant height data was taken on 09/14/2023 but Pusher was run on 09/15/2023 so I changed this date to make them merge-able 
df$Date = replace(df$Date, df$Date == "2024-08-14", "2024-08-16") #Plant height data was taken on 08/14/2024 but Pusher was run on 08/16/2024 so I changed this date to make them merge-able
data = data[,c(2,5,8)]
data = merge(data, df, by = c("ID", "Date"))
data$RatioSD = data$EI / data$StalkDiameter
counts = data %>%
  group_by(ID) %>%
  summarise(count = n())
PlotB = ggplot(data, aes(x=DAP, y=RatioSD, color = Accession)) +
  geom_line(aes(x=DAP, y=RatioSD, group = ID, color = Accession), linewidth=0.25, na.rm = TRUE )+
  geom_point(aes(x=DAP, y=RatioSD), color="black", size=0.1)+
  ylab("Ratio EI/Stalk Diameter")+
  xlab("Days After Planting")+
  scale_color_manual(values = c("pink4", "slategray3", "mistyrose3"))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(Year~Accession)
PlotB
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData","df", "PlotA","PlotB")))
#Pythium#
data0 = read.csv("Pythium_Diameters_2023.csv", header = TRUE, na.strings = "NA")
data=data0
df1 = data[,c(1:8)]
df1 = df1 %>%
  pivot_longer(
    cols = July.4:Aug.28,       # Select the columns to reshape
    names_to = "Date",                 # Name of the new key column
    values_to = "Diameter"                # Name of the new value column
  )
df1 = as.data.frame(df1)
df1$Date = replace(df1$Date, df1$Date == "July.4", "07/04/2023")
df1$Date = replace(df1$Date, df1$Date == "July.10", "07/10/2023")
df1$Date = replace(df1$Date, df1$Date == "July.17", "07/17/2023")
df1$Date = replace(df1$Date, df1$Date == "July.31", "07/31/2023")
df1$Date = replace(df1$Date, df1$Date == "Aug.7", "08/07/2023")
df1$Date = replace(df1$Date, df1$Date == "Aug.21", "08/21/2023")
df1$Date = replace(df1$Date, df1$Date == "Aug.28", "08/28/2023")
df1$Date = as.Date(df1$Date, format="%m/%d/%Y")
df1 = na.omit(df1)
colnames(df1)[1] = "ID"
df1$ID = gsub("-", "_", df1$ID)
data = merge(PythiumPusherData, df1, by = c("ID", "Date"), all = TRUE)
PythiumPusherData = data
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData","df", "PlotA","PlotB")))
data0 = read.csv("Pythium_Phenotype_2023.csv", header = TRUE, na.strings = "NA")
data=data0
data = within(data, {
  original_Grid = Grid
  Grid = ifelse(grepl("4", original_Grid), gsub("4", "3", original_Grid), original_Grid)
  Pair....Rep. = ifelse(grepl("4", original_Grid), Pair....Rep. + 5, Pair....Rep.)
})
original_Grid =  NULL
data = data[,c(2,3,5:16)]
data$ID = paste(data$Grid, data$Pair....Rep., data$Treatment, sep = "_")
data = data[,c(15,3:5,7:14)]
data$Date = "08/28/2023"
data$Date = as.Date(data$Date, format="%m/%d/%Y")
data = merge(PythiumPusherData, data, by = c("ID", "Date"), all = TRUE)
counts = data %>%
  group_by(ID) %>%
  summarise(count = n())
summary_counts = counts %>%
  group_by(count) %>%
  summarise(rows_with_count = n())
summary_counts
data = data %>%
  group_by(ID) %>%
  filter(n_distinct(EI) > 2) %>%
  ungroup()
PlotC = ggplot(data, aes(x=DAP, y=Diameter, color = Treatment)) +
  geom_line(aes(x=DAP, y=Diameter, group = ID, color = Treatment), linewidth=0.25)+
  geom_point(aes(x=DAP, y=Diameter), color="black", size=0.1)+
  ylab("Stalk Diameter (mm)")+
  xlab("Days After Planting")+
  scale_color_manual(values = c("indianred4", "burlywood"))+
  scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(.~Field)
PlotC
data$RatioSD = data$EI / data$Diameter
PlotD = ggplot(data, aes(x=DAP, y=RatioSD, color = Treatment)) +
  geom_line(aes(x=DAP, y=RatioSD, group = ID, color = Treatment), linewidth=0.25)+
  geom_point(aes(x=DAP, y=RatioSD), color="black", size=0.1)+
  ylab("Ratio EI/SD")+
  xlab("Days After Planting")+
  scale_color_manual(values = c("indianred4", "burlywood"))+
  scale_y_continuous(limits=c(0,2), breaks=seq(0,2,0.25))+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(.~Field)
PlotD
ggdraw() +
  draw_plot(PlotA, x = 0, y = 0.3, width = 0.5, height = 0.70) +
  draw_plot(PlotB, x = 0.5, y = 0.3, width = 0.5, height = 0.70) +
  draw_plot(PlotC, x = 0, y = 0, width = 0.5, height = 0.30) +
  draw_plot(PlotD, x = 0.5, y = 0, width = 0.5, height = 0.30) +
  draw_plot_label(label = c("A", "C", "B", "D"), 
                  size = 12,
                  x = c(0,0.5,0,0.5), 
                  y = c(1,1,0.3,0.3))
rm(list = setdiff(ls(), c("PusherData", "ModelData","FullPusherDataFile","PythiumPusherData","PythiumModelData")))
#Figure 3 and S4####
data0 = read.csv("InstronDatabase_BMMB.csv", header = TRUE, na.strings = "NA")
data=data0
data = subset(data, Notes != "node")
data$TestingDate = as.Date(data$TestingDate, format = "%m/%d/%y")
data$TestingDate = replace(data$TestingDate, data$TestingDate == "2020-08-06", "2024-08-06")
data$TestingDate = replace(data$TestingDate, data$TestingDate == "2020-08-10", "2024-08-10")
data$TestingDate = replace(data$TestingDate, data$TestingDate == "2020-08-19", "2024-08-19")
data$TestingDate = replace(data$TestingDate, data$TestingDate == "2020-08-22", "2024-08-22")
data$TestingDate = replace(data$TestingDate, data$TestingDate == "2020-09-06", "2024-09-06")
data$TestingDate = as.Date(data$TestingDate, format = "%m/%d/%y")
colnames(data)[17] = "K"
data$a0 = (data$Vertical.Specimen.Diameter)/2
data$b0 = (data$Horizontal.Specimen.Diameter)/2
data$I = (pi/4)*(((data$a0)^3)*(data$b0))
data$E = data$K * ((200^3)/(48*(data$I)))
data$PlantingDate = data$Year
data$PlantingDate = replace(data$PlantingDate, data$PlantingDate == "2024", "5/20/2024")
data$PlantingDate <- as.Date(data$PlantingDate, format="%m/%d/%Y")
data$DAP = difftime(data$TestingDate, data$PlantingDate, units = "days")
data$DAP = as.numeric(data$DAP)
data$Genotype <- factor(data$Genotype, levels = c("B73","Mo17", "CML258"))
data$DAP2 = data$DAP
data$DAP2 =  as.numeric(as.character(data$DAP2))
head(data)
X = data %>%
  group_by(Genotype) %>%
  summarise(unique_values = list(unique(DAP)))
X = as.data.frame(X)
X
data = data %>%
  mutate(DAP2 = case_when(
    DAP2 %in% c(57, 58) ~ 57,
    DAP2 %in% c(63, 65) ~ 64,
    DAP2 %in% c(91, 94) ~ 92,
    TRUE ~ DAP2 # Keep other values unchanged
  ))
head(data)
colnames(data)[5] = "Accession"
df = ModelData
df = subset(df, Year == "2024")
PlotB= ggplot() +
  geom_line(data = df, aes(x = DAP, y = Pred_EI*20, group = Accession), color = "gray80") +
  geom_boxplot(data = data, aes(x = DAP2, y = E, group = DAP2, fill = Accession))+
  scale_fill_manual(values = c("pink4", "slategray3", "burlywood"))+
  xlab("Days After Planting")+
  scale_y_continuous(name = "E (MPa)",
                     sec.axis = sec_axis(~ . / 20, name = "Predicted EI (Nm2)")) +
  theme_bw()+
  theme(    
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(~Accession)
PlotB
PlotA = ggplot() +
  geom_line(data = df, aes(x = DAP, y = Pred_EI, group = Accession), color = "gray80") +
  geom_boxplot(data = data, aes(x = DAP2, y = K, group = DAP2, fill = Accession))+
  scale_fill_manual(values = c("pink4", "slategray3", "burlywood"))+
  xlab("Days After Planting")+
  scale_y_continuous(name = "K (N/mm)",
                     sec.axis = sec_axis(~ . / 1, name = "Predicted EI (Nm2)")) +
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(~Accession)
PlotA
ggdraw() +
  draw_plot(PlotA, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(PlotB, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), 
                  size = 12,
                  x = c(0,0), 
                  y = c(1,0.5))
FigureS4 = ggplot() +
  geom_boxplot(data = data, aes(x = DAP2, y = I, group = DAP2, fill = Accession))+
  scale_fill_manual(values = c("pink4", "slategray3", "burlywood"))+
  xlab("Days After Planting")+
  ylab("I (mm4)")+
  theme_bw()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(size=10, angle = 60, hjust=1),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_grid(~Accession)
FigureS4

attach(data)
head(data)
lm_x <- lm(I ~ DAP2*Accession)
anova(lm_x)
lm_x_aov=aov(lm_x) 
HSD.test(lm_x_aov, trt = c("Accession","DAP2"), console = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data$I)
hist(data$I)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data to meet assumptions
par(mfrow=c(2,2))
data$tukey <- transformTukey(
  data$I,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data$tukey)
lm_x2 <- lm(data$tukey ~ DAP2*Accession)
anova(lm_x2)
lm_x2_aov2=aov(lm_x2) 
resid = residuals(object = lm_x2_aov2)
shapiro.test(x=resid)
HSD.test(lm_x2_aov2, trt = c("Accession","DAP2"), console = TRUE)
data$tukey = NULL
detach(data)
