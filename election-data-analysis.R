data <- read.csv("All_States_GE.csv")
print(colnames(data))
data_f <- data[data$Assembly_No == 17 & data$Poll_No == 0, ]
data_f$Constituency_Name <- paste(data_f$Constituency_Name, "_", data_f$State_Name, sep = "")
data_w <- data_f[data_f$Position == 1, c('Constituency_Name', 'Vote_Share_Percentage')]
colnames(data_w) <- c('Constituency_Name', 'Vote_Share_Percentage_W')
data_r <- data_f[data_f$Position == 2, c('Constituency_Name', 'Vote_Share_Percentage')]
colnames(data_r) <- c('Constituency_Name', 'Vote_Share_Percentage_RU')
data_f$Vote_Share_Percentage_W <- data_w$Vote_Share_Percentage_W[match(data_f$Constituency_Name, data_w$Constituency_Name)]
data_f$Vote_Share_Percentage_RU <- data_r$Vote_Share_Percentage_RU[match(data_f$Constituency_Name, data_w$Constituency_Name)]
data_f$Win_Margin <- ifelse(data_f$Position == 1, data_f$Vote_Share_Percentage - data_f$Vote_Share_Percentage_RU, data_f$Vote_Share_Percentage - data_f$Vote_Share_Percentage_W)
data_f$Win_Margin <- data_f$Win_Margin/100
print(data_f$Party_ID)
data_f <- data_f[data_f$Party == "BJP", ]
print(data_f$State_Name)

BJP_ruled_UT <- c("Assam", "Bihar", "Goa", "Gujarat", "Haryana", "Himachal_Pradesh", "Jharkhand", "Maharashtra",
               "Manipur", "Nagaland", "Tripura", "Uttar_Pradesh", "Uttarakhand", "Dadra_&_Nagar_Haveli", "Daman_&_Diu", "Andaman_&_Nicobar_Islands", "Lakshadweep", "Chandigarh")

BJP_ruled <- c("Assam", "Bihar", "Goa", "Gujarat", "Haryana", "Himachal_Pradesh", "Jharkhand", "Maharashtra",
               "Manipur", "Nagaland", "Tripura", "Uttar_Pradesh", "Uttarakhand")


data_fb <- data_f[data_f$State_Name %in% BJP_ruled, ]
data_fnb <- data_f[!(data_f$State_Name %in% BJP_ruled), ]

print(data_f)
library(rdd)
DCdensity(data_f$Win_Margin)
DCdensity(data_fb$Win_Margin)
DCdensity(data_fnb$Win_Margin)


