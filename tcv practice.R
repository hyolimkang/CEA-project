# export data file from excel
library(readxl)
tcv<-read_excel("TCV_COI_merged.xlsx")
# export certain variables in tcv data
select(tcv, cd_all, caretype) %>%
View()
