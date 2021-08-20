##### rectal cancer segmentation training/test performance ggplot


library(readxl)
library(ggplot2)
library(dplyr)

df <- read_excel("training_result.xlsx")
vincent <- read_excel("test_vincent.xlsx")

df2 <- subset(df, Flag =="mean")
df3 <- df2[c(1,4,6,8,10,12),]
df3[1,1] <- "Overall"
df3[2,1] <- "Fold 0"
df3[3,1] <- "Fold 1"
df3[4,1] <- "Fold 2"
df3[5,1] <- "Fold 3"
df3[6,1] <- "Fold 4"

df3 <- rename(df3,CV_number = File)

ggplot(df3,aes(x = CV_number, y = Dice, color = CV_number, fill = CV_number, label= round(Dice,digits = 4))) +
  geom_col() +
  geom_text(nudge_y=0.02, size = 5)+
  ggtitle("Perfomance Bar Graph for the Training set") + 
  xlab("CV #") + ylab("DICE") +
  theme_bw() + theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust = 0.5, size=20))+
  theme(axis.title = element_text(size=14))+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size =15)) +
  ylim(0L, 1L)


ggplot(vincent,aes(x = Case, y = Dice, color = Case, fill = Case, label= round(Dice,digits = 4))) +
  geom_col() +
  geom_text(nudge_y=0.02, size = 5)+
  ggtitle("Perfomance Bar Graph for the Test set - Vincent") + 
  xlab("Case #") + ylab("DICE") +
  theme_bw() + theme(legend.position = "right") + 
  theme(plot.title = element_text(hjust = 0.5, size=20))+
  theme(axis.title = element_text(size=14))+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size =15)) +
  ylim(0L, 1L)

