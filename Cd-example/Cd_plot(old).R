#create cadmium plot for report

library (ggplot2)
library(readr)

Cd<- read_csv("C:/R-repositories/SSD-fitting/Cd-example/Cd_BC_Multi.csv")
View(Cd)

#specificy factor for later sorting
Cd$Group<-factor(Cd$Group,levels = c("Amphibian","Fish (non-salmonid)","Fish (salmonid)","Invertebrate","Plant"))

Cdplot<-ggplot(data=Cd,aes(x=Conc,y=Species,shape=Group))+
  geom_point(position="jitter")+
  scale_x_log10(breaks=10^(-2:3))+
  scale_shape_manual(values = c(0,1,2,3,4))+
  geom_point(size=2)+
  xlab("Log Concentration")
  #scale_y_discrete(limits=c("Amphibian","Fish (non-salmonid)","Fish (salmonid)","Invertebrate","Plant"))

Cdplot

ggsave("Cd-example/Plots/Cdplot.png", Cdplot, width = 5.97, height = 5.97)
