my_theme<-theme(plot.title = element_text(hjust = 0.5), #fromPlotsofPlots0000000000
                panel.grid.minor.y = element_blank(),
                plot.margin = margin(5, 10, 5, 5, "pt"),
                axis.text.x = element_text(angle=45, hjust = 1,size=8),axis.title.x = element_blank(),legend.position="none")


p1<-sigMP2 %>%
  ggplot(aes(x = ASV, y = lfc_mp,fill=direction_mp)) + 
  scale_fill_manual(values = c("Positive" = "#00A86B", "Negative" = "#D2C500"))+
  ylab("Log fold change")+ggtitle("MP vs PP")+
  geom_bar(stat = "identity", width = 1.75/length(unique(sigMM2$ASV))) +
  geom_errorbar(aes(ymin = lfc_mp - se_mp, ymax = lfc_mp + se_mp), 
                width = 0.2, position = position_dodge(0.05), color = "black")+
  scale_x_discrete(labels = word(sigMP2$Species,1))+my_theme



p2<-sigMM2 %>%
  ggplot(aes(x = ASV, y = lfc_mm,fill=direction_mm)) + 
  scale_fill_manual(values = c("Positive" = "#002594", "Negative" = "#D2C500"))+
  ylab("Log fold change")+ggtitle("MM vs PP")+
  geom_bar(stat = "identity",width=0.4) +
  geom_errorbar(aes(ymin = lfc_mm - se_mm, ymax = lfc_mm + se_mm), 
                width = 0.2, position = position_dodge2(0.05), color = "black")+ 
  scale_x_discrete(labels = word(sigMM2$Species,1))+my_theme


p3<-sigPM2 %>%
  ggplot(aes(x = ASV, y = lfc_pm,fill=direction_pm)) + 
  scale_fill_manual(values = c("Positive" = "#E0B2CD", "Negative" = "#D2C500"))+
  ylab("Log fold change")+ggtitle("PM vs PP")+
  geom_bar(stat = "identity", width = 1.75/length(unique(sigMM2$ASV))) +
  geom_errorbar(aes(ymin = lfc_pm - se_pm, ymax = lfc_pm + se_pm), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  scale_x_discrete(labels = word(sigPM2$Species,1))+my_theme


p4<-sigMP %>%
  ggplot(aes(x = ASV, y = lfc_mp,fill=direction_mp)) + 
  scale_fill_manual(values = c("Positive" = "#00A86B", "Negative" = "#002594"))+
  ylab("Log fold change")+ggtitle("MP vs MM")+
  geom_bar(stat = "identity", width = 1.75/length(unique(sigMM2$ASV))) +
  geom_errorbar(aes(ymin = lfc_mp - se_mp, ymax = lfc_mp + se_mp), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  scale_x_discrete(labels = word(sigMP$Species,1))+my_theme



p5<-sigPP %>%
  ggplot(aes(x = ASV, y = lfc_pp,fill=direction_pp)) + 
  scale_fill_manual(values = c("Positive" = "#D2C500", "Negative" = "#002594"))+
  ylab("Log fold change")+ggtitle("PP vs MM")+
  geom_bar(stat = "identity", width=1.75/length(unique(sigMM2$ASV))) +
  geom_errorbar(aes(ymin = lfc_pp - se_pp, ymax = lfc_pp + se_pp), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  scale_x_discrete(labels = word(sigPP$Species,1))+my_theme


p6<-sigPM %>%
  ggplot(aes(x = ASV, y = lfc_pm,fill=direction_pm)) + 
  scale_fill_manual(values = c("Positive" = "#E0B2CD", "Negative" = "#002594"))+
  ylab("Log fold change")+ggtitle("PM vs MM")+
  geom_bar(stat = "identity", width = 1.5/length(unique(sigMM2$ASV))) +
  geom_errorbar(aes(ymin = lfc_pm - se_pm, ymax = lfc_pm + se_pm), 
                width = 0.2, position = position_dodge(0.05), color = "black")+ 
  scale_x_discrete(labels = word(sigPM$Species,1))+my_theme

p1c<-p1+coord_cartesian(ylim = c(-6, 7),clip = "off")
p2c<-p2+coord_cartesian(ylim = c(-6, 7),clip = "off")
p3c<-p3+coord_cartesian(ylim = c(-6, 7),clip = "off")
p4c<-p4+coord_cartesian(ylim = c(-6, 7),clip = "off")
p5c<-p5+coord_cartesian(ylim = c(-6, 7),clip = "off")
p6c<-p6+coord_cartesian(ylim = c(-6, 7),clip = "off")


p3c+p1c+p2c+p4c+p6c+p5c+plot_layout(ncol=3)
