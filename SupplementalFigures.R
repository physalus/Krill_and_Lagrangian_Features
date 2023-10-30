#### Supplemental Figure 1 ####
p_ctdDensityAllv2 <- plot_models(D50_GLMM24,D50_GLMM48,
                                 D50_GLMM120,D50_GLMM240,
                                 D40_GLMM24,D40_GLMM48,
                                 D40_GLMM120,D40_GLMM240,
                                 D30_GLMM24,D30_GLMM48,
                                 D30_GLMM120,D30_GLMM240,
                                 D20_GLMM24,D20_GLMM48,
                                 D20_GLMM120,D20_GLMM240,
                                 D10_GLMM24,D10_GLMM48,
                                 D10_GLMM120,D10_GLMM240, 
                                 ssSigmaT_GLMM24,ssSigmaT_GLMM48,
                                 ssSigmaT_GLMM120,ssSigmaT_GLMM240,
                                 vline.color = "red",show.p=TRUE,
                                 show.values=TRUE, spacing = 0.925, 
                                 show.legend = TRUE, dot.size = 5.5, 
                                 value.size = 8, line.size = 2,
                                 axis.labels = c("FTLE (240 Hrs)","FTLE (120 Hrs)",
                                                 "FTLE (48 Hrs)","FTLE (24 Hrs)")) +
  scale_y_continuous(breaks=c(-0.4,0,0.4,0.8), limits = c(-0.4,1.4)) +
  theme_minimal() +
  theme(axis.title = element_text(face="bold", size=30),
        axis.text = element_text( face="bold", size=30),
        plot.title = element_text(face="bold", size=25,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=20,hjust = 0.5),
        plot.caption = element_text( face="bold", size=12,hjust = 0),
        legend.position = c(0.785,0.5),
        legend.direction = "vertical",
        legend.key.size = unit(1.2,"cm"),
        legend.background = element_rect(fill="white",color="black"),
        legend.text = element_text(face="bold", size=28),
        legend.title = element_text(face="bold", size=30, hjust = 0.5))
p_ctdDensityAllv2$guides$colour$title <- "Response Variable"
p_ctdDensityAllv2$data$group <- factor(x = c(rep("Density at 50 meters",4), 
                                             rep("Density at 40 meters",4), 
                                             rep("Density at 30 meters",4), 
                                             rep("Density at 20 meters",4), 
                                             rep("Density at 10 meters",4), 
                                             rep("Density at Surface",4)), 
                                       levels = c("Density at 50 meters",
                                                  "Density at 40 meters",
                                                  "Density at 30 meters",  
                                                  "Density at 20 meters",
                                                  "Density at 10 meters", 
                                                  "Density at Surface") )
p_ctdDensityAllv2 <- p_ctdDensityAllv2 + 
  scale_color_manual(values=c("Density at 50 meters" =  "#0B775E",
                              "Density at 40 meters" =  "#35B779FF",
                              "Density at 30 meters" =  "#46ACC8",
                              "Density at 20 meters" =  "#7294D4",
                              "Density at 10 meters" =  "#31688EFF",
                              "Density at Surface" =  "#440154FF"
  ) 
  )
p_ctdDensityAllv2
ggsave(sprintf("./Output/SupplementalFigure1.jpg"),plot = p_ctdDensityAllv2,
       width=14, height=18, units = "in", dpi=600)

#### Supplemental Figure 2 ####

p_ctdSigmathetaAll <- plot_models(Sigma26_GLMM24,Sigma26_GLMM48,
                                  Sigma26_GLMM120,Sigma26_GLMM240,
                                  Sigma25_5_GLMM24,Sigma25_5_GLMM48,
                                  Sigma25_5_GLMM120,Sigma25_5_GLMM240, 
                                  Sigma25_GLMM24,Sigma25_GLMM48,
                                  Sigma25_GLMM120,Sigma25_GLMM240,
                                  vline.color = "red",
                                  show.p=TRUE,show.values=TRUE, 
                                  spacing = .85, p.shape=FALSE,
                                  show.legend = TRUE, dot.size = 3.5, 
                                  value.size = 4.5, line.size = 1.25,
                                  axis.labels = c("FTLE (240 Hrs)","FTLE (120 Hrs)",
                                                  "FTLE (48 Hrs)","FTLE (24 Hrs)")) +
  scale_y_continuous(breaks=c(-25,0,25), limits = c(-27.5,27.5)) +
  theme_minimal() +
  theme(axis.title = element_text(face="bold", size=30),
        axis.text = element_text( face="bold", size=30),
        plot.title = element_text(face="bold", size=25,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=20,hjust = 0.5),
        plot.caption = element_text( face="bold", size=12,hjust = 0),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.size = unit(1.2,"cm"),
        legend.background = element_rect(fill="white",color="black"),
        legend.text = element_text(face="bold", size=28),
        legend.title = element_text(face="bold", size=30, hjust = 0.5))

p_ctdSigmathetaAll$guides$colour$title <- "Response Variable"
p_ctdSigmathetaAll$data$group <- factor(x = c(rep("Depth at Sigma Theta 26.0",4), 
                                              rep("Depth at Sigma Theta 25.5",4),
                                              rep("Depth at Sigma Theta 25.0",4)),
                                        levels = c("Depth at Sigma Theta 26.0",
                                                   "Depth at Sigma Theta 25.5",
                                                   "Depth at Sigma Theta 25.0") )
p_ctdSigmathetaAll <- p_ctdSigmathetaAll + 
  scale_color_manual(values=c("Depth at Sigma Theta 26.0" =  "#0B775E",
                              "Depth at Sigma Theta 25.5" =  "#46ACC8",
                              "Depth at Sigma Theta 25.0" =  "#440154FF"
  ) 
  )
p_ctdSigmathetaAll
ggsave(sprintf("./Output/SupplementalFigure2.jpg"), plot=p_ctdSigmathetaAll, 
       width=14, height=12, units = "in", dpi=600)


#### Supplemental Figure 3 ####

p_DensityAll <- plot_models(krillDensity_MaxGLMMPQLLogLink24,
                            krillDensity_MaxGLMMPQLLogLink48,
                            krillDensity_MaxGLMMPQLLogLink120,
                            krillDensity_MaxGLMMPQLLogLink240,
                            krillDensity_MNZGLMMPQLLogLink24,
                            krillDensity_MNZGLMMPQLLogLink48,
                            krillDensity_MNZGLMMPQLLogLink120,
                            krillDensity_MNZGLMMPQLLogLink240,
                            krillDensity_gMNZGLMMPQLLogLink24,
                            krillDensity_gMNZGLMMPQLLogLink48,
                            krillDensity_gMNZGLMMPQLLogLink120,
                            krillDensity_gMNZGLMMPQLLogLink240,
                            vline.color = "red",show.p=TRUE,
                            show.values=TRUE, spacing = .75, 
                            #axis.lim = c(0.75,7.5),
                            show.legend = TRUE, dot.size = 6, 
                            value.size = 10, line.size = 2,
                            axis.labels = c("FTLE (240 Hrs)","FTLE (120 Hrs)",
                                            "FTLE (48 Hrs)","FTLE (24 Hrs)")) +
  # ggtitle("Metrics of Krill Density ~ FTLE", 
  # subtitle = "glmmPQL Log-Link with Spatial Autocorrelation Structure")+
  ylim(0.25,10)+
  theme_minimal() +
  theme(axis.title = element_text(face="bold", size=30),
        axis.text = element_text( face="bold", size=30),
        plot.title = element_text(face="bold", size=25,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=20,hjust = 0.5),
        plot.caption = element_text( face="bold", size=12,hjust = 0),
        legend.position = c(0.72,0.5),
        legend.direction = "vertical",
        legend.background = element_rect(fill="white",color="black"),
        legend.text = element_text(face="bold", size=28),
        legend.title = element_text(face="bold", size=30, hjust = 0.5))

p_DensityAll$guides$colour$title <- "Response Variable"
p_DensityAll$data$group <- factor(x = c(rep("Max of Non-Zero Cells",4), 
                                        rep("Mean of Non-Zero Cells",4), 
                                        rep("GeoMean of Non-Zero Cells",4)), 
                                  levels = c( "Max of Non-Zero Cells",
                                              "Mean of Non-Zero Cells",
                                              "GeoMean of Non-Zero Cells") )
p_DensityAll <- p_DensityAll + 
  scale_color_manual(breaks=c("Max of Non-Zero Cells",
                              "Mean of Non-Zero Cells" ,    
                              "GeoMean of Non-Zero Cells"),
                     values=c( "Mean of Non-Zero Cells" =  "#440154FF",    
                               "GeoMean of Non-Zero Cells" =  "#31688EFF",
                               "Max of Non-Zero Cells" =  "#35B779FF"
                     ) )
p_DensityAll

ggsave(sprintf("./Output/SupplementalFigure3.jpg"), width=14, height=14, units = "in", dpi=600) # No gaps


#### Supplemental Figure 4 ####

p_WhaleAll <- plot_models(bw_glmmPQL_24_Hrs,bw_glmmPQL_48_Hrs,bw_glmmPQL_120_Hrs,bw_glmmPQL_240_Hrs,
                          mn_glmmPQL_24_Hrs,mn_glmmPQL_48_Hrs,mn_glmmPQL_120_Hrs,mn_glmmPQL_240_Hrs,
                          vline.color = "red",show.p=TRUE,show.values=TRUE, spacing = .75, 
                          #axis.lim = c(0.75,7.5),
                          show.legend = TRUE, dot.size = 6, value.size = 10, line.size = 2,
                          axis.labels = c("FTLE (240 Hrs)","FTLE (120 Hrs)","FTLE (48 Hrs)","FTLE (24 Hrs)")) +
  scale_y_log10( limits = c(-100,200))+
  theme_minimal() +
  theme(axis.title = element_text(face="bold", size=30),
        axis.text = element_text( face="bold", size=30),
        plot.title = element_text(face="bold", size=25,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=20,hjust = 0.5),
        plot.caption = element_text( face="bold", size=12,hjust = 0),
        legend.position = c(0.8,0.5),
        legend.direction = "vertical",
        legend.background = element_rect(fill="white",color="black"),
        legend.text = element_text(face="bold", size=28),
        legend.title = element_text(face="bold", size=30, hjust = 0.5))

p_WhaleAll$guides$colour$title <- "Response Variable"

p_WhaleAll$data$group <- factor(x = c(rep("Blue Whales",4), rep("Humpback Whales",4)), levels = c("Humpback Whales" , "Blue Whales") )
p_WhaleAll <- p_WhaleAll + scale_color_manual(values=c( "Humpback Whales" =  "#440154FF",    
                                                        "Blue Whales" =  "#31688EFF"
) )
p_WhaleAll
ggsave("./Output/SupplementalFigure4.jpg", width=14, height=18, units = "in", dpi=600) # No gaps
