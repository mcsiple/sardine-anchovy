# Plotting code! I just dumped all this code here for now
#library(beyonce)
#source("/Users/mcsiple/Dropbox/ChapterX-synthesis/Theme_Black.R")
#sa.col <- c("red","darkblue")
#ylabel <- c("Landings","SSB","Recruitment")
# pal <- beyonce_palette(11)
# plot(1:3,sp[,2],xaxt='n',ylim=c(0,1),pch=21,bg = pal[c(1,3,5)],ylab="Prob(WMR < 0.5)", xlab="")
# axis(1, at = c(1,2,3), labels = c("<5 yr","5-10 yr","10+ yr"))
# arrows(x0 = 1:3, x1 = 1:3,y0 = sp[,1], y1 = sp[,3],col = pal[c(1,3,5)],lwd = 1.5,length = 0.03,angle=90,code = 3)

# Start here if not running sims again! -----------------------------------
load("NullModelDistributions_07cutoff_std.RData")
load("TrueDF.RData")
true.df$ou <- NA
true.df$ou[which(true.df$obs>df$X95.)] <- "Asynchronous"
true.df$ou[which(true.df$obs<df$X5.)] <- "Synchronous"

df <- df %>% subset(variable != "rec")
true.df <- true.df %>% subset(variable !="rec")
true.df$ou[is.na(true.df$ou)] <- "not.significant" 
pali <- c(beyonce_palette(48)[3],"grey",
          beyonce_palette(48)[6])

pdf("ExpectationPlot_07cutoff_black.pdf",width=8,height=7,useDingbats = FALSE)
ggplot(df, aes(x=scale,y=X50.)) + 
  #geom_point(size=0.5) + 
  facet_grid(region~variable) + 
  scale_x_discrete(limits=c("ten.plus","five.ten","less.than.5"),labels=c("Long-term","Medium","Short")) +
  coord_flip() +
  geom_linerange(aes(x=scale, ymin = X5.,ymax=X95.),lwd=1.1,colour='darkgrey') +
  geom_linerange(aes(x=scale,ymin=X50.,ymax=X75.),lwd=2.5,colour="darkgrey") +
  #theme_classic(base_size=14) +
  geom_point(data = true.df,
             aes(x=scale,y=obs,colour=ou),size=4) +
  scale_colour_manual("",labels = c("Asynchronous","not significant","Synchronous"),values=pali) +
  theme(strip.background = element_blank()) +
  ylab("Degree of asynchrony") +
  xlab("Time scale") +
  theme_black(base_size=14)
dev.off()

#dev.off()

# 