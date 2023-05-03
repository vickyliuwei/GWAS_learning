load("D:\\svsnp\\snp.RData")
library(qqman)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
don <- uvbsnp %>%  group_by(chr) %>% 
  summarise(chr_len=max(position)) %>%   mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%  left_join(uvbsnp, ., by=c("chr"="chr")) %>%
  arrange(chr, position) %>%
  mutate( BPcum=position+tot)
axisdf <- don %>% group_by(chr) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
png("D:\\u_mahattan_7.png",width = 480*5,height = 480)
ggplot(don, aes(x=BPcum, y=-log10(Pvalue))) +
  scale_color_manual(values = rep(c("#44c1f0", "#44c1f0"), 22 ))+
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=2,shape = 19) +
  scale_x_continuous(name = "Chromosome",labels = as.character(axisdf$chr), breaks=axisdf$center,expand = c(0,0)) +
  geom_point(data=subset(don, is_highlight=="yes"), color="#f88421", size=2,shape = 17)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,10),breaks = seq(0,10,5),labels = c(0,5,10) ) +     # remove space between plot area and x axis
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  theme(axis.line = element_line(size=1, colour = "black")) +
  theme(axis.text.x = element_text(size=15,face = "bold")) +
  theme(axis.title.x=element_text(vjust=2, size=20,face = "bold",margin = margin(1,0,0,0,"cm"))) +
  theme(axis.title.y=element_text(vjust=2, size=20,face = "bold")) +
  theme(axis.text.y=element_text(size=15,face = "bold")) +
  geom_hline(yintercept = c(5.47,7.55),linetype = "longdash", size = 1) +
  geom_vline(xintercept = c(32618760,55814404,81923307,104218866),alpha = 0.8,color = "grey") 
dev.off()

  
