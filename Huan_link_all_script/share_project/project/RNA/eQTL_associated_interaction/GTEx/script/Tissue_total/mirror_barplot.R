library(ggplot2)
library(patchwork)

df <- data.frame(x = paste0("hhhhhh",1:10),
                 y1 = abs(rnorm(10)),
                 y2 = -abs(rnorm(10)))

p1 = ggplot(df) + 
  geom_bar(aes(x = x, y = y1), 
           stat = "identity",  width = 0.4) + 
  labs(x = "", y = "", title = "aaaa") + 
  coord_flip() + 
  scale_y_continuous(limits = c(0,2),breaks = c(0,1,2),labels = c(0,1,2))+
  theme_bw() 

p2 = ggplot(df) + 
  geom_bar(aes(x = x, y = y2), 
           stat = "identity",  width = 0.4) + 
  labs(x = "", y = "", title = "bbbb") + 
  coord_flip() + 
  scale_y_continuous(limits = c(-2,0),breaks = c(-2,-1,0),labels = c(2,1,0))+
  theme_bw()

p2 + p1

p3 = p1 + theme(panel.background = element_rect(fill='transparent', color="white"),
                panel.grid.minor=element_line(colour='transparent', color='white'),
                panel.grid.major=element_line(colour='transparent', color='white'), 
                panel.border = element_rect(fill='transparent', color='white'),
                strip.background = element_rect(color = "black"),
                strip.text = element_text(size = 10),
                plot.title = element_text(size = 12.5, hjust = 0.5), # title
                plot.subtitle = element_text(size = 11.5, hjust = 0), # subtitle
                legend.key = element_rect( fill = "white"),
                axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'), # face = "bold"
                axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'), # face = "bold"
                axis.ticks = element_blank(), 
                axis.text.y=element_blank(),
                axis.line.x = element_line(colour = "black"),
                legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                legend.key.size = unit(0.9, "cm") ) 

p4 = p2 + theme(panel.background = element_rect(fill='transparent', color="white"),
                panel.grid.minor=element_line(colour='transparent', color='white'),
                panel.grid.major=element_line(colour='transparent', color='white'), 
                panel.border = element_rect(fill='transparent', color='white'),
                strip.background = element_rect(color = "black"),
                strip.text = element_text(size = 10),
                plot.title = element_text(size = 12.5, hjust = 0.5), # title
                plot.subtitle = element_text(size = 11.5, hjust = 0), # subtitle
                legend.key = element_rect( fill = "white"),
                axis.ticks = element_blank(), 
                axis.line.x = element_line(colour = "black"),
                axis.line.y = element_line(colour = "black"),
                legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                legend.key.size = unit(0.9, "cm") ) 

p4 + p3 + plot_layout(nrow = 1, width = c(1,1)) 
