rst = read.table("~/Desktop/pongo_repeats/rst", header = TRUE, row.names = 1)
rst_data <- as.data.frame(t(rst))
rst_data = rst_data[,c("P.abelii", "P.p.morio", "P.p.pygmaeus")]
mean(rst_data$P.abelii) # 0.5867372
mean(rst_data$P.p.pygmaeus) # 0.0692665
mean(rst_data$P.p.morio) # 0.2098921

wilcox.test(rst_data$P.abelii, rst_data$P.p.pygmaeus)
wilcox.test(rst_data$P.abelii, rst_data$P.p.morio)
wilcox.test(rst_data$P.p.pygmaeus, rst_data$P.p.morio)

#all less than  2.2e-16

boxplot(rst_data)
install.packages("tidyr")
library(tidyr)
library(ggplot2)
rst_long = gather(rst_data, species, rst, P.abelii:P.p.pygmaeus, factor_key=TRUE)
ggplot(rst_long, aes(y=rst, x=species, color=species)) + geom_violin(aes(fill = rst), show.legend=F) + geom_jitter(height = 0, width = 0.1) + labs(y = "STR dosage") 
ggplot(rst_long, aes(y=rst, x=species)) + geom_boxplot(aes(fill = rst), show.legend=F) 

rst_plot <- ggplot(rst_long, aes(y=rst, x=species)) + 
  geom_violin(aes(fill=species, face="bold"), width=0.5,
              position= position_dodge(width = .5)) + 
  #coord_flip() + 
  #geom_boxplot(width=0.03) + 
  #ggtitle("Heterozygosity for wild orangutan species") + 
  xlab("(Sub-)species") + 
  #labs(x=NULL) +
  ylab("Rst") +   
  theme(axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=9),
        plot.title = element_text(size=20, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"), 
        legend.position="none")

plot(rst_plot)

library(lattice)
library(grid)
library(gridExtra)


grid.arrange(het_plot, rst_plot, heights=c(1,1.2))
