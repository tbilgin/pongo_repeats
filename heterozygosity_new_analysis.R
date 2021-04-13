mean(pongo_abelii_heterozygosity$V4)
t.test(pongo_abelii_heterozygosity$V4)
mean(pongo_pygmaeus_pygmaeus_heterozygosity$V4)
t.test(pongo_pygmaeus_pygmaeus_heterozygosity$V4)
mean(pongo_pygmaeus_morio_heterozygosity$V4)
t.test(pongo_pygmaeus_morio_heterozygosity$V4)

heterozygositiesSpecies = list()
heterozygositiesSpecies[["Pongo abelii"]] =  pongo_abelii_heterozygosity$V4
heterozygositiesSpecies[["Pongo pygmaeus pygmaeus"]] = pongo_pygmaeus_pygmaeus_heterozygosity$V4
heterozygositiesSpecies[["Pongo pygmaeus wurmbii"]] = pongo_pygmaeus_morio_heterozygosity$V4

d2 <- data.frame(x = unlist(heterozygositiesSpecies), 
                 grp = rep(names(heterozygositiesSpecies)[1:length(heterozygositiesSpecies)],
                           times = sapply(heterozygositiesSpecies,length)))

p2 <- ggplot(d2, aes(x = grp, y = x)) + 
  geom_violin(aes(fill=grp, face="bold"), width=0.5,
              position= position_dodge(width = .5)) + 
  coord_flip() + 
  geom_boxplot(width=0.03) + 
  #ggtitle("Heterozygosity for wild orangutan species") + 
  xlab("(Sub-)species") + 
  ylab("Expected heterozygosity") +   
  theme(axis.text.x= element_text(face="bold", size=13),
        axis.text.y= element_text(face="bold", size=12),
        plot.title=element_text(size=20, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"), 
        legend.position="none")
plot(p2)



