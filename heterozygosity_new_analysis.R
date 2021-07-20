mean(pongo_abelii_heterozygosity$V4)
t.test(pongo_abelii_heterozygosity$V4)
mean(pongo_pygmaeus_pygmaeus_heterozygosity$V4)
t.test(pongo_pygmaeus_pygmaeus_heterozygosity$V4)
mean(pongo_pygmaeus_morio_heterozygosity$V4)
t.test(pongo_pygmaeus_morio_heterozygosity$V4)

hetPongoAbelii = read.table("~/Desktop/pongo_repeats/alina.voicu/data-from-analyses/heterozygosity-data/het_pongo_abelii_genome_wide_nei.txt", sep = ',',header=FALSE)
hetPongoPygmaeusMorio = read.table("~/Desktop/pongo_repeats/alina.voicu/data-from-analyses/heterozygosity-data/het_pongo_pygmaeus_morio_genome_wide_nei.txt", sep = ',', header=FALSE)
hetPongoPygmaeusPygmaeus = read.table("~/Desktop/pongo_repeats/alina.voicu/data-from-analyses/heterozygosity-data/het_pongo_pygmaeus_pygmaeus_genome_wide_nei.txt", sep = ',', header=FALSE)



heterozygositiesSpecies = list()
heterozygositiesSpecies[["P.abelii"]] =  hetPongoAbelii$V4
heterozygositiesSpecies[["P.p.morio"]] = hetPongoPygmaeusMorio$V4
heterozygositiesSpecies[["P.p.pygmaeus"]] = hetPongoPygmaeusPygmaeus$V4

d2 <- data.frame(x = unlist(heterozygositiesSpecies), 
                 grp = rep(names(heterozygositiesSpecies)[1:length(heterozygositiesSpecies)],
                           times = sapply(heterozygositiesSpecies,length)))

het_plot <- ggplot(d2, aes(x = grp, y = x)) + 
  geom_violin(aes(fill=grp, face="bold"), width=0.5,
              position= position_dodge(width = .5)) + 
  #coord_flip() + 
  #geom_boxplot(width=0.03) + 
  #ggtitle("Heterozygosity for wild orangutan species") + 
  #xlab("(Sub-)species") +
  labs(x=NULL) +
  ylab("Heterozygosity") +   
  theme(axis.text.x= element_text(face="bold", size=10),
        axis.text.y= element_text(face="bold", size=9),
        plot.title=element_text(size=20, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"), 
        legend.position="none")
plot(het_plot)



