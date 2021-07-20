mydata=read.table("dos_var.data", header = TRUE)
attach(mydata)


mean(NumberVariants[species == "P.p.pygmaeus"])
mean(NumberVariants[species == "P.abelii"])
mean(NumberVariants[species == "P.p.morio"])

cor.test(NumberVariants,Depth, type="Spearman")
cor.test(Depth,Dosage, type="Spearman")
cor.test(NumberVariants,Dosage, type="Spearman")

##### Figures
library(gridExtra)
install.packages("gridExtra")
library(ggplot2)

p1 = ggplot(mydata, aes(x=Depth, y=NumberVariants, color=species)) + geom_point(size = 3, show.legend=F) + labs(x="Genome Depth", y = "Number of genotyped STR variants")
p2=ggplot(mydata, aes(y=Dosage, x=Depth, color=species)) + geom_point(size = 3, show.legend=F) + labs(x="Genome Depth", y = "STR dosage")
ggplot(mydata, aes(y=Dosage, x=NumberVariants, color=species)) + geom_point(size = 3) + labs(x="Number of genotyped STR variants", y = "STR dosage")

grid.arrange(p1, p2,p3, nrow = 1)

ggplot(mydata, aes(y=Dosage, x=species)) + geom_violin(aes(fill = species), show.legend=F) + geom_jitter(height = 0, width = 0.1) + labs(y = "STR dosage", x = "(sub-)species") + scale_x_discrete(limits=c("P.abelii","P.p.pygmaeus","P.p.morio","P.tapanuliensis"))
ggplot(mydata, aes(y=NumberVariants, x=species)) + geom_violin(aes(fill = species),show.legend=F) + geom_jitter(height = 0, width = 0.1) + labs(y = "STR dosage", x = "populations") + scale_x_discrete(limits=c("P.abelii","P.p.pygmaeus","P.p.morio","P.tapanuliensis"))
mean(Dosage[species=="P.abelii"])
stderr(Dosage[species=="P.abelii"])
mean(Dosage[species=="P.p.morio"])
mean(Dosage[species=="P.p.pygmaeus"])
Dosage[species=="P.tapanuliensis"]

###### Exons
install.packages('reshape')
library(reshape)
install.packages('reshape2')
library(reshape2)
exon_data=read.table("ExonicVariationsCount.txt", sep=',', header = TRUE)
exon_NumberVariants=exon_data[,2]
ggplot(mydata, aes(x=Depth, y=exon_NumberVariants, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "Number of genotyped STR variants")


heatmap.exons=read.table("../../data-from-analyses/exons/exons-genome-wide-corr/HeatmapDataExons.txt", header = TRUE, row.names = 1)
ordered_exons=heatmap.exons[as.vector(sort(apply(data.matrix(heatmap.exons),1,sum),index.return=TRUE)$ix),]
melted_heatmap.exons=melt(ordered_exons)
attach(melted_heatmap.exons)
ordered_positions=rep(as.vector(sort(apply(data.matrix(heatmap.exons),1,sum),index.return=TRUE)$ix),28)
ggplot(melted_heatmap.exons, aes(x=reorder(Coord, ordered_positions), y=variable, fill=value)) +geom_tile()


exon_Dosage=as.data.frame(read.table("exons/ExonDosage.heatmap_true"))
ggplot(exon_Dosage, aes(x=V1, y=V2, fill=V3)) +geom_tile() 
exon_Dosage_Matrix=dcast(exon_Dosage, V1 ~ V2)
exon_Dosage_Matrix2 <- exon_Dosage_Matrix[,-1]
rownames(exon_Dosage_Matrix2) <- exon_Dosage_Matrix[,1]
exon_Dosage_Matrix<-exon_Dosage_Matrix2
which.max(apply(data.matrix(exon_Dosage_Matrix),2, sum,na.rm=TRUE))
exon.dosage=apply(data.matrix(exon_Dosage_Matrix),1, mean,na.rm=TRUE)

mean(exon.dosage[species=="P.abelii"])
mean(exon.dosage[species=="P.p.morio"])
mean(exon.dosage[species=="P.p.pygmaeus"])
exon.dosage[species=="P.tapanuliensis"]

cor.test(exon_NumberVariants,Depth, type="Spearman")
cor.test(Depth,exon.dosage, type="Spearman")
cor.test(exon_NumberVariants,exon.dosage, type="Spearman")

ggplot(mydata, aes(x=Depth, y=exon.nv, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "Dosage")

mydata_exon=data.frame(exon.dosage,mydata$species, exon_NumberVariants)
ggplot(mydata_exon, aes(x=Depth, y=exon_NumberVariants, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "Number of genotyped STR variants")
ggplot(mydata_exon, aes(x=Depth, y=exon.dosage, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "STR dosage")
ggplot(mydata_exon, aes(y=exon.dosage, x=mydata.species)) + geom_violin(aes(fill = mydata.species), show.legend=F) + geom_jitter(height = 0, width = 0.1) + labs(y = "STR dosage") + scale_x_discrete(limits=c("P.abelii","P.p.pygmaeus","P.p.morio","P.tapanuliensis"))
ggplot(mydata_exon, aes(y=exon_NumberVariants, x=mydata.species)) + geom_violin(aes(fill = mydata.species), show.legend=F) + geom_jitter(height = 0, width = 0.1) + labs(y = "Number of genotyped STR variants") + scale_x_discrete(limits=c("P.abelii","P.p.pygmaeus","P.p.morio","P.tapanuliensis"))


#### promoters

promoter.variant.data=read.table("../Promoter_Variants", sep=","); promoter_NumberVariants=promoter.variant.data$V2;
promoter.dosage.data=read.table("../Promoter_Dosage", sep=",");promoter.dosage=promoter.dosage.data$V2;

mydata_promoter=data.frame(mydata$species,promoter_NumberVariants, promoter.dosage)
ggplot(mydata_promoter, aes(x=Depth, y=promoter_NumberVariants, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "Number of genotyped STR variants")
ggplot(mydata_promoter, aes(x=Depth, y=promoter.dosage, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "STR dosage")
ggplot(mydata_promoter, aes(x=promoter_NumberVariants, y=promoter.dosage, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "STR dosage")
ggplot(mydata_promoter, aes(y=promoter.dosage, x=mydata.species)) + geom_violin(aes(fill = mydata.species), show.legend=F) + geom_jitter(height = 0, width = 0.1) + labs(y = "STR dosage") + scale_x_discrete(limits=c("P.abelii","P.p.pygmaeus","P.p.morio","P.tapanuliensis"))


### REGIONAL COMPARISON

boxplot(Dosage, promoter.dosage, exon.dosage, notch= TRUE)
wilcox.test(promoter.dosage,exon.dosage)

#genomeDosage=read.table("../genomicSTRdosage.all"); colnames(genomeDosage)=c("sample", "Dosage")
#ggplot(genomeDosage, aes(y=V2, x=V1)) + geom_violin(show.legend=F) 

#genomeDosage_sample=merge(genomeDosage, mydata, by = "sample")
#ggplot(genomeDosage_sample, aes(x=species, y=Dosage.x)) + geom_violin(show.legend=F,outlier.shape = NA) + coord_cartesian(ylim = c(0, 12))

##### notes on bash codes:

###awk '{print $1}' dos_var.data | egrep -v '#' | while read l; do grep $l ../data-from-analyses/variants_10kb/MeanAbsDosageVariants10kb.txt; done > Promoter_Dosage
##awk '{print $1}' ../../../data-from-analyses/exons/exons-genome-wide-corr/HeatmapDataExons.txt | while read l; do grep $l ../../../data-from-analyses/abs-dosage-genome-wide-results/* | cut -d \/ -f 6; done > temp.file
##cat temp.file | while read l; do genome=$(echo $l | cut -d _ -f 1,2); pos=$(echo $l | cut -d : -f 2 | cut -d , -f 1,2); dosage=$(echo $l | awk '{print $2}'); echo $genome $pos $dosage; done > ExonDosage.heatmap


#aggregate(genomeDosage$V2 ~ genomeDosage$V1, genomeDosage, mean)



# STR variants in regions
require(gridExtra)
ggplot(mydata, aes(x=Depth, y=NumberVariants, color=species)) + geom_point(size = 3, show.legend=F) + labs(x="Genome Depth", y = "Number of all STR variants") 
ggplot(mydata_exon, aes(x=Depth, y=exon_NumberVariants, color=species)) + geom_point(size = 3, show.legend=F) + labs(x="Genome Depth", y = "Number of STR variants in exons")
ggplot(mydata_promoter, aes(x=Depth, y=promoter_NumberVariants, color=species)) + geom_point(size = 3) + labs(x="Genome Depth", y = "Number of STR variants in promoters")

grid.arrange(p1, p2, p3, ncol=3)



