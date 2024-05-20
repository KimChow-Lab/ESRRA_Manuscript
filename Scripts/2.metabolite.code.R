library(plyr)
library(dplyr)
library(pheatmap)
setwd("D:/Alzheimer/ROSMAP/syn26401311")

Clinical=read.csv("../ROSMAP_clinical.csv",header=T,row.names=1)
medcor=read.table("metabolomics/rosmap_metabolomics_processed_medcor_data.txt",header=T,row.names=1,check.names=F)

medcor.name=intersect(colnames(medcor),rownames(Clinical))
length(medcor.name)
#500


# Status defination 
medcor.clinical=Clinical[medcor.name,]
medcor.clinical.cluster=medcor.clinical[,c("braaksc","ceradsc","cogdx","dcfdx_lv")]
medcor.clinical.cluster$ceradsc_Modified=5-medcor.clinical.cluster$ceradsc
medcor.clinical.cluster$ceradsc=NULL

t=pheatmap(medcor.clinical.cluster,clustering_method="ward.D2",show_rownames=F)
sampleGroup=cutree(t$tree_row,k=2)
table(sampleGroup)
#  1   2 
#277 223 

medcor.clinical.cov=medcor.clinical[,c("msex","educ","apoe_genotype","age_at_visit_max","pmi")]
medcor.clinical.cov$Gender=ifelse(medcor.clinical.cov$msex==0, "Female", "Male")
medcor.clinical.cov$Age=as.numeric(ifelse(medcor.clinical.cov$age_at_visit_max=="90+","90",medcor.clinical.cov$age_at_visit_max))
medcor.clinical.cov$Apoe=as.character(medcor.clinical.cov$apoe_genotype)
medcor.clinical.cov$msex=NULL
medcor.clinical.cov$age_at_visit_max=NULL
medcor.clinical.cov$apoe_genotype=NULL
medcor.clinical.cov$Status=ifelse(sampleGroup==1,"ND","AD")
ApoeColor=brewer.pal(6, "BuPu")
names(ApoeColor)=sort(unique(medcor.clinical.cov$Apoe))
ann_colors = list(
	Gender = c(Female = "Magenta", Male = "BlueViolet"),
    Age = c("white", "DarkOrange"),
    educ = c("white", "Gold"),
    pmi = c("white", "SandyBrown"),
    Apoe = ApoeColor,
    Status = c(AD ="OrangeRed", ND="SlateBlue")
)
pdf("metabolomics/Heatmap_Status.pdf",height=7,width=6)
pheatmap(medcor.clinical.cluster,clustering_method="ward.D2",color=inferno(10),show_rownames=F,annotation_row=medcor.clinical.cov,annotation_colors = ann_colors)
dev.off()


all(rownames(medcor.clinical.cov)==rownames(medcor.clinical.cluster))
write.table(cbind(medcor.clinical.cov,medcor.clinical.cluster),file="GLM/StatusInfo.txt",sep="\t",quote=F)
all(colnames(medcor)==rownames(medcor.clinical.cov))
medcor.clinical.cov[is.na(medcor.clinical.cov)]=0
medcor.adjust= limma::removeBatchEffect(medcor,covariates=model.matrix(~educ+pmi+Age,data=medcor.clinical.cov),design = model.matrix(~Status+Gender,data=medcor.clinical.cov))
write.table(format(medcor.adjust,digits=3),"D:/Alzheimer/ROSMAP/syn26401311/GLM/rosmap_metabolomics_processed_medcor_data.adjusted.txt",sep="\t",row.names=TRUE,quote=F) # vst transformed and adjusted expression matrix



#### DEM identification ###########
get_effect_pval_Status <- function(metab_id, model_data) {
  readings <- cbind(model_data[,1:ncol(medcor.clinical.cov)], model_data[,metab_id]) #ncol(clinical.female)==ncol(clinical.male)
  colnames(readings)[ncol(readings)] <- "reading"  
  readings <- readings[!is.na(readings$reading),] #remove the sample with missing value
  model <- glm(reading ~ Status+educ+pmi+Age,data = readings, family = gaussian)
  data.frame(
    metab_id = metab_id,
    effect = as.matrix(summary(model)$coefficients)["StatusAD","Estimate"],
    pval = as.matrix(summary(model)$coefficients)["StatusAD","Pr(>|t|)"]
  )
}


medcor.clinical.cov.female=medcor.clinical.cov[medcor.clinical.cov$Gender=="Female",]
medcor.data.female=medcor[,rownames(medcor.clinical.cov.female)]
all(colnames(medcor.data.female)==rownames(medcor.clinical.cov.female))
medcor.data.female.clinical.integrated=cbind(medcor.clinical.cov.female,t(medcor.data.female))
medcor.data.female.clinical.integrated$Status=factor(medcor.data.female.clinical.integrated$Status,levels=c("ND","AD"))
# Create data frame with effects, pvals, and padj for each metabolite
Status_effect_pval_df_female <- ldply(colnames(medcor.data.female.clinical.integrated[,-1:-ncol(medcor.clinical.cov)]), 
                             get_effect_pval_Status,
                             model_data = medcor.data.female.clinical.integrated)
Status_effect_pval_df_female$padj <- p.adjust(Status_effect_pval_df_female$pval, method = "BH")

medcor.clinical.cov.male=medcor.clinical.cov[medcor.clinical.cov$Gender=="Male",]
medcor.data.male=medcor[,rownames(medcor.clinical.cov.male)]
all(colnames(medcor.data.male)==rownames(medcor.clinical.cov.male))
medcor.data.male.clinical.integrated=cbind(medcor.clinical.cov.male,t(medcor.data.male))
medcor.data.male.clinical.integrated$Status=factor(medcor.data.male.clinical.integrated$Status,levels=c("ND","AD"))
# Create data frame with effects, pvals, and padj for each metabolite
Status_effect_pval_df_male <- ldply(colnames(medcor.data.male.clinical.integrated[,-1:-ncol(medcor.clinical.cov)]), 
                             get_effect_pval_Status,
                             model_data = medcor.data.male.clinical.integrated)
Status_effect_pval_df_male$padj <- p.adjust(Status_effect_pval_df_male$pval, method = "BH")

all(Status_effect_pval_df_female$metab_id==Status_effect_pval_df_male$metab_id)
Status_effect_pval_df=cbind(Gender="Female",Status_effect_pval_df_female,Gender="Male",Status_effect_pval_df_male)

write.table(Status_effect_pval_df,file="GLM/StatusScore.txt",sep="\t",row.names=F,quote=F)



#### volcanos for Female DEM ###########
t=ggplot(data = hs_data, aes(x = effect, y = -log10(pval), size=-log10(pval),colour=Pattern, label =name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.01)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-value)",title="DEM In Female") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(hs_data, pval<0.01),
    aes(label = name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
    )
pdf("D:/Alzheimer/ROSMAP/syn26401311/GLM/Female_Status_DEM.pdf",height=7,width=10)
print(t)
dev.off()


#### volcanos for Male DEM ###########
Metabolic_ROSMAP_Male=StatusScoreFromROSMAP[,c(6:10)]
Metabolic_ROSMAP_Male$Pattern=ifelse(Metabolic_ROSMAP_Male$pval>0.01,"NotSig",ifelse(Metabolic_ROSMAP_Male$effect>0,"UpInAD_Male_ROSMAP","DnInAD_Male_ROSMAP"))
hs_data=Metabolic_ROSMAP_Male
hs_data=merge(hs_data,metabolite_info,by.x="metab_id",by.y="metabolite")
write.table(hs_data,file="D:/Alzheimer/ROSMAP/syn26401311/GLM/Male_Status_DEM.info.txt",sep="\t",quote=F,row.names=F)

t=ggplot(data = hs_data, aes(x = effect, y = -log10(pval), size=-log10(pval),colour=Pattern, label =name)) +
  geom_point(alpha=0.4) +
  theme_bw() + 
  #scale_size_manual(values=c(1,rep(3,targetMeSetType)))+
  #scale_color_manual(values=c("black","red")) +
  #scale_color_brewer(palette="Dark2")+
  scale_color_manual(values=c("Navy","LightGrey","Firebrick3"))+
  xlim(c(-1, 1)) + ylim(0,12.5)+
  geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = c(-log10(0.01),-log10(0.05)),lty=4,col="red",lwd=0.8) +
  labs(x="effect",y="-log10 (p-value)",title="DEM In Male") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right") +  
  ggrepel::geom_text_repel(
    data = subset(hs_data, pval<0.01),
    aes(label = name),
    size = 3,
    max.overlaps=10,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
    )
pdf("D:/Alzheimer/ROSMAP/syn26401311/GLM/Male_Status_DEM.pdf",height=7,width=10)
print(t)
dev.off()




#####Expression visualization for target metabolic #######
Clinical=read.csv("../ROSMAP_clinical.csv",header=T,row.names=1)
medcor=read.table("metabolomics/rosmap_metabolomics_processed_medcor_data.txt",header=T,row.names=1,check.names=F)
projectId=intersect(rownames(Clinical),colnames(medcor))
length(projectId)
#500
Clinical.df=Clinical[projectId,]
medcor.df=medcor[,projectId]
Clinical.df$age=ifelse(Clinical.df$age_at_visit_max=="90+",90,Clinical.df$age_at_visit_max)
Clinical.df$Gender=ifelse(Clinical.df$msex==0,"Female","Male")
all(rownames(Clinical.df)==colnames(medcor.df))

medcor.clinical=cbind(Clinical.df,t(medcor.df))
medcor.clinical[1:6,1:20]


medcor.clinical.cogdx=medcor.clinical[medcor.clinical$cogdx %in% c(1,2,4,5),]
target=`N.acetyl.aspartyl.glutamate..NAAG.`

g=ggplot(medcor.clinical.cogdx, aes(x=as.character(cogdx), y=`N.acetyl.aspartyl.glutamate..NAAG.`,color=as.character(cogdx)))+
geom_boxplot()+labs(title="",x="Cogdx score",y="Medications corrected values")+ #ylim(c(7,9))+
#ggpubr::stat_compare_means()+
theme_bw()+
scale_color_brewer(palette="Dark2")+
scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+
labs(title="N.acetyl.aspartyl.glutamate..NAAG.")+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),color="black"),axis.text.y = element_text(size=rel(1.0)))+
facet_grid(~Gender)

pdf("D:/Alzheimer/ROSMAP/syn26401311/Visualization/N.acetyl.aspartyl.glutamate..NAAG.pdf",width=6,height=3)
print(g)
dev.off()
