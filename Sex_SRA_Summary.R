library(extrafont)

## load colors
  bright= c(red=rgb(222,45,38, maxColorValue=255), #de2d26
            pink=rgb( 255, 102, 153, maxColorValue=255), #ff6699
            orange=rgb(232,121,12, maxColorValue=255),   #e8790c
            yellow=rgb(255,222,13, maxColorValue=255), #ffde0d          
            green=rgb(12,189,24, maxColorValue=255),  #0cbd18           
            teal=rgb(59,196,199, maxColorValue=255), #3bc4c7
            blue=rgb(58,158,234, maxColorValue=255), #3a9eea
            purple=rgb(148,12,232, maxColorValue=255)) #940ce8 

## load predicted phenotypes
load('/dcl01/leek/data/sellis/barcoding/output/PredictedPhenotypes_v0.0.03.rda')
df = PredictedPhenotypes #70479


## load SRA metadata
### Load in SRA metadata
load('/dcl01/leek/data/recount-website/metadata/metadata_sra.Rdata')
metadata <- metadata[!is.na(metadata$bigwig_path), ]
sra_meta = metadata
rm(metadata)
pd = read_csv("https://raw.githubusercontent.com/nellore/runs/master/sra/v2/hg38/SraRunInfo.csv")
sra_meta = left_join(as.data.frame(sra_meta),pd,by=c("run"="Run","sample"="Sample"))
colnames(sra_meta)[4] <- "sample_id"


## combine to look at sex across SRA
meta = left_join(sra_meta,df)

## Sex Breakdown w/n SRA
(sex_SRA <- meta %>%  group_by(predicted_sex) %>% select(predicted_sex) %>% summarise(Count = n())) 
# pdf("plots/Sex_SRA.pdf",width=16,family="Roboto Condensed")

pdf("plots/Sex_SRA_Summary.pdf",width=8, height=8,family="Roboto Condensed")
ggplot(data = sex_SRA, aes(x=predicted_sex, y = Count,label = Count)) + 
	labs(y="No. of Samples",x="Data Set",title="Predicted Sex") +
    geom_bar(stat="identity", aes(fill = predicted_sex),position="dodge") +
    geom_text(aes(fill=predicted_sex),size = 6, position = position_dodge(width = 0.9),colour="black") +
    scale_fill_manual(values=c("#940CE8", "#0CBD18", "grey48"))+
    theme_bw()+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))
dev.off()

(proj <- meta %>%  group_by(ProjectID) %>% select(ProjectID)  %>% summarise(n = n())) 
(sex_summ <- meta %>%  group_by(ProjectID,predicted_sex) %>% select(predicted_sex) %>% summarise(Count = n())) 

df = left_join(sex_summ,proj)
df$prop = df$Count/df$n
df$ProjectID = as.character(df$ProjectID)


#give label as to type (male only, female only, Unassigned only, both male and female, all three, UA+m, UA+fm )
one = df[(df$prop==1),]
both = df[df$prop<1,]
df$type <- rep("unassigned",nrow(df))
df$type[(df$prop==1 & df$predicted_sex=="female")] <- "female only"
df$type[(df$prop==1 & df$predicted_sex=="male")] <- "male only"
#get projects with 2 values
#projects with 3 values 
a <- as.data.frame(table(df$ProjectID))
touse<-dplyr::filter(a, a$Freq==3)
df$type[df$ProjectID %in% touse$Var1] <- "all"

touse<-dplyr::filter(a, a$Freq==2)
df$type[df$ProjectID %in% touse$Var1] <- "both"

##figure out which type of both it is
IDs = unique(touse$Var1) 
fx <- function(x){
	a<- df[df$ProjectID==x,]
	if("female" %in% a$predicted_sex){
		if("male" %in% a$predicted_sex){
			out = "female & male"
		}else{
			out = "female & unassigned"
		}
	}else{
		out = "male & unassigned"
	}
}
labels<-sapply(IDs,fx)
names(labels) <- IDs

for(i in 1:length(labels)){
 df$type[df$ProjectID==names(labels[i])] <- labels[i]
}

levels(df$type)

ggplot(df, aes(x=type, y=prop)) + 
  geom_boxplot(outlier.shape=NA) +
  labs(y="Proportion",x="Study Type",title="Sex within Study Type")+ 
	geom_jitter(position=position_jitter(width=.1, height=0),aes(size=n, colour=predicted_sex))+
  scale_colour_manual(values=c("#940CE8", "#0CBD18", "grey48"))+
  theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))

## sample size in each study type
 pdf("plots/Sex_SRA_SampleSize.pdf",width=16,family="Roboto Condensed")
 ggplot(df, aes(x=type, y=n)) + 
  geom_boxplot(outlier.shape=NA)+
  labs(y="Sample Size",x="Study Type",title="Sample Size within Study Type")+ 
	geom_jitter(position=position_jitter(width=.1, height=0),aes(size=prop, colour=predicted_sex))+
  #geom_point(position=position_dodge(width=0.75),aes(size=prop, colour=predicted_sex))+
  scale_colour_manual(values=c("#940CE8", "#0CBD18", "grey48"))+
  theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))
dev.off()

#only graph studies where there's a mix of sex
 df_only <- df %>% dplyr::filter(type!="unassigned" & type!="female only" & type!="male only")
 
  pdf("plots/Sex_SRA_StudyType.pdf",width=16,family="Roboto Condensed")
 ggplot(df_only, aes(x=type, y=prop, fill=predicted_sex)) + 
  geom_boxplot() +
  labs(y="Proportion",x="Study Type",title="Sex within Study Type")+ 
	geom_point(position=position_dodge(width=0.75),aes(size=n, colour=predicted_sex))+
  scale_fill_manual(values=c("#940CE8", "#0CBD18", "grey48"))+
   scale_colour_manual(values=c("black","black","black"))+
  theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))
dev.off()

#sample sizes for all male or all female studies
(allone = one  %>%  group_by(predicted_sex) %>% select(ProjectID,n) )

table(allone$predicted_sex)


## plot samples that are all male or all female
df_mix = both

pdf("plots/Sex_SRA_ProjectType.pdf",width=8, height=8,family="Roboto Condensed")
summ<-as.data.frame(table(allone$predicted_sex))
summ2 <- cbind("mixed",length(unique(both$ProjectID)))
colnames(summ2) <- colnames(summ)
summ <- rbind(summ,summ2)
summ$ProjectType <- c("female only", "male only", "Unassigned only", "mixed")
summ$Freq <- as.numeric(summ$Freq)
ggplot(data = summ, aes(x=ProjectType, y = Freq,label = Freq)) + 
	labs(y="No. of Projects",x="Project Type",title="Project Type Summary")+ 
	geom_bar(stat="identity", aes(fill = ProjectType),position="dodge") +
    geom_text(aes(fill=ProjectType),size = 6, position = position_dodge(width = 0.9),colour="black") +
    scale_fill_manual(values=c("#940CE8", "#0CBD18","#E8790C", "grey48"))+
    theme_bw()+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))
dev.off()
## number of samples in each type
data <- one %>% group_by(predicted_sex) %>% select(predicted_sex,Count)

ggplot(data = one, aes(x=Count)) + 
	labs(y="No. of Samples",x="Predicted Sex",title="Sample Summary")+ 
	stat_density(aes(group=predicted_sex, colour=predicted_sex),geom="line") +
	 scale_colour_manual(values=c("#940CE8", "#0CBD18","#E8790C", "grey48")) +
   	theme_bw()	+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))


## what about in studies with >=50 samples
## number of samples in each type
data2 <- one %>% dplyr::filter(Count>=50) %>% dplyr::filter(predicted_sex=="male"| predicted_sex=="female") %>% group_by(predicted_sex)  %>% select(predicted_sex, Count) 

ggplot(data = data2, aes(x=Count)) + 
	labs(x="No. of Samples",y="",title="Sample Summary: exclusively male or female (n>=50)")+ 
	stat_density(aes(group=predicted_sex, colour=predicted_sex),geom="line") +
	 scale_colour_manual(values=c("#940CE8", "#0CBD18","#E8790C", "grey48")) +
   	theme_bw()	+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))

## what about when both sexes are included
data3 <- both %>% dplyr::filter(Count>=50) %>% dplyr::filter(predicted_sex=="male"| predicted_sex=="female") %>% group_by(predicted_sex)  %>% select(ProjectID,predicted_sex, Count) 

ggplot(data = data3, aes(x=Count)) + 
	labs(x="No. of Samples",y="",title="Sample Summary : both male and female (n>=50)")+ 
	stat_density(aes(group=predicted_sex, colour=predicted_sex),geom="line") +
	 scale_colour_manual(values=c("#940CE8", "#0CBD18","#E8790C", "grey48")) +
   	theme_bw()	+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))

## plot samples that are all male or all female
df_female = df %>% dplyr::filter(n>=mean(proj$n)) %>% dplyr::filter(predicted_sex=="female")

## Sex by project
ggplot(data = df_female, aes(x=reorder(ProjectID,Count), y = prop,label = n)) + 
	labs(y="Proportion",x="Project",title=paste0("Proportion of study : female (n>=",round(mean(proj$n),2),")")) +
    geom_bar(stat="identity", aes(fill = predicted_sex),position="dodge") +
    geom_text(aes(fill=predicted_sex),size = 4, position = position_dodge(width = 0.9),colour="black") +
    scale_fill_manual(values=c("#940CE8"))+
    theme_bw()+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"),axis.text.x=element_blank(),axis.ticks.x=element_blank())

 ##male
df_male = df %>% dplyr::filter(n>=mean(proj$n)) %>% dplyr::filter(predicted_sex=="male")

## Sex by project
ggplot(data = df_male, aes(x=reorder(ProjectID,Count), y = prop,label = n)) + 
	labs(y="Proportion",x="Project",title=paste0("Proportion of study : male (n>=",round(mean(proj$n),2),")")) +
    geom_bar(stat="identity", aes(fill = predicted_sex),position="dodge") +
    geom_text(aes(fill=predicted_sex),size = 4, position = position_dodge(width = 0.9),colour="black") +
    scale_fill_manual(values=c("#0CBD18"))+
    theme_bw()+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"),axis.text.x=element_blank(),axis.ticks.x=element_blank())

# Proportion female vs. study size (n)
df2 = df %>% dplyr::filter(n>=mean(proj$n)) %>% dplyr::filter(predicted_sex=="female") 

ggplot(data = df2, aes(x=n, y = prop)) +
	labs(y="Proportion",x="Study Size (n)",title="") +
    geom_point(colour="#940CE8",aes(size=n))+
    theme_bw()+
    theme_bw()+
    scale_x_continuous()+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))

# Proportion male vs. study size (n)
df2 = df %>% dplyr::filter(n>=mean(proj$n)) %>% dplyr::filter(predicted_sex=="male") 

ggplot(data = df2, aes(x=n, y = prop)) +
	labs(y="Proportion",x="Study Size (n)",title="") +
    geom_point(colour="#0CBD18",aes(size=n, label=n))+
    theme_bw()+
    scale_x_continuous()+
    theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size=20), panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(colour="black"))
dev.off() 

df2 = df %>% dplyr::filter(n>=mean(proj$n))