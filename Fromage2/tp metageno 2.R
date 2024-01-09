rm(list=ls())
setwd("C:/Users/michel/Desktop/MASTER/M2_Bioinfo/metageno/Fromage")   #a changer





################################################################################
PhyChim = read.csv("Data_PhyChim.csv",header = TRUE,sep = 
           ",")
Abondance = read.csv("Abondance_Roquefort.csv",header = TRUE,sep = 
                         ",")
Bois = Abondance$Population == "Bois"
Ensilage = Abondance$Population == "Ensilage"
Roquefort = Abondance$Population == "Roquefort"
non_Roquefort = Abondance$Population == "non-Roquefort"
jour9= Abondance$Jour == '9'
jour20= Abondance$Jour == '20'
Fabrication_A=Abondance$Fabrication == "A"
Fabrication_B=Abondance$Fabrication == "B"
Fabrication_C=Abondance$Fabrication == "C"
replicat_1 = Abondance$Replica == "1"
replicat_2 = Abondance$Replica == "2"
replicat_3 = Abondance$Replica == "3"

#enleve les 4 premieres lignes de abondance
Abondance=Abondance[-1];Abondance=Abondance[-1];Abondance=Abondance[-1];Abondance=Abondance[-1]

x=Abondance


################################################################ 1) DIVERSITE ALPHA


install.packages("vegan")
library("vegan")

data(x)
S <- specnumber(x) # observed number of species
(raremax <- min(rowSums(x)))
#> [1] 340
Srare <- rarefy(x, raremax)
png(file = "image.png", width = 800, height = 700)                                    # ,,,?.

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")    
abline(0, 1)
dev.off()


png(file = "courbe-rarefaction.png", width = 800, height = 700)                                    
rarecurve(x, step = 20, sample = raremax, col = "blue", cex = 0.6)                    #  courbe de rarefaction
dev.off()


shannon_index <- diversity(x, index = "shannon") 
                                     #  Indice de shannon
#  Indice de shannon
print(shannon_index)
hist(shannon_index)
barplot(shannon_index)

mean(shannon_index)
mean(shannon_index[Bois])
mean(shannon_index[Ensilage])
mean(shannon_index[Roquefort])
mean(shannon_index[non_Roquefort])
mean(shannon_index[jour9])
mean(shannon_index[jour20])
mean(shannon_index[replicat_1])
mean(shannon_index[replicat_2])
mean(shannon_index[replicat_3])
mean(shannon_index[Fabrication_A])
mean(shannon_index[Fabrication_B])
mean(shannon_index[Fabrication_C])



################################################################ 2) ACP


#acp_resultat=prcomp(Abondance)
#summary(acp_resultat)
#biplot(acp_resultat)

install.packages("FactoMineR")
library("FactoMineR")
png(file = "ACP_var.png", width = 800, height = 700)                                    #ACP variables
abondance.pca = PCA(Abondance)
dev.off()

png(file = "ACP_individu.png", width = 800, height = 700)                                   #ACP individus
plot(abondance.pca)
dev.off()

summary(abondance.pca)

png(file = "ACP_barplot.png", width = 800, height = 700)                            #Barplot acp represente les dimensions contenant le plus d'infos
barplot(abondance.pca$eig[,1],ylim=c(0,5))
dev.off()

############################################################## 3) ANALYSES DISCRIMINANTES 

DataAb = read.csv("Abondance_Roquefort.csv",header = TRUE,sep = 
                 ",")
MesurePC = read.csv("Data_PhyChim.csv",header = TRUE,sep = 
                     ",")


library(MASS)
LDA<-lda(x=DataAb[,5:18],grouping=DataAb$Population)
png(file = "barplot_analyse_discriminante_coeff_population.png", width = 800, height = 700)                            #coefficient represente taxons qui discirminent le + les fromages
par(mar=c(15,5,5,5))
barplot(LDA$scaling[,1],names.arg = rownames(LDA$scaling),las=2,
        col='black',cex.names = 0.76)
dev.off()

##analyse discriminante pour le parametre "Population" (uniquement "Bois" et "ensilage")

png(file = "analyse_discriminante_ensilage_bois.png", width = 800, height = 700)                            #Population  (Bois Ensilage )
LDA1_Ensilage<-colSums(apply(subset(DataAb,Population=='Ensilage')[,5:18],1,
                             function(x){LDA$scaling[,1]*x}))
LDA1_Bois<-colSums(apply(subset(DataAb,Population=='Bois')[,5:18],1,
                         function(x){LDA$scaling[,1]*x}))
#Représentation sous forme d'histogramme
hist(LDA1_Bois,col='blue',xlab='LD1 means',ylab='Frequency',main='Predicted Values',
     xlim=c(-5,5),breaks=10)
hist(LDA1_Ensilage,col=rgb(1,0,0,0.5),add=T,breaks=10)
legend('topleft',pch=19,col=c('blue','red'),legend = c('Bois','Ensilage'))
dev.off()
##



LDA1_Ensilage<-colSums(apply(subset(DataAb,pH=='Ensilage')[,5:18],1,
                             function(x){LDA$scaling[,1]*x}))

#Représentation sous forme d'histogramme
hist(LDA1_Bois,col='blue',xlab='LD1 means',ylab='Frequency',main='Predicted Values',
     xlim=c(-5,5),breaks=10)
hist(LDA1_Ensilage,col=rgb(1,0,0,0.5),add=T,breaks=10)
legend('topleft',pch=19,col=c('blue','red'),legend = c('Bois','Ensilage'))






##analyse discriminante pour le parametre "Population" (uniquement "Bois" et "ensilage")

png(file = "analyse_discriminante_Population.png", width = 800, height = 700)                            #Population  
LDA1_Ensilage<-colSums(apply(subset(DataAb,Population=='Ensilage')[,5:18],1,
                             function(x){LDA$scaling[,1]*x}))
LDA1_Bois<-colSums(apply(subset(DataAb,Population=='Bois')[,5:18],1,
                         function(x){LDA$scaling[,1]*x}))
LDA1_Non-Roquefort<-colSums(apply(subset(DataAb,Population=='Non-Roquefort')[,5:18],1,
                         function(x){LDA$scaling[,1]*x}))
LDA1_Roquefort<-colSums(apply(subset(DataAb,Population=='Roquefort')[,5:18],1,
                         function(x){LDA$scaling[,1]*x}))
#Représentation sous forme d'histogramme
hist(LDA1_Roquefort,col='blue',xlab='LD1 means',ylab='Frequency',main='Predicted Values',
     xlim=c(-5,5),breaks=10)
hist(LDA1_Non-Roquefort,col=rgb(1,0,0,0.5),add=T,breaks=10)
hist(LDA1_Bois,col=rgb(1,0,0,0.5),add=T,breaks=10)
hist(LDA1_Ensilage,col=rgb(1,0,0,0.5),add=T,breaks=10)
legend('topleft',pch=19,col=c('blue','red'),legend = c('Roquefort','Non-Roquefort','Bois','Ensilage'))
dev.off()
##





##analyse discriminante pour le parametre "Jour" ( "9" et "20")

png(file = "analyse_discriminante__JOUR.png", width = 800, height = 700)                                   #JOUR (9 20)
LDA<-lda(x=DataAb[,5:18],grouping=DataAb$Jour)
LDA1_1<-colSums(apply(subset(DataAb,Jour=='9')[,5:18],1,
                             function(x){LDA$scaling[,1]*x}))
LDA1_2<-colSums(apply(subset(DataAb,Jour=='20')[,5:18],1,
                         function(x){LDA$scaling[,1]*x}))
#Représentation sous forme d'histogramme
hist(LDA1_2,col='blue',xlab='LD1 means',ylab='Frequency',main='Predicted Values',
     xlim=c(-5,5),breaks=10)
hist(LDA1_1,col=rgb(1,0,0,0.5),add=T,breaks=10)
legend('topleft',pch=19,col=c('blue','red'),legend = c('20','9'))
dev.off()
##

png(file = "barplot_analyse_discriminante_coeff_jour.png", width = 800, height = 700)                            #coefficient represente taxons qui discirminent le + les fromages
par(mar=c(15,5,5,5))
barplot(LDA$scaling[,1],names.arg = rownames(LDA$scaling),las=2,
        col='black',cex.names = 0.76)
dev.off()



##analyse discriminante pour le parametre "Fabrication" ( "A" et "B" et "C")


png(file = "analyse_discriminante__fabrication.png", width = 800, height = 700)                            #FABRICATION (A B C)
LDA<-lda(x=DataAb[,5:18],grouping=DataAb$Fabrication)
LDA1_1<-colSums(apply(subset(DataAb,Fabrication=='A')[,5:18],1,
                      function(x){LDA$scaling[,1]*x}))
LDA1_2<-colSums(apply(subset(DataAb,Fabrication=='B')[,5:18],1,
                      function(x){LDA$scaling[,1]*x}))
LDA1_3<-colSums(apply(subset(DataAb,Fabrication=='C')[,5:18],1,
                      function(x){LDA$scaling[,1]*x}))
#Représentation sous forme d'histogramme
hist(LDA1_3,col='blue',xlab='LD1 means',ylab='Frequency',main='Predicted Values',
     xlim=c(-5,5),breaks=10)
hist(LDA1_2,col=rgb(1,0,0,0.5),add=T,breaks=10)
hist(LDA1_1,col=rgb(1,0,0,0.5),add=T,breaks=10)

legend('topleft',pch=19,col=c('blue','red'),legend = c('C','B','A'))
dev.off()

png(file = "barplot_analyse_discriminante_coeff_fabrication.png", width = 800, height = 700)                            #coefficient represente taxons qui discirminent le + les fromages
par(mar=c(15,5,5,5))
barplot(LDA$scaling[,1],names.arg = rownames(LDA$scaling),las=2,
        col='black',cex.names = 0.76)
dev.off()


################################################################### 4) COEFF DISCRIMINENT FROMAGES
#Quels taxons discriminent le + les fromages

png(file = "barplot_analyse_discriminante_coeff.png", width = 800, height = 700)                            #coefficient represente taxons qui discirminent le + les fromages
par(mar=c(15,5,5,5))
barplot(LDA$scaling[,1],names.arg = rownames(LDA$scaling),las=2,
        col='black',cex.names = 0.76)
dev.off()









################################################################ 5) REGRESSION LINEAIRE MULTIPLE
# pour le PH
png(file = "reg_lin_mult_PH.png", width = 800, height = 700)                            #PH
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$pH ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()


# pour humidite
png(file = "reg_lin_mult_HFD.png", width = 800, height = 700)                            #HFD  (humidité ????????????????????)
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$HFD_. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()

# pour lelactate
png(file = "reg_lin_mult_lactate.png", width = 800, height = 700)                            #lactate
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Lactate.g.L.1. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()



# pour jour
png(file = "reg_lin_mult_jour.png", width = 800, height = 700)                            #jour
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Jour ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()

png(file = "reg_lin_mult_population.png", width = 800, height = 700)                            #population
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Population ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()

png(file = "reg_lin_mult_fabrication.png", width = 800, height = 700)                            #lactate
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Fabrication ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()



png(file = "reg_lin_mult_extrait_sec.png", width = 800, height = 700)                            #extrait_sec
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Extrait_Sec._. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()


png(file = "reg_lin_mult_mat_grasse.png", width = 800, height = 700)                            #mat_grasse
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Mat_Grasse.g.L.1. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()

png(file = "reg_lin_mult_gras_sec.png", width = 800, height = 700)                            #gras_sec
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Gras.Sec_. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()


png(file = "reg_lin_mult_acetate.png", width = 800, height = 700)                            #acetate
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Acetate.g.L.1. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()