rm(list=ls())
setwd("C:/Users/michel/Desktop/MASTER/M2_Bioinfo/metageno/TP2")   #a changer





################################################################################
PhyChim = read.csv("Data_PhyChim.csv",header = TRUE,sep = 
           ",")
Abondance = read.csv("Abondance_Roquefort.csv",header = TRUE,sep = 
                         ",")

Abondance=Abondance[-1]  #enleve les 4 premieres lignes de abondance
Abondance=Abondance[-1]
Abondance=Abondance[-1]
Abondance=Abondance[-1]

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


shannon_index <- diversity(x, index = "shannon")                                      #  Indice de shannon
print(shannon_index)
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


# pour le PH
png(file = "reg_lin_mult_HFD.png", width = 800, height = 700)                            #HFD  (humidité ????????????????????)
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$HFD_. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()

# pour le PH
png(file = "reg_lin_mult_lactate.png", width = 800, height = 700)                            #lactate
LM<-apply(DataAb[,5:18],2,function(x){lm(MesurePC$Lactate.g.L.1. ~ x)$coefficients})
## On regarde les coefficients
par(mar=c(15,5,5,5))
barplot(LM[2,],las=2,cex.names = 0.8)
dev.off()


#il n' ya pas de NaCl ???????????????????

