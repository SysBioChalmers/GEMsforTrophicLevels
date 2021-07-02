
library(readr)
library(readxl)
library("RColorBrewer")
library("gplots")
library(drc)
library(BSDA)
library(vioplot)
library(lmodel2)
library(data.table)
library('corrplot')
library(plyr)
library(lmodel2)
library("Hmisc")
library("cluster")
library("factoextra")
library("magrittr")
library("NbClust")
library("rjson")
library(nortest)
library(stringr)
library(rvest)
library(data.table)




###collect fish trophic levels - python - run fishbase_crawler.py

fishbase_TL <- read.csv("Fish trophic levels")
fishbase_TL <- fishbase_TL[,c(2,5)]
fishbase_TL <- na.omit(fishbase_TL)

#load fish data from Vaitla_2018:
Vaitla_2018 <- read_csv("Vaitla 2018 troph prot fat.csv")
Vaitla_2018 <- Vaitla_2018[,c(1,8)]
colnames(Vaitla_2018) <- colnames(fishbase_TL)
for(i in 1:nrow(Vaitla_2018)){
  j <- strsplit(as.character(Vaitla_2018[i,"Scientific"]), "_")[[1]]
  k <- paste(j, collapse=" ")
  Vaitla_2018[i,"Scientific"] <- k
}
Vaitla_2018 <- na.omit(Vaitla_2018)

TL_diet_fish <- c(Vaitla_2018$Scientific, fishbase_TL$Scientific)
TL_diet_fish <- as.data.frame(unique(TL_diet_fish))
colnames(TL_diet_fish) <- "Scientific"
TL_diet_fish <- as.data.frame(TL_diet_fish[order(TL_diet_fish$Scientific),])
colnames(TL_diet_fish) <- "Scientific"


TL_diet_fish$TL <- NA
for(i in 1:nrow(TL_diet_fish)){
  j <- as.character(TL_diet_fish[i,"Scientific"])
  k <- as.numeric(fishbase_TL[fishbase_TL$Scientific==j,"TL"])
  l <- Vaitla_2018[Vaitla_2018$Scientific==j,"TL"]
  if(nrow(l)>1){l <- mean(l$TL, na.rm=T)}else{l <- as.numeric(l)}
  
  TL_diet_fish[i,"TL"] <- mean(c(k,l), na.rm=T)
}
TL_diet_fish <- TL_diet_fish[2:nrow(TL_diet_fish),]
TL_diet_fish <- TL_diet_fish[TL_diet_fish$TL>2,]
rownames(TL_diet_fish) <- 1:nrow(TL_diet_fish)


BirdFuncDat <- read.delim("EltonTraits_BirdFuncDat.txt", header=T)
MamFuncDat <- read.delim("EltonTraits_MamFuncDat.txt", header=T)


#
#collect species for TL calculations: Diet Certainty= ABC and Body Mass = species level data
TL_diet_bird <- BirdFuncDat[BirdFuncDat$Diet.Certainty=="A" | BirdFuncDat$Diet.Certainty=="B" |
                              BirdFuncDat$Diet.Certainty=="C" , ]
TL_diet_bird <- TL_diet_bird[TL_diet_bird$BodyMass.SpecLevel==1, c(1:19,36)]

TL_diet_mamm <- MamFuncDat[MamFuncDat$Diet.Certainty=="ABC",]
#the diet of homo sapiens is D2?? lol add this manually
TL_diet_mamm <- rbind(TL_diet_mamm, MamFuncDat[830,])
#no - do the human stuff separately
TL_diet_mamm <- TL_diet_mamm[TL_diet_mamm$BodyMass.SpecLevel==1, c(1:13,24)]


#TL calc for birds
TL_diet_bird$Trophic.level <- NA
TL_diet_bird$Trophic.level <- 1 + (TL_diet_bird$Diet.Inv/100*2 + 
                                     (TL_diet_bird$Diet.Vend + TL_diet_bird$Diet.Vect + TL_diet_bird$Diet.Vfish + TL_diet_bird$Diet.Vunk)/100*2.5 +
                                     TL_diet_bird$Diet.Scav/100*2 +
                                     (TL_diet_bird$Diet.Fruit + TL_diet_bird$Diet.Nect + TL_diet_bird$Diet.Seed + TL_diet_bird$Diet.PlantO)/100*1
)

#TL calc for mammals
TL_diet_mamm$Trophic.level <- NA
TL_diet_mamm$Trophic.level <- 1 + (TL_diet_mamm$Diet.Inv/100*2 + 
                                     (TL_diet_mamm$Diet.Vend + TL_diet_mamm$Diet.Vect + TL_diet_mamm$Diet.Vfish + TL_diet_mamm$Diet.Vunk)/100*2.5 +
                                     TL_diet_mamm$Diet.Scav/100*2 +
                                     (TL_diet_mamm$Diet.Fruit + TL_diet_mamm$Diet.Nect + TL_diet_mamm$Diet.Seed + TL_diet_mamm$Diet.PlantO)/100*1
)


#load diet composition
diet_composition <- read_csv("diet_composition.csv")


#split Diet-Inv in TL_diet_mamm and TL_diet_bird: load data from ADW
diet_inv_terrestrial <- read_excel("ADW diet.inv.terrestrial animals report-202105041108.xls")
diet_inv_aquatic <- read_excel("ADW diet.inv.aquatic animals report-202105041144.xls")
diet_inv_aq_add <- read_excel("ADW diet.inv.aquatic addendum report-202105051002.xls")

diet_inv_aquatic <- merge(diet_inv_aquatic, diet_inv_aq_add, all=T)
remove(diet_inv_aq_add)

diet_inv_aquatic$genus <- "genus"
for(i in 1:nrow(diet_inv_aquatic)){
  diet_inv_aquatic[i,"genus"] <- strsplit(as.character(diet_inv_aquatic[i,"Species"]), " ")[[1]][1]
}
test <- diet_inv_aquatic[diet_inv_aquatic$Class=="Aves" | diet_inv_aquatic$Class=="Mammalia",]
diet_inv_aquatic <- as.data.frame(unique(test$genus))
colnames(diet_inv_aquatic) <- "genus"

diet_inv_terrestrial$genus <- "genus"
for(i in 1:nrow(diet_inv_terrestrial)){
  diet_inv_terrestrial[i,"genus"] <- strsplit(as.character(diet_inv_terrestrial[i,"Species"]), " ")[[1]][1]
}
test <- diet_inv_terrestrial[diet_inv_terrestrial$Class=="Aves" | diet_inv_terrestrial$Class=="Mammalia",]
diet_inv_terrestrial <- as.data.frame(unique(test$genus))
colnames(diet_inv_terrestrial) <- "genus"

#split TL_diet_mamm 
x <- TL_diet_mamm[TL_diet_mamm$Diet.Inv > 0,c("Scientific", "Diet.Inv")]
x$genus <- "genus"
x$aq <- ""
x$ter <- ""
x$Diet.Inv.aquatic <- 0
x$Diet.Inv.terrestrial <- x$Diet.Inv
for(i in 1:nrow(x)){
  x[i,"genus"] <- strsplit(x[i,"Scientific"], " ")[[1]][1]
  if(x[i,"genus"] %in% diet_inv_aquatic$genus){
    x[i,"aq"] <- "yes"
  }
  if(x[i,"genus"] %in% diet_inv_terrestrial$genus){
    x[i,"ter"] <- "yes"
  }
  
  if(x[i,"aq"]=="yes" && x[i,"ter"]==""){ 
    x[i,"Diet.Inv.aquatic"] <- x[i,"Diet.Inv"]  
    x[i,"Diet.Inv.terrestrial"] <- 0
  } 
  if(x[i,"aq"]=="yes" && x[i,"ter"]=="yes"){ 
    x[i,"Diet.Inv.aquatic"] <- x[i,"Diet.Inv"]/2
    x[i,"Diet.Inv.terrestrial"] <- x[i,"Diet.Inv"]/2
  } 
}

TL_diet_mamm <- merge(TL_diet_mamm, x[,c("Scientific", "Diet.Inv.aquatic", "Diet.Inv.terrestrial")], all=T)
TL_diet_mamm$Diet.Inv <- NULL
TL_diet_mamm[is.na(TL_diet_mamm)] <- 0


#split TL_diet_bird
x <- TL_diet_bird[TL_diet_bird$Diet.Inv > 0,c("Scientific", "Diet.Inv")]
x$genus <- "genus"
x$aq <- ""
x$ter <- ""
x$Diet.Inv.aquatic <- 0
x$Diet.Inv.terrestrial <- x$Diet.Inv
for(i in 1:nrow(x)){
  x[i,"genus"] <- strsplit(x[i,"Scientific"], " ")[[1]][1]
  if(x[i,"genus"] %in% diet_inv_aquatic$genus){
    x[i,"aq"] <- "yes"
  }
  if(x[i,"genus"] %in% diet_inv_terrestrial$genus){
    x[i,"ter"] <- "yes"
  }
  
  if(x[i,"aq"]=="yes" && x[i,"ter"]==""){ 
    x[i,"Diet.Inv.aquatic"] <- x[i,"Diet.Inv"]  
    x[i,"Diet.Inv.terrestrial"] <- 0
  } 
  if(x[i,"aq"]=="yes" && x[i,"ter"]=="yes"){ 
    x[i,"Diet.Inv.aquatic"] <- x[i,"Diet.Inv"]/2
    x[i,"Diet.Inv.terrestrial"] <- x[i,"Diet.Inv"]/2
  } 
}

TL_diet_bird <- merge(TL_diet_bird, x[,c("Scientific", "Diet.Inv.aquatic", "Diet.Inv.terrestrial")], all=T)
TL_diet_bird$Diet.Inv <- NULL
TL_diet_bird[is.na(TL_diet_bird)] <- 0

TL_diet_bird <- TL_diet_bird[,c(1:9,22,21,10:20)]
TL_diet_mamm <- TL_diet_mamm[,c(1:3,16,15,4:14)]

TL_diet_bird$diet_percent_prot <- 0
TL_diet_bird$diet_percent_carb <- 0
TL_diet_bird$diet_percent_fat <- 0

for(i in 1:11){
  k <- TL_diet_bird[,i+9] * as.numeric(diet_composition[i,2]) /100
  TL_diet_bird$diet_percent_prot <- TL_diet_bird$diet_percent_prot + k
  
  k <- TL_diet_bird[,i+9] * as.numeric(diet_composition[i,3]) /100
  TL_diet_bird$diet_percent_carb <- TL_diet_bird$diet_percent_carb + k
  
  k <- TL_diet_bird[,i+9] * as.numeric(diet_composition[i,4]) /100
  TL_diet_bird$diet_percent_fat <- TL_diet_bird$diet_percent_fat + k
}


TL_diet_mamm$diet_percent_prot <- 0
TL_diet_mamm$diet_percent_carb <- 0
TL_diet_mamm$diet_percent_fat <- 0

for(i in 1:11){
  k <- TL_diet_mamm[,i+3] * as.numeric(diet_composition[i,2]) /100
  TL_diet_mamm$diet_percent_prot <- TL_diet_mamm$diet_percent_prot + k
  
  k <- TL_diet_mamm[,i+3] * as.numeric(diet_composition[i,3]) /100
  TL_diet_mamm$diet_percent_carb <- TL_diet_mamm$diet_percent_carb + k
  
  k <- TL_diet_mamm[,i+3] * as.numeric(diet_composition[i,4]) /100
  TL_diet_mamm$diet_percent_fat <- TL_diet_mamm$diet_percent_fat + k
}


#fish data: 
#estimate diet categories
TL_diet_fish <- TL_diet_fish[order(TL_diet_fish$TL),]
TL_diet_fish$Diet.Algae <- 0
TL_diet_fish$Diet.Inv.aquatic <- 0
TL_diet_fish$Diet.Vfish <- 0

for(i in 1:nrow(TL_diet_fish)){
  if(as.numeric(TL_diet_fish[i,"TL"]) == 2){
    TL_diet_fish[i,"Diet.Algae"] <- 100
  }
  
  if(as.numeric(TL_diet_fish[i,"TL"]) > 2 & as.numeric(TL_diet_fish[i,"TL"]) < 3 ){
    j <- as.numeric(TL_diet_fish[i,"TL"]) - 1
    TL_diet_fish[i,"Diet.Inv.aquatic"] <- (j - 1)*100
    TL_diet_fish[i,"Diet.Algae"] <- (2 - j)*100
  }
  
  if(as.numeric(TL_diet_fish[i,"TL"] == 3)){
    TL_diet_fish[i,"Diet.Inv.aquatic"] <- 100
  }
  
  if(as.numeric(TL_diet_fish[i,"TL"]) > 3 & as.numeric(TL_diet_fish[i,"TL"]) < 4 ){
    j <- as.numeric(TL_diet_fish[i,"TL"]) - 1
    TL_diet_fish[i,"Diet.Vfish"] <- (j - 2)*100
    TL_diet_fish[i,"Diet.Inv.aquatic"] <- (3 - j)*100
  }
  
  if(as.numeric(TL_diet_fish[i,"TL"]) >= 4){
    TL_diet_fish[i,"Diet.Vfish"] <- 100
  }
}

# calculate diet composition 
TL_diet_fish$diet_percent_prot <- 0
TL_diet_fish$diet_percent_carb <- 0
TL_diet_fish$diet_percent_fat <- 0

for(i in 1:3){
  x <- 3:5 #diet columns in TL_diet_fish
  y <- c(12,2,5) #diet rows in diet_compositoin
  
  k <- unlist(TL_diet_fish[,x[i]]) * as.numeric(diet_composition[y[i],2]) /100
  TL_diet_fish$diet_percent_prot <- as.numeric(TL_diet_fish$diet_percent_prot + k)
  
  k <- unlist(TL_diet_fish[,x[i]]) * as.numeric(diet_composition[y[i],3]) /100
  TL_diet_fish$diet_percent_carb <- as.numeric(TL_diet_fish$diet_percent_carb + k)
  
  k <- unlist(TL_diet_fish[,x[i]]) * as.numeric(diet_composition[y[i],4]) /100
  TL_diet_fish$diet_percent_fat <- as.numeric(TL_diet_fish$diet_percent_fat + k)
}

colnames(TL_diet_fish)[2] <- "Trophic.level"


# diet composition plots:
plot(TL_diet_fish$Trophic.level, TL_diet_fish$diet_percent_fat, axes=F, xlab="", ylab="", col="skyblue", pch="f",
     xlim=c(1,5), ylim=c(0,20))
par(new=T)
plot(TL_diet_mamm$Trophic.level, TL_diet_mamm$diet_percent_fat, xlab="Trophic level", ylab="% fat", main="fat in diet", col="orange2", pch="m",
     xlim=c(1,5), ylim=c(0,20))
par(new=T)
plot(TL_diet_bird$Trophic.level, TL_diet_bird$diet_percent_fat,xlab="", ylab="", pch="b",
     xlim=c(1,5), ylim=c(0,20), col="purple", axes=F)


plot(TL_diet_fish$Trophic.level, TL_diet_fish$diet_percent_prot, axes=F, xlab="", ylab="", col="skyblue", pch="f",
     xlim=c(1,5), ylim=c(0,25))
par(new=T)
plot(TL_diet_mamm$Trophic.level, TL_diet_mamm$diet_percent_prot, 
     xlab="Trophic level", ylab="% protein", main="protein in diet", col="orange2", pch="m",
     xlim=c(1,5), ylim=c(0,25))
par(new=T)
plot(TL_diet_bird$Trophic.level, TL_diet_bird$diet_percent_prot, xlab="", ylab="", pch="b",
     xlim=c(1,5), ylim=c(0,25), col="purple", axes=F)


plot(TL_diet_fish$Trophic.level, TL_diet_fish$diet_percent_carb, axes=F, xlab="", ylab="", col="skyblue", pch="f",
     xlim=c(1,5), ylim=c(0,25))
par(new=T)
plot(TL_diet_mamm$Trophic.level, TL_diet_mamm$diet_percent_carb, 
     xlab="Trophic level", ylab="% carbohydrates", main="carbohydrates in diet", col="orange2", pch="m",
     xlim=c(1,5), ylim=c(0,25))
par(new=T)
plot(TL_diet_bird$Trophic.level, TL_diet_bird$diet_percent_carb, xlab="", ylab="", pch="b",
     xlim=c(1,5), ylim=c(0,25), col="purple", axes=F)



speciesFullGEMs_20210407 <- read_excel("speciesFullGEMs_20210407.xlsx")

sppCLPmodels <- merge(speciesFullGEMs_20210407, TL_diet_bird[,c("Scientific","Trophic.level", 
                                                                "diet_percent_prot", "diet_percent_carb","diet_percent_fat")], 
                      by.x="ScientificName", by.y="Scientific")

x <- TL_diet_mamm[,c("Scientific","Trophic.level", 
                      "diet_percent_prot", "diet_percent_carb","diet_percent_fat")]
colnames(x)[1] <- "ScientificName"
y <- merge(speciesFullGEMs_20210407, x)
sppCLPmodels <- rbind(sppCLPmodels, y)



z <- TL_diet_fish[,c("Scientific","Trophic.level", 
                     "diet_percent_prot", "diet_percent_carb","diet_percent_fat")]
colnames(z)[1] <- "ScientificName"
y <- merge(speciesFullGEMs_20210407, z)
#this only pulls out 13 fishes
#from Fishbase many entries only have genus info - go by genus instead

x <- z[grep("spp", z$ScientificName),]
x$ScientificName <- gsub(",", " ", x$ScientificName)
x$ScientificName <- gsub("spp", "", x$ScientificName)
x$ScientificName <- gsub("Ex ", "", x$ScientificName)
x$ScientificName <- gsub("\\(=", " ", x$ScientificName)
x$ScientificName <- gsub("\\(", " ", x$ScientificName)
x$ScientificName <- gsub("\\.", " ", x$ScientificName)
x$ScientificName <- gsub(")", "", x$ScientificName)
x$ScientificName <- gsub("\\s+", " ", x$ScientificName) #one or more spaces --> one space

for(i in 1:nrow(x)){
  j <- trimws(x[i,"ScientificName"])
  if(grepl(" ",j)){
    k <- (strsplit(j, " ")[[1]])
    for(l in 1:length(k)){
      if(length(grep(k[l], speciesFullGEMs_20210407$ScientificName))>0){
        m <- grep(k[l], speciesFullGEMs_20210407$ScientificName)
        for(n in 1:length(m)){
          if(length(grep(speciesFullGEMs_20210407[m[n],"ScientificName"], y$ScientificName))==0){
            
            y[nrow(y)+1,] <- c(unlist(speciesFullGEMs_20210407[m[n],1:4]), 
                               unlist(x[i,2:5]))
          }
        }
      }
    }
  }else{
    if(length(grep(j, speciesFullGEMs_20210407$ScientificName))>0){
      m <- grep(j, speciesFullGEMs_20210407$ScientificName)
      for(n in 1:length(m)){
        if(length(grep(speciesFullGEMs_20210407[m[n], "ScientificName"], y$ScientificName))==0){
          y[nrow(y)+1,] <- c(unlist(speciesFullGEMs_20210407[m[n],1:4]), 
                             unlist(x[i,2:5]))
        }
      }
    }
  }
}


sppCLPmodels <- rbind(sppCLPmodels, y)

#duplicated entries --> take average
test <- as.data.frame(table(sppCLPmodels$ScientificName))
test2 <- as.character(test[test$Freq>1,1])

sppCLPmodels$Trophic.level <- as.numeric(sppCLPmodels$Trophic.level)
sppCLPmodels$diet_percent_prot <- as.numeric(sppCLPmodels$diet_percent_prot)
sppCLPmodels$diet_percent_carb <- as.numeric(sppCLPmodels$diet_percent_carb)
sppCLPmodels$diet_percent_fat <- as.numeric(sppCLPmodels$diet_percent_fat)
for(i in 1:length(test2)){
  j <- grep(test2[i], sppCLPmodels$ScientificName)
  sppCLPmodels[j[1],c("Trophic.level","diet_percent_prot","diet_percent_carb","diet_percent_fat")] <- 
    colMeans(sppCLPmodels[j,c("Trophic.level","diet_percent_prot","diet_percent_carb","diet_percent_fat")])
  sppCLPmodels <- sppCLPmodels[-j[2:length(j)],]
}


test <- as.data.frame(table(speciesFullGEMs_20210407$Category))
test2 <- as.data.frame(table(sppCLPmodels$Category))
colnames(test2)[2] <- "Freq2"
test3 <- merge(test, test2, all=T)

barplot(test3$Freq, las=2, ylim=c(0,50),main="84 full GEMs", names.arg = test3$Var1)
barplot(test3$Freq2, las=2, ylim=c(0,50), names.arg=test3$Var1,
        main=paste(nrow(sppCLPmodels), "with diet data", sep=" "))

remove(test, test2, test3, x,y,z,i,j,k,l,m,n)

sppCLPmodels$pch <- "o"
sppCLPmodels$color <- "black"
for(i in 1:nrow(sppCLPmodels)){
  if(sppCLPmodels[i,"Category"]=="Birds"){sppCLPmodels[i,c("pch", "color")] <- c("b", "purple")}
  if(sppCLPmodels[i,"Category"]=="Mammals"){sppCLPmodels[i,c("pch", "color")] <- c("m", "orange2")}
  if(sppCLPmodels[i,"Category"]=="Fishes"){sppCLPmodels[i,c("pch", "color")] <- c("f", "skyblue")}
  if(sppCLPmodels[i,"Category"]=="Cartilaginous fishes"){sppCLPmodels[i,c("pch", "color")] <- c("f", "skyblue")}
}

plot(sppCLPmodels$Trophic.level, sppCLPmodels$diet_percent_fat, xlim=c(1,5), ylim=c(0,10), main="fat in diet",
     pch=sppCLPmodels$pch, col=sppCLPmodels$color)
plot(sppCLPmodels$Trophic.level, sppCLPmodels$diet_percent_prot, xlim=c(1,5), ylim=c(0,20), main="protein in diet",
     pch=sppCLPmodels$pch, col=sppCLPmodels$color)
plot(sppCLPmodels$Trophic.level, sppCLPmodels$diet_percent_carb, xlim=c(1,5), ylim=c(0,15), main="carbs in diet",
     pch=sppCLPmodels$pch, col=sppCLPmodels$color)


x <- sppCLPmodels[sppCLPmodels$Trophic.level>3 & sppCLPmodels$diet_percent_fat<5 &
                    sppCLPmodels$Category=="Mammals",]



Human_data <- TL_diet_mamm[grep("Homo sapiens", TL_diet_mamm$Scientific),
                           c("Scientific","Trophic.level","diet_percent_prot","diet_percent_carb","diet_percent_fat"  )]
Human_data$taxonID <- 9606
Human_data$CommonName <- "human"
Human_data$Category <- "Mammals"
Human_data$pch <- "m"
Human_data$color <- "orange2"
colnames(Human_data)[1] <- "ScientificName"
Human_data <- merge(Human_data,sppCLPmodels, all.x=T)

spp_GEMs_for_sims <- join(Human_data, sppCLPmodels, type="full")
spp_GEMs_for_sims <- spp_GEMs_for_sims[,colnames(sppCLPmodels)]

x <- spp_GEMs_for_sims[,1:8] #exclude pch and color


write.csv(x, file = "spp_GEMs_for_sims.csv", row.names = FALSE, quote=F)

write.csv(Human_data, file = "Human_data.csv", row.names = FALSE)



### GEM construction and modeling - MATLAB - run CLPmodel_FBA_randsamp_main.m top to bottom


### GEM data analyses & figures

CLP_maxATP <- read_csv("CLPmodel_maxATP_with_AMP_allowance.txt")
colnames(CLP_maxATP)[1] <- "taxonID"
CLP_maxATP <- merge(spp_GEMs_for_sims, CLP_maxATP[,c(1,2,5)])
#remove platypus
CLP_maxATP[grep("Platypus", CLP_maxATP$CommonName),] <- NA
CLP_maxATP <- na.omit(CLP_maxATP)

CLP_maxATP[grep("Whale", CLP_maxATP$CommonName), "color"] <- "skyblue"

Table_S1 <- CLP_maxATP[,c(1:5,7,8,6,11:12,10)]
colnames(Table_S1)[7] <- "diet_percent_lipids"

write.csv(Table_S1, file = "Table_S1.csv", row.names = FALSE)



#Fig 1

plot(Table_S1$Trophic.level, Table_S1$diet_percent_carb, col=Table_S1$color,
     xlim=c(1,5), ylim=c(0,12), cex=1.2, las=1,
     xlab="trophic level", ylab="% (g/g diet)", main="carbohydrates")

plot(Table_S1$Trophic.level, Table_S1$diet_percent_lipids, col=Table_S1$color,
     xlim=c(1,5), ylim=c(0,10), cex=1.2, las=1,
     xlab="trophic level", ylab="% (g/g diet)", main="lipids")

plot(Table_S1$Trophic.level, Table_S1$diet_percent_prot, col=Table_S1$color,
     xlim=c(1,5), ylim=c(0,20), cex=1.2, las=1,
     xlab="trophic level", ylab="% (g/g diet)", main="proteins")


#Fig 2A

plot(Table_S1$Trophic.level, Table_S1$CLPmodel_max_ATP, col=Table_S1$color,
     xlim=c(1,5), ylim=c(0,50), cex=1.2, las=1,
     xlab="trophic level", ylab="max ATP (mmol)")

plot(Table_S1$Trophic.level, Table_S1$CLPmodel_maxATP_AMP10g, pch=21, bg=Table_S1$color, 
     xlim=c(1,5), ylim=c(0,80), cex=1.2, las=1,
     xlab="trophic level", ylab="max ATP (mmol)")
par(new=T)
plot(Table_S1$Trophic.level, Table_S1$CLPmodel_max_ATP, col=Table_S1$color,
     xlim=c(1,5), ylim=c(0,80), cex=1.2, axes=F, xlab="", ylab="")


plot(c(1.3,1.3), c(40,32), col=c("orange2", "skyblue"),
     xlim=c(1,5), ylim=c(0,50), cex=1.2, las=1,
     xlab="trophic level", ylab="max ATP (mmol)")
text(1.7,40, "terrestrial ", adj=0 )
text(1.7,32, "aquatic ", adj=0)


# load random sampling data:
#creat file list
filelist = list.files(pattern = "^randomSample")

library(data.table)
# Read all the files and create a FileName column to store filenames
DT_1 <- rbindlist( sapply(filelist, fread, simplify = FALSE),
                 use.names = TRUE, idcol = "FileName" )
DT_1$randomSample_subsystem <- NULL


# grab all sample names
CLPtaxonID <- unique(DT$FileName)
for(i in 1:length(CLPtaxonID)){
  CLPtaxonID[i] <- strsplit(CLPtaxonID[i],"_")[[1]][4]
  CLPtaxonID[i] <- sub("*.txt$", "", CLPtaxonID[i])
}

#build dataframe
randomSample <- as.data.frame(unique(DT$randomSample_rxns))
colnames(randomSample) <- "randomSample_rxns"
for(i in 1:length(CLPtaxonID)){
  dframe <- (DT[grep(CLPtaxonID[i], DT$FileName), c("randomSample_rxns", "randomSample_flux")])
  randomSample <- merge(randomSample, dframe, all.x=T)
  colnames(randomSample)[i+1] <- paste("mean_flux", CLPtaxonID[i], sep="_") 
}

remove(filelist, DT, DT_1, DT_2)


rSamp_clean <- randomSample[,2:34]
rSamp_clean[rSamp_clean < 1e-6 & rSamp_clean > -1e-6] <- NA
rSamp_clean$randomSample_rxns <- randomSample$randomSample_rxns
  
rSamp_clean$keep <- "yes"
for(i in 1:nrow(rSamp_clean)){
  if(all(is.na(rSamp_clean[i,1:33]))){rSamp_clean[i,"keep"] <- "no"}
}

rSamp_clean <- rSamp_clean[rSamp_clean$keep=="yes",]
randomSample <- rSamp_clean[,c(34, 1:33)]

HumanGEM_1_6_0_rxns <- read_excel("HumanGEM_1_6_0_rxns.xlsx")
randomSample <- merge(HumanGEM_1_6_0_rxns, randomSample, by.x="ID", by.y="randomSample_rxns", all.y=T)


#grab trophic levels

CLPtroph <- CLPtaxonID
for(i in 1:length(CLPtaxonID)){
  CLPtroph[i] <- Table_S1[grep(CLPtaxonID[i], Table_S1$taxonID), "Trophic.level"]
}
CLPtroph <- as.numeric(CLPtroph)


randomSample$Spearman_rho <- NA
randomSample$number_data_points <- NA
for(i in 1:nrow(randomSample)){
  randomSample[i,"number_data_points"] <- length(na.omit(unlist(randomSample[i,5:37])))
  randomSample[i,"Spearman_rho"] <- cor(CLPtroph, unlist(randomSample[i,5:37]), use="complete.obs", method="spearman")
  
}

randomSample <- randomSample[order(randomSample$SUBSYSTEM),]

remove(rSamp_clean)

Table_S2 <- randomSample
write.csv(Table_S2, file = "Table_S2.csv", row.names = FALSE)



randSamp_pos <- Table_S2[Table_S2$Spearman_rho > 0.75 & Table_S2$number_data_points > 17,]
randSamp_neg <- Table_S2[Table_S2$Spearman_rho < -0.75 & Table_S2$number_data_points > 17,]
Table_S2_plots <- rbind(randSamp_pos, randSamp_neg)
Table_S2_plots <- Table_S2_plots[order(Table_S2_plots$SUBSYSTEM),]

remove(randSamp_pos, randSamp_neg)


Fig_3 <- as.data.frame(table(Table_S2_plots$SUBSYSTEM))
colnames(Fig_3) <- c("subsystem", "n in Table_S2_plots")
Fig_3 <- Fig_3[Fig_3$`n in Table_S2_plots`>2,]

x <- as.data.frame(table(HumanGEM_1_6_0_rxns$SUBSYSTEM))
colnames(x) <- c("subsystem", "n in HumanGEM_1_6_0")

Fig_3 <- merge(Fig_3, x)

Fig_3$percent <- Fig_3$`n in Table_S2_plots` / Fig_3$`n in HumanGEM_1_6_0` * 100

Fig_cat <- read_csv("subsystem categories.csv")

Fig_3 <- merge(Fig_3, Fig_cat)
Fig_3 <- Fig_3[Fig_3$category>2,]
Fig_3 <- Fig_3[order(Fig_3$category, Fig_3$percent),]
x <- sub("Tricarboxylic acid cycle", "TCA cycle", Fig_3$subsystem)
Fig_3$subsystem <- x


#Fig 3: 600x500
par(oma=c(0,17,0,0))
barplot(Fig_3$percent, col=Fig_3$col, horiz = T, xlim=c(0,30), axes=F,
        names.arg = Fig_3$subsystem, las=1, xlab="% rxns correlated with TL")
axis(side=1, at=c(0,10,20,30), labels=c(0,10,20,30))
par(oma=c(0,0,0,0))


#glycolysis
glycolysis <- Table_S2_plots[grep("Glycolysis", Table_S2_plots$SUBSYSTEM),]
order <- c("HMR_4375", "HMR_4391", "HMR_4373", "HMR_4368", "HMR_4365", "HMR_4363")
order <- as.data.frame(order)
colnames(order) <- "ID"
glyc_order <- join(order, glycolysis)

glyc_hm <- glyc_order[,5:37]
glyc_hm <- glyc_hm[,order(CLPtroph)]
glyc_hm[c(1,3,5),] <- -glyc_hm[c(1,3,5),]



# Fig 4 glycolysis: 600x550
par(oma=c(0,0,0,5))
heatmap.2(data.matrix(glyc_hm),
          labRow=c("HMR_4375r", "HMR_4391", "HMR_4373r", "HMR_4368", "HMR_4365r", "HMR_4363"),
          cexCol=1, labCol="",#labCol=sort(CLPtroph), #
          scale="none", dendrogram="none", Colv=F, Rowv=F, 
          trace="none", margins=c(6,6), col=rev(brewer.pal(10,"RdBu")),
          #col=colorRampPalette(c("ghostwhite","ghostwhite","ghostwhite", "salmon", "red4"))(n=10),
          key = F, cexRow = 0.9
          #density.info="none", keysize="2", key.title=NA, lhei=c(1.5,7), key.xlab="flux"
) 
par(oma=c(0,0,0,0))
#plot(1:32, sort(CLPtroph), ylim=c(2,4.5), type="o")


#TCA cycle
TCA_cycle <- Table_S2_plots[grep("Tricarboxylic acid cycle", Table_S2_plots$SUBSYSTEM),]
order <- c("HMR_4145", "HMR_3787", "HMR_4652", "HMR_4410", "HMR_4408")
order <- as.data.frame(order)
colnames(order) <- "ID"
TCA_order <- join(order, TCA_cycle)

TCA_hm <- TCA_order[,5:37]
TCA_hm <- TCA_hm[,order(CLPtroph)]
TCA_hm[c(1,3),] <- -TCA_hm[c(1,3),]

# Fig 4 TCA cycle: 625x560
heatmap.2(data.matrix(TCA_hm),
          #labRow=TCA_cycle$ID, 
          labRow=c("HMR_4145r", "HMR_3787", "HMR_4652r", "HMR_4410", "HMR_4408"),
          cexCol=1, labCol="",#labCol=sort(CLPtroph), #
          scale="none", dendrogram="none", Colv=F, Rowv=F, 
          trace="none", margins=c(6,6), col=(brewer.pal(9,"Reds")),
          #col=colorRampPalette(c("ghostwhite","ghostwhite","ghostwhite", "salmon", "red4"))(n=10),
          #key = F, cexRow = 0.9
          density.info="none", keysize="2", key.title=NA, lhei=c(1.5,7), key.xlab="flux"
) 


#Lys_degradation
Lys_degradation <- Table_S2_plots[grep("Lysine metabolism", Table_S2_plots$SUBSYSTEM),]
Lys_degradation <- rbind(Lys_degradation, Table_S2_plots[grep("Tryptophan metabolism", Table_S2_plots$SUBSYSTEM),])
order <- c("HMR_4288", "HMR_4239", "HMR_4599", "HMR_6415", "HMR_3164", "HMR_3166")
order <- as.data.frame(order)
colnames(order) <- "ID"
Lys_order <- join(order, Lys_degradation)

Lys_hm <- Lys_order[,5:37]
Lys_hm <- Lys_hm[,order(CLPtroph)]

# Fig 5 Lysine: 600x550
par(oma=c(0,0,0,5))
heatmap.2(data.matrix(Lys_hm),
          #labRow=TCA_cycle$ID, 
          labRow=c("HMR_4288", "HMR_4239", "HMR_4599", "HMR_6415", "HMR_3164", "HMR_3166"),
          cexCol=1, labCol="",#labCol=sort(CLPtroph), #
          scale="none", dendrogram="none", Colv=F, Rowv=F, 
          trace="none", margins=c(6,6), col=(brewer.pal(9,"Reds")),
          #col=colorRampPalette(c("ghostwhite","ghostwhite","ghostwhite", "salmon", "red4"))(n=10),
          #key = F, cexRow = 0.9
          density.info="none", keysize="2", key.title=NA, lhei=c(1.5,7), key.xlab="flux"
)
par(oma=c(0,0,0,0))


#his --> IMP
His_IMP <- c("HMR_4437",
             "HMR_4712",
             "HMR_4660",
             "HMR_4658",
             "HMR_3929",
             "HMR_4503" )
His_IMP <- as.data.frame(His_IMP)
colnames(His_IMP) <- "ID"
His_IMP <- join(His_IMP, Table_S2_plots)

#Thr/Asn --> IMP
AA_IMP <- c("HMR_4284", "HMR_3890", "HMR_4172")
AA_IMP <- as.data.frame(AA_IMP)
colnames(AA_IMP) <- "ID"
AA_IMP <- join(AA_IMP, Table_S2_plots)

#PRPP --> IMP
PRPP_IMP <- c("HMR_4406",
              "HMR_4799",
              "HMR_4806",
              "HMR_4808",
              "HMR_4802",
              "HMR_4804",
              "HMR_4810",
              "HMR_4812",
              "HMR_4814",
              "HMR_4038",
              "HMR_4042",
              "HMR_4412" )
PRPP_IMP <- as.data.frame(PRPP_IMP)
colnames(PRPP_IMP) <- "ID"
PRPP_IMP <- join(PRPP_IMP, Table_S2_plots)


#
x <- as.matrix(t(diet_composition[,2:4]))

colnames(x) <- diet_composition$`EltronTriats Diet`

#272x620
par(oma=c(0,3,0,0))
barplot(x[,c(9,8,11,10,1,#7,
             4,3,
             12,2,5)], beside=T, horiz=T,las=2, cex.names=0.9,
        names.arg = c( "Nectar","Fruit","Vegetation", "Seed", "Terrestrial \ninvertebrates", #"Carrion",
                      "Vertebrate \nectoderms", "Vertebrate \nendoderms",
                      "Algae", "Aquatic \ninvertebrates", "Fish"),
        legend=c("proteins", "carbohydrates", "lipids"),
        col=c("papayawhip", "orangered", "darkred"),
        xlim=c(0,25), xlab="nutrient composition (%)")

par(oma=c(0,0,0,0))





