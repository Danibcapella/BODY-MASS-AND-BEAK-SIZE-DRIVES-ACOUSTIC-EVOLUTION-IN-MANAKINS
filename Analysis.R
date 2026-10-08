###### PACOTES E DIRETÓRIO ######

setwd("~/Mestrado/Projeto")

library(phytools)
library(phylosignal)
library(geiger)
library(nlme)
library(vegan)
library(raster)
library(geodata)
library(ggmap)
library(ggplot2)
library(lmodel2)
library(MuMIn)
library(sjPlot)
library(caper)
library(tibble)
library(phylobase)
library(DHARMa)
library(MASS)
library(lme4)

load(".RData")

###### IMPORTAR DADOS ######

##importar arvóre filogenética
Pipridae.Tree <- read.tree("mafft-edge-uce-95-raxml-pfinder.bipartitions.tre")


##corrigir nomes das spp

Pipridae.Tree$tip.label <- gsub('(.*)_\\w+', '\\1', Pipridae.Tree$tip.label)
Pipridae.Tree$tip.label[Pipridae.Tree$tip.label == "Dixiphia_pipra"] <- "Pseudopipra_pipra"  
Pipridae.Tree$tip.label[Pipridae.Tree$tip.label == "Illicura_militaris"] <- "Ilicura_militaris" 
Pipridae.Tree$tip.label[Pipridae.Tree$tip.label == "Xenopipo_atronitrens"] <- "Xenopipo_atronitens" 

###para objetos multi.phylo
for(i in 1:length(Pipridae.Tree)) {
  Pipridae.Tree[[i]]$tip.label <- gsub('(.*)_\\w+', '\\1', Pipridae.Tree[[i]]$tip.label)
  Pipridae.Tree[[i]]$tip.label[Pipridae.Tree[[i]]$tip.label == "Dixiphia_pipra"] <- "Pseudopipra_pipra"   
  
  
}
consensus.pipridae <- consensus(Pipridae.Tree)


## enraizar árvore com grupo externo

Pipridae.Tree <- root(Pipridae.Tree, 1:3, resolve.root = T)
plot.phylo(Pipridae.Tree)

##/////IMPORTAR DATAFRAME\\\\\##
Pipridae.Traits <- read.csv2("Pipridae.csv", stringsAsFactors = T)

#cortar arvore
pruned.tree<-drop.tip(Pipridae.Tree, setdiff(Pipridae.Tree$tip.label, Pipridae.Traits$especie))
pruned.tree$node.label <- NULL
plot.phylo(pruned.tree)
is.rooted(pruned.tree)

##/////IMPORTAR DATA FRAME COM MÉDIAS/////##
Traits.mean <- read.csv2("Pipridae_mean.csv", stringsAsFactors = T, row.names = 1)

###### OBTENÇÃO DE DADOS DO WORLDCLIM E INCORPORAÇÃO DESTES NO DATAFRAME######
bioclim <- worldclim_global(var = "bio", res = 10, path = tempdir())

coords <- data.frame(x=Pipridae.Traits$longitude,y=Pipridae.Traits$latitude)

points <- SpatialPoints(coords)

values <- extract(bioclim,coords)

df <- cbind.data.frame(coordinates(points),values)
Pipridae.Traits[,44:62] <- df[,4:22]


## matriz de correlação das variaveis bioclimaticas de temperatura
#cor(Pipridae.Traits[,44:54], method = "pearson", use = "complete.obs")

## matriz de correlação das variaveis bioclimaticas de precipitação
#cor(Pipridae.Traits[,55:62], method = "pearson", use = "complete.obs")

## remover variaveis bioclimaticas correlacionadas (de temperatura com bio01 e de precipitação com bio12)
# <- Pipridae.Traits[,-c(48, 49, 51, 52, 53, 54)]


##fazer o mapa
#mapa <- borders("world", regions = c("Brazil", "Uruguay", "Argentina", "French Guiana", "Suriname", "Colombia", "Venezuela",
                                     #"Bolivia", "Ecuador", "Chile", "Paraguay", "Peru", "Guyana", "Panama", "Costa Rica", 
                                     #"Nicaragua", "Honduras", "El Salvador", "Belize", "Guatemala", "Mexico", "Trinidad and Tobago",
                                     #"Caribe", "Puerto Rico", "Dominican Republic", "Haiti", "Jamaica", "Cuba", "Bahamas", "Antiles",
                                     #"Dominica", "Saba"), 
               # fill = "grey70", colour = "black")

#mapa.plot <- ggplot() + mapa + theme_bw() + xlab("") + ylab("") + 
  #theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "grey80"), panel.grid.minor = element_blank())
  
#mapa.plot + 
  #geom_point(data = Pipridae.Traits, aes(x = x, y = y),col="red", size=1)+
  #theme_bw()
 

##médias das variáveis bioclimáticas com detour pra remover NAs
bioclim.mean <- Pipridae.Traits[-214, c(2, 44:62)]
bioclim.mean <- aggregate(bioclim.mean[,2:20], list(bioclim.mean$especie), FUN=mean)

cor(bioclim.mean[,2:20], method = "pearson", use = "complete.obs")

bioclim.cor <- as.data.frame(c(bioclim.mean[,1:5],bioclim.mean[,13:20]))
bioclim.cor$bio7 <- bioclim.mean$bio7

##Padronização dos dados bioclimaticos
bioclim.std <- decostand(bioclim.cor[2:14], method="standardize")

##PCA das variáveis bioclimaticas
bioclim.pca <- rda(bioclim.std, scale=T)
#summary(bioclim.pca)

##Percentual de explicação dos eixos da PPCA
explic.bioclim <- (bioclim.pca$CA$eig)/(sum(bioclim.pca$CA$eig))*100 
explic.bioclim

##scores da pca bioclimátca
bioclim.score <- scores(bioclim.pca, display="wa", choices=1:2)
bioclim.score <- as.data.frame(bioclim.score)
#bioclim.score
scores(bioclim.pca, display="sp", choices=1:2)




###### ADICIONAR COLUNAS COM VALORES EM LOG######
Pipridae.stats$logFmin <- log(Pipridae.stats$FMIN)
Pipridae.stats$logFdom <- log(Pipridae.stats$FDOM)
Pipridae.stats$logMass <- log(Pipridae.stats$Mass)
Pipridae.stats$logDur <- log(Pipridae.stats$DUR)
Pipridae.stats$logBan <- log(Pipridae.stats$BAN)
Pipridae.stats$logFmax <- log(Pipridae.stats$FMAX)
Pipridae.stats$logPace <- log(Pipridae.stats$Pace)


## matriz de correlação das variaveis de tamanho do bico
cor(Traits.mean[,8:11], method = "pearson", use = "complete.obs")

##Separar atributos de tamanho do bico em um dataframe
Beak.mean <- Traits.mean[,8:11]

##Padronização dos dados do bico
Beak.std <- decostand(Beak.mean, method="standardize")

##PPCA das variáveis de bico utlizando modelo de evolução Browniana
Beak.ppca <- phyl.pca(pruned.tree, Beak.std, method = "BM", mode = "cov")
#Beak.ppca

##Percentual de explicação dos eixos da PPCA
#explic.beak <- (Beak.ppca$Eval)/(sum(Beak.ppca$Eval))*100 
#explic.beak

##Scores da PPCA
Beak.scores.ppca <- Beak.ppca$S
#Beak.scores.ppca
Beak.scores.ppca <- as.data.frame(Beak.scores.ppca)


#name.check(Pipridae.Tree, Traits.mean)
#name.check(scores.ppca, Traits.mean)

##Juntar escores da ppca no dataframe para analise
Pipridae.stats$BeakPC1 <- Beak.scores.ppca$PC1

Pipridae.stats$BeakPC1.cor <- Pipridae.stats$BeakPC1 * -1



## Juntar escores da PCA bioclimática
Pipridae.stats$BioclimPC1 <- bioclim.score$PC1

Pipridae.stats$BioclimPC2 <- bioclim.score$PC2


###### CONSTRUÇÃO DOS MODELOS ######

## juntar filogenia e data frame para analises de pgls
Pipridae.stats <- rownames_to_column(Pipridae.stats, "species")

Pipridae.stats[,31:49] <- bioclim.mean[,2:20]

Pipridae.caper <- comparative.data(phy = pruned.tree, data = Pipridae.stats, names.col = species,  vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)




## calcular todos AICCs e juntar em um dataframe
AICC_brown <- numeric()
AICC_lambda <- numeric()
AICC_delta <- numeric()
vec_resp <- factor(c("NNOT","NNOT","NNOT","NNOT", "NNOT", "Fmin","Fmin","Fmin","Fmin", "Fmin",
                     "Fdom","Fdom","Fdom","Fdom", "Fdom","Dur", "Dur",
                     "Dur","Dur","Dur", "Ban","Ban","Ban","Ban", "Ban", "Fmax",
                     "Fmax","Fmax","Fmax", "Fmax", "Pace", "Pace", "Pace", "Pace", "Pace"))
vec_exp <- factor(rep(c("dicromatismo", "display", "bioclim", "mass", "beak"), 7))
for (i in c(6, 23, 24, 26:29)) {
  for (j in c(19:21, 25, 30)) {
    x <- pgls(Pipridae.caper[[2]][,i] ~ Pipridae.caper[[2]][,j], data = Pipridae.caper)
    AICC_brown <- c(AICC_brown, x$aicc)
    y <- pgls(Pipridae.caper[[2]][,i] ~ Pipridae.caper[[2]][,j], data = Pipridae.caper, lambda = "ML")
    AICC_lambda <- c(AICC_lambda, y$aicc)
    z <- pgls(Pipridae.caper[[2]][,i] ~ Pipridae.caper[[2]][,j], data = Pipridae.caper, delta = "ML")
    AICC_delta <- c(AICC_delta, z$aicc)
  }
}
AICC_df <- data.frame(vec_resp, vec_exp, AICC_brown, AICC_lambda, AICC_delta)

write.csv2(AICC_df, "tabela_AICC.csv", row.names = F)


##pgls  - FDOM x CLIMA

# modelo maximum likelihood delta
pglsModel_FDom_bioclim <- pgls(logFdom ~ BioclimPC1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_FDom_bioclim)
par(mfrow = c(2,2))
#plot.pgls(pglsModel_FDom_bioclim)


##pgls FMIN x CLIMA

# delta
pglsModel_FMin_bioclim <- pgls(logFmin ~ BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_FMin_bioclim)
#plot.pgls(pglsModel_FMin_bioclim)


## pgls FMAX x CLIMA

# delta
pglsModel_FMax_bioclim <- pgls(logFmax ~ BioclimPC1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_FMax_bioclim)
plot.pgls(pglsModel_FMax_bioclim)


## pgls BAN x CLIMA

# delta
pglsModel_Ban_bioclim <- pgls(logBan ~ BioclimPC1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Ban_bioclim)
#plot.pgls(pglsModel_Ban_bioclim)


##pgls  - DUR x CLIMA

# delta
pglsModel_Dur_bioclim <- pgls(logDur ~ BioclimPC1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Dur_bioclim)
#plot.pgls(pglsModel_Dur_bioclim)


## pgls NNOT x CLIMA

# lambda
pglsModel_Not_bioclim <- pgls(NNOT ~ BioclimPC1, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Not_bioclim)
#plot.pgls(pglsModel_Not_bioclim)


## pgls PACE x CLIMA

# lambda
pglsModel_Pace_bioclim<- pgls(logPace~ BioclimPC1, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Pace_bioclim)
#plot.pgls(pglsModel_Pace_bioclim)


## pgls FMIN x DICROMATISMO

# delta
pglsModel_Fmin_dicrom <- pgls(logFmin ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fmin_dicrom)
#plot.pgls(pglsModel_Fmin_dicrom)



## pgls FDOM x DICROMATISMO

# delta
pglsModel_Fdom_dicrom <- pgls(logFdom ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fdom_dicrom)
#plot.pgls(pglsModel_Fdom_dicrom)



## pgls FMAX x DICROMATISMO

# delta
pglsModel_Fmax_dicrom <- pgls(logFmax ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fmax_dicrom)
#plot.pgls(pglsModel_Fmax_dicrom)



##pgls DUR x DICROMATISMO

# delta
pglsModel_Dur_dicrom <- pgls(logDur ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Dur_dicrom)
#plot.pgls(pglsModel_Dur_dicrom)



## pgls NNOT x DICROMATISMO

# delta
pglsModel_nnot_dicrom <- pgls(NNOT ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_nnot_dicrom)
#plot.pgls(pglsModel_nnot_dicrom)



## pgls BAN x DICROMATISMO

# delta
pglsModel_Ban_dicrom <- pgls(logBan ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Ban_dicrom)
#plot.pgls(pglsModel_Ban_dicrom)


## pgls PACE x DICROMATISMO

# delta
pglsModel_Pace_dicrom <- pgls(logPace ~ Dicromatismo, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Pace_dicrom)
#plot.pgls(pglsModel_Pace_dicrom)


##pgls FMIN X DISPLAY

# delta
pglsModel_Fmin_display <- pgls(logFmin ~ Display.Complex, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fmin_display)
#plot.pgls(pglsModel_Fmin_display)


##pgls FDOM X DISPLAY

# delta
pglsModel_Fdom_display <- pgls(logFdom ~ Display.Complex, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fdom_display)
#plot.pgls(pglsModel_Fdom_display)


##pgls FMAX X DISPLAY

# lambda
pglsModel_Fmax_display <- pgls(logFmax ~ Display.Complex, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Fmax_display)
#plot.pgls(pglsModel_Fmax_display)


##pgls DUR X DISPLAY

# lambda
pglsModel_Dur_display <- pgls(logDur ~ Display.Complex, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Dur_display)
#plot.pgls(pglsModel_Dur_display)


## pgls BAN x DISPLAY

# lambda
pglsModel_Ban_display <- pgls(logBan ~ Display.Complex, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Ban_display)
#plot.pgls(pglsModel_Ban_display)


## pgls PACE x DISPLAY

# lambda
pglsModel_Pace_display_lambda <- pgls(logPace ~ Display.Complex, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Pace_display_lambda)
#plot.pgls(pglsModel_Pace_display_lambda)


## pgls NNOT x DISPLAY

# lambda
pglsModel_Nnot_display <- pgls(NNOT ~ Display.Complex, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Nnot_display)
#plot.pgls(pglsModel_Nnot_display)


## pgls NNOT x MASSA

# lambda
pglsModel_Nnot_mass <- pgls(NNOT ~ logMass, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Nnot_mass)
#plot.pgls(pglsModel_Nnot_mass)


## pgls FMIN x MASSA

# lambda
pglsModel_Fmin_mass <- pgls(logFmin ~ logMass, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Fmin_mass)
#plot.pgls(pglsModel_Fmin_mass)


## pgls FDOM x MASSA

# delta
pglsModel_Fdom_mass <- pgls(logFdom ~ logMass, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fdom_mass)
#plot.pgls(pglsModel_Fdom_mass)


## pgls DUR x MASSA

# delta
pglsModel_Dur_mass <- pgls(logDur ~ logMass, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Dur_mass)
#plot.pgls(pglsModel_Dur_mass)


## pgls BAN x MASSA

# lambda
pglsModel_Ban_mass <- pgls(logBan ~ logMass, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Ban_mass)
#plot.pgls(pglsModel_Ban_mass)


## pgls FMAX x MASSA

# delta
pglsModel_Fmax_mass <- pgls(logFmax ~ logMass, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fmax_mass)
#plot.pgls(pglsModel_Fmax_mass)


## pgls PACE x MASSA

# lambda
pglsModel_Pace_mass <- pgls(Pace ~ logMass, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Pace_mass)
#plot.pgls(pglsModel_Pace_mass)



## pgls NNOT x BICO

# lambda
pglsModel_Nnot_beak <- pgls(NNOT ~ BeakPC1, data = Pipridae.caper, lambda = "ML")
summary(pglsModel_Nnot_beak)
#plot.pgls(pglsModel_Nnot_beak)


## pgls FMIN x BICO

# lambda
pglsModel_Fmin_beak <- pgls(logFmin ~ BeakPC1, data = Pipridae.caper, lambda = "ML")
summary(pglsModel_Fmin_beak)
#plot.pgls(pglsModel_Fmin_beak)


## pgls FDOM x BICO

# delta
pglsModel_Fdom_beak <- pgls(logFdom ~ BeakPC1, data = Pipridae.caper, delta = "ML")
summary(pglsModel_Fdom_beak)
#plot.pgls(pglsModel_Fdom_beak)


## pgls DUR x BICO

# delta
pglsModel_Dur_beak <- pgls(logDur ~ BeakPC1, data = Pipridae.caper, delta = "ML")
summary(pglsModel_Dur_beak)
#plot.pgls(pglsModel_Dur_beak)


## pgls BAN x BICO

# lambda
pglsModel_Ban_beak <- pgls(logBan ~ BeakPC1, data = Pipridae.caper, lambda = "ML")
summary(pglsModel_Ban_beak)
#plot.pgls(pglsModel_Ban_beak)


## pgls FMAX x BICO

# delta
pglsModel_Fmax_beak <- pgls(logFmax ~ BeakPC1, data = Pipridae.caper, delta = "ML")
summary(pglsModel_Fmax_beak)
#plot.pgls(pglsModel_Fmax_beak)
par(mfrow=c(2,2))
plot.pgls(pglsModel_Pace_null)

## pgls PACE x BICO

# lambda
pglsModel_Pace_beak <- pgls(Pace ~ BeakPC1, data = Pipridae.caper, lambda = "ML")
summary(pglsModel_Pace_beak)
#plot.pgls(pglsModel_Pace_beak)

pglsModel_Pace_null <- pgls(Pace ~ 1, data = Pipridae.caper)
summary(pglsModel_Pace_null)

http://127.0.0.1:29123/graphics/plot_zoom_png?width=1765&height=541
###### MODELOS ADITIVOS ######

## calcular todos AICCs e juntar em um dataframe
AICC_brown <- numeric()
AICC_lambda <- numeric()
AICC_delta <- numeric()
vec_resp <- factor(c("NNOT","NNOT","NNOT","NNOT", "Fmin","Fmin","Fmin","Fmin",
                     "Fdom","Fdom","Fdom","Fdom","Dur", "Dur",
                     "Dur","Dur", "Ban","Ban","Ban","Ban", "Fmax",
                     "Fmax","Fmax","Fmax", "Pace", "Pace", "Pace", "Pace"))
vec_exp <- factor(rep(c("dicromatismo", "display", "mass", "beak"), 7))
for (i in c(6, 23, 24, 26:29)) {
  for (j in c(19:20, 25, 30)) {
    x <- pgls(Pipridae.caper[[2]][,i] ~ Pipridae.caper[[2]][,j] + Pipridae.caper[[2]][,21], data = Pipridae.caper)
    AICC_brown <- c(AICC_brown, x$aicc)
    y <- pgls(Pipridae.caper[[2]][,i] ~ Pipridae.caper[[2]][,j] + Pipridae.caper[[2]][,21], data = Pipridae.caper, lambda = "ML")
    AICC_lambda <- c(AICC_lambda, y$aicc)
    z <- pgls(Pipridae.caper[[2]][,i] ~ Pipridae.caper[[2]][,j] + Pipridae.caper[[2]][,21], data = Pipridae.caper, delta = "ML")
    AICC_delta <- c(AICC_delta, z$aicc)
  }
}
AICC_df_add <- data.frame(vec_resp, vec_exp, AICC_brown, AICC_lambda, AICC_delta)

write.csv2(AICC_df_add, "tabela_AICC_add.csv", row.names = F)



## matriz de correlação das variaveis explicativas
cor(Pipridae.stats[,c(20:22, 26, 31)], method = "pearson", use = "complete.obs")


##pgls NNOT x DICROMATISMO + CLIMA

# delta
pglsModel_Nnot_dicrom_bioclim <- pgls(NNOT ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Nnot_dicrom_bioclim)
#plot.pgls(pglsModel_Nnot_dicrom_bioclim)

##pgls NNOT x DISPLAY + CLIMA

# lambda
pglsModel_Nnot_display_bioclim <- pgls(NNOT ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Nnot_display_bioclim)
#plot.pgls(pglsModel_Nnot_display_bioclim)

##pgls NNOT x MASSA + CLIMA

# lambda
pglsModel_Nnot_mass_bioclim <- pgls(NNOT ~ logMass + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Nnot_mass_bioclim)
#plot.pgls(pglsModel_Nnot_mass_bioclim)

##pgls NNOT x BICO + CLIMA

# lambda
pglsModel_Nnot_beak_bioclim <- pgls(NNOT ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Nnot_beak_bioclim)
#plot.pgls(pglsModel_Nnot_beak_bioclim)


##pgls FMIN x DICROMATISMO + CLIMA

# delta
pglsModel_Fmin_dicrom_bioclim <- pgls(logFmin ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fmin_dicrom_bioclim)
#plot.pgls(pglsModel_Fmin_dicrom_bioclim)

##pgls FMIN x DISPLAY + CLIMA

# lambda
pglsModel_Fmin_display_bioclim <- pgls(logFmin ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Fmin_display_bioclim)
#plot.pgls(pglsModel_Fmin_display_bioclim)

##pgls FMIN x MASSA + CLIMA

# delta
pglsModel_Fmin_mass_bioclim <- pgls(logFmin ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fmin_mass_bioclim)
#plot.pgls(pglsModel_Fmin_mass_bioclim)

##pgls FMIN x BICO + CLIMA

# delta
pglsModel_Fmin_beak_bioclim <- pgls(logFmin ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fmin_beak_bioclim)
#plot.pgls(pglsModel_Fmin_beak_bioclim)


##pgls FDOM x DICROMATISMO + CLIMA

# delta
pglsModel_Fdom_dicrom_bioclim <- pgls(logFdom ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fdom_dicrom_bioclim)
#plot.pgls(pglsModel_Fdom_dicrom_bioclim)

##pgls FDOM x DISPLAY + CLIMA

# lambda
pglsModel_Fdom_display_bioclim <- pgls(logFdom ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Fdom_display_bioclim)
#plot.pgls(pglsModel_Fdom_display_bioclim)

##pgls FDOM x MASSA + CLIMA

# delta
pglsModel_Fdom_mass_bioclim <- pgls(logFdom ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fdom_mass_bioclim)
#plot.pgls(pglsModel_Fdom_mass_bioclim)

##pgls FDOM x BICO + CLIMA

# delta
pglsModel_Fdom_beak_bioclim <- pgls(logFdom ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fdom_beak_bioclim)
#plot.pgls(pglsModel_Fdom_beak_bioclim)


##pgls DUR x DICROMATISMO + CLIMA

# delta
pglsModel_Dur_dicrom_bioclim <- pgls(logDur ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Dur_dicrom_bioclim)
#plot.pgls(pglsModel_Dur_dicrom_bioclim)

##pgls DUR x DISPLAY + CLIMA

# lambda
pglsModel_Dur_display_bioclim <- pgls(logDur ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Dur_display_bioclim)
#plot.pgls(pglsModel_Dur_display_bioclim)

##pgls DUR x MASSA + CLIMA

# delta
pglsModel_Dur_mass_bioclim <- pgls(logDur ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Dur_mass_bioclim)
#plot.pgls(pglsModel_Dur_mass_bioclim)

##pgls DUR x BICO + CLIMA

# delta
pglsModel_Dur_beak_bioclim <- pgls(logDur ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Dur_beak_bioclim)
#plot.pgls(pglsModel_Dur_beak_bioclim)


##pgls BAN x DICROMATISMO + CLIMA

# delta
pglsModel_Ban_dicrom_bioclim <- pgls(logBan ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Ban_dicrom_bioclim)
#plot.pgls(pglsModel_Ban_dicrom_bioclim)

##pgls BAN x DISPLAY + CLIMA

# lambda
pglsModel_Ban_display_bioclim <- pgls(logBan ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Ban_display_bioclim)
#plot.pgls(pglsModel_Ban_display_bioclim)

##pgls BAN x MASSA + CLIMA

# lambda
pglsModel_Ban_mass_bioclim <- pgls(logBan ~ logMass + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Ban_mass_bioclim)
#plot.pgls(pglsModel_Ban_mass_bioclim)

##pgls BAN x BICO + CLIMA

# lambda
pglsModel_Ban_beak_bioclim <- pgls(logBan ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Ban_beak_bioclim)
#plot.pgls(pglsModel_Ban_beak_bioclim)


##pgls FMAX x DICROMATISMO + CLIMA

# delta
pglsModel_Fmax_dicrom_bioclim <- pgls(logFmax ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fmax_dicrom_bioclim)
#plot.pgls(pglsModel_Fmax_dicrom_bioclim)

##pgls FMAX x DISPLAY + CLIMA

# lambda
pglsModel_Fmax_display_bioclim <- pgls(logFmax ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Fmax_display_bioclim)
#plot.pgls(pglsModel_Fmax_display_bioclim)

##pgls FMAX x MASSA + CLIMA

# delta
pglsModel_Fmax_mass_bioclim <- pgls(logFmax ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fmax_mass_bioclim)
#plot.pgls(pglsModel_Fmax_mass_bioclim)

##pgls FMAX x BICO + CLIMA

# delta
pglsModel_Fmax_beak_bioclim <- pgls(logFmax ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Fmax_beak_bioclim)
#plot.pgls(pglsModel_Fmax_beak_bioclim)


##pgls PACE x DICROMATISMO + CLIMA

# delta
pglsModel_Pace_dicrom_bioclim <- pgls(logPace ~ Dicromatismo + BioclimPC1, data = Pipridae.caper, delta= "ML")
#summary(pglsModel_Pace_dicrom_bioclim)
#plot.pgls(pglsModel_Pace_dicrom_bioclim)

##pgls PACE x DISPLAY + CLIMA

# lambda
pglsModel_Pace_display_bioclim <- pgls(logPace ~ Display.Complex + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Pace_display_bioclim)
#plot.pgls(pglsModel_Pace_display_bioclim

##pgls PACE x MASSA + CLIMA

# lambda
pglsModel_Pace_mass_bioclim <- pgls(logPace ~ logMass + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Pace_mass_bioclim)
#plot.pgls(pglsModel_Pace_mass_bioclim)

##pgls PACE x BICO + CLIMA

# lambda
pglsModel_Pace_beak_bioclim <- pgls(logPace ~ BeakPC1 + BioclimPC1, data = Pipridae.caper, lambda= "ML")
#summary(pglsModel_Pace_beak_bioclim)
#plot.pgls(pglsModel_Pace_beak_bioclim)

# lambda
pglsModel_Pace_beak_bioclim_corrected <- pgls(logPace ~ BeakPC1.cor + BioclimPC1, data = Pipridae.caper, lambda= "ML")


###### MODELOS NULOS ######

## calcular todos AICCs e juntar em um dataframe
AICC_brown <- numeric()
AICC_lambda <- numeric()
AICC_delta <- numeric()
vec_resp <- factor(c("NNOT","Fmin","Fdom","Dur","Ban","Fmax","Pace"))
for (i in c(6, 23, 24, 26:29)) {
    x <- pgls(Pipridae.caper[[2]][,i] ~ 1, data = Pipridae.caper)
    AICC_brown <- c(AICC_brown, x$aicc)
    y <- pgls(Pipridae.caper[[2]][,i] ~ 1, data = Pipridae.caper, lambda = "ML")
    AICC_lambda <- c(AICC_lambda, y$aicc)
    z <- pgls(Pipridae.caper[[2]][,i] ~ 1, data = Pipridae.caper, delta = "ML")
    AICC_delta <- c(AICC_delta, z$aicc)
  }
AICC_df_null <- data.frame(vec_resp, AICC_brown, AICC_lambda, AICC_delta)

write.csv2(AICC_df_null, "tabela_AICC_null.csv", row.names = F)


## pgls NNOT NULO

# lambda
pglsModel_NNot_null <- pgls(NNOT ~ 1, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_NNot_null)
#plot.pgls(pglsModel_NNot_null)


## pgls FMIN NULO

# lambda
pglsModel_Fmin_null <- pgls(logFmin ~ 1, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Fmin_null)
#plot.pgls(pglsModel_Fmin_null)


## pgls FDOM NULO

# delta
pglsModel_Fdom_null <- pgls(logFdom ~ 1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fdom_null)
#plot.pgls(pglsModel_Fdom_null)


## pgls DUR NULO

# delta
pglsModel_Dur_null <- pgls(logDur ~ 1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Dur_null)
#plot.pgls(pglsModel_Dur_null)


## pgls BAN NULO

# lambda
pglsModel_Ban_null <- pgls(logBan ~ 1, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Ban_null)
#plot.pgls(pglsModel_Ban_null)


## pgls FMAX NULO

# delta
pglsModel_Fmax_null <- pgls(logFmax ~ 1, data = Pipridae.caper, delta = "ML")
#summary(pglsModel_Fmax_null)
#plot.pgls(pglsModel_Fmax_null)


## pgls PACE NULO

# lambda
pglsModel_Pace_null <- pgls(logPace ~ 1, data = Pipridae.caper, lambda = "ML")
#summary(pglsModel_Pace_null)
#plot.pgls(pglsModel_Pace_null)



###### RMA - Constantes alométricas ######

windows()
par(mfrow = c(1,3))

#RMA_fmin_mass <- lmodel2(logFmin ~ logMass, data = Pipridae.stats, "relative", "relative", 100)
#RMA_fmin_mass
plot(RMA_fmin_mass, "RMA", pch = 19, cex = 1.2, ylab = "log Fmin", xlab = "log Mass")

#RMA_fdom_mass <- lmodel2(logFdom ~ logMass, data = Pipridae.stats, "relative", "relative", 100)
#RMA_fdom_mass
plot(RMA_fdom_mass, "RMA", pch = 19, cex = 1.2, ylab = "log Fdom", xlab = "log Mass")

#RMA_fmax_mass <- lmodel2(logFmax ~ logMass, data= Pipridae.stats, "relative", "relative", 100)
#RMA_fmax_mass
plot(RMA_fmax_mass, "RMA", pch = 19, cex = 1.2, ylab = "log Fmax", xlab = "log Mass")


###### PLOTAR GRÁFICOS ######

windows()
par(mfrow = c(2,2))

## Fmin ~ display + bioclim (r² = 0,874)

plot(logFmin ~ Display.Complex, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Fmin_display_bioclim, col = "gray70", lwd = 3)
#summary(pglsModel_Fmin_display_bioclim)

## Fdom ~ display + bioclim (r² = 0,818)

plot(logFdom ~ Display.Complex, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Fdom_display_bioclim, col = "gray70", lwd = 3)
#summary(pglsModel_Fdom_display_bioclim)

## Fmax ~ display + bioclim (r² = 0,818)

plot(logFmax ~ Display.Complex, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Fmax_display_bioclim, col = "gray70", lwd = 3)
#summary(pglsModel_Fmax_display_bioclim)

## Fmax ~ display (r² = 0,727)

plot(logFmax ~ Display.Complex, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Fmax_display, col = "gray70", lwd = 3)
#summary(pglsModel_Fmax_display)

windows()
par(mfrow = c(2,2))

## Pace ~ beak + bioclim (r² = 0,157)

plot(logPace ~ BeakPC1.cor, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Pace_beak_bioclim_corrected, col = "gray70", lwd = 3)
#summary(pglsModel_Pace_beak_bioclim)

## Ban ~ mass + bioclim (r² = 0,156)

plot(logBan ~ logMass, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Ban_mass_bioclim, col = "gray70", lwd = 3)
#summary(pglsModel_Ban_mass_bioclim)

## Dur ~ massa (r² = 0,104)

plot(logDur ~ logMass, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_Dur_mass, col = "gray70", lwd = 3)
#summary(pglsModel_Dur_mass)

## FMax ~ Bioclim (r² = 0,098)

plot(logFmax ~ BioclimPC1, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pglsModel_FMax_bioclim, col = "gray70", lwd = 3)
#summary(pglsModel_FMax_bioclim)

###### PHYLO4D E ANALISES RELACIONADAS ######
##/////MONTAR PHYLO4D\\\\\##
Pipridae.p4d <- Pipridae.stats[-33,]

Pipridae.p4d <- phylo4d(pruned.tree, Pipridae.p4d, match.data=F)

##gráficos mostrando valores dos atributos em meio a filogenia

barplot.phylo4d(Pipridae.p4d, trait = "logFdom")


fdom.cg <- phyloCorrelogram(Pipridae.p4d, trait = "logFdom")
plot(fdom.cg)

##cálculo de sinal filogenético dos atributos
signal <- phyloSignal(Pipridae.p4d[,-1], reps = 1000)
signal

library(openxlsx)

df.signal <- as.data.frame(signal$pvalue)
df.signal$parameter <- rownames(df.signal)
write.xlsx(df.signal, "table 5.xlsx")

##correlograma de sinal filogenético

par(mfrow = c(2,2))

fdom.cg <- phyloCorrelogram(Pipridae.p4d, trait = "logFdom")
plot(fdom.cg, main = "Peak Frequency")

ban.cg <- phyloCorrelogram(Pipridae.p4d, trait = "logBan")
plot(ban.cg, main = "Frequency bandwidth")

dur.cg <- phyloCorrelogram(Pipridae.p4d, trait = "logDur")
plot(dur.cg, main = "Call duration")

fmax.cg <- phyloCorrelogram(Pipridae.p4d, trait = "logFmax")
plot(fmax.cg, main = "Maximum frequency")

##encontrar hotspots de autocorrelação filogenética

fdom.lipa <- lipaMoran(Pipridae.p4d, trait = "logFdom", prox.phylo = "nNodes", as.p4d = TRUE)

points.fdom <- lipaMoran(Pipridae.p4d, trait = "logFdom", prox.phylo = "nNodes")$p.value

points.fdom <- ifelse(points.fdom < 0.05, "red", "black")
dotplot.phylo4d(fdom.lipa, dot.col = points.fdom)



ban.lipa <- lipaMoran(Pipridae.p4d, trait = "logBan", prox.phylo = "nNodes", as.p4d = TRUE)

points.ban <- lipaMoran(Pipridae.p4d, trait = "logBan", prox.phylo = "nNodes")$p.value

points.ban <- ifelse(points.ban < 0.05, "red", "black")
dotplot.phylo4d(ban.lipa, dot.col = points.ban)



dur.lipa <- lipaMoran(Pipridae.p4d, trait = "logDur", prox.phylo = "nNodes", as.p4d = TRUE)

points.dur <- lipaMoran(Pipridae.p4d, trait = "logDur", prox.phylo = "nNodes")$p.value

points.dur <- ifelse(points.dur < 0.05, "red", "black")
dotplot.phylo4d(dur.lipa, dot.col = points.dur)



fmax.lipa <- lipaMoran(Pipridae.p4d, trait = "logFmax", prox.phylo = "nNodes", as.p4d = TRUE)

points.fmax <- lipaMoran(Pipridae.p4d, trait = "logFmax", prox.phylo = "nNodes")$p.value

points.fmax <- ifelse(points.fmax < 0.05, "red", "black")
dotplot.phylo4d(fmax.lipa, dot.col = points.fmax)


###### NOVAS ANALISES ######

##/////IMPORTAR DATAFRAME\\\\\##
Pipridae.New <- read.csv2("Pipridae.csv", stringsAsFactors = T)
Pipridae.New <- cbind.data.frame(Pipridae.New, df)

str(Pipridae.New)

Pipridae.New$lat <- abs(Pipridae.New$y)

Pipridae.New <- subset(Pipridae.New, select = -c(amostra, localidade, estado, pais, longitude, latitude, altitude,
                                                 ano, mes, dia, ocorrencia, registro_de_analise, Inference, Traits.inferred,
                                                 Reference.species, Habitat, HabitaDensity, Migration, Trophic.Level, Trophic.Niche,
                                                 Primary.Lifestyle, SNR, x, y))
library(dplyr)
Pipridae.New <- Pipridae.New %>% mutate_if(is.integer, as.numeric)

boxplot(Pipridae.New)


Pipridae.Std <- decostand(x = Pipridae.New[,-1], method = "standardize")

Pipridae.Std$especie <- Pipridae.New$especie

boxplot(Pipridae.Std)


###????
glm_fdom1 <- glmer(FDOM ~ bio5 + bio6 + lat + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom1)

glm_fdom2 <- glmer(FDOM ~ bio5 + bio6 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom2)

glm_fdom3 <- glmer(FDOM ~ bio5 + lat + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom3)

glm_fdom4 <- glmer(FDOM ~ bio6 + lat + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom4)

glm_fdom5 <- glmer(FDOM ~ lat + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom5)

glm_fdom6 <- glmer(FDOM ~ bio6 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom6)

glm_fdom7 <- glmer(FDOM ~ bio5 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fdom7)

model.sel(glm_fdom1, glm_fdom2, glm_fdom3, glm_fdom4, glm_fdom5, glm_fdom6, glm_fdom7)


plot_model(glm_fdom4, type = "pred")



glm_fmin1 <- glmer(FMIN ~ bio5 + bio6 + lat + (1|especie),  family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmin1)

glm_fmin2 <- glmer(FMIN~ bio5 + bio6 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmin2)

glm_fmin3 <- glmer(FMIN ~ bio5 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmin3)

glm_fmin4 <- glmer(FMIN ~ bio6 + (1|especie), family =  family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmin4)

model.sel(glm_fmin1, glm_fmin2, glm_fmin3)

plot_model(glm_fmin1, type = "slope")



glm_fmax1 <- glmer(FMAX ~ bio5 + bio6 + lat + (1|especie),  family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmax1)

glm_fmax2 <- glmer(FMAX ~ bio5 + bio6 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmax2)

glm_fmax3 <- glmer(FMAX ~ bio5 + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmax3)

glm_fmax4 <- glmer(FMAX ~ bio6 + (1|especie),  family = gaussian(link = "inverse"), data = Pipridae.Std, 
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_fmax4)

model.sel(glm_fmax1, glm_fmax2, glm_fmax3, glm_fmax4)



### Duração com lei de Bergman

glm_dur1 <- glmer(DUR ~ bio5 + bio6 + lat + (1|especie), family = gaussian(link = "inverse"), data = Pipridae.Std, 
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm_dur1)



# Calcular os residuos
simu <- simulateResiduals(fittedModel = glm_fdom5, plot = T)
testOutliers(simu, type = "bootstrap")
#### Residuos muito dispersos, preciso encontrar um modelo melhor

plot_model(glm_dur1, type = "eff")

plot(FDOM ~ lat, data = Pipridae.New)

###### MODELOS NOVOS ADITIVOS ######

cor(Traits.mean[,8:17], method = "pearson", use = "complete.obs")

# Modelo frequência dominante #

pgls_fdom1 <- pgls(logFdom ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom1)
pgls_fdom1$aic

pgls_fdom2 <- pgls(logFdom ~ logMass + bio5 + bio6 + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom2)
pgls_fdom2$aic

pgls_fdom3 <- pgls(logFdom ~ logMass + bio6 + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom3)
pgls_fdom3$aic

pgls_fdom4 <- pgls(logFdom ~ logMass + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom4)
pgls_fdom4$aic

pgls_fdom5 <- pgls(logFdom ~ Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom5)
pgls_fdom5$aic

##*
pgls_fdom6 <- pgls(logFdom ~ Beak.Length_Nares + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom6)
pgls_fdom6$aic

pgls_fdom7 <- pgls(logFdom ~ Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom7)
pgls_fdom7$aic

pgls_fdom0 <-pgls(logFdom ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom0)
pgls_fdom0$aic

x11()
par(mfrow=c(2,2))
plot(pgls_fdom6)

#O melhor modelo pelo AICc foi o modelo 6

plot(logFdom ~ logMass, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pgls_fdom6, col = "gray70", lwd = 3)

# Modelo frequência máxima #

pgls_fmax1 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax1)
pgls_fmax1$aic

pgls_fmax2 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax2)
pgls_fmax2$aic

pgls_fmax3 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Nares + Beak.Depth, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax3)
pgls_fmax3$aic

pgls_fmax4 <- pgls(logFmax ~ logMass + bio6 + Beak.Length_Nares + Beak.Depth, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax4)
pgls_fmax4$aic

pgls_fmax5 <- pgls(logFmax ~ logMass + bio6 + Beak.Depth, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax5)
pgls_fmax5$aic

pgls_fmax6 <- pgls(logFmax ~ logMass + bio6, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax6)
pgls_fmax6$aic

pgls_fmax7 <- pgls(logFmax ~ bio6, data = Pipridae.caper, delta= "ML")
summary(pgls_fmax7)
pgls_fmax7$aic

pgls_fmax0 <-pgls(logFmax ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax0)
pgls_fmax0$aic

#Nenhum modelo foi adequado


#Modelo frequência mínima #

pgls_fmin1 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin1)
pgls_fmin1$aic

pgls_fmin2 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin2)
pgls_fmin2$aic

pgls_fmin3 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin3)
pgls_fmin3$aic

pgls_fmin4 <- pgls(logFmin ~ logMass + bio6 + Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin4)
pgls_fmin4$aic

pgls_fmin5 <- pgls(logFmin ~ logMass + bio6 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin5)
pgls_fmin5$aic

pgls_fmin6 <- pgls(logFmin ~ bio6 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin6)
pgls_fmin6$aic

pgls_fmin7 <- pgls(logFmin ~ Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin7)
pgls_fmin7$aic

pgls_fmin0 <-pgls(logFmin ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin0)
pgls_fmin0$aic

#Modelo pace #

pgls_pace1 <- pgls(logPace ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace1)
pgls_pace1$aic

pgls_pace2 <- pgls(logPace ~ logMass + bio5 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace2)
pgls_pace2$aic

pgls_pace3 <- pgls(logPace ~ logMass + bio5 + Beak.Length_Culmen + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace3)
pgls_pace3$aic

pgls_pace4 <- pgls(logPace ~ logMass + Beak.Length_Culmen + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace4)
pgls_pace4$aic

##*
pgls_pace5 <- pgls(logPace ~ logMass + Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace5)
pgls_pace5$aic

pgls_pace6 <- pgls(formula = logPace ~ Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, lambda = "ML", delta = "ML")
summary(pgls_pace6)
pgls_pace6$aic

pgls_pace7 <- pgls(logPace ~ Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace7)
pgls_pace7$aic

pgls_pace0 <- pgls(logPace ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace0)
pgls_pace0$aic


x11()
par(mfrow=c(2,2))
plot(pgls_pace7)

plot(logPace ~ logMass, data = Pipridae.stats, pch = 19, cex = 1.2)
abline(pgls_pace5, col = "gray70", lwd = 3)

x11()
par(mfrow=c(2,2))
plot(pgls_pace7)


#Modelo Banda #

pgls_ban1 <- pgls(logBan ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban1)
pgls_ban1$aic

pgls_ban2 <- pgls(logBan ~ logMass + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban2)
pgls_ban2$aic

pgls_ban3 <- pgls(logBan ~ logMass + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban3)
pgls_ban3$aic

pgls_ban4 <- pgls(logBan ~ bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban4)
pgls_ban4$aic

pgls_ban5 <- pgls(logBan ~ bio6 + Beak.Length_Culmen + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban5)
pgls_ban5$aic

pgls_ban6 <- pgls(logBan ~ bio6 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban6)
pgls_ban6$aic

pgls_ban7 <- pgls(logBan ~ Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban7)
pgls_ban7$aic

pgls_ban0 <- pgls(logBan ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban0)
pgls_ban0$aic

#Modelo Duração #

pgls_dur1 <- pgls(logDur ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur1)
pgls_dur1$aic

pgls_dur2 <- pgls(logDur ~ logMass + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur2)
pgls_dur2$aic

pgls_dur3 <- pgls(logDur ~ logMass + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur3)
pgls_dur3$aic

pgls_dur4 <- pgls(logDur ~ logMass + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur4)
pgls_dur4$aic

##* 
pgls_dur5 <- pgls(logDur ~ logMass + Beak.Length_Culmen + Beak.Length_Nares, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur5)
pgls_dur5$aic

pgls_dur6 <- pgls(logDur ~ logMass + Beak.Length_Nares, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur6)
pgls_dur6$aic

pgls_dur7 <- pgls(logDur ~ Beak.Length_Nares, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur7)
pgls_dur7$aic

pgls_dur0 <- pgls(logDur ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur0)
pgls_dur0$aic

x11()
par(mfrow=c(2,2))
plot(pgls_dur5)


#Modelo Número de Notas #

pgls_not1 <- pgls(NNOT ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not1)
pgls_not1$aic

pgls_not2 <- pgls(NNOT ~ logMass + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not2)
pgls_not2$aic

pgls_not3 <- pgls(NNOT ~ logMass + bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not3)
pgls_not3$aic

pgls_not4 <- pgls(NNOT ~ logMass + bio6 + Beak.Length_Culmen + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not4)
pgls_not4$aic

pgls_not5 <- pgls(NNOT ~ logMass + bio6 + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not5)
pgls_not5$aic

pgls_not6 <- pgls(NNOT ~ logMass + bio6, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not6)
pgls_not6$aic

pgls_not7 <- pgls(NNOT ~ logMass, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not7)
pgls_not7$aic

pgls_not0 <- pgls(NNOT ~ 1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not0)
pgls_not0$aic

###### MODELOS NOVOS COM DICROMATISMO ######

#Modelo frequência dominante

pgls_fdom_dicro1 <- pgls(logFdom ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_dicro1)
pgls_fdom_dicro1$aic

pgls_fdom_dicro2 <- pgls(logFdom ~ logMass + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_dicro2)
pgls_fdom_dicro2$aic

pgls_fdom_dicro3 <- pgls(logFdom ~ logMass + bio6 + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_dicro3)
pgls_fdom_dicro3$aic

pgls_fdom_dicro4 <- pgls(logFdom ~ logMass + bio6 + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_dicro4)
pgls_fdom_dicro4$aic

pgls_fdom_dicro5 <- pgls(logFdom ~ logMass + bio6 + Beak.Width + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_dicro5)
pgls_fdom_dicro5$aic

###////////////////////

#Modelo frequência máxima #

pgls_fmax_dicro1 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_dicro1)
pgls_fmax_dicro1$aic

pgls_fmax_dicro2 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_dicro2)
pgls_fmax_dicro2$aic

###///////////////////////


#Modelo frequência mínima #

pgls_fmin_dicro1 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_dicro1)
pgls_fmin_dicro1$aic

pgls_fmin_dicro2 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_dicro2)
pgls_fmin_dicro2$aic

pgls_fmin_dicro3 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_dicro3)
pgls_fmin_dicro3$aic

pgls_fmin_dicro4 <- pgls(logFmin ~ logMass + bio5 + Beak.Length_Culmen + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_dicro4)
pgls_fmin_dicro4$aic

pgls_fmin_dicro5 <- pgls(logFmin ~ logMass + bio5 + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_dicro5)
pgls_fmin_dicro5$aic

pgls_fmin_dicro6 <- pgls(logFmin ~ bio5 + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_dicro6)
pgls_fmin_dicro6$aic

###////////////////

##Modelo banda #

pgls_ban_dicro1 <- pgls(logBan ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_dicro1)
pgls_ban_dicro1$aic

pgls_ban_dicro2 <- pgls(logBan ~ bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_dicro2)
pgls_ban_dicro2$aic

###//////////////


##Modelo pace #

pgls_pace_dicro1 <- pgls(logPace ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_dicro1)
pgls_pace_dicro1$aic

pgls_pace_dicro2 <- pgls(logPace ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_dicro2)
pgls_pace_dicro2$aic

###///////////////


##Modelo duração

pgls_dur_dicro1 <- pgls(logDur ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_dicro1)
pgls_dur_dicro1$aic

pgls_dur_dicro2 <- pgls(logDur ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_dicro2)
pgls_dur_dicro2$aic

pgls_dur_dicro3 <- pgls(logDur ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Width + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_dicro3)
pgls_dur_dicro3$aic

pgls_dur_dicro4 <- pgls(logDur ~ logMass + bio5 + bio6+ Beak.Width + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_dicro4)
pgls_dur_dicro4$aic

pgls_dur_dicro4 <- pgls(logDur ~ logMass + bio6 + Beak.Width + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_dicro4)
pgls_dur_dicro4$aic

pgls_dur_dicro5 <- pgls(logDur ~ logMass + bio6 + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_dicro5)
pgls_dur_dicro5$aic

###////////////



##Modelo número de notas #

pgls_not_dicro1 <- pgls(NNOT ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_dicro1)
pgls_not_dicro1$aic

pgls_not_dicro2 <- pgls(NNOT ~ bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_dicro2)
pgls_not_dicro2$aic

pgls_not_dicro3 <- pgls(NNOT ~ bio5 + bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_dicro3)
pgls_not_dicro3$aic

pgls_not_dicro4 <- pgls(NNOT ~ bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_dicro4)
pgls_not_dicro4$aic

pgls_not_dicro5 <- pgls(NNOT ~ bio6 + Beak.Width + Beak.Depth + Dicromatismo, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_dicro5)
pgls_not_dicro5$aic

###/////////////


###### MODELOS NOVOS COM DISPLAY ######

#Modelo frequência dominante

pgls_fdom_disp1 <- pgls(logFdom ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp1)
pgls_fdom_disp1$aic

pgls_fdom_disp2 <- pgls(logFdom ~ bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp2)
pgls_fdom_disp2$aic

pgls_fdom_disp3 <- pgls(logFdom ~ bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp3)
pgls_fdom_disp3$aic

pgls_fdom_disp4 <- pgls(logFdom ~ + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp4)
pgls_fdom_disp4$aic

pgls_fdom_disp5 <- pgls(logFdom ~ + bio6 + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp5)
pgls_fdom_disp5$aic

pgls_fdom_disp6 <- pgls(logFdom ~ + bio6 + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp6)
pgls_fdom_disp6$aic

pgls_fdom_disp7 <- pgls(logFdom ~ + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_disp7)
pgls_fdom_disp7$aic

###/////////////////

#Modelo frequência máxima

pgls_fmax_disp1 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_disp1)
pgls_fmax_disp1$aic

pgls_fmax_disp2 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_disp2)
pgls_fmax_disp2$aic

pgls_fmax_disp3 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_disp3)
pgls_fmax_disp3$aic

pgls_fmax_disp4 <- pgls(logFmax ~ logMass + bio5 + bio6 + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_disp4)
pgls_fmax_disp4$aic


###/////////////////

#Modelo frequência mínima

pgls_fmin_disp1 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp1)
pgls_fmin_disp1$aic

pgls_fmin_disp2 <- pgls(logFmin ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp2)
pgls_fmin_disp2$aic

pgls_fmin_disp3 <- pgls(logFmin ~ bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp3)
pgls_fmin_disp3$aic

pgls_fmin_disp4 <- pgls(logFmin ~ bio5 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp4)
pgls_fmin_disp4$aic

pgls_fmin_disp5 <- pgls(logFmin ~ bio5 + Beak.Length_Nares + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp5)
pgls_fmin_disp5$aic

pgls_fmin_disp6 <- pgls(logFmin ~ bio5 + Beak.Length_Nares + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp6)
pgls_fmin_disp6$aic

pgls_fmin_disp7 <- pgls(logFmin ~ Beak.Length_Nares + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_disp7)
pgls_fmin_disp7$aic

###////////////////////

#Modelo pace

pgls_pace_disp1 <- pgls(logPace ~ logMass + bio5 + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_disp1)
pgls_pace_disp1$aic

pgls_pace_disp2 <- pgls(logPace ~ logMass + bio6 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_disp2)
pgls_pace_disp2$aic

pgls_pace_disp3 <- pgls(logPace ~ logMass + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_disp3)
pgls_pace_disp3$aic

pgls_pace_disp4 <- pgls(logPace ~ Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth + Display.Complex, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_disp4)
pgls_pace_disp4$aic

###/////////////////////
###### MODELOS UTILIZANDO APENAS VARIAVEIS COM EXPLICAÇÃO BIOLÓGICA ######

pgls_biodur1 <- pgls(logDur ~ logMass + bio5 + bio6, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_biodur1)
pgls_biodur1$aic

pgls_biodur2 <- pgls(logDur ~ logMass + bio6, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_biodur2)
pgls_biodur2$aic

pgls_biodur3 <- pgls(logDur ~ logMass, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_biodur3)
pgls_biodur3$aic

Pipridae.stats %>% 
  ggplot(aes(x = logMass, y = logDur)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "log(Mass)", y = "log(Dur)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())



###### MODELOS COM PCA BIOCLIMATICA ######


### FDOM
pgls_fdom_pca1 <- pgls(logFdom ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_pca1)
pgls_fdom_pca1$aic

pgls_fdom_pca2 <- pgls(logFdom ~ logMass + BioclimPC1 + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_pca2)
pgls_fdom_pca2$aic

pgls_fdom_pca3 <- pgls(logFdom ~ BioclimPC1 + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_pca3)
pgls_fdom_pca3$aic

pgls_fdom_pca4 <- pgls(logFdom ~ BioclimPC1 + Beak.Length_Nares + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fdom_pca4)
pgls_fdom_pca4$aic


### FMAX
pgls_fmax_pca1 <- pgls(logFmax ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_pca1)
pgls_fmax_pca1$aic

pgls_fmax_pca2 <- pgls(logFmax ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_pca2)
pgls_fmax_pca2$aic

pgls_fmax_pca3 <- pgls(logFmax ~ logMass + BioclimPC1 + Beak.Length_Nares + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_pca3)
pgls_fmax_pca3$aic

pgls_fmax_pca4 <- pgls(logFmax ~ logMass + BioclimPC1 + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_pca4)
pgls_fmax_pca4$aic

pgls_fmax_pca5 <- pgls(logFmax ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_pca5)
pgls_fmax_pca5$aic

pgls_fmax_pca6 <- pgls(logFmax ~ BioclimPC1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmax_pca6)
pgls_fmax_pca6$aic



### FMIN
pgls_fmin_pca1 <- pgls(logFmin ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_pca1)
pgls_fmin_pca1$aic

pgls_fmin_pca2 <- pgls(logFmin ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_pca2)
pgls_fmin_pca2$aic

pgls_fmin_pca3 <- pgls(logFmin ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_pca3)
pgls_fmin_pca3$aic

pgls_fmin_pca4 <- pgls(logFmin ~ logMass + BioclimPC1 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_pca4)
pgls_fmin_pca4$aic

pgls_fmin_pca5 <- pgls(logFmin ~ BioclimPC1 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_pca5)
pgls_fmin_pca5$aic

pgls_fmin_pca6 <- pgls(logFmin ~ BioclimPC1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_fmin_pca6)
pgls_fmin_pca6$aic


### BAN
pgls_ban_pca1 <- pgls(logBan ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_pca1)
pgls_ban_pca1$aic

pgls_ban_pca2 <- pgls(logBan ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_pca2)
pgls_ban_pca2$aic

pgls_ban_pca3 <- pgls(logBan ~ BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_pca3)
pgls_ban_pca3$aic

pgls_ban_pca4 <- pgls(logBan ~ BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_pca4)
pgls_ban_pca4$aic

pgls_ban_pca5 <- pgls(logBan ~ BioclimPC1 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_ban_pca5)
pgls_ban_pca5$aic


### DUR
pgls_dur_pca1 <- pgls(logDur ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_dur_pca1)
pgls_dur_pca1$aic


### PACE
pgls_pace_pca1 <- pgls(Pace ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_pca1)
pgls_pace_pca1$aic

pgls_pace_pca2 <- pgls(Pace ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_pca2)
pgls_pace_pca2$aic

pgls_pace_pca3 <- pgls(Pace ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_pca3)
pgls_pace_pca3$aic

pgls_pace_pca4 <- pgls(Pace ~ BioclimPC1 + Beak.Length_Culmen + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_pca4)
pgls_pace_pca4$aic

pgls_pace_pca5 <- pgls(Pace ~ BioclimPC1 + Beak.Length_Culmen, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_pace_pca5)
pgls_pace_pca5$aic


### NNOT
pgls_not_pca1 <- pgls(NNOT ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Length_Nares + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_pca1)
pgls_not_pca1$aic

pgls_not_pca2 <- pgls(NNOT ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Width + Beak.Depth, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_pca2)
pgls_not_pca2$aic

pgls_not_pca3 <- pgls(NNOT ~ logMass + BioclimPC1 + Beak.Length_Culmen + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_pca3)
pgls_not_pca3$aic

pgls_not_pca4 <- pgls(NNOT ~ logMass + BioclimPC1 + Beak.Width, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_pca4)
pgls_not_pca4$aic

pgls_not_pca5 <- pgls(NNOT ~ logMass + BioclimPC1, data = Pipridae.caper, delta= "ML", lambda = "ML")
summary(pgls_not_pca5)
pgls_not_pca5$aic

###### Modelos significativos ######


## Fdom ~ nares + width, R² = 0,1
summary(pgls_fdom6)
pgls_fdom6$aic
## Coefficients:
##Estimate Std. Error t value Pr(>|t|)    
##(Intercept)        8.382230   0.478134 17.5311  < 2e-16 ***
##  Beak.Length_Nares  0.143139   0.065393  2.1889  0.03435 *  
##  Beak.Width        -0.301139   0.123190 -2.4445  0.01889 * 


Pipridae.stats %>% 
  ggplot(aes(x = Beak.Width, y = logFdom)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "Beak width", y = "log(FDom)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())

Pipridae.stats %>% 
  ggplot(aes(x = Beak.Depth, y = logFdom)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "Beak depth", y = "log(FDom)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())



## Pace ~ mass + culmen + depth, R² = 0,16
summary(pgls_pace5)
pgls_pace5$aic
##Coefficients:
##Estimate Std. Error t value  Pr(>|t|)    
##(Intercept)         3.000023   0.618938  4.8470 2.395e-05 ***
##  logMass            -0.217290   0.394074 -0.5514    0.5848    
##Beak.Length_Culmen -0.152380   0.072931 -2.0894    0.0438 *  
##  Beak.Depth          0.227251   0.219146  1.0370    0.3067  


Pipridae.stats %>% 
  ggplot(aes(x = Beak.Length_Culmen, y = logPace)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "Beak Length", y = "log(Pace)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())



## Dur ~ mass + bio6 (tmin), R² = 0,19
summary(pgls_biodur2)
pgls_biodur2$aic
##Coefficients:
##  Estimate Std. Error t value  Pr(>|t|)    
##(Intercept) -3.7073864  1.0240632 -3.6203 0.0008014 ***
##  logMass      0.8253748  0.3446334  2.3949 0.0212759 *  
##  bio6         0.0062535  0.0022308  2.8032 0.0076938 ** 


Pipridae.stats %>% 
  ggplot(aes(x = bio6, y = logDur)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "Bio6", y = "log(Dur)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())

Pipridae.stats %>% 
  ggplot(aes(x = logMass, y = logDur)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "log(Mass)", y = "log(Dur)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())


## Dur ~ mass + culmen + nares, R² = 0,39
summary(pgls_dur5)
pgls_dur5$aic
##Coefficients:
##  Estimate Std. Error t value  Pr(>|t|)    
##(Intercept)        -3.82035    0.65562 -5.8270 1.183e-06 ***
##  logMass             0.93894    0.39213  2.3945  0.021975 *  
##  Beak.Length_Culmen -0.23630    0.13704 -1.7244  0.093222 .  
##Beak.Length_Nares   0.51236    0.16439  3.1168  0.003583 ** 


Pipridae.stats %>% 
  ggplot(aes(x = Beak.Length_Nares, y = logDur)) + 
  geom_point(alpha = 0.5, color = "darkgray") + 
  labs(x = "Beak lenght", y = "log(Dur)") + 
  geom_smooth(method = "lm", color = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())


plot(predict(pgls_fdom6, se.fit = TRUE))
predict_response(pgls_fdom6)
