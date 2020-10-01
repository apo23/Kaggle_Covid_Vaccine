setwd("Documents/M2/ML/Projet/")


library(jsonlite)
require(MASS)

out=lapply(readLines("test.json"),fromJSON)


df <- data.frame(Reduce ( rbind, out))
head(df)
dim(df)

#Conversion du df pour le rendre plus lisible
conversion <- function(df) {
  seq_score <- rep(0, dim(df)[1])
  seq_len <- rep(0, dim(df)[1])
  sequence <- rep(0, dim(df)[1])
  structure <- rep(0, dim(df)[1])
  predicted_loop_type <- rep(0, dim(df)[1])
  for(i in 1:dim(df)[1]) {
    seq_score[i] <- df$seq_scored[[i]]
    seq_len[i] <- df$seq_length[[i]]
    structure[i] <- as.character(df$structure[[i]])
    sequence[i] <- as.character(df$sequence[[i]])
    predicted_loop_type[i] <- as.character(df$predicted_loop_type[[i]])
  }
  rdf=as.data.frame(cbind(sequence, structure, predicted_loop_type, seq_len, seq_score))
  rdf$sequence <- as.character(rdf$sequence)
  rdf$structure <- as.character(rdf$structure)
  rdf$predicted_loop_type <- as.character(rdf$predicted_loop_type)
  return(rdf)
}


df2 <- conversion(df)


#Analyse de la taille et score des séquences
seq_score <- as.integer(as.vector(df2$seq_score))
seq_len <- as.integer(as.vector(df2$seq_len))
sum(seq_score == min(seq_score))
sum(seq_len == min(seq_len))
par(mfrow=c(1,2))
truehist(seq_score, col="#a182a3", main = "Distribution score")
truehist(seq_len, col="#f9ec57", main = "Distribution longeur")
#Distribution équivalente entre la taille et le score de laséqeunce donc une des deux variables inutile
#VAriable a 2 classe


##Etude de la répartition des acides aminées :

aa_detail <- function(AA_chain) {
  result = rep(0,4)
  #Pass la chaine en un vecteur de taillelongeur de chaine avec chaque élément un AA
  vec = substring(AA_chain, seq(1, nchar(AA_chain), 1), seq(1, nchar(AA_chain), 1))
  #Comptabilise le nombre de chaque AA dans un vecteur
  result[1] = sum(vec == "A")
  result[2] = sum(vec == "U")
  result[3] = sum(vec == "G")
  result[4] = sum(vec == "C")
  return(result)
}

#Creation d'une matrice contenant 4 colonne pour les AA et a chaque ligne
# Le nombre de chaque AA présent dans une chaine
tmpvec <- rep(0,4)
AA_visu = matrix(data=NA, nrow=0, ncol=4)
colnames(AA_visu) = c("A", "U", "G", "C")

for(i in 1:dim(df2)[1]){
  AA_visu = rbind(AA_visu, aa_detail(df2$sequence[i]))
}
rownames(AA_visu)=1:dim(df)[1]


par(mfrow=c(2,2))
truehist(AA_visu[,1], col="steelblue", main="Histogramme des A", xlim=c(0,max(AA_visu)),
         xlab = paste("Mean  = ",round(mean(AA_visu[,1]),2)))
abline(v=mean(AA_visu[,1]), col="red", lty=5, lwd=3)
truehist(AA_visu[,2], col="#ff8e32", main="Histogramme des U", xlim=c(0,max(AA_visu)),
         xlab = paste("Mean  = ",round(mean(AA_visu[,2]),2)))
abline(v=mean(AA_visu[,2]), col="red", lty=5, lwd=2)
truehist(AA_visu[,3], col="#a182a3", main="Histogramme des G", xlim=c(0,max(AA_visu)),
         xlab = paste("Mean  = ",round(mean(AA_visu[,3]),2)))
abline(v=mean(AA_visu[,3]), col="red", lty=5, lwd=2)
truehist(AA_visu[,4], col="#f9ec57", main="Histogramme des C", xlim=c(0,max(AA_visu)),
         xlab = paste("Mean  = ",round(mean(AA_visu[,4]),2)))
abline(v=mean(AA_visu[,4]), col="red", lty=5, lwd=2)


#Etude des structures secondaires

ss_detail <- function(AA_chain) {
  result = rep(0,5)
  #Pass la chaine en un vecteur de taillelongeur de chaine avec chaque élément est la SS d'un AA
  vec = substring(AA_chain, seq(1, nchar(AA_chain), 1), seq(1, nchar(AA_chain), 1))
  #Comptabilise le nombre de chaque AA dans un vecteur
  result[1] = sum(vec == "E")
  result[2] = sum(vec == "S")
  result[3] = sum(vec == "B")
  result[4] = sum(vec == "H")
  result[5] = sum(vec == "X")
  return(result)
}


tmpvec <- rep(0,5)
SS_visu = matrix(data=NA, nrow=0, ncol=5)
colnames(SS_visu) = c("E", "S", "B", "H", "X")

for(i in 1:dim(df2)[1]){
  SS_visu = rbind(SS_visu, ss_detail(df2$predicted_loop_type[i]))
}
rownames(SS_visu)=1:dim(df)[1]
SS_visu




par(mfrow=c(1,2))
truehist(SS_visu[,1], col="steelblue", main="Histogramme des E", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,1]),2)), prob = FALSE)
abline(v=mean(SS_visu[,1]), col="red", lty=5, lwd=3)
truehist(SS_visu[,2], col="#ff8e32", main="Histogramme des S", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,2]),2)), prob = FALSE)
abline(v=mean(SS_visu[,2]), col="red", lty=5, lwd=2)
truehist(SS_visu[,3], col="#a182a3", main="Histogramme des B", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,3]),2)), prob = FALSE)
abline(v=mean(SS_visu[,3]), col="red", lty=5, lwd=2)
truehist(SS_visu[,4], col="#f9ec57", main="Histogramme des H", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,4]),2)), prob = FALSE)
abline(v=mean(SS_visu[,4]), col="red", lty=5, lwd=2)
truehist(SS_visu[,5], col="violetred", main="Histogramme des X", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,5]),2)), prob = FALSE)
abline(v=mean(SS_visu[,5]), col="red", lty=5, lwd=2)

