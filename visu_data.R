setwd("Documents/M2/ML/Projet/")


library(jsonlite)
require(MASS)

out=lapply(readLines("train.json"),fromJSON)


df <- data.frame(Reduce ( rbind, out))
head(df)
dim(df)
dim(na.omit(df))
#Pas de NA fichier lisible et analysable


#Conversion du df pour le rendre plus lisible
conversion <- function(df) {
  nligne = dim(df)[1]
  seq_score <- rep(0, nligne)
  seq_len <- rep(0, nligne)
  sequence <- rep(0, nligne)
  structure <- rep(0, nligne)
  predicted_loop_type <- rep(0, nligne)
  reactivity <- rep(0, nligne)
  deg_Mg_pH10 <- rep(0, nligne)
  deg_ph10 <- rep(0, nligne)
  deg_Mg_50C <- rep(0, nligne)
  deg_50C <- rep(0, nligne)
  for(i in 1:nligne) {
      seq_score[i] <- df$seq_scored[[i]]
      seq_len[i] <- df$seq_length[[i]]
      structure[i] <- as.character(df$structure[[i]])
      sequence[i] <- as.character(df$sequence[[i]])
      predicted_loop_type[i] <- as.character(df$predicted_loop_type[[i]])
      reactivity[i] <- df$reactivity[i]
      deg_Mg_pH10[i] <- df$deg_Mg_pH10[i]
      deg_ph10[i] <- df$deg_pH10[i]
      deg_Mg_50C[i] <- df$deg_Mg_50C[i]
      deg_50C[i] <- df$deg_50C[i]
  }
  rdf=as.data.frame(cbind(sequence, structure, predicted_loop_type, seq_len, seq_score,
                          reactivity, deg_Mg_pH10, deg_ph10, deg_Mg_50C, deg_50C))
  rdf$sequence <- as.character(rdf$sequence)
  rdf$structure <- as.character(rdf$structure)
  rdf$predicted_loop_type <- as.character(rdf$predicted_loop_type)
  
  #Garde que les colonnes avec des noises < 1
  rdf = rdf[which(df$SN_filter == 0, arr.ind = TRUE),]
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

aa_count <- function(AA_chain) {
  result = rep(0,5)
  #Pass la chaine en un vecteur de taillelongeur de chaine avec chaque élément un AA
  vec = substring(AA_chain, seq(1, nchar(AA_chain), 1), seq(1, nchar(AA_chain), 1))
  #Comptabilise le nombre de chaque AA dans un vecteur
  result[1] = sum(vec == "A")
  result[2] = sum(vec == "U")
  result[3] = sum(vec == "G")
  result[4] = sum(vec == "C")
  result[5] = sum(tail(vec, n=20) == "A")
  return(result)
}

#Creation d'une matrice contenant 4 colonne pour les AA et a chaque ligne
# Le nombre de chaque AA présent dans une chaine
tmpvec <- rep(0,5)
AA_visu = matrix(data=NA, nrow=0, ncol=5)
colnames(AA_visu) = c("A", "U", "G", "C", "A tail")

for(i in 1:dim(df2)[1]){
  AA_visu = rbind(AA_visu, aa_count(df2$sequence[i]))
}
rownames(AA_visu)=1:dim(df2)[1]


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

#Moyen des AA A 2x supéreures aux autres. Due a des queues polyA pour stab l'ARN ?
par(mfrow = c(1,1))
truehist(AA_visu[,5], col="steelblue", main="Histogramme des A en queue (20 dernières pos)", 
         xlim=c(0,max(AA_visu)), xlab = paste("Mean  = ",round(mean(AA_visu[,5]),2)))
#On a 14 A en queue sur les 20 dernères position de chaque ARN
#Bug ?



#Etude des structures secondaires

ss_count <- function(AA_chain) {
  result = rep(0,7)
  #Passe la chaine en un vecteur de taillelongeur de chaine avec chaque élément est la SS d'un AA
  vec = substring(AA_chain, seq(1, nchar(AA_chain), 1), seq(1, nchar(AA_chain), 1))
  #Comptabilise le nombre de chaque AA dans un vecteur
  result[1] = sum(vec == "E")
  result[2] = sum(vec == "S")
  result[3] = sum(vec == "B")
  result[4] = sum(vec == "H")
  result[5] = sum(vec == "X")
  result[6] = sum(vec == "M")
  result[7] = sum(vec == "I")
  return(result)
}


SS_visu = matrix(data=NA, nrow=0, ncol=7)
colnames(SS_visu) = c("E", "S", "B", "H", "X", "M", "I")

for(i in 1:dim(df2)[1]){
  SS_visu = rbind(SS_visu, ss_count(df2$predicted_loop_type[i]))
}
rownames(SS_visu)=1:dim(df2)[1]



par(mfrow=c(2,2))
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
truehist(SS_visu[,6], col="#9b0909", main="Histogramme des HM", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,6]),2)), prob = FALSE)
abline(v=mean(SS_visu[,6]), col="red", lty=5, lwd=2)
truehist(SS_visu[,7], col="#2b2760", main="Histogramme des I", xlim=c(0,max(SS_visu)),
         xlab = paste("Mean  = ",round(mean(SS_visu[,7]),2)), prob = FALSE)
abline(v=mean(SS_visu[,7]), col="red", lty=5, lwd=2)

#Beucoup de "E" -> dangling end -> fin incertaines donc aa pas sur de la SS ?
