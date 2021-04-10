# Importacion de librerias
library('stringr')
library(dplyr)
library(cluster)
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ape))

################################################ FUNCIONES ################################################

# Correr blastp
blasteo <- function(filename, pathPan){
  query <- paste('-query', paste(pathPan, filename, sep = '/'))
  subject <- paste('-subject', paste(pathPan, filename, sep = '/'))
  cleanFileName <- gsub(".faa", '', filename)
  fileOut <- paste(cleanFileName, '-blastp.out', sep= '')
  out <- paste('-out', fileOut)
  blastp <- paste('blastp', query, subject, '-num_threads 5 -outfmt 7 -max_hsps 1 -use_sw_tback', out)
  system(blastp)
  return(0)
}

#Solo se mantienen aquellas lineas que no sean comentarios
quitComments <- function(line){
  if (startsWith(x = line, prefix= '#') == FALSE) {
    return(line)}
  else {
    return(NA)}
}

#Convertir a un vector una linea separada por tabs
splitLines <- function(line) {
  splitedLine <- as.vector(str_split(line, pattern= '\t', simplify =TRUE))
  return (splitedLine)
}

#Convertir un archivo output de blastp (tabular con comentarios) a un data frame sin comentarios
fileToDataFrame <- function(filename){
  fileContent <- readLines(filename)
  cleanFileContent <- as.vector(sapply(fileContent, quitComments))
  cleanestFileContent <- cleanFileContent[!is.na(cleanFileContent)]
  fragmentedFileContent <- sapply(cleanestFileContent, splitLines)
  preDataframeFileContent <- as.data.frame(fragmentedFileContent)
  if (dim(preDataframeFileContent)[1] > 0){
    names(preDataframeFileContent) <- NULL
    dataFrameFileContent <- as.data.frame(t(preDataframeFileContent))
    names(dataFrameFileContent) <- c('query_acc.ver', 'subject_acc.ver', '%_identity', 'alignment_length', 'mismatches', 'gap_opens', 'q._start', 'q._end', 's._start', 's. end', 'evalue', 'bit_score')
    return(dataFrameFileContent)}
  else {
    return(data.frame())}
}


################################################ CODIGO PRINCIPAL ################################################

# Correr blastp para ABC.faa (en Buluc)
blasteo('ABC.faa', '/export/storage/users/aescobed/Bioinfo_LCG/Modulo-4/Proyecto')

# Lectura de la salida de blastp en un dataframe
data <- fileToDataFrame('ABC-blastp.out')

# Creacion de matriz de con los bitscore
disMatrix <- matrix(rep(0,10000), nrow= 100, ncol=100)
protsNames <- unique(data$query_acc.ver)
colnames(disMatrix) <- protsNames
row.names(disMatrix) <- protsNames
for (protA in protsNames) {
  for (protB in protsNames) {
    if (protA != protB){
      filterData <- filter(data, query_acc.ver== protA, subject_acc.ver == protB)
      if (dim(filterData)[1] != 0){
        bitScore <- filterData$bit_score[1]
        disMatrix[protA, protB] <- as.numeric(bitScore)}
      else{
        disMatrix[protA, protB] <- as.numeric(0)}
    }
    else{
      disMatrix[protA, protB] <- as.numeric(1)}
  }
}

# Normalizacion de matriz en bitscore en una matriz de disimilitud
normDisMatrix <-matrix(rep(0,10000), nrow= 100, ncol=100)
colnames(normDisMatrix) <- protsNames
row.names(normDisMatrix) <- protsNames
maxBitScore <-  max(apply(disMatrix, 1, max))
for (protA in protsNames) {
  for (protB in protsNames) {
    if (protA != protB){
      norm <- disMatrix[protA, protB]/maxBitScore
      dis <- 1- norm
      normDisMatrix[protA, protB] <- dis
    }
  }
}

# Creacion de arbol con funciones tradicionales (con metodo single)
csin <- hclust(dist(normDisMatrix, method= 'euclidean'), method = "single")
singleTree <- as.phylo(csin)
write.tree(phy=singleTree, file="single.tree")

# Creacion de arbol con funciones tradicionales (con metodo average)
cave <- hclust(dist(normDisMatrix, method= 'euclidean'), method = "average")
averageTree <- as.phylo(cave)
write.tree(phy= averageTree, file= "average.tree")

# Creacion de arbol con funciones tradicionales (con metodo complete)
ccom <- hclust(dist(normDisMatrix, method= 'euclidean'), method = "complete")
completeTree <- as.phylo(ccom)
write.tree(phy= completeTree, file= "complete.tree")

# Creacion de arbol con funciones tradicionales (con metodo ward)
cwar <- hclust(dist(normDisMatrix, method= 'euclidean'), method = "ward.D2")
wardTree <- as.phylo(cwar)
write.tree(phy= wardTree , file= "ward.tree")


