annovar = read.csv("epi-2022-CCP-260-t/epi-2022-CCP-260-t.annovar.splice.hgvs.tsv",sep = "\t",na.strings = ".")

annovar = data.table::fread("epi-2022-CCP-260-t/epi-2022-CCP-260-t.annovar.splice.tsv",header = F,dec = ".",na.strings = ".")
colnames.annovar = colnames(read.table("annovar_colnames.txt", sep= "\t",header = T))
annovar = as.data.frame(annovar)
colnames(annovar) = colnames.annovar

rare.annovar = subset(annovar, gnomAD_exome_ALL < 0.05)

annovar = data.table::fread("epi-2022-CCP-260-t/epi-2022-CCP-260-t.annovar.splice.tsv",header = F,dec = ".",na.strings = ".",sep = "\t")
annovar[,15] > 0.05
library(dplyr)
annovar = read.csv2("epi-2022-CCP-260-t/epi-2022-CCP-260-t.annovar.splice.tsv",header = F,sep = "\t",na.strings = ".", stringsAsFactors=FALSE,dec = ".")
column.names = as.vector(as.matrix(annovar[1,]))
annovar = as.data.frame(data.table::fread("epi-2022-CCP-260-t/epi-2022-CCP-260-t.annovar.splice.tsv",header = F,dec = ".",na.strings = ".",sep = "\t"))
colnames(annovar)[1:151] = column.names[1:151]
rare.annovar = annovar[annovar["gnomAD_exome_ALL"] < 0.05 | is.na(annovar["gnomAD_exome_ALL"]),]

getRareVariants = function(filename){
    
    #В первый раз из таблицы, забираем строку заголовка, которая по длине <= количество столбцов
    annovar = read.csv2(filename, header = F,sep = "\t",na.strings = ".", stringsAsFactors=FALSE,dec = ".")
    column.names = as.vector(as.matrix(annovar[1,]))
    
    #Считываем заново, т.к. read.csv не парсит адекватно числа
    annovar = as.data.frame(data.table::fread(filename,header = F,dec = ".",na.strings = ".",sep = "\t"))
    colnames(annovar)[1:151] = column.names[1:151]# Берем только ту часть, которая точно всегда есть
    
    annovar[annovar["gnomAD_exome_ALL"] < 0.05 | is.na(annovar["gnomAD_exome_ALL"]),]
}

runs = c("wbc.4739"="epi-2022-BRCA-4739-wbc/epi-2022-BRCA-4739-wbc.annovar.splice.tsv", "wbc.4740"="epi-2022-BRCA-4740-wbc/epi-2022-BRCA-4740-wbc.annovar.splice.tsv")
mut.cols = c("Chr","Start","End","Ref","Alt")
rares = lapply(runs, getRareVariants)
rares = lapply(rares, function(x){x[,mut.cols]})
intersection = lapply(rares,dplyr::intersect)
m = matrix(nrow = length(rares),ncol = length(rares))
num.intersect.SNP = unlist(lapply(seq_along(rares),
                                  function(i){
                                    lapply(seq_along(rares),
                                           function(j, i) {
                                                nrow(dplyr::intersect(rares[[i]], rares[[j]]))   
                                           }
                                           ,i)
                                  }
                                )
                           )
m.inter.SNP = matrix(num.intersect.SNP,nrow = length(runs),byrow = T)
