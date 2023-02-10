getRareVariants = function(filename){
    
    #В первый раз из таблицы, забираем строку заголовка, которая по длине <= количество столбцов
    annovar = read.csv2(filename, header = F,sep = "\t",na.strings = ".", stringsAsFactors=FALSE,dec = ".")
    column.names = as.vector(as.matrix(annovar[1,]))
    
    #Считываем заново, т.к. read.csv не парсит адекватно числа
    annovar = as.data.frame(data.table::fread(filename,header = F,dec = ".",na.strings = ".",sep = "\t"))
    colnames(annovar)[1:151] = column.names[1:151]# Берем только ту часть, которая точно всегда есть
    
    annovar[annovar["gnomAD_exome_ALL"] < 0.05 | is.na(annovar["gnomAD_exome_ALL"]),]
}

runs = c("wbc.4739"="epi-2022-BRCA-4739-wbc/epi-2022-BRCA-4739-wbc.annovar.tsv", 
         "wbc.4740"="epi-2022-BRCA-4740-wbc/epi-2022-BRCA-4740-wbc.annovar.tsv")
mut.cols = c("Chr","Start","End","Ref","Alt")
rares = lapply(runs, getRareVariants)
rares = lapply(rares, function(x){x[,mut.cols]})
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
print(m.inter.SNP)