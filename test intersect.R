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

runs.pcblt = c("pcblt.4738"="epi-2022-BRCA-4738-pcblt/epi-2022-BRCA-4738-pcblt.annovar.tsv",
               "pcblt.4741"="epi-2022-BRCA-4741-pcblt/epi-2022-BRCA-4741-pcblt.annovar.tsv",
               "pcblt.4742"="epi-2022-BRCA-4742-pcblt/epi-2022-BRCA-4742-pcblt.annovar.tsv")
mut.cols = c("Chr","Start","End","Ref","Alt")
rares = lapply(runs.pcblt, getRareVariants)
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
m.inter.SNP = matrix(num.intersect.SNP,nrow = length(runs.pcblt),byrow = T)
colnames(m.inter.SNP) = names(runs.pcblt)
row.names(m.inter.SNP) = names(runs.pcblt)
print(m.inter.SNP)
xlsx::write.xlsx(m.inter.SNP,"test_pcblt.xlsx")

intersect.SNP = lapply(seq_along(rares),
                              function(i){
                                  lapply(seq_along(rares),
                                         function(j, i) {
                                             dplyr::intersect(rares[[i]], rares[[j]])
                                         }
                                         ,i)
                              }
)
for (i in names(rares)) {
    for (j in names(rares)) {
        path = paste0("intersection/",i,"-",j,".csv")
        path1 = paste0("intersection/",j,"-",i,".csv")
        if(!file.exists(path1)){
            write.csv2(dplyr::intersect(rares[[i]],rares[[j]]),path)
        }
    }
}
