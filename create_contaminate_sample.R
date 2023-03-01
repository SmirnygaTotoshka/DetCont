dirs = list.dirs(recursive = F)
sample.dirs = stringr::str_starts(dirs,"./epi")
sample.dirs = dirs[sample.dirs]
names = paste0(stringr::str_sub(sample.dirs,start = 3),".fastq")
fastq.paths = paste(sample.dirs,names,sep = "/")

set.seed(100)
brca.4738 = ShortRead::readFastq(fastq.paths[1])
print(brca.4738[1:10])
ShortRead::sread(brca.4738[1:10])

brca = fastq.paths[grepl("BRCA",fastq.paths)]
ccp = fastq.paths[grepl("CCP",fastq.paths)]
cont.brca.comb = as.data.frame(expand.grid(brca,brca))
colnames(cont.brca.comb) = c("Sample","Cont_source")
cont.brca.comb = cont.brca.comb[cont.brca.comb$Sample != cont.brca.comb$Cont_source,]

cont.brca.comb$Sample = as.character(cont.brca.comb$Sample)
cont.brca.comb$Cont_source = as.character(cont.brca.comb$Cont_source)
cont.brca.comb$Sample.id = stringr::str_remove(stringr::str_split(cont.brca.comb$Sample,"/",simplify = T)[,2],"epi-2022-")
cont.brca.comb$Cont_source.id = stringr::str_remove(stringr::str_split(cont.brca.comb$Cont_source,"/",simplify = T)[,2],"epi-2022-")
procent = c(0.01,0.05,0.1,0.25,0.5)
for (i in seq_len(nrow(cont.brca.comb))) {
    for(p in procent){
       name = paste0(cont.brca.comb[i,"Sample.id"],"-",cont.brca.comb[i,"Cont_source.id"],"-",p,".fastq") 
       sample = ShortRead::readFastq(cont.brca.comb[i,"Sample"])
       cont.source = ShortRead::readFastq(cont.brca.comb[i,"Cont_source"])
       l = min(length(sample),length(cont.source))
       ind = sample(1:l, size = p * l)
       sample[ind] = cont.source[ind]
       ShortRead::writeFastq(sample, paste0("contaminate_fastq/", name), compress = F)
    }
}
