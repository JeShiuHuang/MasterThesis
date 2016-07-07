table = read.table("/Users/josh/Desktop/1000_genome_country/country",header = FALSE)
TW_row = data.frame(c(rep("TW",496)))
#colnames(TW_row) = "V1"


for (i in 1:22){
  onethousandgenome = paste("sort_chr",i,"_sum",sep = "")
  root = paste("/Users/josh/Desktop/each_chr_of_samples_mutate/",onethousandgenome,sep = "")
  global = read.table(root,header = FALSE,sep = "\t")
  table = cbind(table,global[,2])
}

for(i in 1:22){
  taionethousand = paste("chr",i,"_sum",sep= "")
  root2 = paste("/Users/josh/Desktop/tai_mutate_count/",taionethousand,sep = "")
  tai = read.table(root2,header = FALSE,sep = "\t")
  TW_row = cbind(TW_row,tai[,2])
}

colnames(table) = c("country","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
colnames(TW_row) = c("country","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
table = rbind(table,TW_row)


#different country people_count
AMR_people = 347
EUR_people = 503
AFR_people = 661
SAS_people = 489
EAS_people = 504
TW_people = 496

#Normalized variant by people_count
''
# table[table$country=="AMR",2:length(table)] = table[table$country=="AMR",2:length(table)]/AMR_people
# table[table$country=="EUR",2:length(table)] = table[table$country=="EUR",2:length(table)]/EUR_people
# table[table$country=="AFR",2:length(table)] = table[table$country=="AFR",2:length(table)]/AFR_people
# table[table$country=="SAS",2:length(table)] = table[table$country=="SAS",2:length(table)]/SAS_people
# table[table$country=="EAS",2:length(table)] = table[table$country=="EAS",2:length(table)]/EAS_people
# table[table$country=="TW",2:length(table)] = table[table$country=="TW",2:length(table)]/TW_people
''

#change here -------------------------------->
column = 23

#and here ------------------------------------>
interval = seq(40000,80000,5000)

#chart main topic
main_name = colnames(table[column])

#find the count between each interval
AMR_table = c(table(cut(table[table$country=="AMR",column],breaks = interval)))
EUR_table = c(table(cut(table[table$country=="EUR",column],breaks = interval)))
AFR_table = c(table(cut(table[table$country=="AFR",column],breaks = interval)))
SAS_table = c(table(cut(table[table$country=="SAS",column],breaks = interval)))
EAS_table = c(table(cut(table[table$country=="EAS",column],breaks = interval)))
TW_table = c(table(cut(table[table$country=="TW",column],breaks = interval)))


#creat a vector to store each count
AMR_column = c()
EUR_column = c()
AFR_column = c()
SAS_column = c()
EAS_column = c()
TW_column = c()

for (x in AMR_table){
  AMR_column = append(AMR_column,x/AMR_people)
  #AMR_column = append(AMR_column,x)
}
for (x in EUR_table){
  EUR_column = append(EUR_column,x/EUR_people)
  #EUR_column = append(EUR_column,x)
}
for (x in AFR_table){
  AFR_column = append(AFR_column,x/AFR_people)
  #AFR_column = append(AFR_column,x)
}
for (x in SAS_table){
  SAS_column = append(SAS_column,x/SAS_people)
  #SAS_column = append(SAS_column,x)
}
for (x in EAS_table){
  EAS_column = append(EAS_column,x/EAS_people)
  #EAS_column = append(EAS_column,x)
}
for (x in TW_table){
  TW_column = append(TW_column,x/TW_people)
  #TW_column = append(TW_column,x)
}


proportional_table = cbind(AMR_column,EUR_column,AFR_column,SAS_column,EAS_column,TW_column)
colnames(proportional_table) = c("AMR","EUR","AFR","SAS","EAS","TW")

y = interval[1:length(interval)-1]

matplot(y,proportional_table,pch=18,type = "l",lty = 1,col = 1:6,ylab = "proportional",xlab = "variant count")
legend("topleft",c("AMR","EUR","AFR","SAS","EAS","TW"),lty = 1,col = 1:6,cex = 0.8)
title(main = main_name,cex.main = 2)