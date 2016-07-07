setwd("/Users/josh/Desktop/TW_1000_variant_in_CDS/allele/merge/TW_1000_same_pos_ref_alt/full_compare_line/")
y = 8
daa_path = paste("chr",y,"_table",sep = "")
daa = read.table(daa_path,header = F,sep = "\t")


AMR_people = 347*2
EUR_people = 503*2
AFR_people = 661*2  
SAS_people = 489*2
EAS_people = 504*2
TW_people = 496*2
Total_people = 3000*2

#chr = paste("chr",y,sep = "")
#daa = cbind(c(rep(chr,dim(daa)[1])),daa)

more = c()
less = c()

for (i in 1:dim(daa)[1]){
  TW_AMR = prop.test(c(daa[i,12],daa[i,7]),c(TW_people,AMR_people))
  TW_EUR = prop.test(c(daa[i,12],daa[i,8]) ,c(TW_people,EUR_people))
  TW_AFR = prop.test(c(daa[i,12],daa[i,9]) ,c(TW_people,AFR_people))
  TW_SAS = prop.test(c(daa[i,12],daa[i,10]) ,c(TW_people,SAS_people))
  TW_EAS = prop.test(c(daa[i,12],daa[i,11]),c(TW_people,EAS_people))
  TW_Total = prop.test(c(daa[i,12],daa[i,13]) ,c(TW_people,Total_people))

  if(is.nan(TW_AMR$p.value) || is.nan(TW_EUR$p.value) || is.nan(TW_AFR$p.value) || is.nan(TW_SAS$p.value) || is.nan(TW_Total$p.value) || is.nan(TW_EAS$p.value)){
    next()    
  }
  else if(TW_Total$p.value < 0.05 & TW_EUR$p.value < 0.05 & TW_AFR$p.value < 0.05 & TW_SAS$p.value < 0.05 & TW_EAS$p.value < 0.05 &TW_Total$p.value < 0.05){
    TW_Total_greater = prop.test(c(daa[i,12],daa[i,13]) ,c(TW_people,Total_people),alternative="greater")  
    TW_Total_less = prop.test(c(daa[i,12],daa[i,13]) ,c(TW_people,Total_people),alternative="less")
    TW_AMR_greater = prop.test(c(daa[i,12],daa[i,7]),c(TW_people,AMR_people),alternative = "greater")
    TW_EUR_greater = prop.test(c(daa[i,12],daa[i,8]) ,c(TW_people,EUR_people),alternative = "greater")
    TW_AFR_greater = prop.test(c(daa[i,12],daa[i,9]) ,c(TW_people,AFR_people),alternative = "greater")
    TW_SAS_greater = prop.test(c(daa[i,12],daa[i,10]) ,c(TW_people,SAS_people),alternative = "greater")
    TW_EAS_greater = prop.test(c(daa[i,12],daa[i,11]),c(TW_people,EAS_people),alternative = "greater")
    TW_AMR_less = prop.test(c(daa[i,12],daa[i,7]),c(TW_people,AMR_people),alternative = "less")
    TW_EUR_less = prop.test(c(daa[i,12],daa[i,8]) ,c(TW_people,EUR_people),alternative = "less")
    TW_AFR_less = prop.test(c(daa[i,12],daa[i,9]) ,c(TW_people,AFR_people),alternative = "less")
    TW_SAS_less = prop.test(c(daa[i,12],daa[i,10]) ,c(TW_people,SAS_people),alternative = "less")
    TW_EAS_less = prop.test(c(daa[i,12],daa[i,11]),c(TW_people,EAS_people),alternative = "less")
    if(TW_Total_greater$p.value < 0.05 & TW_Total_less$p.value > 0.05 & TW_EUR_greater$p.value < 0.05 & TW_EUR_less$p.value > 0.05 
       & TW_AMR_greater$p.value < 0.05 & TW_AMR_less$p.value > 0.05 & TW_AFR_greater$p.value < 0.05 & TW_AFR_less$p.value > 0.05
       & TW_SAS_greater$p.value < 0.05 & TW_SAS_less$p.value > 0.05 & TW_EAS_greater$p.value < 0.05 & TW_EAS_less$p.value > 0.05){
      TW_AMR_pvalue = -1*log(TW_AMR_greater$p.value)
      TW_EUR_pvalue = -1*log(TW_EUR_greater$p.value)
      TW_AFR_pvalue = -1*log(TW_AFR_greater$p.value)
      TW_SAS_pvalue = -1*log(TW_SAS_greater$p.value)
      TW_EAS_pvalue = -1*log(TW_EAS_greater$p.value)
      TW_Total_pvalue = -1*log(TW_Total_greater$p.value)
      a = cbind(daa[i,1:6],TW_AMR_pvalue,TW_EUR_pvalue,TW_AFR_pvalue,TW_SAS_pvalue,TW_EAS_pvalue,TW_Total_pvalue,daa[i,12]/TW_people)
      more = rbind(more,a)
    }
    else if(TW_Total_greater$p.value > 0.05 & TW_Total_less$p.value < 0.05 & TW_EUR_greater$p.value > 0.05 & TW_EUR_less$p.value < 0.05 
            & TW_AMR_greater$p.value > 0.05 & TW_AMR_less$p.value < 0.05 & TW_AFR_greater$p.value > 0.05 & TW_AFR_less$p.value < 0.05
            & TW_SAS_greater$p.value > 0.05 & TW_SAS_less$p.value < 0.05 & TW_EAS_greater$p.value > 0.05 & TW_EAS_less$p.value < 0.05){
      TW_AMR_pvalue = log(TW_AMR_less$p.value)
      TW_EUR_pvalue = log(TW_EUR_less$p.value)
      TW_AFR_pvalue = log(TW_AFR_less$p.value)
      TW_SAS_pvalue = log(TW_SAS_less$p.value)
      TW_EAS_pvalue = log(TW_EAS_less$p.value)
      TW_Total_pvalue = log(TW_Total_less$p.value)
      b = cbind(daa[i,1:6],TW_AMR_pvalue,TW_EUR_pvalue,TW_AFR_pvalue,TW_SAS_pvalue,TW_EAS_pvalue,TW_Total_pvalue,daa[i,12]/TW_people)
      less = rbind(less,b)
    }
  }  
}

more[more=="Inf"] <- -1000
more[more==-1000] = max(more[,7:12])



less[less=="-Inf"] <- 1000
less[less== 1000] = min(less[,7:12])

#more[,7] = more[,7]/TW_people
#less[,7] = less[,7]/TW_people

more_path = paste("/Users/josh/Desktop/TW_more_in_CDS/chr",y,sep = "")
write.table(more,file = more_path,sep = "\t",col.names = F,row.names = F,quote = F)
less_path = paste("/Users/josh/Desktop/TW_less_in_CDS/chr",y,sep = "")
write.table(less,file = less_path,sep = "\t",col.names = F,row.names = F,quote = F)
