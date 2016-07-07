setwd("/Users/josh/Desktop/TW_1000_variant_in_CDS/allele/merge/TW_1000_same_pos_ref_alt/full_compare_line/")
y = 1
daa_path = paste("chr",y,"_table",sep = "")
daa = read.table(daa_path,header = F,sep = "\t")


AMR_people = 347*2
EUR_people = 503*2
AFR_people = 661*2  
SAS_people = 489*2
EAS_people = 504*2
TW_people = 496*2
Total_people = 3000*2

data = daa[,7:dim(daa)[2]]
AMR = data[,1]
EUR = data[,2]
AFR = data[,3]
SAS = data[,4]
EAS = data[,5]
TW = data[,6]
Total = data[,7]

judge_pvalue = function(a,b,a_p,b_p){
  a_b = prop.test(c(a,b) ,c(a_p,b_p))
  if(is.nan(a_b$p.value)){
    return(0)
  }
  else{
    if(a_b$p.value < 0.05){
      a_b_greater = prop.test(c(a,b),c(a_p,b_p),alternative="greater")  
      a_b_less = prop.test(c(a,b),c(a_p,b_p),alternative="less")
      if(a_b_greater$p.value < 0.05 & a_b_less$p.value > 0.05){
        return(-1*log(a_b_greater$p.value))
      }
      else if(a_b_greater$p.value > 0.05 & a_b_less$p.value < 0.05){
        return(log(a_b_less$p.value))
      }
    }
    else if(a_b$p.value == 0.05){
      if(a/a_p > b/b_p){
        return(-1*log(a_b$p.value))
      }
      else if(a/a_p < b/b_p){
        return(log(a_b$p.value))
      }
    }
    else{
      if(a/a_p > b/b_p){
        return(-1*log(a_b$p.value))
      }
      else if(a/a_p < b/b_p){
        return(log(a_b$p.value))
      }
    }
  }
}

Total_tt = c()
for (i in 1:dim(daa)[1]){
  Total_AMR = judge_pvalue(AMR[i],Total[i],AMR_people,Total_people)
  Total_EUR = judge_pvalue(EUR[i],Total[i],EUR_people,Total_people)
  Total_AFR = judge_pvalue(AFR[i],Total[i],AFR_people,Total_people)
  Total_SAS = judge_pvalue(SAS[i],Total[i],SAS_people,Total_people)
  Total_EAS = judge_pvalue(EAS[i],Total[i],EAS_people,Total_people)
  Total_TW = judge_pvalue(TW[i],Total[i],TW_people,Total_people)
  Total_t = cbind(Total_AMR,Total_EUR,Total_AFR,Total_SAS,Total_EAS,Total_TW)
  Total_tt = rbind(Total_tt,Total_t)
}
Total_tt[Total_tt=="-Inf"] = min(Total_tt[Total_tt < 0 & Total_tt != "-Inf"])
Total_tt[Total_tt=="Inf"] = max(Total_tt[Total_tt > 0 & Total_tt != "Inf"])
Total_tt = data.matrix(Total_tt)
colnames(Total_tt) = c("AMR","EUR","AFR","SAS","EAS","TW")
Total_path = paste("/Users/josh/Desktop/variant_ratio_in_group_show_pvalue/allele/In_CDS/Total/chr",y,".png",sep = "")
png(file=Total_path,width = 750,height = 1000)
aheatmap(Total_tt)
dev.off()
