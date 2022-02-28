


library(data.table)
library(dplyr)
S7<-fread("Data/S7.csv")
S7_names<-paste(S7$`Gene symbol`,S7$`Chromosome band`,sep="_")

ss<-expand.grid(S7_names,S7_names)
ss<-ss[apply(ss,1,function(x)length(unique(x)))==2,]

ss2<-expand.grid(S7$`Gene symbol`,S7$`Gene symbol`)
ss2<-ss2[apply(ss2,1,function(x)length(unique(x)))==2,]

m<-readRDS("raw/CUP_amplified.rds")
m<-as.data.frame(m)
m$rn<-rownames(m)
setDT(m)
setkey(m,"rn")
df_list<-list()
for (i in 1:NROW(ss2)) {
  print(i)
  
  x<-unlist(ss2[i,])%>%as.character()
  #x<-c("MET","MYC")
  print(x)
  a<-m[,.SD,.SDcols=x]
  #fisher two sided test:
  pval=fisher.test(a,
                   workspace = 2000000,
                   simulate.p.value=F,
                   alternative = "two.sided",
                   conf.int=F,
                   )$p.value
  
  #Amplified cohort:Gene2 amps from Gene1 amps  
  n1=sum(a[,1])
  n2=sum(rowSums(a)==2)
  n2/n1
  ac<-paste0(round(n2/n1,2),"(",n2,"/",n1,")")
  
  #Non-amplified cohort: Gene 2 amps from not gene1 amp 
  n1_B=nrow(a) - n1
  n2_B=sum(a[,2]==1)-n2
  
  Nac<-paste0(round(n2_B/n1_B,2),"(",n2_B,"/",n1_B,")")
  df<-data.frame(pval=pval,
                 amplified_cohort=ac,
                 Non_Amplified_cohort=Nac)
  df_list[[i]]<-df
}

S8<-do.call("rbind",df_list)
S8_all<-cbind(ss,S8)
writexl::write_xlsx(S8_all,"S8.xlsx")