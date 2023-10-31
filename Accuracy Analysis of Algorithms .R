
library(GEOquery)
gds=getGEO("GDS4794")

eset=GDS2eSet(gds,do.log2 = TRUE)
print(eset)

colnames(pData(eset))
durum=pData(eset)$disease.state

levels(durum)[levels(durum)=="small cell lung cancer"]="SCLC"

kayip=nrow(which(is.na(exprs(eset)),arr.ind = TRUE))

#kayip verilerin yerine yeni degerler tahmin ederek yerlestirmek icin inpute paketinin 
#inpute.knn fonksiyonu kullanilir.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
veri=t(exprs(eset))

library(impute)
veri.imputed<-impute.knn(as.matrix(veri))
kayipyok=veri.imputed$data
dim(eset)

#filtreleme islemini yerine getirdikten sonra veri kumesi boyutlari kuculur

library(genefilter)
varFiltered=varFilter(t(kayipyok),var.cutoff = 0.9)
dim(varFiltered)

sonveri=data.frame(t(varFiltered))

#verinin bolunmesi %70 egitim %30 test
set.seed(123)
sinir=floor(.7*length(durum))
ind=sample.int(n=length(durum),size=sinir,replace = F)

veriegitim=sonveri[ind,]
sinifegitim=durum[ind]
siniftest=durum[-ind]
veritest=sonveri[-ind,]

#ayirma islemi sonucunda 45 gozlemden olusan egitim, 20 gozlemden olusan test 
#veri kumeleri olusturulmustur

dim(veriegitim)
dim(veritest)
names(train) 

learn_nb <- naiveBayes(train[,-c(1,2)], train$diagnosis) #Egitim
pre_nb <- predict(learn_nb, test[,-c(1,2)]) #Test 
cm_nb <- confusionMatrix(pre_nb, test$diagnosis) #Confusion matris 
cm_nb


learn_svm <- svm(diagnosis~., data=train[,-1])
pre_svm <- predict(learn_svm, test[,-c(1,2)])
cm_svm <- confusionMatrix(pre_svm, test$diagnosis)
cm_svm

#LOGISTIC REGRESSION 

log_model <- glm(sinifegitim ~ ., data = veriegitim, family = binomial())
summary(log_model)
show(log_model)

ongoru=predict(log_model,veritest)
con=table(ongoru,siniftest)
dogruluk=sum(diag(con))/sum(con)


#ANN

library(neuralnet)
inputs <- data[,1:4]
outputs <- data[,5]

# Verileri normalize etme
inputs <- scale(inputs)

# Yapay sinir ann modelini olu�turma
model <- neuralnet(outputs ~ inputs, data = iris, hidden = c(5, 3), linear.output = FALSE)

# Tahminleri alma
predictions <- compute(model, inputs)

# Tahminlerin sinifa donusmesi
predicted.classes <- max.col(predictions$net.result)
table(predicted.classes, outputs)

library(nnet)
ann_model <- nnet(veriegitim, sinifegitim, size=c(0.1,0.5),maxit=1)
ann_tahmin <- predict(ann_model, veritest,type="class")
ann_tahmin_sinif <- ifelse(ann_tahmin >= 0.5, "SCLC","Normal")
