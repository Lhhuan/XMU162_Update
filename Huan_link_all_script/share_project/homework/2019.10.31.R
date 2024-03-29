x<-c(0.45,0.5,1.5,0.4,0.9,0.8,1.8,0.6,1.7,0.65,0.4,2,2,2.4,3,1,2.8,1.45,1.5,1.5,0.9,0.65,1.83,2)
y<-c(0.066,0.076,0.001,0.17,0.156,0.12,0.04,0.12,0.1,0.129,0.135,0.099,0.005,0.011,0.003,0.14,0.039,0.059,0.087,0.039,0.222,0.145,0.029,0.099)
xy<-lm(y~x)
plot(x,y,xlab="风速（m/s)",ylab="NO浓度（*10^6)")
abline(xy)
summary(xy)
y.res<-residuals(xy)
shapiro.test(y.res)
x0<-data.frame(x=1.5)
lm.pred<-predict(xy,x0,interval="prediction",level=0.95)
lm.con<-predict(xy,x0,interval="confidence",level=0.95)