L=function(x){
  q=200
  
  if(x < 0.99*q) return( 0.9*(0.99*q-x))
  if(abs(x-q) <= 0.01*q) return( 0)
  if(x > 1.01*q) return( 0.1*(x-1.01
                              *q))
  
}

x=seq(180,220,length=1000)
out <-x 
for(  i in 1:1000) out[i] = L(x[i])
plot(x,out,type="l")

library(ggplot2)
# Basic line plot with points
p=ggplot(data=data.frame(x=x,
                       out=out), aes(x=x, y=out, group=1)) +
  geom_line()+theme(
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20))+
  geom_point()+xlab( expression(hat(theta)))+ylab( expression(L(theta,hat(theta))))

ggsave(p, filename =paste0("../Figures/asym_loss.pdf"), height = 7, width = 7, bg = "transparent" )



ests<-read.csv("NBE/intermediates/estimates/estimates.csv")
p2=ggplot(data=data.frame(x=ests$truth,
                          out=ests$estimate), aes(x=x, y=out, group=1)) +theme(
                            axis.title.x = element_text(size = 20),
                            axis.text.y = element_text(size = 18),
                            axis.text.x = element_text(size = 18),
                            axis.title.y = element_text(size = 20))+
  geom_point()+xlab( "True quantile")+ylab( "Estimated quantile")+  geom_abline(intercept = 0, slope = 1)  + coord_fixed(ratio=20)

p2

ests<-read.csv("NBE/intermediates/estimates/estimates.csv")
p2=ggplot(data=data.frame(x=ests$truth,
                          out=ests$estimate), aes(x=x, y=out, group=1)) +theme(
                            axis.title.x = element_text(size = 20),
                            axis.text.y = element_text(size = 18),
                            axis.text.x = element_text(size = 18),
                            axis.title.y = element_text(size = 20))+
  geom_point()+ xlim(min(ests$truth), max(ests$truth))+ ylim(200,203)+
xlab( "True quantile")+ylab( "Estimated quantile")+  geom_abline(intercept = 0, slope = 1)  + coord_fixed(ratio=10)

p2

L=function(x,q){
 
  if(x < 0.99*q) return( 0.9*(0.99*q-x))
  if(abs(x-q) <= 0.01*q) return( 0)
  if(x > 1.01*q) return( 0.1*(x-1.01
                              *q))
  
}

temp=apply(cbind(ests$truth,ests$estimate),1,function(x) L(x[2],x[1]))
temp[temp>0]=1
p2=ggplot(data=data.frame(x=ests$truth,
                          out=ests$estimate, L=temp), aes(x=x, y=out,color=L,group=1)) +
        theme(
                            axis.title.x = element_text(size = 20),
                            axis.text.y = element_text(size = 18),
                            axis.text.x = element_text(size = 18),
                            axis.title.y = element_text(size = 20))+
  geom_point()+xlim(min(ests$truth), max(ests$truth))+ ylim(200,203)+
  xlab( "True quantile")+ylab( "Estimated quantile")+  geom_abline(intercept = 0, slope = 1)  + coord_fixed(ratio=10)+ guides(alpha = "none", color = "none")

p2
ggsave(p2, filename =paste0("../Figures/NBE_test.pdf"), height = 7, width = 7, bg = "transparent" )

boot.ests<-read.csv("NBE/intermediates/estimates/boot_estimates.csv")
quantile(as.numeric(boot.ests),probs=c(0.025,0.975))
